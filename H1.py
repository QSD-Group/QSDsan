# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
        
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.

'''
import qsdsan as qs
from qsdsan import (
    WasteStream,
    processes as pc,
    sanunits as su,
    )
from qsdsan.utils import ospath, time_printer, load_data, get_SRT
from exposan.werf import data_path#, default_fctss_init

__all__ = ('create_h1_system',)

ID = 'H1'
#%%
dfs = load_data(
    ospath.join(data_path, 'initial_conditions.xlsx'), 
    sheet=None,
    )
asinit = dfs['H']
fcinit = asinit.iloc[-1].to_dict()
default_fctss_init = [8, 9, 9, 30, 350, 350, 350, 350, 350, 9000]
adinit = dfs['adm'].loc['H'].to_dict()

MGD2cmd = 3785.412
Temp = 273.15+20 # temperature [K]
T_ad = 273.15+35

def create_h1_system(flowsheet=None, default_init_conds=True):
    flowsheet = flowsheet or qs.Flowsheet(ID)
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    pc.create_masm2d_cmps()
    asm = pc.mASM2d(electron_acceptor_dependent_decay=True)
    rww = pc.create_masm2d_inf(
        'RWW', 10, 'MGD', T=Temp, 
        COD=358, NH4_N=25.91, PO4_P=5,
        fr_SI=0.05, fr_SF=0.16, fr_SA=0.024, fr_XI=0.2,
        )
    carb = WasteStream(
        'carbon', T=Temp, units='kg/hr', 
        S_A=80,
        # S_A=24.5,      # insufficient for denitrification
        )
    thermo_asm = qs.get_thermo()
    
    MD = su.MetalDosage.from_mASM2d(
        'MD', ins=[rww, 'reject'], model=asm,
        metal_ID='X_AlOH', metal_dosage=10
        )
    PC = su.PrimaryClarifier(
        'PC', ins=MD-0, 
        outs=('PE', 'PS'),
        isdyanmic=True,
        sludge_flow_rate=0.074*MGD2cmd,
        solids_removal_efficiency=0.6
        )
    
    GT = su.IdealClarifier(
        'GT', PC-1, outs=['', 'thickened_PS'],
        sludge_flow_rate=0.026*MGD2cmd,
        solids_removal_efficiency=0.9
        )
    
    n_zones = 6
    V_tot = 3.66 * MGD2cmd
    fr_V = [0.115, 0.115, 0.314, 0.314, 0.115, 0.027]
    
    # ae_kwargs = dict(aeration=2.0, DO_ID='S_O2', suspended_growth_model=asm, gas_stripping=True)
    # an_kwargs = dict(aeration=None, DO_ID='S_O2', suspended_growth_model=asm, gas_stripping=True)
        
    # A1 = su.CSTR('A1', [rww, 'RAS', carb, 'intr'], V_max=fr_V[0]*V_tot, **an_kwargs)
    # A2 = su.CSTR('A2', A1-0, V_max=fr_V[1]*V_tot, **an_kwargs)
    # O3 = su.CSTR('O3', A2-0, V_max=fr_V[2]*V_tot, **ae_kwargs)
    # O4 = su.CSTR('O4', O3-0, ('', 3-A1), split=[30, 14], V_max=fr_V[3]*V_tot, **ae_kwargs)
    # A5 = su.CSTR('A5', O4-0, V_max=fr_V[4]*V_tot, **an_kwargs)
    # O6 = su.CSTR('A6', A5-0, 'treated', V_max=fr_V[5]*V_tot, **ae_kwargs)
    
    ASR = su.PFR(
        'ASR', ins=[PC-0, 'RAS', carb], outs='treated',
        N_tanks_in_series=n_zones,
        V_tanks=[f*V_tot for f in fr_V],
        influent_fractions=[
            [1,0,0,0,0,0],          # PE
            [1,0,0,0,0,0],          # RAS
            [1,0,0,0,0,0],          # carb
            ],
        internal_recycles=[(3,0,30*MGD2cmd)], DO_ID='S_O2',
        kLa=[0, 0, 50, 50, 0, 50], 
        DO_setpoints=[0, 0, 2.0, 2.0, 0, 2.0],
        suspended_growth_model=asm,
        gas_stripping=True
        )
    
    FC = su.FlatBottomCircularClarifier(
        'FC', ins=ASR-0, outs=['SE', 1-ASR, 'WAS'],
        # 'FC', ins=O6-0, outs=['SE', 1-A1, 'WAS'],
        underflow=0.4*10*MGD2cmd, wastage=0.1*MGD2cmd,
        surface_area=1579.352, height=3.6576, N_layer=10, feed_layer=5,
        X_threshold=3000, v_max=410, v_max_practical=274,
        rh=4e-4, rp=0.1, fns=0.01, 
        maximum_nonsettleable_solids=8.0
        )

    MT = su.IdealClarifier(
        'MT', FC-2, outs=['', 'thickened_WAS'],
        sludge_flow_rate=0.017*MGD2cmd,
        solids_removal_efficiency=0.95
        )
    M1 = su.Mixer('M1', ins=[GT-1, MT-1])

    pc.create_adm1p_cmps()
    thermo_adm = qs.get_thermo()
    adm = pc.ADM1p(kLa=10.0)
    
    J1 = su.mASM2dtoADM1p('J1', upstream=M1-0, thermo=thermo_adm, isdynamic=True, 
                          adm1_model=adm, asm2d_model=asm)
    AD = su.AnaerobicCSTR(
        'AD', ins=J1-0, outs=('biogas', 'digestate'), 
        V_liq=1.1*MGD2cmd, V_gas=0.11*MGD2cmd, 
        # V_liq=0.85*MGD2cmd, V_gas=0.085*MGD2cmd, 
        fixed_headspace_P=False, fraction_retain=0,
        T=T_ad, model=adm,
        pH_ctrl=7.0,
        )
    AD.algebraic_h2 = False
    J2 = su.ADM1ptomASM2d('J2', upstream=AD-1, thermo=thermo_asm, isdynamic=True, 
                          adm1_model=adm, asm2d_model=asm)
    qs.set_thermo(thermo_asm)
     
    DW = su.IdealClarifier(
        'DW', J2-0, outs=('', 'cake'),
        sludge_flow_rate=0.00593*MGD2cmd,   # aim for 18% TS
        solids_removal_efficiency=0.9
        )
    M2 = su.Mixer('M2', ins=[GT-0, MT-0, DW-0])
    
    HD = su.HydraulicDelay('HD', ins=M2-0, outs=1-MD)
    
    if default_init_conds:
        # asdct = asinit.to_dict('index')
        # for i in (A1, A2, O3, O4, A5, O6):
        #     i.set_init_conc(**asdct[i.ID])
        ASR.set_init_conc(concentrations=asinit)
        FC.set_init_solubles(**fcinit)
        FC.set_init_sludge_solids(**fcinit)
        FC.set_init_TSS(default_fctss_init)
        AD.set_init_conc(**adinit)

    sys = qs.System(
        ID, 
        path=(MD, PC, GT, ASR, FC, MT, M1, J1, AD, J2, DW, M2, HD),
        # path=(MD, PC, GT, A1, A2, O3, O4, A5, O6, FC, 
        #       MT, M1, J1, AD, J2, DW, M2, HD),
        recycle=(FC-1, HD-0)
        )

    sys.set_dynamic_tracker(FC-0, AD)

    return sys

#%%
@time_printer
def run(sys, t, t_step, method=None, **kwargs):
    msg = f'Method {method}'
    print(f'\n{msg}\n{"-"*len(msg)}')
    print(f'Time span 0-{t}d \n')
    
    sys.simulate(
        # state_reset_hook='reset_cache',
        t_span=(0,t),
        # t_eval=np.arange(0, t+t_step, t_step),
        method=method,
        print_t=True,
        # rtol=1e-2,
        # atol=1e-3,
        # export_state_to=f'results/sol_{t}d_{method}.xlsx',
        **kwargs)
    

#%%
if __name__ == '__main__':
    sys = create_h1_system()
    dct = globals()
    dct.update(sys.flowsheet.to_dict())
    
    t = 300
    t_step = 1
    # method = 'RK45'
    # method = 'RK23'
    # method = 'DOP853'
    # method = 'Radau'
    method = 'BDF'
    # method = 'LSODA'
    
    run(sys, t, t_step, method=method)
    # biomass_IDs = ('X_H', 'X_PAO', 'X_AUT')
    # srt = get_SRT(sys, biomass_IDs,
    #               wastage=[WAS],
    #               # active_unit_IDs=('A1', 'A2', 'O3', 'O4', 'A5', 'O6'))
    #               active_unit_IDs=('ASR'))
    # if srt: print(f'Estimated SRT assuming at steady state is {round(srt, 2)} days')
    
    # from exposan.werf import figures_path
    # sys.diagram(format='png', file=ospath.join(figures_path, f'{ID}'))