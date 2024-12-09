import qsdsan as qs
import numpy as np
from qsdsan import processes as pc, sanunits as su

from qsdsan.utils import (
    ospath, 
    time_printer, 
    load_data, 
    get_P_blower, get_power_utility, 
    get_cost_sludge_disposal,
    get_normalized_energy, 
    get_daily_operational_cost, 
    get_total_operational_cost,
    get_GHG_emissions_sec_treatment,
    get_GHG_emissions_discharge,
    get_GHG_emissions_electricity,
    get_GHG_emissions_sludge_disposal,
    get_CO2_eq_WRRF,
    get_total_CO2_eq,
    
    )

# from exposan.bsm1 import data_path
# from . import data_path
folder = ospath.dirname(__file__)

__all__ = (
    'biomass_IDs',
    'create_system',
    'default_asm2d_kwargs', 
    'steady_state_ad_init_conds',
    'domestic_ww', 'Q_domestic', 'Q_brewery',
    'Q_ras', 'Q_was', 'Temp', 'V_ae',
    )
# %%

# =============================================================================
# Parameters and util functions
# =============================================================================

Q_domestic = 38000      # influent flowrate [m3/d] = 10 MGD (WERF report)                        
Temp = 273.15+20 # temperature [K]
Q_was = 378.54 # [m3/day] 0.1 MGD (WERF report: L1) 
# Q_ras = 0.40*(39147) # [m3/day] 40% of effluent from PC (WERF report: L1)
# V_ae = 4982 # WERF report: L1

biomass_IDs = ('X_H', 'X_AUT', 'X_PAO')

domestic_ww = {
   'S_I': 20,
   'X_I': 40,
   'S_F': 45,
   'S_A': 63,
   'X_S': 160,
   'S_NH4': 25,
   'S_PO4': 4.5,
   'X_PP': 0,
   'X_PHA': 10,
   'X_H': 10,
   'X_AUT': 10, 
   'X_PAO': 5, 
   'X_MeOH': 32, 
   'S_ALK':7*12,
    }

default_asm2d_kwargs = dict(iN_SI=0.01, iN_SF=0.03, iN_XI=0.02, iN_XS=0.04, iN_BM=0.07,
            iP_SI=0.0, iP_SF=0.01, iP_XI=0.01, iP_XS=0.01, iP_BM=0.02,
            iTSS_XI=0.75, iTSS_XS=0.75, iTSS_BM=0.9,
            f_SI=0.0, Y_H=0.625, f_XI_H=0.1,
            Y_PAO=0.625, Y_PO4=0.4, Y_PHA=0.2, f_XI_PAO=0.1,
            Y_A=0.24, f_XI_AUT=0.1,
            K_h=3.0, eta_NO3=0.6, eta_fe=0.4, K_O2=0.2, K_NO3=0.5, K_X=0.1,
            mu_H=6.0, q_fe=3.0, eta_NO3_H=0.8, b_H=0.4, K_O2_H=0.2, K_F=4.0,
            K_fe=4.0, K_A_H=4.0, K_NO3_H=0.5, K_NH4_H=0.05, K_P_H=0.01, K_ALK_H=0.1,
            q_PHA=3.0, q_PP=1.5, mu_PAO=1.0, eta_NO3_PAO=0.6, b_PAO=0.2, b_PP=0.2,
            b_PHA=0.2, K_O2_PAO=0.2, K_NO3_PAO=0.5, K_A_PAO=4.0, K_NH4_PAO=0.05,
            K_PS=0.2, K_P_PAO=0.01, K_ALK_PAO=0.1,
            K_PP=0.01, K_MAX=0.34, K_IPP=0.02, K_PHA=0.01,
            mu_AUT=1.0, b_AUT=0.15, K_O2_AUT=0.5, K_NH4_AUT=1.0, K_ALK_AUT=0.5, K_P_AUT=0.01,
            k_PRE=1.0, k_RED=0.6, K_ALK_PRE=0.5, 
            # path=os.path.join(data_path, '_asm2d.tsv'),
            )

# Using effluent DG sludge concentrations from steady state sys consisting of only AD

steady_state_ad_init_conds = {
      'S_su': 1.094e+01,
       # 'S_aa': 4.713e+00, # original ss
      # 'S_aa': 4.713e+01,
       'S_aa': 4.713e+03,
      'S_fa': 9.150e+01,
      'S_va': 9.512e+00,
      'S_bu': 1.255e+01,
      'S_pro': 1.486e+01,
      'S_ac': 3.749e+01,
       'S_h2': 2.193e-04,
      'S_ch4': 7.692e+01,
      'S_IC': 6.708e+02,
      'S_IN': 1.035e+03,
      # --
      'S_IP': 1.035e+00, #random
      # --
      'S_I': 2.325e+02,
      # --
      # 'X_c': 2.171e+02,
      # --
       # 'X_ch': 6.432e+01, # original ss
      # 'X_ch': 6.432e+00, 
       'X_ch': 6.432e-01,
      'X_pr': 6.719e+01,
       # 'X_li': 1.436e+02, # original ss
      # 'X_li': 1.436e+01,
       'X_li': 1.436e+00,
      'X_su': 1.116e+03,
      'X_aa': 8.755e+02,
      'X_fa': 1.274e+03,
      'X_c4': 3.730e+02,
      'X_pro': 1.751e+02,
      'X_ac': 1.405e+03,
      'X_h2': 6.820e+02,
      'X_I': 2.413e+04,
      # --
      'X_PHA': 10.000e-01, 
      'X_PP':  10.000e-01,  
      
      'X_PAO': 10.000e+00, 
      
      'S_K': 10.000e-001, 
      'S_Mg': 10.000e-001, 
      'X_MeOH': 10.000e-001, 
      'X_MeP': 10.000e-001, 
      # --
       'S_cat': 0.000e+00, 
       'S_an': 4.772e+00
      }

def batch_init(sys, path, sheet):
    # isa = isinstance
    df = load_data(path, sheet)
    dct = df.to_dict('index')
    u = sys.flowsheet.unit # unit registry
    u.DG.set_init_conc(**steady_state_ad_init_conds)
#%%

def create_components():
    cmps = pc.create_asm2d_cmps(False)
    # cmps.X_I.i_N = 0.0600327162
    # cmps.X_I.i_P = 0.01
    # cmps.S_I.i_N = 0.0600327162
    # cmps.S_I.i_P = 0.01
    # cmps.X_PAO.i_N = 0.07
    # cmps.refresh_constants()
    cmps.compile()
    return cmps

def create_system(flowsheet=None, inf_kwargs={}, asm_kwargs={}, init_conds={},
                  aeration_processes=()):
    # Create ASM2d components and set initial thermo
    cmps_asm2d = create_components()
    qs.set_thermo(cmps_asm2d)
    thermo_asm2d = qs.get_thermo()
    
    # Create streams with ASM2d thermo
    global dom_ww, effluent, sludge_PC, sludge_DG, gas, PC, J1, DG
    dom_ww = qs.WasteStream('domestic_wastewater', T=Temp)
    dom_ww.set_flow_by_concentration(Q_domestic, 
                                   concentrations=domestic_ww, 
                                   units=('m3/d', 'mg/L'))
    
    effluent = qs.WasteStream('effluent', T=Temp)
    sludge_PC = qs.WasteStream('sludge_PC', T=Temp)
    
    # Create Primary Clarifier with ASM2d thermo
    asm_kwargs = asm_kwargs or default_asm2d_kwargs
    asm2d = pc.ASM2d()
    
    PC = su.PrimaryClarifier('PC', ins=[dom_ww],
                            outs=(effluent, sludge_PC), 
                            isdynamic=True,
                            init_with='WasteStream', 
                            thermo=thermo_asm2d,
                            sludge_flow_rate=280, 
                            solids_removal_efficiency=0.6)
    
    # Create ADM1 components and streams
    cmps_adm1 = qs.processes.create_adm1_p_extension_cmps()
    cmps_adm1.X_PAO.i_N = 0.07
    cmps_adm1.X_PAO.i_P = 0.02
    cmps_adm1.refresh_constants()
    qs.set_thermo(cmps_adm1)
    thermo_adm1 = qs.get_thermo()
    
    # Create ADM1 streams with ADM1 thermo
    sludge_DG = qs.WasteStream('sludge_DG', T=Temp, thermo=thermo_adm1)
    gas = qs.WasteStream('gas', thermo=thermo_adm1)
    
    adm1 = qs.processes.ADM1_p_extension()
    
    # Create junction and digester with ADM1 thermo
    J1 = su.ASM2dtomADM1('J1', upstream=[sludge_PC], 
                         thermo=thermo_adm1, 
                         isdynamic=True, 
                         adm1_model=adm1, 
                         asm2d_model=asm2d)
    
    DG = su.AnaerobicCSTR(ID='DG', 
                          ins=J1.outs[0], 
                          outs=[gas, sludge_DG], 
                          model=adm1, 
                          thermo=thermo_adm1, 
                          V_liq=4542, 
                          V_gas=454)
    DG.algebraic_h2 = True
    
    # Create system
    sys = qs.System('L1_WERF', path=(PC, J1, DG))
    sys.set_tolerance(rmol=1e-6)
    sys.maxiter = 500

    return sys

@time_printer
def run(sys, t, method=None, **kwargs):
    
    batch_init(sys, "C:/Users/Aravi/Downloads/initial_conditions_ASM2d.xlsx", sheet='L1')
    
    sys.set_dynamic_tracker(*sys.products)
    
    msg = f'Method {method}'
    print(f'\n{msg}\n{"-"*len(msg)}') # long live OCD!
    print(f'Time span 0-{t}d \n')    
    
    sys.simulate(
        state_reset_hook='reset_cache',
        t_span=(0,t),
        method=method,
        # print_t=True,
        **kwargs)
    
    # print('-------------------Influent-------------------')
    # dom_ww.show()

    # print('-------------------Effluent-------------------')
    # effluent.show()

    # print('-------------------Clarifier Sludge-------------------')
    # sludge_PC.show()

    print('--------------------Intermediate--------------------')
    J1.outs[0].show()
    J1.outs[1].show()

    # print('-------------------Digester Sludge-------------------')
    # sludge_DG.show()

    # print('-------------------Biogas-------------------')
    # gas.show()

    # sys.diagram()
    # return sys
    
if __name__ == '__main__':
    t = 0.001
    # method = 'RK45'
    method = 'RK23'
    # method = 'DOP853'
    # method = 'Radau'
    # method = 'BDF'
    # method = 'LSODA'

    sys = create_system()
    fs = sys.flowsheet.stream
    fu = sys.flowsheet.unit
    
    # sys._setup()
    # sys.converge()
    run(sys, t, method=method)  