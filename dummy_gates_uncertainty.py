# -*- coding: utf-8 -*-
'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    
    Saumitra Rai <raisaumitra9@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

import os, qsdsan as qs
import numpy as np
from qsdsan import processes as pc, sanunits as su, Model as mod
from chaospy import distributions as shape #chaospy is a library for performing uncertainty quantification

from qsdsan.utils import (ospath, time_printer, estimate_ww_treatment_energy_demand,
                          estimate_N_removal_energy_demand, estimate_P_removal_energy_demand)
from exposan.bsm1 import data_path
folder = ospath.dirname(__file__)


__all__ = ('default_asm1_kwargs', 'domestic_ww', 'Q_domestic', 'Temp',)

# %%

Q_domestic = 38000      # influent flowrate [m3/d] = 10 MGD (WERF report)

Temp = 273.15+20 # temperature [K]

# Based on WERF report
domestic_ww = {
  'S_I': 10, 
  'S_S': 68,  
  'X_I': 28, 
  'X_S': 242,
  'X_BH': 5, 
  'X_BA': 5, 
  'S_NH': 25,  
  'S_ND': 12,
  'S_ALK':7*12,
    }

default_asm1_kwargs = dict(
    Y_A=0.24, Y_H=0.67, f_P=0.08, i_XB=0.08, i_XP=0.06,
    mu_H=4.0, K_S=10.0, K_O_H=0.2, K_NO=0.5, b_H=0.3,
    eta_g=0.8, eta_h=0.8, k_h=3.0, K_X=0.1, mu_A=0.5,
    K_NH=1.0, b_A=0.05, K_O_A=0.4, k_a=0.05, fr_SS_COD=0.75,
    path=os.path.join(data_path, '_asm1.tsv'),
    )

#%%

def create_components():
     asm1_cmps = pc.create_asm1_cmps(False)
     cmps = qs.Components([*asm1_cmps])
     cmps.compile()
     return cmps

def create_system(flowsheet=None, inf_kwargs={}, asm_kwargs={}, init_conds={}, lifetime=100, discount_rate=0.1,
                  aeration_processes=()):

    # Components and stream
    cmps = create_components()
    qs.set_thermo(cmps)
    thermo_asm1 = qs.get_thermo()
    
    dom_ww = qs.WasteStream('domestic_wastewater', T=Temp)
    
    # inf_kwargs = inf_kwargs or default_inf_kwargs
    dom_ww.set_flow_by_concentration(Q_domestic, 
                                     concentrations=domestic_ww, 
                                     units=('m3/d', 'mg/L'))

    effluent = qs.WasteStream('effluent', T=Temp)
    
    # asm_kwargs = asm_kwargs or default_asm1_kwargs
    # asm1 = pc.ASM1(**asm_kwargs)
    
    Mixer = su.Mixer('Mixer', dom_ww, effluent)
    
    Mixer.energy_consumption = 49631
    Mixer.ww_per_capita = 0.175
    Mixer.protein_intake = 68.6
    Mixer.N_pro = 0.13
    Mixer.N_excreted = 1
    
    Mixer.ani_protein_intake = 12.39,
    Mixer.plant_protein_intake = 40.29, 
    Mixer.P_ani_pro = 0.011
    Mixer.P_plant_pro = 0.022
    Mixer.P_excreted = 1
    
    sys = qs.System('Gates_WERF', path=(Mixer, ))
        
    return sys
#%%

@time_printer
def run(t, method=None, **kwargs):
    sys = create_system()    
    
    # sys.simulate(
    #     state_reset_hook='reset_cache',
    #     t_span=(0,t),
    #     method=method,
    #     # print_t=True,
    #     **kwargs)
    
    sys._setup()
    sys.converge() 
    
    # sys.diagram()
    return sys
    
if __name__ == '__main__':
    t = 5
    method = 'RK23' 
    msg = f'Method {method}'
    print(f'\n{msg}\n{"-"*len(msg)}') # long live OCD!
    print(f'Time span 0-{t}d \n')
    sys = run(t, method=method)
    
    sys.diagram()
    fs = sys.flowsheet.stream
    fu = sys.flowsheet.unit
    
# ------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------
    Gates_model = mod(sys)
    
    param = Gates_model.parameter
    
    # daily_energy_demand : float
    #     Energy consumption at a centralized WRRF. (kWh/day)
    
    # BOD removal: min = 4200, max = 6737, median = 4323
    # dist_ed = shape.Triangle(lower= 4200, midpoint= 4323, upper= 6737)
    
    # Nitrification: min = 7036, max = 9439, median = 7697
    # dist_ed = shape.Triangle(lower= 7036, midpoint= 7697, upper= 9439)
    
    # BNR: min = 14538, max = 35738, median = 19834
    # dist_ed = shape.Triangle(lower= 14538, midpoint= 19834, upper= 35738)
    
    # ENR: min = 32125, max = 44581, median = 38353
    dist_ed = shape.Uniform(lower= 32125, upper= 44581)
    
    @param(name = 'Energy consumption', 
            element = fu.Mixer, 
            kind = 'coupled', 
            units = 'kWH/day', 
            baseline = fu.Mixer.energy_consumption, # mean of upper and lower limits
            distribution = dist_ed )
    def set_ed(i):
        fu.Mixer.energy_consumption = i
        
    # From "Jones, E. R., Van Vliet, M. T., Qadir, M., & Bierkens, M. F. (2021). Country-level and 
    # gridded estimates of wastewater production, collection, treatment and reuse. Earth System Science Data, 
    # 13(2), 237-254."
    dist_ww_per_capita = shape.Triangle(lower= 0.03, midpoint= 0.175, upper= 0.5739)
    @param(name = 'Wastewater generated pcpd', 
            element = fu.Mixer, 
            kind = 'coupled', 
            units = 'm3/cap/day', 
            baseline = fu.Mixer.ww_per_capita, # mean of upper and lower limits
            distribution = dist_ww_per_capita)
    def set_ww_pcpd(i):
        fu.Mixer.ww_per_capita = i
    
    # per_capita_protein_intake : float
    #     Per capita protein intake. (g/cap/day)
    dist_protein_intake = shape.Uniform(lower = 45.2, upper = 92)
    @param(name = 'Protein intake', 
            element = fu.Mixer, 
            kind = 'coupled', 
            units = 'g/cap/day', 
            baseline = fu.Mixer.protein_intake, # mean of upper and lower limits
            distribution = dist_protein_intake)
    def set_pi(i):
        fu.Mixer.protein_intake = i
        
    
    # N_in_pro : float
    #     % of N in protein. (%)
    dist_N_in_pro = shape.Uniform(lower = 0.13, upper = 0.19)
    @param(name = 'N content in protein', 
            element = fu.Mixer, 
            kind = 'coupled', 
            units = 'percentage', 
            baseline = fu.Mixer.N_pro, # mean of upper and lower limits
            distribution = dist_N_in_pro)
    def set_N_pro(i):
        fu.Mixer.N_pro = i
    
    # N_excreted : float
    #     % of N intake that is excreted. (%)
    dist_N_excreta = shape.Uniform(lower = 0.99, upper = 1)
    @param(name = 'N excreted', 
            element = fu.Mixer, 
            kind = 'coupled', 
            units = 'percentage', 
            baseline = fu.Mixer.N_excreted, # mean of upper and lower limits
            distribution = dist_N_excreta)
    def set_N_excreted(i):
        fu.Mixer.N_excreted = i
    
    # 36.26 - 44.32 = vegetable
    # 11.15 - 13.63 = animal
    
    # pc_ani_protein_intake : float, optional
    #     DESCRIPTION. The default is 12.39.
    dist_pc_ani_protein_intake = shape.Uniform(lower = 11.15, upper = 13.63)
    @param(name = 'Animal protein intake', 
            element = fu.Mixer, 
            kind = 'coupled', 
            units = 'g/cap/day', 
            baseline = fu.Mixer.ani_protein_intake, # mean of upper and lower limits
            distribution = dist_pc_ani_protein_intake)
    def set_pc_ani_protein_intake(i):
        fu.Mixer.ani_protein_intake = i
    
    # pc_plant_protein_intake : TYPE, optional
    #     DESCRIPTION. The default is 40.29.
    dist_pc_plant_protein_intake = shape.Uniform(lower = 36.26, upper = 44.32)
    @param(name = 'Plant protein intake', 
            element = fu.Mixer, 
            kind = 'coupled', 
            units = 'g/cap/day', 
            baseline = fu.Mixer.plant_protein_intake, # mean of upper and lower limits
            distribution = dist_pc_plant_protein_intake)
    def set_pc_plant_protein_intake(i):
        fu.Mixer.plant_protein_intake = i
    
    # P_ani_pro : TYPE, optional
    #     DESCRIPTION. The default is 0.011.
    dist_P_ani_pro = shape.Triangle(lower= 0.002, midpoint= 0.011, upper=0.032)
    @param(name = 'P content in animal protein', 
            element = fu.Mixer, 
            kind = 'coupled', 
            units = 'percentage', 
            baseline = fu.Mixer.P_ani_pro, # mean of upper and lower limits
            distribution = dist_P_ani_pro)
    def set_P_ani_pro(i):
        fu.Mixer.P_ani_pro = i
    
    # P_plant_pro : TYPE, optional
    #     DESCRIPTION. The default is 0.022.
    dist_P_plant_pro = shape.Triangle(lower= 0.004, midpoint= 0.022, upper=0.048)
    @param(name = 'P content in plant protein', 
            element = fu.Mixer, 
            kind = 'coupled', 
            units = 'percentage', 
            baseline = fu.Mixer.P_plant_pro, # mean of upper and lower limits
            distribution = dist_P_plant_pro)
    def set_P_plant_pro(i):
        fu.Mixer.P_plant_pro = i
    
    # P_excreted : TYPE, optional
    #     DESCRIPTION. The default is 1.
    dist_P_excreta = shape.Uniform(lower = 0.99, upper = 1)
    @param(name = 'P excreted', 
            element = fu.Mixer, 
            kind = 'coupled', 
            units = 'percentage', 
            baseline = fu.Mixer.P_excreted, # mean of upper and lower limits
            distribution = dist_P_excreta)
    def set_P_excreted(i):
        fu.Mixer.P_excreted = i
    
    metric = Gates_model.metric
    
    @metric(name='Energy consumed in ww treatment', units='kWh/day', element= fu.Mixer)
    def energy_ww_treatment():
        return estimate_ww_treatment_energy_demand(daily_energy_demand = fu.Mixer.energy_consumption,  
                                                    daily_flow = 37854, 
                                                    ww_pcpd = fu.Mixer.ww_per_capita)
    
    @metric(name='Energy consumed in N removal', units='kWh/day', element= fu.Mixer)
    def energy_N_removal():
        return estimate_N_removal_energy_demand(daily_energy_demand = fu.Mixer.energy_consumption, 
                                          effluent_N_conc = 5, 
                                          daily_flow = 37854, 
                                          influent_N_conc = 40, 
                                          per_capita_protein_intake = fu.Mixer.protein_intake, 
                                          N_in_pro = fu.Mixer.N_pro, 
                                          N_excreted = fu.Mixer.N_excreted)
    
    @metric(name='Energy consumed in P removal', units='kWh/day', element= fu.Mixer)
    def energy_P_removal():
        return estimate_P_removal_energy_demand(daily_energy_demand = fu.Mixer.energy_consumption, 
                                                effluent_TP_conc = 1, 
                                                daily_flow = 37854, 
                                                influent_TP_conc = 7, 
                                                
                                                pc_ani_protein_intake = fu.Mixer.ani_protein_intake, 
                                                pc_plant_protein_intake = fu.Mixer.plant_protein_intake, 
                                                P_ani_pro = fu.Mixer.P_ani_pro, 
                                                P_plant_pro = fu.Mixer.P_plant_pro, 
                                                P_excreted = fu.Mixer.P_excreted
                                                )
        				
    np.random.seed(3221) # setting the seed ensures you getting the same sample
    
    #seed fixes the sample
    samples = Gates_model.sample(N=5000, rule='L')
    Gates_model.load_samples(samples)
    
    sys.isdynamic = False
    Gates_model.evaluate()
    
    all_outputs = []
    for sample in samples:
        for p, v in zip(Gates_model.parameters, sample):
            p.setter(v)
    
        out = [m() for m in Gates_model.metrics]
        all_outputs.append(out)
    
    Gates_model.table.iloc[:,-3:] = all_outputs
    
    Gates_model.table

    import pandas as pd
    file_path = r"C:\Users\Aravi\OneDrive\Desktop\workspace\Book1.xlsx"
    sheet_name = "ENR new"  
    with pd.ExcelWriter(file_path, engine='openpyxl', mode='a') as writer:
        Gates_model.table.to_excel(writer, sheet_name=sheet_name)
    
    # fig, ax = qs.stats.plot_uncertainties(Gates_model, x_axis=Gates_model.metrics[0], y_axis=Gates_model.metrics[1],
    #                                   kind='kde-kde', center_kws={'fill': True})
    # fig

    r_df, p_df = qs.stats.get_correlations(Gates_model, kind='Spearman')
    print(r_df, p_df)
    fig, ax = qs.stats.plot_correlations(r_df)
    fig