#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 24 11:14:48 2024

@author: saumitrarai
"""

import qsdsan as qs
import numpy as np
from qsdsan import processes as pc, sanunits as su

from qsdsan.utils import (
    time_printer
    )

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

biomass_IDs = ('X_H', 'X_AUT', 'X_PAO')

domestic_ww = {
   'S_I': 20,
   'X_I': 120,
   'S_F': 45,
   'S_A': 63,
   'X_S': 480,
   'S_NH4': 25,
   'S_PO4': 4.5,
   'X_PP': 100,
   'X_PHA': 100,
   'X_H': 100,
   'X_AUT': 500, 
   'X_PAO': 500, 
   'X_MeOH': 320, 
   'S_ALK':7*12,
    }

cmps = pc.create_asm2d_cmps(False)
cmps.compile()
qs.set_thermo(cmps)
thermo_asm2d = qs.get_thermo()

dom_ww = qs.WasteStream('domestic_wastewater', T=Temp)
dom_ww.set_flow_by_concentration(Q_domestic, 
                                concentrations=domestic_ww, 
                                units=('m3/d', 'mg/L'))
    
effluent = qs.WasteStream('effluent', T=Temp)
WAS = qs.WasteStream('WAS', T=Temp)
RAS = qs.WasteStream('RAS', T=Temp)

C1 = su.FlatBottomCircularClarifier('C1', ins=dom_ww, outs=[effluent, RAS, WAS],
                                    N_layer=10, 
                                    feed_layer=5, thermo = thermo_asm2d)
    
sys = qs.System('L1_WERF', path=(C1,), 
                recycle = [RAS])
sys.set_tolerance(rmol=1e-6)
sys.maxiter = 500

sys.set_dynamic_tracker(*sys.products)

if __name__ == '__main__':
    t = 0.001
    method = 'RK23'

# msg = f'Method {method}'
# print(f'\n{msg}\n{"-"*len(msg)}') # long live OCD!
# print(f'Time span 0-{t}d \n')    

sys.simulate(
        state_reset_hook='reset_cache',
        t_span=(0,t),
        method=method,
        )

# print("\n")    
# print("--------------------Influent wastewater properties---------------------")
# dom_ww.show()

# print("\n")    
# print("--------------------Effluent wastewater properties---------------------")
# effluent.show()

# print("\n")    
# print("---------------------------WAS properties------------------------------")
# WAS.show()

# print("\n")    
# print("---------------------------RAS properties------------------------------")
# RAS.show()

# sys.diagram()

# #%%

# # def create_components():
# #     cmps = pc.create_asm2d_cmps(False)
# #     cmps.compile()
# #     return cmps

# def create_system():
#     # Components and stream
#     # cmps = create_components()
#     # qs.set_thermo(cmps)
#     # thermo_asm2d = qs.get_thermo()
    
#     dom_ww = qs.WasteStream('domestic_wastewater', T=Temp)
#     dom_ww.set_flow_by_concentration(Q_domestic, 
#                                      concentrations=domestic_ww, 
#                                      units=('m3/d', 'mg/L'))
    
#     effluent = qs.WasteStream('effluent', T=Temp)
#     WAS = qs.WasteStream('WAS', T=Temp)
#     RAS = qs.WasteStream('RAS', T=Temp)

#     C1 = su.FlatBottomCircularClarifier('C1', dom_ww, [effluent, RAS, WAS],
#                                         N_layer=10, 
#                                         feed_layer=5, thermo = thermo_asm2d)
    
#     sys = qs.System('L1_WERF', path=(C1,), 
#                     recycle = [RAS])
#     sys.set_tolerance(rmol=1e-6)
#     sys.maxiter = 500
#     dom_ww.show()
#     effluent.show()
#     WAS.show()
#     RAS.show()
#     return sys
# #%%

# @time_printer
# def run(sys, t, method=None, **kwargs):
    
#     sys.set_dynamic_tracker(*sys.products)
    
#     msg = f'Method {method}'
#     print(f'\n{msg}\n{"-"*len(msg)}') # long live OCD!
#     print(f'Time span 0-{t}d \n')    
    
#     sys.simulate(
#         state_reset_hook='reset_cache',
#         t_span=(0,t),
#         method=method,
#         # print_t=True,
#         **kwargs)
    

#     sys.diagram()
#     # return sys
    
# if __name__ == '__main__':
#     t = 0.001
#     # method = 'RK45'
#     method = 'RK23'
#     # method = 'DOP853'
#     # method = 'Radau'
#     # method = 'BDF'
#     # method = 'LSODA'

#     sys = create_system()
#     run(sys, t, method=method)
    
#     sys.diagram()
    
# ------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------

    
#     # ACCOUNTING FOR PUMPING POWER (AND OTHER ELECTRICITY COSTS SUCH AS MOTOR IN CENTRIFUGE)
#     act_power_units = [u.ID for u in sys.units if \
#     isinstance(u, (su.FlatBottomCircularClarifier))]
        
#     required_pumping_power = get_power_utility(sys, active_unit_IDs=act_power_units)
#     print(f'Required pumping (and other equipment) power at steady state is {required_pumping_power:.2f} kW\n')
    
#     disposed_sludge = (fs.sludge_DU, )
#     sludge_disposal_costs = get_cost_sludge_disposal(sludge = disposed_sludge, unit_weight_disposal_cost = 375)
#     print(f'Sludge disposal cost = {sludge_disposal_costs} USD/day\n')

# # ------------------------------------------------------------------------------------------------------------------
# # ------------------------------------------------------------------------------------------------------------------
#     influent_ST = (fs.effluent_PC,)

#     effluent_ST = (fs.effluent, fs.WAS)
#     GHG_ST = get_GHG_emissions_sec_treatment(influent = influent_ST, effluent = effluent_ST)
#     print(f'CH4 and N2O emissions during secondary treatment equals {GHG_ST[0]} kg CH4/day and\
#       {GHG_ST[1]} kg N2O-N/day respectively\n')
    
#     effluent = (fs.effluent, )
#     GHG_discharge = get_GHG_emissions_discharge(effluent = effluent)
#     print(f'CH4 and N2O emissions at discharge equals {GHG_discharge[0]} kg CH4/day and \
#       {GHG_discharge[1]} kg N2O-N/day respectively\n')
     
#     GHG_electricity = get_GHG_emissions_electricity(system=sys, power_blower=power_blower, 
#                                                     power_pump=required_pumping_power)
#     print(f'CO2 emissions due to electricity consumption equals {GHG_electricity} kg-CO2-eq/day\n')
    
    
#     GHG_sludge_disposal = get_GHG_emissions_sludge_disposal(sludge = disposed_sludge)
#     print(f'CH4 emissions due to sludge_disposal equals {GHG_sludge_disposal} kg-CH4/day\n')
    
# # # ------------------------------------------------------------------------------------------------------------------
# # # ------------------------------------------------------------------------------------------------------------------
    
#     normalized_energy = get_normalized_energy(system = sys, aeration_power = power_blower, \
#                                               pumping_power = required_pumping_power, \
#                                               miscellaneous_power = 0)
#     print(f'Normalized energy = {normalized_energy} kWh/m3\n')
    
#     daily_operational_cost = get_daily_operational_cost(system = sys, aeration_power = power_blower, 
#                             pumping_power = required_pumping_power, miscellaneous_power = 0, 
#                             sludge_disposal_cost = sludge_disposal_costs)
#     print(f'Daily operational costs = {daily_operational_cost} USD/day\n')
    
#     total_daily_operational_cost = get_total_operational_cost(q_air = airflow, sludge = disposed_sludge, 
#                                               system = sys, active_unit_IDs= act_power_units)
#     print(f'Total daily operational costs = {total_daily_operational_cost} USD/day\n')
    
#     CO2_eq_WRRF = get_CO2_eq_WRRF(system = sys, GHG_treatment = GHG_ST, GHG_discharge = GHG_discharge, \
#                                     GHG_electricity = GHG_electricity, GHG_sludge_disposal = GHG_sludge_disposal)
#     print(f'GHG emissions = {CO2_eq_WRRF} kg CO2 eq/m3\n')
    
#     total_CO2_eq_WRRF = get_total_CO2_eq(system = sys, q_air = airflow, 
#                                     influent_sc = influent_ST, effluent_sc = effluent_ST, effluent_sys = effluent, 
#                                     active_unit_IDs=act_power_units, sludge= disposed_sludge)
#     print(f'Total GHG emissions = {total_CO2_eq_WRRF} kg CO2 eq/m3\n')
    
    
#     print(f'Airflow = {airflow} m3/min\n')
    
#     sludge_prod = np.array([sludge.composite('solids', True, particle_size='x', unit='ton/d') \
#                             for sludge in disposed_sludge]) # in ton/day
#     sludge_prod = np.sum(sludge_prod)
#     print(f'Total sludge produced = {sludge_prod} ton/day\n')
    
    
#     effluent_COD = fs.effluent.COD
#     print(f'Effluent COD = {effluent_COD} mg/L\n')
    
#     effluent_TN = fs.effluent.TN
#     print(f'Effluent TN = {effluent_TN} mg/L\n')
    
#     mass_degradable_sludge = np.array([slg.composite("C", flow=True, exclude_gas=True, subgroup=None, particle_size=None,
#                   degradability="b", organic=True, volatile=None, specification=None, unit="kg/day") for slg in disposed_sludge])
#     mass_degradable_sludge = np.sum(mass_degradable_sludge)
#     print(f'Mass degradable sludge = {mass_degradable_sludge} kg/day\n')
    
#     mass_N_sludge = np.array([slg.composite("N", flow=True, exclude_gas=True, subgroup=None, particle_size=None,
#                   degradability="b", organic=True, volatile=None, specification=None, unit="kg/day") for slg in disposed_sludge])
#     mass_N_sludge = np.sum(mass_N_sludge)
#     print(f'Mass N sludge = {mass_N_sludge} kg/day\n')
    
#     fig, axis = fs.effluent.scope.plot_time_series(('S_A', 'S_F', 'X_S', 'S_NH4', 'X_I', 'S_I', 'S_N2')) 
#     fig
    
#     # _dstate of units should be close to zero to imply steady state 
    
# # # ------------------------------------------------------------------------------------------------------------------
# # # ------------------------------------------------------------------------------------------------------------------