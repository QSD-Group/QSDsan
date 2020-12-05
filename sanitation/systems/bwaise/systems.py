#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Sanitation Explorer: Sustainable design of non-sewered sanitation technologies
Copyright (C) 2020, Sanitation Explorer Development Group

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the UIUC open-source license. Please refer to 
https://github.com/QSD-for-WaSH/sanitation/blob/master/LICENSE.txt
for license details.

Ref:
    [1] Trimmer et al., Navigating Multidimensional Social–Ecological System
        Trade-Offs across Sanitation Alternatives in an Urban Informal Settlement.
        Environ. Sci. Technol. 2020, 54 (19), 12641–12653.
        https://doi.org/10.1021/acs.est.0c03296.

TODO:
    [1] Recheck unit consistency, power law degradation, truck

'''


# %%

import numpy as np
import biosteam as bst
import sanitation # import sanplorer as sp
# from sklearn.linear_model import LinearRegression
from sanitation import units, WasteStream, ImpactItem, StreamImpactItem, \
    SimpleTEA, LCA
import bwaise
cmps = bwaise._cmps.cmps

bst.settings.set_thermo(cmps)
e = bst.PowerUtility
e.price = 0.17
currency = sanitation.currency = 'USD'
sanitation.CEPCI = sanitation.CEPCI_by_year[2018]
items = ImpactItem._items
GWP = sanitation.ImpactIndicator._indicators['GWP']


# %%

# =============================================================================
# Assumptions
# =============================================================================

toilet_user = 16 # four people per household, four households per toilet
ppl_existing = 4e4 # number of people served by the existing plant (SceA and SceC)
ppl_alternative = 5e4 # number of people served by the alternative plant (SceB)

exchange_rate = 3700 # UGX per USD, triangular of 3600, 3700, 3900
discount_rate = 0.05 # uniform of 0.03-0.06

# Time take for full degradation, [yr]
tau_deg = 2
# Log reduction at full degradation
log_deg = 3
# Get reduction rate constant k for COD and N, use a function so that k can be
# changed during uncertainty analysis
def get_decay_k(tau_deg=2, log_deg=3):
    k = (-1/tau_deg)*np.log(10**-log_deg)
    return k

max_CH4_emission = 0.25

#!!! How this was selected? 
truck_cost = {
    'TankerTruck1': 8e4/exchange_rate*1.15, # 15% additional, per m3
    'TankerTruck2': 12e4/exchange_rate*1.15, # 15% additional, per m3
    'TankerTruck3': 2e5/exchange_rate*1.15, # 15% additional, per m3
    'TankerTruck4': 25e4/exchange_rate*1.15, # 15% additional, per m3
    'HandCart': 0.01, # per cap/d
    'CBSTruck': 23e3/exchange_rate # per m3
    }

# Assume density is 1 tonne/m3 (as water)
V = (3, 4.5, 8, 15, 1, 1)
truck_V = dict.fromkeys(truck_cost.keys())
for i, j in zip(truck_V.keys(), V):
    truck_V[i] = j

# Trucking = items['Trucking']
# for i, j in truck_cost.items():
#     new = Trucking.copy()
#     new.price = j
#     new.functional_unit = 'm3'
#     setattr(new, 'ID', i)

items['Concrete'].price = 194
items['Steel'].price = 2.665

N_AD_rx = 3


# Nutrient loss during applciation
app_loss = dict.fromkeys(('NH3', 'NonNH3', 'P', 'K', 'Mg', 'Ca'), 0.02)
app_loss['NH3'] = 0.05





# %%

# =============================================================================
# Scenario A: existing system
# =============================================================================

fugitive_CH4 = WasteStream('fugitive_CH4', phase='g')
fugitive_N2O = WasteStream('fugitive_N2O', phase='g')
# The product is actually liquid fertilizer, but separated here to be solids
# just to avoid showing WasteStream related properties
liq_N = WasteStream('liq_N', phase='s', price=1.507)
sol_N = WasteStream('sol_N', phase='s', price=1.507)
liq_P = WasteStream('liq_P', phase='s', price=3.983)
sol_P = WasteStream('sol_P', phase='s', price=3.983)
liq_K = WasteStream('liq_K', phase='s', price=1.333)
sol_K = WasteStream('sol_K', phase='s', price=1.333)

print('\n----------Scenario A----------\n')
A1 = units.Excretion('A1', outs=('urine', 'feces'))

A2 = units.PitLatrine('A2', ins=(A1-0, A1-1,
                                  'toilet_paper', 'flushing_water',
                                  'cleansing_water', 'desiccant'),
                      outs=('mixed_waste', 'leachate', '', ''),
                      N_user=toilet_user, N_toilet=ppl_existing/toilet_user,
                      OPEX_over_CAPEX=0.05,
                      decay_k_COD=get_decay_k(tau_deg, log_deg),
                      decay_k_N=get_decay_k(tau_deg, log_deg),
                      max_CH4_emission=max_CH4_emission)

truck = 'TankerTruck1' # assumed
interval = (A2.emptying_period*365*truck_V[truck])/A2.pit_V
A3 = units.Trucking('A3', ins=A2-0, outs=('transported', 'loss'),
                    load_type='mass', load=truck_V[truck], load_unit='tonne',
                    distance=5, distance_unit='km',
                    interval=interval, interval_unit='day',
                    fee=truck_cost[truck],
                    loss_ratio=0.02)

A4 = units.SedimentationTank('A4', ins=A3-0,
                              outs=('liq', 'sol', '', ''),
                              decay_k_COD=get_decay_k(tau_deg, log_deg),
                              decay_k_N=get_decay_k(tau_deg, log_deg),
                              max_CH4_emission=max_CH4_emission)

A5 = units.Lagoon('A5', ins=A4-0, outs=('anaerobic_treated', '', ''),
                  design_type='anaerobic',
                  decay_k_N=get_decay_k(tau_deg, log_deg),
                  max_CH4_emission=max_CH4_emission)

A6 = units.Lagoon('A6', ins=A5-0, outs=('facultative_treated', '', ''),
                  design_type='facultative',
                  decay_k_N=get_decay_k(tau_deg, log_deg),
                  max_CH4_emission=max_CH4_emission)

A7 = units.DryingBed('A7', ins=A4-1, outs=('dried_sludge', 'evaporated', '', ''),
                     design_type='unplanted',
                     decay_k_COD=get_decay_k(tau_deg, log_deg),
                     decay_k_N=get_decay_k(tau_deg, log_deg),
                     max_CH4_emission=max_CH4_emission)

A8 = units.CropApplication('A8', ins=A6-0, outs=('liquid_fertilizer', 'loss'),
                           loss_ratio=app_loss)

A9 = units.Mixer('A9', ins=(A2-2, A4-2, A5-1, A6-1, A7-2),
                 outs=fugitive_CH4)
A9.line = 'CH4 mixer'

A10 = units.Mixer('A10', ins=(A2-3, A4-3, A5-2, A6-2, A7-3),
                 outs=fugitive_N2O)
A9.line = 'N2O mixer'


A11 = units.ComponentSplitter('A11', ins=A7-0,
                              outs=(sol_N, sol_P, sol_K, 'sol_non_fertilizers'),
                              splits=(('NH3', 'NonNH3'), 'P', 'K'))

A12 = units.ComponentSplitter('A12', ins=A8-0,
                              outs=(liq_N, liq_P, liq_K, 'liq_non_fertilizers'),
                              splits=(('NH3', 'NonNH3'), 'P', 'K'))


def adjust_NH3_loss():
    A8._run()
    # Assume the slight higher loss of NH3 does not affect COD,
    # does not matter much since COD not considered in crop application
    A8.outs[0]._COD = A8.outs[1]._COD = A8.ins[0]._COD
A8.specification = adjust_NH3_loss

SceA = bst.System('SceA', path=(A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11, A12))
SceA.simulate()

# Emissions and product credits
#!! The linked_ws will be flushed out after simulation
CH4_item = StreamImpactItem(fugitive_CH4, GWP=28)
N2O_item = StreamImpactItem(fugitive_N2O, GWP=265)
liq_N_item = StreamImpactItem(liq_N, GWP=-5.4)
sol_N_item = StreamImpactItem(sol_N, GWP=-5.4)
liq_P_item = StreamImpactItem(liq_P, GWP=-4.9)
sol_P_item = StreamImpactItem(sol_P, GWP=-4.9)
liq_K_item = StreamImpactItem(liq_K, GWP=-1.5)
sol_K_item = StreamImpactItem(sol_K, GWP=-1.5)

SceA_tea = SimpleTEA(system=SceA, discount_rate=discount_rate, start_year=2018,
                     life_time=8, uptime_ratio=1, CAPEX=18606700, lang_factor=None,
                     annual_maintenance=0, annual_labor=12*3e6*12/exchange_rate,
                     system_add_OPEX=57120*e.price,
                     construction_schedule=None)
SceA_tea.show()
print('\n')

e_item = ImpactItem(ID='e_item', functional_unit='kWh', GWP=0.15)
get_e_price = lambda: e.price
# 57120 is the annual electricity usage for the whole treatment plant
get_annual_e = lambda: A2.add_OPEX/get_e_price()+57120

SceA_lca = LCA(system=SceA, life_time=8, life_time_unit='yr', uptime_ratio=1,
               # assuming all additional OPEX from electricity
               e_item=get_annual_e()*8)
SceA_lca.show()
print('\n')


# %%

# =============================================================================
# Summaries
# =============================================================================

get_AOC_cap = lambda: SceA_tea.AOC/ppl_existing
get_EAC_cap = lambda: SceA_tea.EAC/ppl_existing
print(f'Without CAPEX, the net cost is {get_AOC_cap():.1f} {currency}/cap/yr.')
print(f'With CAPEX, the net cost is {get_EAC_cap():.1f} {currency}/cap/yr.')

get_GWP = lambda: SceA_lca.total_impacts['GlobalWarming']/8/ppl_existing
print(f'Net emission is {get_GWP():.1f} {GWP.unit}/cap/yr.')

get_total_N = lambda: \
    (A1.outs[0].imass['NH3', 'NonNH3']+A1.outs[1].imass['NH3', 'NonNH3']).sum()
get_liq_N_recovery = lambda: liq_N.F_mass/ppl_existing/get_total_N()
get_sol_N_recovery = lambda: sol_N.F_mass/ppl_existing/get_total_N()
get_N_recovery = lambda: get_liq_N_recovery()+get_sol_N_recovery()
print(f'Total N recovery is {get_N_recovery():.1%}, '
      f'{get_liq_N_recovery():.1%} in liquid, '
      f'{get_sol_N_recovery():.1%} in solid.')

get_total_P = lambda: \
    (A1.outs[0].imass['P']+A1.outs[1].imass['P']).sum()
get_liq_P_recovery = lambda: liq_P.F_mass/ppl_existing/get_total_P()
get_sol_P_recovery = lambda: sol_P.F_mass/ppl_existing/get_total_P()
get_P_recovery = lambda: get_liq_P_recovery()+get_sol_P_recovery()
print(f'Total P recovery is {get_P_recovery():.1%}, '
      f'{get_liq_P_recovery():.1%} in liquid, '
      f'{get_sol_P_recovery():.1%} in solid.')

get_total_K = lambda: \
    (A1.outs[0].imass['K']+A1.outs[1].imass['K']).sum()
get_liq_K_recovery = lambda: liq_K.F_mass/ppl_existing/get_total_K()
get_sol_K_recovery = lambda: sol_K.F_mass/ppl_existing/get_total_K()
get_K_recovery = lambda: get_liq_K_recovery()+get_sol_K_recovery()
print(f'Total K recovery is {get_K_recovery():.1%}, '
      f'{get_liq_K_recovery():.1%} in liquid, '
      f'{get_sol_K_recovery():.1%} in solid.')

get_COD = lambda stream: stream.COD*stream.F_vol/1e3
get_total_COD = lambda: get_COD(A1.outs[0])+get_COD(A1.outs[1])
get_liq_COD_recovery = lambda: get_COD(A12.ins[0])/ppl_existing/get_total_COD()
get_sol_COD_recovery = lambda: get_COD(A11.ins[0])/ppl_existing/get_total_COD()
get_COD_recovery = lambda: get_liq_COD_recovery()+get_sol_COD_recovery()
print(f'Total COD recovery is {get_COD_recovery():.1%}, '
      f'{get_liq_COD_recovery():.1%} in liquid, '
      f'{get_sol_COD_recovery():.1%} in solid.')



# %%

# ws1 = A2.outs[0].copy('ws1')
# A4 = units.AnaerobicDigestion('A4', ins=ws1,
#                               outs=('treated', 'CH4', 'N2O'),
#                               # tau_previous=A2.emptying_period*365,
#                               decay_k_N=get_decay_k(tau_deg, log_deg),
#                                 max_CH4_emission=max_CH4_emission)
# A4.simulate()
# # A4.show()

# ws2 = A2.outs[0].copy('ws2')
# A5 = units.SludgeSeparator('A5', ins=ws2, outs=('liq', 'sol'))

# A5.simulate()
# # A5.show()




# ws5 = A2.outs[0].copy('ws5')
# A8 = units.DryingBed('A8', ins=ws5, outs=('solid', 'evaporated', 'CH4', 'N2O'),
#                      design_type='unplanted',
#                      decay_k_COD=get_decay_k(tau_deg, log_deg),
#                      decay_k_N=get_decay_k(tau_deg, log_deg),
#                      max_CH4_emission=max_CH4_emission)

# A8.simulate()
# # A8.show()

# ws6 = A2.outs[0].copy('ws6')
# A9 = units.AnaerobicBaffledReactor('A9', ins=ws6, outs=('treated', 'CH4', 'N2O'),
#                                    decay_k_COD=get_decay_k(tau_deg, log_deg),
#                                    max_CH4_emission=max_CH4_emission)

# A9.simulate()
# # A9.show()

# ws7 = A2.outs[0].copy('ws7')
# A10 = units.LiquidTreatmentBed('A10', ins=ws7, outs=('treated', 'CH4', 'N2O'),
#                                decay_k_COD=get_decay_k(tau_deg, log_deg),
#                                decay_k_N=get_decay_k(tau_deg, log_deg),
#                                max_CH4_emission=max_CH4_emission)

# A10.simulate()
# A10.show()








# biogas = WasteStream('biogas', CH4=1)
# AX2 = units.BiogasCombustion('AX2', ins=(biogas, 'air'),
#                               outs=('used', 'lost', 'wasted'),
#                               if_combustion=True,
#                               biogas_loss=0.1, biogas_eff=0.55)



# SceA = bst.System('SceA', path=(A1, A2, A3, AX, AX2))

# SceA.simulate()


# e_item = ImpactItem(ID='e_item', functional_unit='kWh', GWP=20)

# SceA_lca = LCA(SceA, life_time=8, life_time_unit='yr',
#                e_item=(5000, 'Wh'))



# # %%

# # =============================================================================
# # Scenario B: anaerobic treatment with existing latrines and conveyance
# # =============================================================================


# # %%

# # =============================================================================
# # Scenario C: containaer-based sanitation with existing treatment
# # =============================================================================

# print('\n----------Scenario C----------\n')
# C1 = units.Excretion('C1', outs=('urine', 'feces'), N_user=toilet_user)
# C2 = units.UDDT('C2', ins=(C1-0, C1-1,
#                             'toilet_paper', 'flushing_water',
#                             'cleaning_water', 'desiccant'),
#                 outs=('liquid_waste', 'solid_waste',
#                       'struvite', 'HAP', 'CH4', 'N2O'),
#                 N_user=toilet_user, OPEX_over_CAPEX=0.1,
#                 decay_k_COD=get_decay_k(tau_deg, log_deg),
#                 decay_k_N=get_decay_k(tau_deg, log_deg),
#                 max_CH4_emission=max_CH4_emission)
# C1.simulate()
# C2.simulate()

# truck = 'HandcartAndTruck'
# # Liquid waste
# interval = (C2.collection_period*truck_V[truck])/C2.tank_V
# C3 = units.Trucking('C3', ins=C2-0, outs=('transported_l', 'loss_l'),
#                     load_type='mass', load=truck_V[truck], load_unit='tonne',
#                     distance=5, distance_unit='km',
#                     interval=interval, interval_unit='day',
#                     fee=truck_cost[truck]/truck_V[truck]+0.01*ppl_alternative*interval,
#                     loss_ratio=0.02)

# # Solid waste
# interval = (C2.collection_period*truck_V[truck])/C2.tank_V
# C4 = units.Trucking('C4', ins=C2-1, outs=('transported_s', 'loss_s'),
#                     load_type='mass', load=truck_V[truck], load_unit='tonne',
#                     distance=5, distance_unit='km',
#                     interval=interval, interval_unit='day',
#                     fee=truck_cost[truck]/truck_V[truck]+0.01*ppl_alternative*interval,
#                     loss_ratio=0.02)


# CX = units.CropApplication('CX', ins=WasteStream(), loss_ratio=app_loss)
# def adjust_NH3_loss():
#     CX._run()
#     CX.outs[0]._COD = CX.outs[1]._COD = CX.ins[0]._COD
# CX.specification = adjust_NH3_loss


# SceC = bst.System('SceC', path=(C1, C2, C3, C4, CX))

# SceC.simulate()



# %%

# # =============================================================================
# # Try out impact-related classes
# # =============================================================================

# import biosteam as bst
# from sanitation import *

# H2SO4 = Component('H2SO4', search_ID='H2SO4', phase='l',
#                   particle_size='Soluble', degradability='Undegradable', organic=False)
# H2O = Component('H2O', search_ID='H2O', phase='l',
#                 particle_size='Soluble', degradability='Undegradable', organic=False)

# cmps2 = Components((H2SO4, H2O))
# cmps2.compile()
# cmps2.set_synonym('H2SO4', 'SulfuricAcid')
# cmps2.set_synonym('H2O', 'Water')

# bst.settings.set_thermo(cmps2)

# sulfuric_acid = WasteStream('sulfuric_acid', H2SO4=10, H2O=1000, units='kg/hr')
# sulfuric_acid2 = WasteStream('sulfuric_acid2', H2SO4=20, H2O=2000, units='kg/hr')

# CEDf = ImpactIndicator(ID='CEDf', unit='MJ', method='Cumulative energy demand',
#                         category='Fossil')





