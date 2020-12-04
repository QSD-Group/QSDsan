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
    [1] Use consistent units for retention time, concentration, etc.
    [2] Lang factor for TEA
    [3] Add AOC-product for TEA
    [4] Streamline electricity

'''


# %%

import numpy as np
import biosteam as bst
# from sklearn.linear_model import LinearRegression
from sanitation import units, WasteStream, ImpactItem, StreamImpactItem, \
    SimpleTEA, LCA
# from sanitation import *
import bwaise
cmps = bwaise._cmps.cmps

bst.settings.set_thermo(cmps)
bst.PowerUtility.price = 0.17
items = ImpactItem._items


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
N_fertilizer = WasteStream('N_fertilizer', price=1.507)
P_fertilizer = WasteStream('P_fertilizer', price=3.983)
K_fertilizer = WasteStream('K_fertilizer', price=1.333)
N_item = StreamImpactItem(N_fertilizer, GWP=-5.4)
P_item = StreamImpactItem(P_fertilizer, GWP=-4.9)
K_item = StreamImpactItem(K_fertilizer, GWP=-1.5)

print('\n----------Scenario A----------\n')
A1 = units.Excretion('A1', outs=('urine', 'feces'))

A2 = units.PitLatrine('A2', ins=(A1-0, A1-1,
                                  'toilet_paper', 'flushing_water',
                                  'cleansing_water', 'desiccant'),
                      outs=('mixed_waste', 'leachate', 'CH4', 'N2O'),
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
                              outs=('liq', 'sol', 'CH4', 'N2O'),
                              decay_k_COD=get_decay_k(tau_deg, log_deg),
                              decay_k_N=get_decay_k(tau_deg, log_deg),
                              max_CH4_emission=max_CH4_emission)

A5 = units.Lagoon('A5', ins=A4-0, outs=('anaerobic_treated', 'CH4', 'N2O'),
                  design_type='anaerobic',
                  decay_k_N=get_decay_k(tau_deg, log_deg),
                  max_CH4_emission=max_CH4_emission)

A6 = units.Lagoon('A6', ins=A5-0, outs=('facultative_treated', 'CH4', 'N2O'),
                  design_type='facultative',
                  decay_k_N=get_decay_k(tau_deg, log_deg),
                  max_CH4_emission=max_CH4_emission)

A7 = units.DryingBed('A7', ins=A4-1, outs=('dried_sludge', 'evaporated', 'CH4', 'N2O'),
                     design_type='unplanted',
                     decay_k_COD=get_decay_k(tau_deg, log_deg),
                     decay_k_N=get_decay_k(tau_deg, log_deg),
                     max_CH4_emission=max_CH4_emission)

A8 = units.CropApplication('A8', ins=A6-0, outs=('liquid_fertilizer', 'loss'),
                           loss_ratio=app_loss)

A9 = units.Mixer('A9', ins=(A2-2, A4-2, A5-1, A6-1, A7-2),
                 outs=fugitive_CH4)

A10 = units.Mixer('A10', ins=(A2-3, A4-3, A5-2, A6-2, A7-3),
                 outs=fugitive_N2O)

A11 = units.ComponentSplitter('A11', ins=A8-0,
                              outs=(N_fertilizer, P_fertilizer, K_fertilizer,
                                    'non_fertilizers'),
                              splits=(('NH3', 'NonNH3'), 'P', 'K'))


def adjust_NH3_loss():
    A8._run()
    # Assume the slight higher loss of NH3 does not affect COD,
    # does not matter much since COD not considered in crop application
    A8.outs[0]._COD = A8.outs[1]._COD = A8.ins[0]._COD
A8.specification = adjust_NH3_loss

SceA = bst.System('SceA', path=(A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11))
SceA.simulate()


SceA_tea = SimpleTEA(system=SceA, discount_rate=discount_rate, start_year=2020,
                     life_time=8, uptime_ratio=1, CAPEX=18606700, lang_factor=None,
                     annual_maintenance=0, annual_labor=12*3e6*12/exchange_rate,
                     construction_schedule=None)


e_item = ImpactItem(ID='e_item', functional_unit='kWh', GWP=0.15)

SceA_lca = LCA(system=SceA, life_time=8, life_time_unit='yr',
               # assuming all additional OPEX from electricity
               e_item=(SceA_tea.add_OPEX/bst.PowerUtility.price+57120)*8)

# ADD ELECTRICITY COST TO A STANDALONE UNIT?













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





