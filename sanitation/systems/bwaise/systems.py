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


'''

import numpy as np
import biosteam as bst
from sanitation import units, WasteStream, LCA
import bwaise
cmps = bwaise._cmps.cmps

bst.settings.set_thermo(cmps)


# %%

# =============================================================================
# Assumptions
# =============================================================================

N_user = 16 # four people per household, four households per toilet
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

# Nutrient loss during applciation
app_loss = dict.fromkeys(('NH3', 'NonNH3', 'P', 'K', 'Mg', 'Ca'), 0.02)
app_loss['NH3'] = 0.05

#!!! Power law originally used in costing
truck_cost = {
    'TankerTruck1': 8e4/exchange_rate*1.15, # additional fee for tanker trucks
    'TankerTruck2': 12e4/exchange_rate*1.15,
    'TankerTruck3': 2e5/exchange_rate*1.15,
    'TankerTruck4': 25e4/exchange_rate*1.15,
    'HandcartAndTruck': 23e3/exchange_rate # per m3
    }
V = (3, 4.5, 8, 15, 1)
truck_V = dict.fromkeys(truck_cost.keys())
for i, j in zip(truck_V.keys(), V):
    truck_V[i] = j

bst.PowerUtility.price = 0.17


N_AD_rx = 3





# %%

# =============================================================================
# Scenario A: existing system
# =============================================================================
print('\n----------Scenario A----------\n')
A1 = units.Excretion('A1', outs=('urine', 'feces'), N_user=N_user)

A2 = units.PitLatrine('A2', ins=(A1-0, A1-1,
                                  'toilet_paper', 'flushing_water',
                                  'cleansing_water', 'desiccant'),
                      outs=('mixed_waste', 'leachate', 'CH4', 'N2O'),
                      N_user=N_user, OPEX_over_CAPEX=0.05,
                      decay_k_COD=get_decay_k(tau_deg, log_deg),
                      decay_k_N=get_decay_k(tau_deg, log_deg),
                      max_CH4_emission=max_CH4_emission)


#!!! RECHECK ALL TRANSPORTATION
truck = 'TankerTruck1' # assumed
interval = (A2.emptying_period*365*truck_V[truck])/A2.pit_V
A3 = units.Trucking('A3', ins=A2-0, outs=('transported', 'loss'),
                    load_type='mass', load=truck_V[truck], load_unit='tonne',
                    distance=5, distance_unit='km',
                    interval=interval, interval_unit='day',
                    fee=truck_cost[truck]/interval,
                    loss_ratio=0.02)


#!!! Note that before crop application, the COD is set, but after, is calculated
# based on the assumed ratio
AX = units.CropApplication('AX', ins=WasteStream(), loss_ratio=app_loss)
def adjust_NH3_loss():
    AX._run()
    # Assume the slight higher loss of NH3 does not affect COD,
    # does not matter much since COD not considered in crop application
    AX.outs[0]._COD = AX.outs[1]._COD = AX.ins[0]._COD
AX.specification = adjust_NH3_loss



biogas = WasteStream('biogas', CH4=1)
AX2 = units.BiogasCombustion('AX2', ins=(biogas, 'air'),
                              outs=('used', 'lost', 'wasted'),
                              if_combustion=True,
                              biogas_loss=0.1, biogas_eff=0.55)


SceA = bst.System('SceA', path=(A1, A2, A3, AX, AX2))

SceA.simulate()

SceA_lca = LCA(SceA, life_time=8, life_time_unit='yr')


# %%

ws1 = A2.outs[0].copy('ws1')
A4 = units.AnaerobicDigestion('A4', ins=ws1,
                              outs=('treated', 'CH4', 'N2O'),
                              # tau_previous=A2.emptying_period*365,
                              decay_k_N=get_decay_k(tau_deg, log_deg),
                                max_CH4_emission=max_CH4_emission)
A4.simulate()
# A4.show()

ws2 = A2.outs[0].copy('ws2')
A5 = units.SludgeSeparator('A5', ins=ws2, outs=('liq', 'sol'))

A5.simulate()
# A5.show()

ws3 = A2.outs[0].copy('ws3')
A6 = units.SedimentationTank('A6', ins=ws3,
                              outs=('liq', 'sol', 'CH4', 'N2O'),
                              # tau_previous=A2.emptying_period*365,
                              decay_k_COD=get_decay_k(tau_deg, log_deg),
                              decay_k_N=get_decay_k(tau_deg, log_deg),
                              max_CH4_emission=max_CH4_emission)

A6.simulate()
# A6.show()

ws4 = A2.outs[0].copy('ws4')
A7 = units.Lagoon('A7', ins=ws4, outs=('treated', 'CH4', 'N2O'),
                  design_type='anaerobic',
                  decay_k_N=get_decay_k(tau_deg, log_deg),
                  max_CH4_emission=max_CH4_emission)

A7.simulate()

ws5 = A2.outs[0].copy('ws5')
A8 = units.DryingBed('A8', ins=ws5, outs=('solid', 'evaporated', 'CH4', 'N2O'),
                     design_type='unplanted',
                     decay_k_COD=get_decay_k(tau_deg, log_deg),
                     decay_k_N=get_decay_k(tau_deg, log_deg),
                     max_CH4_emission=max_CH4_emission)

A8.simulate()
# A8.show()

ws6 = A2.outs[0].copy('ws6')
A9 = units.AnaerobicBaffledReactor('A9', ins=ws6, outs=('treated', 'CH4', 'N2O'),
                                   decay_k_COD=get_decay_k(tau_deg, log_deg),
                                   max_CH4_emission=max_CH4_emission)

A9.simulate()
# A9.show()

ws7 = A2.outs[0].copy('ws7')
A10 = units.LiquidTreatmentBed('A10', ins=ws7, outs=('treated', 'CH4', 'N2O'),
                               decay_k_COD=get_decay_k(tau_deg, log_deg),
                               decay_k_N=get_decay_k(tau_deg, log_deg),
                               max_CH4_emission=max_CH4_emission)

A10.simulate()
A10.show()






# %%

# =============================================================================
# Scenario B: anaerobic treatment with existing latrines and conveyance
# =============================================================================


# %%

# =============================================================================
# Scenario C: containaer-based sanitation with existing treatment
# =============================================================================

print('\n----------Scenario C----------\n')
C1 = units.Excretion('C1', outs=('urine', 'feces'), N_user=N_user)
C2 = units.UDDT('C2', ins=(C1-0, C1-1,
                            'toilet_paper', 'flushing_water',
                            'cleaning_water', 'desiccant'),
                outs=('liquid_waste', 'solid_waste',
                      'struvite', 'HAP', 'CH4', 'N2O'),
                N_user=N_user, OPEX_over_CAPEX=0.1,
                decay_k_COD=get_decay_k(tau_deg, log_deg),
                decay_k_N=get_decay_k(tau_deg, log_deg),
                max_CH4_emission=max_CH4_emission)
C1.simulate()
C2.simulate()

truck = 'HandcartAndTruck'
# Liquid waste
interval = (C2.collection_period*truck_V[truck])/C2.tank_V
C3 = units.Trucking('C3', ins=C2-0, outs=('transported_l', 'loss_l'),
                    load_type='mass', load=truck_V[truck], load_unit='tonne',
                    distance=5, distance_unit='km',
                    interval=interval, interval_unit='day',
                    fee=truck_cost[truck]/5/truck_V[truck]+0.01*N_user,
                    loss_ratio=0.02)

# Solid waste
interval = (C2.collection_period*truck_V[truck])/C2.tank_V
C4 = units.Trucking('C4', ins=C2-1, outs=('transported_s', 'loss_s'),
                    load_type='mass', load=truck_V[truck], load_unit='tonne',
                    distance=5, distance_unit='km',
                    interval=interval, interval_unit='day',
                    fee=truck_cost[truck]/5/truck_V[truck]+0.01*N_user,
                    loss_ratio=0.02)


CX = units.CropApplication('CX', ins=WasteStream(), loss_ratio=app_loss)
def adjust_NH3_loss():
    CX._run()
    CX.outs[0]._COD = CX.outs[1]._COD = CX.ins[0]._COD
CX.specification = adjust_NH3_loss


SceC = bst.System('SceC', path=(C1, C2, C3, C4, CX))

SceC.simulate()







