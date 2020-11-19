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
'''

import math
import biosteam as bst
from sanitation import WasteStream as WS

import bwaise
cmps = bwaise._cmps.cmps
units = bwaise._units

bst.settings.set_thermo(cmps)


# %%

# =============================================================================
# Assumptions
# =============================================================================

N_user = 16 # four people per household, four households per toilet

# Time take for full degradation, [yr]
tau_deg = 2
# Log reduction at full degradation
log_deg = 3
# Get reduction rate constant k for COD and N, use a function so that k can be
# changed during uncertainty analysis
def get_decay_k(tau_deg=2, log_deg=3):
    k = (-1/tau_deg)*math.log(10**-log_deg)
    return k

max_CH4_emission = 0.25

# =============================================================================
# Scenario A: existing system
# =============================================================================
print('\n----------Scenario A----------\n')
A1 = units.Excretion('A1', outs=('urine', 'feces'), N_user=N_user)

A2 = units.PitLatrine('A2', ins=(A1-0, A1-1, 'toilet_paper', 'flushing_water',
                                 'cleaning_water', 'desiccant'),
                      outs=('mixed_waste', 'leachate', 'CH4', 'N2O'),
                      N_user=N_user, OPEX_over_CAPEX=0.05,
                      decay_k=get_decay_k(tau_deg, log_deg),
                      max_CH4_emission=max_CH4_emission)



SceA = bst.System('SceA', path=(A1, A2))

SceA.simulate()

A2.show()


# =============================================================================
# Scenario B: anaerobic treatment with existing latrines and conveyance
# =============================================================================




# =============================================================================
# Scenario C: containaer-based sanitation with existing treatment
# =============================================================================

print('\n----------Scenario C----------\n')
C1 = units.Excretion('C1', outs=('urine', 'feces'), N_user=N_user)
C2 = units.UDDT('C2', ins=(C1-0, C1-1, 'toilet_paper', 'flushing_water',
                                 'cleaning_water', 'desiccant'),
                      outs=('liquid_waste', 'solid_waste',
                            'struvite', 'HAP', 'CH4', 'N2O'),
                      N_user=N_user, OPEX_over_CAPEX=0.1,
                      decay_k=get_decay_k(tau_deg, log_deg),
                      max_CH4_emission=max_CH4_emission)

SceC = bst.System('SceC', path=(C1, C2))

SceC.simulate()

C2.show()








