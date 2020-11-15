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

U1 = units.Excretion('U1', outs=('urine', 'feces'), N_user=N_user)

U2 = units.PitLatrine('U2', ins=(U1-0, U1-1, 'toilet_paper', 'flushing_water',
                                 'cleaning_water', 'desiccant'),
                      outs=('mixed_waste', 'leachate', 'CH4', 'N2O'),
                      N_user=N_user, OPEX_over_CAPEX=0.05,
                      decay_k=get_decay_k(tau_deg, log_deg),
                      max_CH4_emission=max_CH4_emission)



SceA = bst.System('SceA', path=(U1, U2))

SceA.simulate()

U2.show()










