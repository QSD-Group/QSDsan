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

import biosteam as bst
from sanitation import WasteStream as WS

import bwaise
cmps = bwaise._cmps.cmps
units = bwaise._units

bst.settings.set_thermo(cmps)


# %%

# =============================================================================
# Human inputs, user-interface, and storage
# =============================================================================

N_user = 16 # four people per household, four households per toilet


U1 = units.Excretion('U1', outs=('urine', 'feces'), N_user=N_user)

U2 = units.PitLatrine('U2', ins=(U1-0, U1-1, 'toilet_paper', 'flushing_water',
                                 'cleaning_water'),
                      outs='mixed_waste',  N_user=N_user)



SceA = bst.System('SceA', path=(U1, U2))

SceA.simulate()












