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
# from sanitation import WasteStream
from sanitation.systems import bwaise
cmps = bwaise._cmps.cmps
units = bwaise._units

bst.settings.set_thermo(cmps)


U1 = units.Excretion('U1', outs=('urine', 'feces'))
U1.simulate()
















