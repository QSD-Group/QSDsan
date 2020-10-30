#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Sanitation Explorer: Sustainable design of non-sewered sanitation technologies
Copyright (C) 2020, Sanitation Explorer Development Group

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>
    Joy Cheung

This module is under the UIUC open-source license. See 
https://github.com/QSD-for-WaSH/sanitation/blob/master/LICENSE.txt
for license details.
'''

import os
os.chdir("/Users/yalinli_cabbi/OneDrive/Coding/sanitation_root/")

import biosteam as bst
from sanitation import Components, WasteStream, units

components = Components.load_default()
components.compile()
bst.settings.set_thermo(components)

ins1 = WasteStream('ins1', SAc=5, H2O=1000, units='kg/hr')
ins2 = WasteStream('ins2', SF=10, H2O=1000, units='kg/hr')


ws1 = WasteStream.from_composite_measures('ws1', 1000)
ws2 = WasteStream.from_composite_measures('ws2', 1000, SNO3=3)

new_r = {'fSF_TotCOD': .3, 'fSNH4_STKN': .85}
ws2 = WasteStream.from_composite_measures('ws2', 10, ratios=new_r)

M1 = units.Mixer('M1', ins=(ins1, ins2, ''), outs='mixture')
M1.simulate()
M1.show()
M1.diagram()

S1 = units.Splitter('S1', ins=M1-0, outs=('', ''), split=0.2)

ins3 = WasteStream('ins3', SNO3=7, H2O=1000, units='kg/hr')
P1 = units.Pump('P1', ins=ins3)

M2 = units.MixTank('M2', ins=(S1-0, P1-0), tau=2)
M2-0-2-M1

System = bst.System('System', path=(M1, S1, P1, M2), recycle=M2-0)
System.show()
System.simulate()
System.show()
System.diagram()

M2.show()
M2.results()



