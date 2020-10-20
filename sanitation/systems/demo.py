#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 18 16:52:05 2020

@author: yalinli_cabbi
"""



import biosteam as bst
import thermosteam as tmo
from sanitation import Component, WasteStream
from sanitation.utils import load_components_from_excel

components = load_components_from_excel(
    '/Users/yalinli_cabbi/OneDrive/Coding/sanitation/sanitation/utils/default_components.xlsx')

H2O = Component.from_chemical('H2O', tmo.Chemical('H2O'),
                              i_charge=0, particle_size='Soluble',
                              degradability='Undegradable', organic=False)
components.append(H2O)

for i in components:
    i.default()
    i.copy_models_from(components.H2O, ['sigma', 'epsilon', 'kappa', 'V', 'Cn', 'mu'])

components.compile()
bst.settings.set_thermo(components)


ws1 = WasteStream('ws1', SAc=5, H2O=1000, units='kg/hr')
ws1.show()
ws2 = WasteStream('ws2', SHAc=10, H2O=1000, units='kg/hr')
ws2.show(flow='kg/hr')

M1 = bst.units.Mixer('M1', ins=(ws1, ws2, ''), outs='mixture')
M1.simulate()
M1.show()
M1.diagram()

S1 = bst.units.Splitter('S1', ins=M1-0, outs=('ws3', 'ws4'), split=0.2)
S1.simulate()
S1.show()

ws5 = WasteStream('ws5', SHCO3=7, H2O=1000, units='kg/hr')
M2 = bst.units.MixTank('M2', ins=(S1-0, ws5), tau=2)
M2-0-2-M1

System = bst.System('System', path=(M1, S1, M2), recycle=M2-0)
System.show()
System.simulate()
System.show()
System.diagram()

M2.show()
M2.results()



