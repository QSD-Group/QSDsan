#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 18 16:52:05 2020

@author: yalinli_cabbi
"""



import biosteam as bst
import thermosteam as tmo
from sanitation import Component, WasteStream, units
from sanitation.utils import load_components_from_excel

components = load_components_from_excel(
    # '/Users/yalinli_cabbi/OneDrive/Coding/sanitation/sanitation/utils/default_components.xlsx'
    )

H2O = Component.from_chemical('H2O', tmo.Chemical('H2O'),
                              i_charge=0, particle_size='Soluble',
                              degradability='Undegradable', organic=False)
components.append(H2O)

for i in components:
    i.default()
    i.copy_models_from(components.H2O, ['sigma', 'epsilon', 'kappa', 'V', 'Cn', 'mu'])

components.compile()
bst.settings.set_thermo(components)


ins1 = WasteStream('ins1', SAc=5, H2O=1000, units='kg/hr')
ins2 = WasteStream('ins2', SHAc=10, H2O=1000, units='kg/hr')

M1 = units.Mixer('M1', ins=(ins1, ins2, ''), outs='mixture')
M1.simulate()
M1.show()
M1.diagram()

S1 = units.Splitter('S1', ins=M1-0, outs=('', ''), split=0.2)

ins3 = WasteStream('ins3', SHCO3=7, H2O=1000, units='kg/hr')
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



