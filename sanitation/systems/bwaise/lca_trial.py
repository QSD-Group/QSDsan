#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 10:33:17 2020

@author: yalinli_cabbi
"""

# %%

# =============================================================================
# Try out impact-related classes
# =============================================================================

import biosteam as bst
from sanitation import Component, Components, WasteStream, ImpactIndicator, \
    ImpactItem, StreamImpactItem

H2SO4 = Component('H2SO4', search_ID='H2SO4', phase='l',
                  particle_size='Soluble', degradability='Undegradable', organic=False)
H2O = Component('H2O', search_ID='H2O', phase='l',
                particle_size='Soluble', degradability='Undegradable', organic=False)

cmps2 = Components((H2SO4, H2O))
cmps2.compile()
cmps2.set_synonym('H2SO4', 'SulfuricAcid')
cmps2.set_synonym('H2O', 'Water')

bst.settings.set_thermo(cmps2)

sulfuric_acid = WasteStream('sulfuric_acid', H2SO4=10, H2O=1000, units='kg/hr')

CEDf = ImpactIndicator(ID='CEDf', unit='MJ', method='Cumulative energy demand',
                       category='Fossil')

item_sulfuric_acid = StreamImpactItem(sulfuric_acid, CEDf=5, GWP=100, Eutrophication=1)



