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


# %%

#!!! Make sure all units are defined in a single registry (same as thermosteam)
#!!! Move Component/WasteStream ones here as well

# # Set up pint units
# from pint import UnitRegistry
# ureg = UnitRegistry()

from .utils.loading import data_path
data_path += 'units_of_measure.txt'

from thermosteam import units_of_measure as uom
ureg = uom.ureg
ureg.load_definitions(data_path)

AbsoluteUnitsOfMeasure = uom.AbsoluteUnitsOfMeasure
construction_units_of_measure = {
    'mass': AbsoluteUnitsOfMeasure('kg'),
    'volume': AbsoluteUnitsOfMeasure('m3'),
    }


definitions = uom.definitions
definitions['GWP'] = 'Global warming potential'