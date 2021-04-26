#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

'''
TODO:
    Move Component/WasteStream ones here as well
'''



# %%

__all__ = ('ureg', 'auom', 'ruom', 'parse_unit')

import os
from thermosteam.units_of_measure import (
    ureg,
    AbsoluteUnitsOfMeasure as auom,
    RelativeUnitsOfMeasure as ruom
    )

path = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                    'units_definition.txt')
ureg.load_definitions(path)

def parse_unit(value):
    str_list = value.split(' ') # for something like 'kg CO2-eq'
    if len(str_list) > 1:
        unit = str_list[0]
        others = ' '.join(str_list.pop(0))
        try: return auom(unit), others
        except: pass
    str_list = value.split('-') # for something like 'MJ-eq'
    if len(str_list) > 1:
        unit = str_list[0]
        others = '-'.join(str_list.pop(0))
        try: return auom(unit), others
        except: pass
    
    # For something like 'MJ' or 'tonne*km',
    # not doing this earlier as something like 'kg N' will be misinterpreted
    try: return auom(value), ''
    except: pass
    
    return None, value



