#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems
Copyright (C) 2020, Quantitative Sustainable Design Group

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the UIUC open-source license. Please refer to 
https://github.com/QSD-Group/QSDsan/blob/master/LICENSE.txt
for license details.
'''

# This tiered importing is because some modules in utils need to be imported
# before the the main modules (e.g., _component) since the main modules depend
# on them, while other modules in utils depend on the main modules

from . import (
    checkers,
    descriptors,
    )

__all__ = (
    *checkers.__all__,
    *descriptors.__all__,
    )

def secondary_importing():
    global __all__
    from . import (
        piping,
        loading,
        formatting,
        setters,
        )
    
    __all__ = (
        *__all__,
        *piping.__all__,
        *loading.__all__,
        *formatting.__all__,
        *setters.__all__,
                )
