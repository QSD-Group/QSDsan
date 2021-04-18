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

# This tiered importing is because some modules in utils need to be imported
# before the the main modules (e.g., _component) since the main modules depend
# on them, while other modules in utils depend on the main modules.

from . import (
    checkers,
    cod,
    descriptors,
    loading,
    parse,
    )

__all__ = (
    *checkers.__all__,
    *cod.__all__,
    *descriptors.__all__,    
    *loading.__all__,
    *parse.__all__,
    )

def _secondary_importing():
    global __all__
    from . import (
        decorators,
        formatting,
        getters,
        setters,
        )
    
    __all__ = (
        *__all__,
        *decorators.__all__,
        *formatting.__all__,
        *getters.__all__,
        *setters.__all__,
                )
