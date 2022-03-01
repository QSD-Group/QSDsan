#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Joy Zhang <joycheung1994@gmail.com>
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

# Add a trailing "_" to differentiate the module from
# the functions within the module

from . import (
    # QSDsan modules, alphabetically
    test_bst_units_,
    test_component_,
    test_dyn_sys_,
    test_process_,
    test_sanunit_,
    test_waste_stream_,
    # EXPOsan systems
    test_exposan_,
    )


from .test_bst_units_ import *
from .test_component_ import *
from .test_dyn_sys_ import *
from .test_process_ import *
from .test_sanunit_ import *
from .test_waste_stream_ import *

from .test_exposan_ import *


__all__ = (
    *test_bst_units_.__all__,
    *test_component_.__all__,
    *test_dyn_sys_.__all__,
    *test_process_.__all__,
    *test_sanunit_.__all__,
    *test_waste_stream_.__all__,

    *test_exposan_.__all__,
    )