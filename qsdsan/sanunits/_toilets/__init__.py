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

from ._toilet import *
from ._pit_latrine import *
from ._uddt import *


from . import (
    _toilet,
    _pit_latrine,
    _uddt,
    )


__all__ = (
    *_toilet.__all__,
    *_pit_latrine.__all__,
    *_uddt.__all__,
           )