#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:

    Yalin Li <mailto.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

import biosteam as bst, qsdsan as qs

__all__ = (
    'ProcessWaterCenter',
    )

_SANUNIT_ADDON_KEYS = (
    'include_construction', 'construction', 'transportation', 'equipment',
    'add_OPEX', 'lifetime', 'F_BM_default', 'isdynamic', 'exogenous_vars',
)


class ProcessWaterCenter(bst.facilities.ProcessWaterCenter, qs.SanUnit):
    '''
    biosteam.facilities.ProcessWaterCenter with QSDsan properties.

    See Also
    --------
    `biosteam.facilities.ProcessWaterCenter <https://biosteam.readthedocs.io/en/latest/API/facilities/ProcessWaterCenter.html>`_
    '''

    def __init__(self, *args, **kwargs):
        addons = {k: kwargs.pop(k) for k in tuple(_SANUNIT_ADDON_KEYS) if k in kwargs}
        super().__init__(*args, **kwargs)
        self._init_sanunit_addons(**addons)