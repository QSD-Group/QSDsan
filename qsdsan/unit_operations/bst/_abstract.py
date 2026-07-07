#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:

    Yalin Li <mailto.yalin.li@gmail.com>

    Joy Zhang <joycheung1994@gmail.com>

    Jianan Feng <jiananf2@illinois.edu>

Part of this module is based on the biosteam package:
https://github.com/BioSTEAMDevelopmentGroup/biosteam

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.

Note: `unit_operations/static/_abstract.py` and `unit_operations/dynamic/_abstract.py`
are unrelated modules that happen to share this filename; each holds a different
set of unit classes.
'''

from biosteam.units import (
    FakeSplitter as BSTFakeSplitter,
    ReversedSplitter as BSTReversedSplitter,
    )
from ... import SanUnit

__all__ = (
    'FakeSplitter',
    'ReversedSplitter',
    )


# %%

class FakeSplitter(SanUnit, BSTFakeSplitter):
    '''
    Similar to :class:`biosteam.units.FakeSplitter`,
    but can be initialized with :class:`qsdsan.SanStream` and :class:`qsdsan.WasteStream`.

    See Also
    --------
    `biosteam.units.FakeSplitter <https://biosteam.readthedocs.io/en/latest/units/splitting.html>`_
    '''


class ReversedSplitter(SanUnit, BSTReversedSplitter):
    '''
    Similar to :class:`biosteam.units.ReversedSplitter`,
    but can be initialized with :class:`qsdsan.SanStream` and :class:`qsdsan.WasteStream`.

    See Also
    --------
    `biosteam.units.ReversedSplitter <https://biosteam.readthedocs.io/en/latest/units/splitting.html>`_
    '''
