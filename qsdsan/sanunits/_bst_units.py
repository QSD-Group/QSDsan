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

import biosteam as bst
from .. import SanUnit

__all__ = (
    'Mixer',
    'Splitter', 'FakeSplitter', 'ReversedSplitter',
    'Pump',
    'Tank', 'StorageTank', 'MixTank',
    )

# =============================================================================
# Units subclassed from biosteam
# =============================================================================

class Mixer(SanUnit, bst.units.Mixer):
    '''
    Similar to the Mixer unit in biosteam, but takes WasteStreams.
    
    biosteam document
    -----------------
    `biosteam.units.Mixer <https://biosteam.readthedocs.io/en/latest/units/Mixer.html>`_

    '''
    __init__ = bst.units.Mixer.__init__


class Splitter(SanUnit, bst.units.Splitter):
    '''
    Similar to the Splitter unit in biosteam, but takes WasteStreams.

    biosteam document
    -----------------
    `biosteam.units.Splitter <https://biosteam.readthedocs.io/en/latest/units/Splitter.html>`_

    '''
    __init__ = bst.units.Splitter.__init__


class FakeSplitter(SanUnit, bst.units.FakeSplitter):
    '''
    Similar to the FakeSplitter unit in biosteam, but takes WasteStreams.
    
    biosteam document
    -----------------
    biosteam.units.FakeSplitter

    '''
    __init__ = bst.units.FakeSplitter.__init__


class ReversedSplitter(SanUnit, bst.units.ReversedSplitter):
    '''
    Similar to the ReversedSplitter unit in biosteam, but takes WasteStreams.
    
    biosteam document
    -----------------
    biosteam.units.FakeSplitter

    '''
    __init__ = bst.units.ReversedSplitter.__init__

    
class Pump(SanUnit, bst.units.Pump):
    '''
    Similar to the Pump unit in biosteam, but takes WasteStreams.
    
    biosteam document
    -----------------
    `biosteam.units.Pump <https://biosteam.readthedocs.io/en/latest/units/Pump.html>`_

    '''
    __init__ = bst.units.Pump.__init__

    
class Tank(SanUnit, bst.units.Tank, isabstract=True):
    '''
    Similar to the Tank unit in biosteam, but takes WasteStreams.
    
    biosteam document
    -----------------
    `biosteam.units.Tank <https://biosteam.readthedocs.io/en/latest/units/Tank.html>`_

    '''
    __init__ = bst.units.Tank.__init__

    
class StorageTank(SanUnit, bst.units.StorageTank):
    '''
    Similar to the StorageTank unit in biosteam, but takes WasteStreams.
    
    biosteam document
    -----------------
    `biosteam.units.StorageTank <https://biosteam.readthedocs.io/en/latest/units/StorageTank.html>`_

    '''
    __init__ = bst.units.StorageTank.__init__

    
class MixTank(SanUnit, bst.units.MixTank):
    '''
    Similar to the MixTank unit in biosteam, but takes WasteStreams.
    
    biosteam document
    -----------------
    `biosteam.units.MixTank <https://biosteam.readthedocs.io/en/latest/units/MixTank.html>`_

    '''
    __init__ = bst.units.MixTank.__init__
    _run = Mixer._run
    
    
    
    
    
    
    