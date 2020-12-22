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

    '''
    __init__ = bst.units.Mixer.__init__
    __doc__ += bst.units.Mixer.__doc__
    __init__.__doc__ = __doc__

class Splitter(SanUnit, bst.units.Splitter):
    '''
    Similar to the Splitter unit in biosteam, but takes WasteStreams.

    biosteam document
    -----------------

    '''
    __init__ = bst.units.Splitter.__init__
    __doc__ += bst.units.Splitter.__doc__
    __init__.__doc__ = __doc__

class FakeSplitter(SanUnit, bst.units.FakeSplitter):
    '''
    Similar to the FakeSplitter unit in biosteam, but takes WasteStreams.
    
    biosteam document
    -----------------

    '''
    __init__ = bst.units.FakeSplitter.__init__
    __doc__ += bst.units.FakeSplitter.__doc__
    __init__.__doc__ = __doc__

class ReversedSplitter(SanUnit, bst.units.ReversedSplitter):
    '''
    Similar to the ReversedSplitter unit in biosteam, but takes WasteStreams.
    
    biosteam document
    -----------------

    '''
    __init__ = bst.units.ReversedSplitter.__init__
    __doc__ += bst.units.ReversedSplitter.__doc__
    __init__.__doc__ = __doc__
    
class Pump(SanUnit, bst.units.Pump):
    '''
    Similar to the Pump unit in biosteam, but takes WasteStreams.
    
    biosteam document
    -----------------

    '''
    __init__ = bst.units.Pump.__init__
    __doc__ += bst.units.Pump.__doc__
    __init__.__doc__ = __doc__
    
class Tank(SanUnit, bst.units.Tank, isabstract=True):
    '''
    Similar to the Tank unit in biosteam, but takes WasteStreams.
    
    biosteam document
    -----------------

    '''
    __init__ = bst.units.Tank.__init__
    __doc__ += bst.units.Tank.__doc__
    __init__.__doc__ = __doc__
    
class StorageTank(SanUnit, bst.units.StorageTank):
    '''
    Similar to the MixTank unit in biosteam, but takes WasteStreams.
    
    biosteam document
    -----------------

    '''
    __init__ = bst.units.StorageTank.__init__
    __doc__ += bst.units.StorageTank.__doc__
    __init__.__doc__ = __doc__
    
class MixTank(SanUnit, bst.units.MixTank):
    '''
    Similar to the MixTank unit in biosteam, but takes WasteStreams.
    
    biosteam document
    -----------------

    '''
    __init__ = bst.units.MixTank.__init__
    __doc__ += bst.units.MixTank.__doc__
    __init__.__doc__ = __doc__
    _run = Mixer._run
    
    
    
    
    
    
    