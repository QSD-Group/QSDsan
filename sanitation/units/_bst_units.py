#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Sanitation Explorer: Sustainable design of non-sewered sanitation technologies
Copyright (C) 2020, Sanitation Explorer Development Group

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the UIUC open-source license. See 
https://github.com/QSD-for-WaSH/sanitation/blob/master/LICENSE.txt
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

class Mixer(SanUnit, bst.Mixer):
    '''Similar to the Mixer unit in biosteam, but takes WasteStreams'''
    __init__ = bst.Mixer.__init__
    
class Splitter(SanUnit, bst.Splitter):
    '''Similar to the Splitter unit in biosteam, but takes WasteStreams'''
    __init__ = bst.Splitter.__init__

class FakeSplitter(SanUnit, bst.FakeSplitter):
    '''Similar to the FakeSplitter unit in biosteam, but takes WasteStreams'''
    __init__ = bst.FakeSplitter.__init__

class ReversedSplitter(SanUnit, bst.ReversedSplitter):
    '''Similar to the ReversedSplitter unit in biosteam, but takes WasteStreams'''
    __init__ = bst.ReversedSplitter.__init__

class Pump(SanUnit, bst.Pump):
    '''Similar to the Pump unit in biosteam, but takes WasteStreams'''
    __init__ = bst.Pump.__init__
    
class Tank(SanUnit, bst.Tank, isabstract=True):
    '''Similar to the Tank unit in biosteam, but takes WasteStreams'''
    __init__ = bst.Tank.__init__
    
class StorageTank(SanUnit, bst.StorageTank):
    '''Similar to the MixTank unit in biosteam, but takes WasteStreams'''
    __init__ = bst.StorageTank.__init__
    
class MixTank(SanUnit, bst.MixTank):
    '''Similar to the MixTank unit in biosteam, but takes WasteStreams'''
    __init__ = bst.MixTank.__init__
    _run = Mixer._run
    
    
    
    
    
    
    