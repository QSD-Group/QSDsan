#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 21 11:26:19 2020

@author: yalinli_cabbi
"""

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
    
    
    
    
    
    
    