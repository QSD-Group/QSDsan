#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/master/LICENSE.txt
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
    Similar to :class:`biosteam.units.Mixer`, but takes :class:`qsdsan.WasteStream` objects.
    
    See Also
    --------
    `biosteam.units.Mixer <https://biosteam.readthedocs.io/en/latest/units/mixing.html>`_

    '''
    __init__ = bst.units.Mixer.__init__


class Splitter(SanUnit, bst.units.Splitter):
    '''
    Similar to :class:`biosteam.units.Splitter`, but takes :class:`qsdsan.WasteStream` objects.

    See Also
    --------
    `biosteam.units.Splitter <https://biosteam.readthedocs.io/en/latest/units/splitting.html>`_

    '''
    __init__ = bst.units.Splitter.__init__


class FakeSplitter(SanUnit, bst.units.FakeSplitter):
    '''
    Similar to :class:`biosteam.units.FakeSplitter`, but takes `qsdsan.WasteStream` objects.
    
    See Also
    --------
    `biosteam.units.FakeSplitter <https://biosteam.readthedocs.io/en/latest/units/splitting.html>`_

    '''
    __init__ = bst.units.FakeSplitter.__init__


class ReversedSplitter(SanUnit, bst.units.ReversedSplitter):
    '''
    Similar to :class:`biosteam.units.ReversedSplitter`, but takes `qsdsan.WasteStream` objects.
    
    See Also
    --------
    `biosteam.units.ReversedSplitter <https://biosteam.readthedocs.io/en/latest/units/splitting.html>`_

    '''
    __init__ = bst.units.ReversedSplitter.__init__

    
class Pump(SanUnit, bst.units.Pump):
    '''
    Similar to the :class:`biosteam.units.Pump`, but takes `qsdsan.WasteStream` objects.
    
    See Also
    --------
    `biosteam.units.Pump <https://biosteam.readthedocs.io/en/latest/units/Pump.html>`_

    '''
    __init__ = bst.units.Pump.__init__

    
class Tank(SanUnit, bst.units.Tank, isabstract=True):
    '''
    Similar to the :class:`biosteam.units.Tank`, but takes `qsdsan.WasteStream` objects.
    
    See Also
    --------
    `biosteam.units.Tank <https://biosteam.readthedocs.io/en/latest/units/Tank.html>`_

    '''
    __init__ = bst.units.Tank.__init__

    
class StorageTank(SanUnit, bst.units.StorageTank):
    '''
    Similar to the :class:`biosteam.units.StorageTank`, but takes `qsdsan.WasteStream` objects.
    
    See Also
    --------
    `biosteam.units.StorageTank <https://biosteam.readthedocs.io/en/latest/units/Tank.html>`_

    '''
    __init__ = bst.units.StorageTank.__init__

    
class MixTank(SanUnit, bst.units.MixTank):
    '''
    Similar to the :class:`biosteam.units.MixTank`, but takes `qsdsan.WasteStream` objects.
    
    See Also
    --------
    `biosteam.units.MixTank <https://biosteam.readthedocs.io/en/latest/units/Tank.html>`_

    '''
    __init__ = bst.units.MixTank.__init__
    _run = Mixer._run
    
    
    
    
    
    
    