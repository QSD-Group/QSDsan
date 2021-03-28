#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>

Part of this module is based on the biosteam package:
https://github.com/BioSTEAMDevelopmentGroup/biosteam

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
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
    Similar to :class:`biosteam.units.Mixer`, but can be initilized with :class:`qsdsan.SanStream` and :class:`qsdsan.WasteStream`.
    
    See Also
    --------    
    `biosteam.units.Mixer <https://biosteam.readthedocs.io/en/latest/units/mixing.html>`_

    '''


class Splitter(SanUnit, bst.units.Splitter):
    '''
    Similar to :class:`biosteam.units.Splitter`, but can be initilized with :class:`qsdsan.SanStream` and :class:`qsdsan.WasteStream`.

    See Also
    --------
    `biosteam.units.Splitter <https://biosteam.readthedocs.io/en/latest/units/splitting.html>`_

    '''

    def __init__(self, ID='', ins=None, outs=(), thermo=None, *, split, order=None,
                 init_with='WasteStream'):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self._isplit = self.thermo.chemicals.isplit(split, order)


class FakeSplitter(SanUnit, bst.units.FakeSplitter):
    '''
    Similar to :class:`biosteam.units.FakeSplitter`, but can be initilized with :class:`qsdsan.SanStream` and :class:`qsdsan.WasteStream`.
    
    See Also
    --------
    `biosteam.units.FakeSplitter <https://biosteam.readthedocs.io/en/latest/units/splitting.html>`_

    '''


class ReversedSplitter(SanUnit, bst.units.ReversedSplitter):
    '''
    Similar to :class:`biosteam.units.ReversedSplitter`, but can be initilized with :class:`qsdsan.SanStream` and :class:`qsdsan.WasteStream`.
    
    See Also
    --------
    `biosteam.units.ReversedSplitter <https://biosteam.readthedocs.io/en/latest/units/splitting.html>`_

    '''

    
class Pump(SanUnit, bst.units.Pump):
    '''
    Similar to the :class:`biosteam.units.Pump`, but can be initilized with :class:`qsdsan.SanStream` and :class:`qsdsan.WasteStream`.
    
    See Also
    --------
    `biosteam.units.Pump <https://biosteam.readthedocs.io/en/latest/units/Pump.html>`_

    '''
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                 P=101325, pump_type='Default', material='Cast iron',
                 dP_design=405300, ignore_NPSH=True, init_with='WasteStream'):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.P = P
        self.pump_type = pump_type
        self.material = material
        self.dP_design = dP_design
        self.ignore_NPSH = ignore_NPSH

    
class Tank(SanUnit, bst.units.Tank, isabstract=True):
    '''
    Similar to the :class:`biosteam.units.Tank`, but can be initilized with :class:`qsdsan.SanStream` and :class:`qsdsan.WasteStream`.
    
    See Also
    --------
    `biosteam.units.Tank <https://biosteam.readthedocs.io/en/latest/units/Tank.html>`_

    '''
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                 vessel_type=None, tau=None, V_wf=None, 
                 vessel_material=None, kW_per_m3=0., init_with='WasteStream'):

        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)

        self.vessel_type = vessel_type or self._default_vessel_type
        self.tau = tau or self._default_tau
        self.V_wf = V_wf or self._default_V_wf
        self.vessel_material = vessel_material or self._default_vessel_material
        self.kW_per_m3 = kW_per_m3 or self._default_kW_per_m3

    
class StorageTank(Tank, bst.units.StorageTank):
    '''
    Similar to the :class:`biosteam.units.StorageTank`, but can be initilized with :class:`qsdsan.SanStream` and :class:`qsdsan.WasteStream`.
    
    See Also
    --------
    `biosteam.units.StorageTank <https://biosteam.readthedocs.io/en/latest/units/Tank.html>`_

    '''

    
class MixTank(Tank, bst.units.MixTank):
    '''
    Similar to the :class:`biosteam.units.MixTank`, but can be initilized with :class:`qsdsan.SanStream` and :class:`qsdsan.WasteStream`.
    
    See Also
    --------
    `biosteam.units.MixTank <https://biosteam.readthedocs.io/en/latest/units/Tank.html>`_

    '''
    
    
    
    
    
    
    