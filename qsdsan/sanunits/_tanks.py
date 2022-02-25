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


from biosteam.units import Tank, MixTank, StorageTank
from .. import SanUnit


__all__ = ('Tank', 'MixTank', 'StorageTank', )


class Tank(SanUnit, Tank, isabstract=True):
    '''
    Similar to the :class:`biosteam.units.Tank`,
    but can be initialized with :class:`qsdsan.SanStream` and :class:`qsdsan.WasteStream`.

    See Also
    --------
    `biosteam.units.Tank <https://biosteam.readthedocs.io/en/latest/units/Tank.html>`_
    '''

    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                  vessel_type=None, tau=None, V_wf=None,
                  vessel_material=None, kW_per_m3=0.,
                  init_with='WasteStream', F_BM_default=None):

        SanUnit.__init__(self, ID, ins, outs, thermo,
                         init_with=init_with, F_BM_default=F_BM_default)

        self.vessel_type = vessel_type or self._default_vessel_type
        self.tau = tau or self._default_tau
        self.V_wf = V_wf or self._default_V_wf
        self.vessel_material = vessel_material or self._default_vessel_material
        self.kW_per_m3 = kW_per_m3 or self._default_kW_per_m3


class MixTank(Tank, MixTank):
    '''
    Similar to the :class:`biosteam.units.MixTank`,
    but can be initialized with :class:`qsdsan.SanStream` and :class:`qsdsan.WasteStream`.

    .. note::

        For dynamic simulation, CSTR should be used.

    See Also
    --------
    `biosteam.units.MixTank <https://biosteam.readthedocs.io/en/latest/units/Tank.html>`_
    '''


class StorageTank(Tank, StorageTank):
    '''
    Similar to the :class:`biosteam.units.StorageTank`,
    but can be initialized with :class:`qsdsan.SanStream` and :class:`qsdsan.WasteStream`.

    See Also
    --------
    `biosteam.units.StorageTank <https://biosteam.readthedocs.io/en/latest/units/Tank.html>`_
    '''