#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

__all__ = ('get_SRT',)

def get_SRT(system, biomass_IDs, active_unit_IDs=None):
    """
    Estimate sludge residence time (SRT) of an activated sludge system.

    Parameters
    ----------
    system : obj
        The system whose SRT will be calculated for.
    biomass_IDs : tuple[str]
        Component IDs of active biomass
    active_unit_IDs : tuple[str], optional
        IDs of activated sludge units. The default is None, meaning to include
        all units in the system.

    Returns
    -------
    [float] Estimated sludge residence time in days.

    .. note::

        [1] This function uses component flowrates of the system's product `WasteStream`
        for calculation, which does not carry time-dependent information. So, it
        should not be used for a system with dynamic influents.
        [2] The units included in calculation must all have :func:`get_retained_mass`
        to calculate the retained biomass.

    Examples
    --------
    `bsm1 system <https://github.com/QSD-Group/EXPOsan/blob/main/exposan/bsm1/system.py>`_
    """
    waste = sum([ws.composite('solids', subgroup=biomass_IDs)*ws.F_vol*24 for ws in system.products])
    units = system.units if active_unit_IDs is None \
        else [u for u in system.units if u.ID in active_unit_IDs]
    retain = sum([u.get_retained_mass(biomass_IDs) for u in units if u.isdynamic])
    return retain/waste