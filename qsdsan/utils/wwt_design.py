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

import numpy as np

__all__ = ('get_SRT', 'get_oxygen_heterotrophs', 'get_oxygen_autotrophs', 'get_airflow')

def get_SRT(system, biomass_IDs, wastage=None, active_unit_IDs=None):
    """
    Estimate sludge residence time (SRT) of an activated sludge system.

    Parameters
    ----------
    system : obj
        The system whose SRT will be calculated for.
    biomass_IDs : tuple[str]
        Component IDs of active biomass.
    wastage : iterable[:class:`WasteStream`]
        Streams with wasted biomass.
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
    if wastage is None: wastage = [ws for ws in system.products if ws.phase in ('l','s')]
    waste = sum([ws.composite('solids', subgroup=biomass_IDs)*ws.F_vol*24 \
                 for ws in wastage])
    units = system.units if active_unit_IDs is None \
        else [u for u in system.units if u.ID in active_unit_IDs]
    retain = sum([u.get_retained_mass(biomass_IDs) for u in units if u.isdynamic])
    return retain/waste

def get_oxygen_heterotrophs(system, influent=None, eff_COD_soluble = None, f_d = 0.1, b_H = 0.4, SRT = 10, Y_H = 0.625):
    """
    Parameters
    ----------
    system : obj
        The system whose airflow will be calculated.
    influent : iterable[:class:`WasteStream`]
        Streams incoming to the process for which required oxygen needs to be calculated. The default is None.
    eff_COD_soluble : TYPE, optional
        DESCRIPTION. The default is None.
    f_d : float
        fraction of biomass that remains as cell debris. Default value is 0.1 gCOD/gCOD based on ASM2d 'f_XI_H'. 
    b_H : float
        Decay of heterotrophs [d^-1]. The default is 0.4 based on ASM2d.
    SRT : float
        Estimated sludge retention time of the system. Default is 10 days. 
    Y_H : float
        Yield of heterotrophs [gCOD/gCOD]. The default is 0.625.

    Returns
    -------
    float
        Oxygen requirement for heterotrophs in kg/day.

    """
    if influent is None:
        influent = [inf for inf in system.feeds if inf.phase in ('l')]
    
    influent_flow = np.array([inf.F_vol*24 for inf in influent]) # in m3/day
    influent_COD = np.array([inf.COD for inf in influent]) # in mg/L
    
    if eff_COD_soluble is None:
        eff_COD_soluble = np.array([eff.composite('COD', particle_size='s', unit = 'mg/L') for eff in system.products if eff.phase in 'l'])
        
    mass_influent_COD = np.sum(influent_flow*influent_COD/1000) # in kg/day
    mass_effluent_COD = np.sum(influent_flow*eff_COD_soluble/1000) # in kg/day
    
    mass_COD_treated =  mass_influent_COD - mass_effluent_COD # kg/day
    aeration_factor = 1 - (1 + f_d*b_H*SRT)*Y_H/(1 + b_H*SRT)
    
    return mass_COD_treated*aeration_factor


def get_oxygen_autotrophs(system, influent=None, eff_COD_soluble = None, f_d = 0.1, b_H = 0.4, b_AUT = 0.15, SRT = 10, 
                           Y_H = 0.625, Y_AUT = 0.24, K_NH = 1, U_AUT = 1, SF_DO = 1.375, ammonia_component_ID = 'S_NH4'):
    """
    
    Parameters
    ----------
    system : TYPE
        DESCRIPTION.
    influent : TYPE, optional
        Streams incoming to activated sludge process. Default is None. 
    eff_COD_soluble : TYPE, optional
        Maximum effluent soluble COD concentration [mg/L]. Default value is None.
    f_d : float
        fraction of biomass that remains as cell debris. Default value is 0.1 gCOD/gCOD based on ASM2d 'f_XI_A'. 
    b_H : float
        Decay of heterotrophs [d^-1]. The default is 0.4 based on ASM2d.
    b_AUT : float
        Decay of autotrophs [d^-1]. The default is 0.15 based on ASM2d.
    SRT : float
        Estimated sludge retention time of the system. Default is 10 days. 
    Y_H : float
        Yield of heterotrophs [gCOD/gCOD]. The default is 0.625.
    Y_AUT : float
        Yield of autotrophs [g COD/g N]. The default is 0.24.
    K_NH: float
        Substrate (Ammonium (nutrient)) half saturation coefficient for autotrophs, in [g N/m^3].
        The default is 1.0.
    U_AUT: float
        Autotrophic maximum specific growth rate, in [d^(-1)]. The default is 1.0.
    SF_DO: float
        Safety factor for dissolved oxygen. The default is 1.375. 
    ammonia_component_ID : tuple
        Component ID for ammonia.. The default is 'S_NH4'.

    Returns
    -------
    float
        Oxygen requirement for heterotrophs in kg/day.

    """

    if influent is None:
        influent = [inf for inf in system.feeds if inf.phase in ('l')]
        
    influent_flow = np.array([inf.F_vol*24 for inf in influent]) # in m3/day
    influent_COD = np.array([inf.COD for inf in influent]) # in mg/L
    
    if eff_COD_soluble is None:
        eff_COD_soluble = np.array([eff.composite('COD', particle_size='s', unit = 'mg/L') for eff in system.products if eff.phase in 'l'])
    
    NR = 0.087*(1 + f_d*b_H*SRT)*Y_H/(1 + b_H*SRT)
    TKN = np.array([inf.TKN for inf in influent])
    S_N_a = np.array([TKN]) - np.array([NR*(influent_COD - eff_COD_soluble)])
    S_NH = K_NH*(1/SRT  + b_AUT)/(U_AUT/SF_DO - (1 + b_AUT/SRT))
    aeration_factor = 4.57 - (1 + f_d*b_AUT*SRT)*Y_AUT/(1 + b_AUT*SRT)
    mass_N_removed = np.sum(influent_flow*(S_N_a - S_NH)/1000) # kg/day
    
    return mass_N_removed*aeration_factor

def get_airflow(oxygen_heterotrophs, oxygen_autotrophs, oxygen_transfer_efficiency = 15):
    """
    
    Parameters
    ----------
    oxygen_heterotrophs : float
        In kg/day.
    oxygen_autotrophs : float
        In kg/day.
    oxygen_transfer_efficiency : float
        Field oxygen transfer efficiency is percentage. The default is 15. 

    Returns
    -------
    Airflow in m3/min. 

    """
    
    required_oxygen = (oxygen_heterotrophs + oxygen_autotrophs)/24 # in kg/hr
    Q_air = 6*required_oxygen/oxygen_transfer_efficiency
    
    return Q_air
    