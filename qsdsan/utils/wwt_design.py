#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Joy Zhang <joycheung1994@gmail.com>
    Saumitra Rai <raisaumitra9@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

import numpy as np

__all__ = ('get_SRT', 
           'get_oxygen_heterotrophs', 
           'get_oxygen_autotrophs', 
           'get_airflow', 
           'get_P_blower')

#%%

def get_SRT(system, biomass_IDs, wastage=None, active_unit_IDs=None):
    """
    Estimate sludge residence time (SRT) of an activated sludge system.

    Parameters
    ----------
    system : :class:`biosteam.System`
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

# def get_oxygen_heterotrophs(system, influent=None, eff_COD_soluble = None, f_d = 0.1, b_H = 0.4, SRT = 10, Y_H = 0.625):
def get_oxygen_heterotrophs(flow, influent_COD, eff_COD_soluble, 
                            f_d=0.1, b_H=0.4, SRT=10, Y_H=0.625):
    """
    Parameters
    ----------
    # system : :class:`biosteam.System`
    #     The system whose aeration airflow will be calculated.
    flow : float
        Volumetric flow rate through the system [m3/d].
    influent_COD : float
        Influent COD concentration [mg/L]. The default is None. 
    eff_COD_soluble : float
        Maximum effluent soluble COD concentration [mg/L]. 
    f_d : float, optional
        Fraction of biomass that remains as cell debris after decay. Default value is 
        0.1 gCOD/gCOD based on ASM2d 'f_XI_H'. 
    b_H : float, optional
        Decay rate constant of heterotrophs [d^-1]. The default is 0.4 based on ASM2d.
    SRT : float, optional
        Estimated sludge retention time of the system. Default is 10 days. 
    Y_H : float, optional
        Yield of heterotrophs [gCOD/gCOD]. The default is 0.625.

    Returns
    -------
    float
        Oxygen requirement for heterotrophs in kg/day.
        
    References
    ----------
    [1] Adapted from Equation 10.10 in GDLF [Grady Jr, C.P.L.; Daigger, G.T.; Love, N.G.; Filipe, C.D.M. Biological 
    Wastewater Treatment. 3rd ed. Boca Raton, FL: CRC Press 2011.]

    """
    # if influent is None:
    #     influent = [inf for inf in system.feeds if inf.phase == 'l']
    
    # influent_flow = np.array([inf.F_vol*24 for inf in influent]) # in m3/day
    # influent_COD = np.array([inf.COD for inf in influent]) # in mg/L
    
    # if eff_COD_soluble is None:
    #     eff_COD_soluble = np.array([eff.composite('COD', particle_size='s', unit = 'mg/L') \
    #                                 for eff in system.products if eff.phase == 'l'])
        
    # mass_influent_COD = np.sum(influent_flow*influent_COD/1000) # in kg/day
    # mass_effluent_COD = np.sum(influent_flow*eff_COD_soluble/1000) # in kg/day
    
    # mass_COD_treated =  mass_influent_COD - mass_effluent_COD # kg/day
    
    mass_COD_treated = flow * (influent_COD - eff_COD_soluble) * 1e-3 # kg/day
    aeration_factor = 1 - (1 + f_d*b_H*SRT)*Y_H/(1 + b_H*SRT)
    
    return mass_COD_treated*aeration_factor

# def get_oxygen_autotrophs(system, influent=None, eff_COD_soluble = None, f_d = 0.1, b_H = 0.4, b_AUT = 0.15, SRT = 10, 
#                            Y_H = 0.625, Y_AUT = 0.24, K_NH = 1, U_AUT = 1, SF_DO = 1.375, ammonia_component_ID = 'S_NH4'):
def get_oxygen_autotrophs(flow, influent_COD, eff_COD_soluble, influent_TKN,
                          f_d=0.1, b_H=0.4, b_AUT=0.15, SRT=10, 
                          Y_H=0.625, Y_AUT=0.24, K_NH=1, U_AUT=1, SF_DO=1.375):
    """
    
    Parameters
    ----------
    # system : :class:`biosteam.System`
    #     The system whose aeration airflow will be calculated.
    # influent : iterable[:class:`WasteStream`], optional
    #     Streams incoming to activated sludge process. Default is None. 
    flow : float
        Volumetric flow rate through the system [m3/d].
    influent_COD : float
        Influent COD concentration [mg/L].
    eff_COD_soluble : float
        Maximum effluent soluble COD concentration [mg/L].
    influent_TKN : float
        Influent TKN concentration [mg-N/L]. The default is None. 
    f_d : float
        fraction of biomass that remains as cell debris. Default value is 
        0.1 gCOD/gCOD based on ASM2d 'f_XI_A'. 
    b_H : float
        Decay rate constant of heterotrophs [d^-1]. The default is 0.4 based on ASM2d.
    b_AUT : float
        Decay rate constant of autotrophs [d^-1]. The default is 0.15 based on ASM2d.
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
    ammonia_component_ID : str
        Component ID for ammonia. The default is 'S_NH4'.

    Returns
    -------
    float
        Oxygen requirement for autotrophs in kg/day.
    
    References
    ----------
    [1] Adapted from Equation 11.16-11.19 in GDLF [Grady Jr, C.P.L.; Daigger, G.T.; Love, N.G.; Filipe, C.D.M. 
    Biological Wastewater Treatment. 3rd ed. Boca Raton, FL: CRC Press 2011.]

    """

    # if influent is None:
    #     influent = [inf for inf in system.feeds if inf.phase == 'l']
        
    # influent_flow = np.array([inf.F_vol*24 for inf in influent]) # in m3/day
    # influent_COD = np.array([inf.COD for inf in influent]) # in mg/L
    
    # if eff_COD_soluble is None:
    #     eff_COD_soluble = np.array([eff.composite('COD', particle_size='s', unit='mg/L') \
    #                                 for eff in system.products if eff.phase == 'l'])
    
    NR = 0.087*(1 + f_d*b_H*SRT)*Y_H/(1 + b_H*SRT)
    # TKN = np.array([inf.composite('N', organic=True) + 
    #                 inf.composite('N', subgroup=(ammonia_component_ID,)) \
    #                     for inf in influent])
    TKN = influent_TKN
    S_N_a = TKN - NR*(influent_COD - eff_COD_soluble)
    S_NH = K_NH*(1/SRT  + b_AUT)/(U_AUT/SF_DO - (1 + b_AUT/SRT))
    aeration_factor = 4.57 - (1 + f_d*b_AUT*SRT)*Y_AUT/(1 + b_AUT*SRT)
    # mass_N_removed = np.sum(influent_flow*(S_N_a - S_NH)/1000) # kg/day
    mass_N_removed = flow*(S_N_a - S_NH)/1000 # kg/day
    
    return mass_N_removed*aeration_factor

def get_airflow(oxygen_heterotrophs, oxygen_autotrophs, oxygen_transfer_efficiency=12):
    """
    
    Parameters
    ----------
    oxygen_heterotrophs : float
        In kg/day.
    oxygen_autotrophs : float
        In kg/day.
    oxygen_transfer_efficiency : float
        Field oxygen transfer efficiency in percentage. The default is 12. 

    Returns
    -------
    Airflow in m3/min. 
    
    References
    ----------
    [1] Adapted from Equation 11.2 in GDLF [Grady Jr, C.P.L.; Daigger, G.T.; Love, N.G.; Filipe, C.D.M. 
    Biological Wastewater Treatment. 3rd ed. Boca Raton, FL: CRC Press 2011.]

    """
    
    required_oxygen = (oxygen_heterotrophs + oxygen_autotrophs)/24 # in kg/hr
    Q_air = 6*required_oxygen/oxygen_transfer_efficiency
    
    return Q_air


def get_P_blower(q_air, T=20, p_atm=101.325, P_inlet_loss=1, 
                 P_diffuser_loss=7, h_submergance=5.18, efficiency=0.7,
                 K=0.283):
    """
    
    Parameters
    ----------
    q_air : float
        Air volumetric flow rate [m3/min].
    T : float
        Air temperature [degree Celsius].
    p_atm : float
        Atmostpheric pressure [kPa]
    P_inlet_loss : float
        Head loss at inlet [kPa].
    P_diffuser_loss : float
        Head loss due to piping and diffuser [kPa].
    h_submergance : float
        Diffuser submergance depth in m. The default is 17 feet (5.18 m)
    efficiency : float
        Blower efficiency. Default is 0.8. 
    K : float, optional
        Equal to (kappa - 1)/kappa, where kappa is the adiabatic exponent of air.
        Default is 0.283, equivalent to kappa = 1.3947.

    Returns
    -------
    Power of blower [kW].
    
    References
    ----------
    [1] Eq.(5-77) in Metcalf & Eddy, Wastewater Engineering: Treatment and 
        Resource Recovery. 5th ed. New York: McGraw-Hill Education. 2014.
    [2] Eq.(4.27) in Mueller, J.; William C.B.; and Popel, H.J. Aeration: 
        Principles and Practice, Volume 11. CRC press, 2002.
    """
    
    T += 273.15
    air_molar_vol = 22.4e-3 * (T/273.15)/(p_atm/101.325)   # m3/mol
    R = 8.31446261815324 # J/mol/K, universal gas constant
    
    p_in = p_atm - P_inlet_loss
    p_out = p_atm + 9.81*h_submergance + P_diffuser_loss
    
    # Q_air = q_air*(24*60) # m3/min to m3/day
    # P_blower = 1.4161e-5*(T + 273.15)*Q_air*((p_out/p_in)**0.283 - 1)/efficiency
    # return P_blower
    
    Q_air = q_air/60 # m3/min to m3/s
    WP = (Q_air / air_molar_vol) * R * T / K * ((p_out/p_in)**K - 1) / efficiency  # in W
    return WP/1000
    
