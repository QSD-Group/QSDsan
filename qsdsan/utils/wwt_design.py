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
from warnings import warn

__all__ = ('get_SRT', 
           'get_oxygen_heterotrophs', 
           'get_oxygen_autotrophs', 
           'get_airflow', 
           'get_P_blower',
           'get_power_utility',
           'get_cost_sludge_disposal',
           'get_normalized_energy', 
           'get_daily_operational_cost',
           'get_total_operational_cost',
           'get_GHG_emissions_sec_treatment',
           'get_GHG_emissions_discharge',
           'get_GHG_emissions_electricity',
           'get_GHG_emissions_sludge_disposal',
           'get_CO2_eq_WRRF',
           'get_total_CO2_eq',
           
           'get_aeration_cost',
           'get_pumping_cost',
           'get_sludge_disposal_costs',
           'get_CH4_CO2_eq_treatment',
           'get_N2O_CO2_eq_treatment',
           'get_CH4_CO2_eq_discharge',
           'get_N2O_CO2_eq_discharge',
           'get_CH4_emitted_during_pl',
           'get_CH4_emitted_after_pl',
           'get_CO2_eq_electricity',
           'get_eq_natural_gas_price',
           
           # Function for gates work (to be removed at a later date)
           
           'estimate_ww_treatment_energy_demand',
           'estimate_N_removal_energy_demand',
           'estimate_P_removal_energy_demand')

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
    if waste == 0:
        warn('Wasted biomass calculated to be 0.')
        return
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
        Blower efficiency. Default is 0.7. 
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

def get_power_utility(system, active_unit_IDs=None):
    '''
    Parameters
    ----------
    system : :class:`biosteam.System`
        The system whose power will be calculated.
    active_unit_IDs : tuple[str], optional
        IDs of all units whose power needs to be accounted for. The default is None.

    Returns
    -------
    Cumulative power of sludge pumps [kW].
    '''
    
    power_consumption = 0
        
    for y in system.flowsheet.unit:
        if y.ID in active_unit_IDs:
            power_consumption += y.power_utility.power
        
    return power_consumption

def get_cost_sludge_disposal(sludge, unit_weight_disposal_cost = 375):
    '''
    Parameters
    ----------
    sludge : : iterable[:class:`WasteStream`], optional
        Effluent sludge from the system for which treatment and disposal costs are being calculated.
        The default is None.
    unit_weight_disposal_cost : float
        The sludge treatment and disposal cost per unit weight (USD/ton).
        Feasible range for this value lies between 100-800 USD/ton [1]. 
        
    Land application: 300 - 800 USD/US ton. [2]
    Landfill: 100 - 650 USD/US ton. [2]
    Incineration: 300 - 500 USD/US ton. [2]
    
    The default is 375 USD/US ton, which is the close to average of lanfill. 

    Returns
    -------
    Cost of sludge treatment and disposal (USD/day).
    
    References
    -------
    [1] Feng, J., Li, Y., Strathmann, T. J., & Guest, J. S. (2024b). 
    Characterizing the opportunity space for sustainable hydrothermal valorization 
    of Wet Organic Wastes. Environmental Science &; Technology. 
    https://doi.org/10.1021/acs.est.3c07394 
    
    [2] Peccia, J., & Westerhoff, P. (2015). We should expect more out of our sewage 
    sludge. Environmental Science & Technology, 49(14), 8271–8276. 
    https://doi.org/10.1021/acs.est.5b01931 

    '''
    
    sludge_prod = np.array([sludge.composite('solids', True, particle_size='x', unit='ton/d') \
                            for sludge in sludge]) # in ton/day
        
    cost_sludge_disposal = np.sum(sludge_prod)*unit_weight_disposal_cost   #in USD/day
    
    return cost_sludge_disposal

def get_normalized_energy(system, aeration_power, pumping_power, miscellaneous_power):
    '''
    Parameters
    ----------
    system : :class:`biosteam.System`
        The system for which normalized energy consumption is being determined.
    aeration_power : float, optional
        Power of blower [kW].
    pumping_power : float, optional
        Power rquired for sludge pumping and other equipments [kW].
    miscellaneous_power : float, optional
        Any other power requirement in the system [kW].

    Returns
    -------
    Normalized energy consumption associated with WRRF (kWh/m3). [numpy array] 

    '''
    
    normalized_aeration_energy = aeration_power/sum([s.F_vol for s in system.feeds])
    normalized_pumping_energy = pumping_power/sum([s.F_vol for s in system.feeds])
    normalized_miscellaneous_energy = miscellaneous_power/sum([s.F_vol for s in system.feeds])

    normalized_energy_WRRF = np.array([normalized_aeration_energy, normalized_pumping_energy, \
                                       normalized_miscellaneous_energy])

    return normalized_energy_WRRF

def get_daily_operational_cost(system, aeration_power, pumping_power, miscellaneous_power, \
                                    sludge_disposal_cost, unit_electricity_cost = 0.161):
    '''
    Parameters
    ----------
    system : :class:`biosteam.System`
        The system for which normalized energy consumption is being determined.
    aeration_power : float, optional
        Power of blower [kW].
    pumping_power : float, optional
        Power rquired for sludge pumping and other equipments [kW].
    miscellaneous_power : float, optional
        Any other power requirement in the system [kW].
    sludge_disposal_cost : float
        Cost of sludge treatment and disposal (USD/day).
    unit_electricity_cost : float
        Unit cost of electricity. Default value is taken as $0.161/kWh [1]. 

    Returns
    -------
    Normalized operational cost associated with WRRF (USD/day). [numpy array]
    
    [1] https://www.bls.gov/regions/midwest/data/averageenergyprices_selectedareas_table.htm
    
    '''
    aeration_cost = aeration_power*24*unit_electricity_cost # in (kWh/day)*(USD/kWh) = USD/day
    pumping_cost = pumping_power*24*unit_electricity_cost # in (kWh/day)*(USD/kWh) = USD/day
    sludge_disposal_costs = sludge_disposal_cost
    miscellaneous_cost = miscellaneous_power*24*unit_electricity_cost # in (kWh/day)*(USD/kWh) = USD/day
    
    operational_costs_WRRF = np.array([aeration_cost, pumping_cost, sludge_disposal_costs, miscellaneous_cost])
    
    operational_costs_WRRF = operational_costs_WRRF/sum([s.F_vol*24 for s in system.feeds])
    
    return operational_costs_WRRF

def get_aeration_cost(q_air, # aeration (blower) power 
                      system, # sludge pumping power 
                      T=20, p_atm=101.325, P_inlet_loss=1, P_diffuser_loss=7, 
                      h_submergance=5.18, efficiency=0.7, K=0.283, # aeration (blower) power
                      unit_electricity_cost = 0.161): 
    '''
    Parameters
    ----------
    
    q_air : float
        Air volumetric flow rate [m3/min].
    T : float
        Air temperature [degree Celsius].
    p_atm : float
        Atmostpheric pressure [kPa]
    P_inlet_loss : float
        Head loss at inlet [kPa]. The default is 1 kPa. 
    P_diffuser_loss : float
        Head loss due to piping and diffuser [kPa]. The default is 7 kPa. 
    h_submergance : float
        Diffuser submergance depth in m. The default is 17 feet (5.18 m)
    efficiency : float
        Blower efficiency. Default is 0.7. 
    K : float, optional
        Equal to (kappa - 1)/kappa, where kappa is the adiabatic exponent of air.
        Default is 0.283, equivalent to kappa = 1.3947.  
        
    ------------------------------------------------------------------------------ 
    
    system : :class:`biosteam.System`
        The system whose power will be calculated.
        
    unit_electricity_cost : float
        Unit cost of electricity. Default value is taken as $0.161/kWh [1]. 

    Returns
    -------
    Normalized cost associated with aeration (USD/m3). [int]
    
    '''
    
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
    
    aeration_power = WP/1000 # in kW
    
    aeration_cost = aeration_power*24*unit_electricity_cost # in (kWh/day)*(USD/kWh) = USD/day
    
    operational_costs_WRRF = aeration_cost/sum([s.F_vol*24 for s in system.feeds])

    return operational_costs_WRRF

def get_pumping_cost(system, active_unit_IDs=None,  # sludge pumping power 
                     unit_electricity_cost = 0.161): 
    '''
    Parameters
    ----------
    system : :class:`biosteam.System`
        The system whose power will be calculated.
    active_unit_IDs : tuple[str], optional
        IDs of all units whose power needs to be accounted for. The default is None.
        
    ------------------------------------------------------------------------------ 
        
    unit_electricity_cost : float
        Unit cost of electricity. Default value is taken as $0.161/kWh [1]. 

    Returns
    -------
    Normalized operational cost associated with sludge pumping (USD/m3). [int]
    
    '''
    power_consumption = 0
        
    for y in system.flowsheet.unit:
        if y.ID in active_unit_IDs:
            power_consumption += y.power_utility.power
            
    pumping_power = power_consumption
    
    pumping_cost = pumping_power*24*unit_electricity_cost # in (kWh/day)*(USD/kWh) = USD/day

    operational_costs_WRRF = pumping_cost/sum([s.F_vol*24 for s in system.feeds])
    
    return operational_costs_WRRF

def get_sludge_disposal_costs(sludge,  # sludge disposal costs 
                              system, unit_weight_disposal_cost = 350, # sludge disposal costs 
                              ): 
    '''
    Parameters
    ----------
    
    system : :class:`biosteam.System`
        The system whose power will be calculated.
        
    ------------------------------------------------------------------------------ 
    
    sludge : : iterable[:class:`WasteStream`], optional
        Effluent sludge from the system for which treatment and disposal costs are being calculated.
        The default is None.
    unit_weight_disposal_cost : float
        The sludge treatment and disposal cost per unit weight (USD/ton).
        Feasible range for this value lies between 100-800 USD/ton [1]. 
        
    Land application: 300 - 800 USD/US ton. [2]
    Landfill: 100 - 650 USD/US ton. [2]
    Incineration: 300 - 500 USD/US ton. [2]
    
    The default is 375 USD/US ton, which is the close to average of lanfill. 

    Returns
    -------
    Normalized operational cost associated with WRRF (USD/m3). [int]
    
    '''
    
    sludge_prod = np.array([sludge.composite('solids', True, particle_size='x', unit='ton/d') \
                            for sludge in sludge]) # in ton/day
    sludge_disposal_costs = np.sum(sludge_prod)*unit_weight_disposal_cost   #in USD/day
    operational_costs_WRRF = sludge_disposal_costs/sum([s.F_vol*24 for s in system.feeds])
    
    return operational_costs_WRRF

def get_eq_natural_gas_price(system, gas, natural_gas_price = 0.0041, CH4_nat_gas = 0.95): 
    '''
    Parameters
    ----------
    
    system : :class:`biosteam.System`
        The system whose power will be calculated.
        
    ------------------------------------------------------------------------------ 
    
    gas : : iterable[:class:`WasteStream`], optional
        Effluent gas from anaerobic digestor. 
        The default is None.
    natural_gas_price : float
        Price of Natural gas in USD/m3 (eia.gov). 
        The default is 0.0041 USD/m3. 
    CH4_nat_gas: Percentage of natural gas that is methane. 
        The default is 0.95.

    Returns
    -------
    Normalized natural gas equivalent cost (USD/m3). [int]
    
    '''
    
    CH4_produced = gas.F_vol*24 # m3/day
    
    eq_nat_gas_produced = CH4_produced/CH4_nat_gas # m3/day
    
    price_nat_gas = eq_nat_gas_produced*natural_gas_price # USD/day
        
    operational_costs_WRRF = price_nat_gas/sum([s.F_vol*24 for s in system.feeds])
    
    return operational_costs_WRRF

def get_total_operational_cost(q_air, # aeration (blower) power 
                                     sludge,  # sludge disposal costs 
                                     system, active_unit_IDs=None,  # sludge pumping power 
                                     T=20, p_atm=101.325, P_inlet_loss=1, P_diffuser_loss=7, 
                                     h_submergance=5.18, efficiency=0.7, K=0.283, # aeration (blower) power 
                                     miscellaneous_power = 0, 
                                     unit_weight_disposal_cost = 350, # sludge disposal costs 
                                     unit_electricity_cost = 0.161): 
    '''
    Parameters
    ----------
    
    q_air : float
        Air volumetric flow rate [m3/min].
    T : float
        Air temperature [degree Celsius].
    p_atm : float
        Atmostpheric pressure [kPa]
    P_inlet_loss : float
        Head loss at inlet [kPa]. The default is 1 kPa. 
    P_diffuser_loss : float
        Head loss due to piping and diffuser [kPa]. The default is 7 kPa. 
    h_submergance : float
        Diffuser submergance depth in m. The default is 17 feet (5.18 m)
    efficiency : float
        Blower efficiency. Default is 0.7. 
    K : float, optional
        Equal to (kappa - 1)/kappa, where kappa is the adiabatic exponent of air.
        Default is 0.283, equivalent to kappa = 1.3947.  
        
    ------------------------------------------------------------------------------ 
    
    system : :class:`biosteam.System`
        The system whose power will be calculated.
    active_unit_IDs : tuple[str], optional
        IDs of all units whose power needs to be accounted for. The default is None.
        
    ------------------------------------------------------------------------------ 
    
    sludge : : iterable[:class:`WasteStream`], optional
        Effluent sludge from the system for which treatment and disposal costs are being calculated.
        The default is None.
    unit_weight_disposal_cost : float
        The sludge treatment and disposal cost per unit weight (USD/ton).
        Feasible range for this value lies between 100-800 USD/ton [1]. 
        
    Land application: 300 - 800 USD/US ton. [2]
    Landfill: 100 - 650 USD/US ton. [2]
    Incineration: 300 - 500 USD/US ton. [2]
    
    The default is 375 USD/US ton, which is the close to average of lanfill. 
    
    ------------------------------------------------------------------------------ 
        
    miscellaneous_power : float, optional
        Any other power requirement in the system [kW].
        
    unit_electricity_cost : float
        Unit cost of electricity. Default value is taken as $0.161/kWh [1]. 

    Returns
    -------
    Normalized operational cost associated with WRRF (USD/m3). [int]
    
    '''
    
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
    
    aeration_power = WP/1000 # in kW
    
    power_consumption = 0
        
    for y in system.flowsheet.unit:
        if y.ID in active_unit_IDs:
            power_consumption += y.power_utility.power
            
    pumping_power = power_consumption
    
    aeration_cost = aeration_power*24*unit_electricity_cost # in (kWh/day)*(USD/kWh) = USD/day
    pumping_cost = pumping_power*24*unit_electricity_cost # in (kWh/day)*(USD/kWh) = USD/day
    
    sludge_prod = np.array([sludge.composite('solids', True, particle_size='x', unit='ton/d') \
                            for sludge in sludge]) # in ton/day
    sludge_disposal_costs = np.sum(sludge_prod)*unit_weight_disposal_cost   #in USD/day
    
    miscellaneous_cost = miscellaneous_power*24*unit_electricity_cost # in (kWh/day)*(USD/kWh) = USD/day
    
    operational_costs_WRRF = np.array([aeration_cost, pumping_cost, sludge_disposal_costs, miscellaneous_cost]) 

    operational_costs_WRRF = operational_costs_WRRF/sum([s.F_vol*24 for s in system.feeds])

    total_operational_cost = np.sum(operational_costs_WRRF)             #5
    
    return total_operational_cost 

def get_GHG_emissions_sec_treatment(system = None, influent=None, effluent = None,
                                    CH4_EF=0.0075, N2O_EF=0.016):
    '''    
    Parameters
    ----------
    system : :class:`biosteam.System`
        The system for which emissions during secondary treatment are being calculated. The default is None.
    influent : : iterable[:class:`WasteStream`], optional
        Influent wastestreams to the system whose wastewater composition determine the potential for GHG emissions. The default is None.
    effluent  : : iterable[:class:`WasteStream`], optional
        The wastestreams which represents the effluent from the treatment process. The default is None.
    CH4_EF :  float, optional.
        The emission factor used to calculate methane emissions in secondary treatment. The default is 0.0075 kg CH4/ kg rCOD. [1]
    N2O_EF : float, optional
        The emission factor used to calculate nitrous oxide emissions in secondary treatment. The default is 0.016 kg N2O-N/ kg N. [1]

    Returns
    -------
    CH4_emitted : float
        The amount of methane emitted during secondary treatment (kg/day).
    N2O_emitted : float
        The amount of nitrous oxide emitted during secondary treatment (kg/day).
        
    References
    ----------
    [1] Chapter - 6, IPCC. (2019). In 2019 Refinement to the 2006 IPCC Guidelines for National Greenhouse Gas Inventories.

    '''
    
    if influent is None:
        influent = [inf for inf in system.feeds if inf.phase == 'l']

    influent_flow = np.array([inf.F_vol*24 for inf in influent]) # in m3/day
    influent_COD = np.array([inf.COD for inf in influent]) # in mg/L
    mass_influent_COD = np.sum(influent_flow*influent_COD/1000) # in kg/day
    
    effluent_flow = np.array([eff.F_vol*24 for eff in effluent])
    effluent_COD =  np.array([eff.COD for eff in effluent]) # in mg/L
    mass_effluent_COD = np.sum(effluent_flow*effluent_COD/1000) # in kg/day
    
    mass_removed_COD = mass_influent_COD - mass_effluent_COD
    CH4_emitted = CH4_EF*mass_removed_COD
    
    influent_N = np.array([inf.TN for inf in influent]) # in mg/L
    mass_influent_N = np.sum(influent_flow*influent_N/1000) # in kg/day
    
    N2O_emitted = N2O_EF*mass_influent_N
    
    return CH4_emitted, N2O_emitted

def get_GHG_emissions_discharge(effluent=None, CH4_EF=0.009, N2O_EF=0.005):
    '''    
    Parameters
    ----------
    effluent : : iterable[:class:`WasteStream`], optional
        Effluent wastestreams from the system whose wastewater composition determine the potential for GHG emissions at discharge. The default is None.
    CH4_EF_discharge :  float, optional.
        The emission factor used to calculate methane emissions in discharge. The default is 0.009 kg CH4/ kg effluent COD. [1]
    N2O_EF_discharge : float, optional
        The emission factor used to calculate nitrous oxide emissions in discharge. The default is 0.005 kg N2O-N/ kg effluent N. [1]

    Returns
    -------
    CH4_emitted : float
        The amount of methane emitted at discharge (kg/day).
    N2O_emitted : float
        The amount of nitrous oxide emitted at discharge (kg/day).
        
    References
    ----------
    [1] Chapter - 6, IPCC. (2019). In 2019 Refinement to the 2006 IPCC Guidelines 
    for National Greenhouse Gas Inventories.

    '''

    effluent_flow = np.array([eff.F_vol*24 for eff in effluent]) # in m3/day
    effluent_COD = np.array([eff.COD for eff in effluent]) # in mg/L
    mass_effluent_COD = np.sum(effluent_flow*effluent_COD/1000) # in kg/day
    
    CH4_emitted = CH4_EF*mass_effluent_COD
    
    effluent_N = np.array([eff.TN for eff in effluent]) # in mg/L
    mass_effluent_N = np.sum(effluent_flow*effluent_N/1000) # in kg/day
    
    N2O_emitted = N2O_EF*mass_effluent_N
    
    return CH4_emitted, N2O_emitted
    
def get_GHG_emissions_electricity(system, power_blower, power_pump, CO2_EF=0.675):
    '''
    Parameters
    ----------
    system : :class:`biosteam.System`
        The system for which tier-2 GHG emissions due to electricity consumption are being calculated. 
    power_blower : float
        Power of blower [kW].
    power_pump : float
        Power required for pumping and other utilities at treatment facility [kW].
    CO2_EF : float
        The emission factor used to calculate tier-2 CO2 emissions due to electricity consumption. 
        The default is 0.675 kg-CO2-Eq/kWh. [1]
        The emission factor is dependent on the region, and is as follows for USA:
            
            {SERC Reliability Corporation (SERC): 0.618 kg-CO2-Eq/kWh
            ReliabilityFirst (RFC): 0.619 kg-CO2-Eq/kWh
            Western Electricity Coordinating Council (WECC): 0.436 kg-CO2-Eq/kWh
            Texas Reliability Entity (TRE): 0.574 kg-CO2-Eq/kWh
            Southwest Power Pool (SPP): 0.733 kg-CO2-Eq/kWh
            Midwest Reliability Organization (MRO): 0.675 kg-CO2-Eq/kWh
            Florida Reliability Coordinating Council (FRCC): 0.531 kg-CO2-Eq/kWh
            Northeast Power Coordinating Council (NPCC): 0.244 kg-CO2-Eq/kWh}
            (HICC): 0.814 kg-CO2-Eq/kWh}
            (ASCC): 0.599 kg-CO2-Eq/kWh}
            
            
    Returns
    -------
    CO2_emitted : float
        The amount of eq. CO2 emitted due to electrity consumption (kg-CO2-Eq/day).

    '''
    
    total_energy_consumed = (power_blower + power_pump)*24 # in kWh/day
    CO2_emissions = total_energy_consumed*CO2_EF # in kg-CO2-Eq/day
    
    return CO2_emissions

def get_GHG_emissions_sludge_disposal(sludge=None, DOC_f = 0.45, MCF = 0.8, k = 0.06, F=0.5, pl=30):
    '''
    Parameters
    ----------
    sludge : : iterable[:class:`WasteStream`], optional
        Effluent sludge from the system for which GHG emissions are being calculated. The default is None.
    DOC_f : float, optional
        fraction of DOC that can decompose (fraction). The default value is 0.5.
    MCF : float, optional
        CH4 correction factor for aerobic decomposition in the year of deposition (fraction). The default is 0.8.
    k : TYPE, optional
        Methane generation rate (k). The default is 0.185. (1/year)
        
        The decomposition of carbon is assumed to follow first-order kinetics 
        (with rate constant k), and methane generation is dependent on the amount of remaining decomposable carbon 
        in the waste. For North America (boreal and temperate climate) the default values for wet and dry climate 
        are:
            k (dry climate) = 0.06
            k (wet climate) = 0.185
    F : float, optional
        Fraction of methane in generated landfill gas (volume fraction).  The default is 0.5.
    pl : float, optional
        The project lifetime over which methane emissions would be calculated. (years)
        The default is 30 years.

    Returns
    -------
    CH4_emitted : float
        The average amount of methane emitted during sludge disposal (kg/day).
        
     References
     ----------
     [1] Chapter - 3: Solid Waste Disposal, IPCC. (2019). In 2019 Refinement to
     the 2006 IPCC Guidelines for National Greenhouse Gas Inventories.

    '''
    # sludge_flow = np.array([slg.F_vol*24 for slg in sludge]) # in m3/day
    # sludge_COD = np.array([slg.COD for slg in sludge]) # in mg/L
    DOC_mass_flow = 0
    
    for slg in sludge:
        DOC_mass_flow += slg.composite("C", flow=True, exclude_gas=True, 
                      subgroup=None, particle_size=None,
                      degradability="b", organic=True, volatile=None,
                      specification=None, unit="kg/day")

    annual_DOC_mass = 365*DOC_mass_flow # in kg/year

    annual_DDOC = annual_DOC_mass*DOC_f*MCF
    
    # decomposed_DOC = 0
    
    # # sum of sum of geometric series
    # for t in range(pl + 1):
    #     # sum of a geometric series where acc_DOC 
    #     acc_DOC = annual_DDOC * (1 - np.exp(-1 * k * t)) / (1 - np.exp(-1 * k)) 
    #     # all the acc_DOC at the start of the year is the only one contributing 
    #     # to decomposition in that one year
    #     decomposed_DOC += acc_DOC*(1 - np.exp(-1*k)) 
        
    #     # make annumpy  array from 0 to pl + 1 outside the for loop
        
    #     # replace t with DOC_ARRAY 
        
    t_vary = np.arange(pl)
    decomposed_DOC = annual_DDOC * (1 - np.exp(-1 * k * t_vary))
    total_decomposed_DOC = np.sum(decomposed_DOC)
    CH4_emitted_during_pl = total_decomposed_DOC*F*16/12
    
    accumulated_DOC_at_pl = annual_DDOC* (1 - np.exp(-1 * k * (pl-1))) / (1 - np.exp(-1 * k)) 
    CH4_emitted_after_pl = accumulated_DOC_at_pl*F*16/12
    
    days_in_year = 365

    return CH4_emitted_during_pl/(pl*days_in_year), CH4_emitted_after_pl/(pl*days_in_year)

def get_CH4_CO2_eq_treatment(system, influent_sc =None, effluent_sc = None,
                       CH4_CO2eq=29.8, CH4_EF_sc =0.0075):
    
    '''
    Parameters
    ----------
    system : :class:`biosteam.System`
        The system for which normalized GHG emission is being determined.
        
    ----Secondary treatment----
     
     influent_sc : : iterable[:class:`WasteStream`], optional
         Influent wastestreams to secondary treatment whose wastewater composition determine the potential for GHG emissions. The default is None.
     effluent_sc  : : iterable[:class:`WasteStream`], optional
         The wastestreams which represents the effluent from the secondary treatment process. The default is None.
     CH4_EF_sc :  float, optional.
         The emission factor used to calculate methane emissions in secondary treatment. The default is 0.0075 kg CH4/ kg rCOD. [1]
        
    --------- Eq CO2 EFs --------- 
   
    CH4_CO2eq : TYPE, optional
        DESCRIPTION. The default is 29.8 kg CO2eq/kg CH4 [1].

    Returns
    -------
    Total Normalized CH4 emissions from secondary treatment (kg CO2 eq./m3) [int].
    '''
    
    # source 1 (on-site)
    if influent_sc is None:
        influent_sc = [inf for inf in system.feeds if inf.phase == 'l']

    influent_sc_flow = np.array([inf.F_vol*24 for inf in influent_sc]) # in m3/day
    influent_sc_COD = np.array([inf.COD for inf in influent_sc]) # in mg/L
    mass_influent_COD_sc = np.sum(influent_sc_flow*influent_sc_COD/1000) # in kg/day
    
    effluent_sc_flow = np.array([eff.F_vol*24 for eff in effluent_sc])
    effluent_sc_COD =  np.array([eff.COD for eff in effluent_sc]) # in mg/L
    mass_effluent_COD_sc = np.sum(effluent_sc_flow*effluent_sc_COD/1000) # in kg/day
    
    mass_removed_COD_sc = mass_influent_COD_sc - mass_effluent_COD_sc
    CH4_emitted_sc = CH4_EF_sc*mass_removed_COD_sc
    
    # source 1 (on-site)
    CH4_CO2_eq_treatment = CH4_emitted_sc*CH4_CO2eq
    
    normalized_CO2_eq_WRRF = CH4_CO2_eq_treatment/sum([24*s.F_vol for s in system.feeds])
    
    return normalized_CO2_eq_WRRF

def get_N2O_CO2_eq_treatment(system, influent_sc =None, N2O_CO2eq=273, F=0.5,  N2O_EF_sc =0.016):
    
    '''
    Parameters
    ----------
    system : :class:`biosteam.System`
        The system for which normalized GHG emission is being determined.
        
    ----Secondary treatment----
     
     influent_sc : : iterable[:class:`WasteStream`], optional
         Influent wastestreams to secondary treatment whose wastewater composition determine the potential for GHG emissions. The default is None.
     N2O_EF_sc : float, optional
         The emission factor used to calculate nitrous oxide emissions in secondary treatment. The default is 0.016 kg N2O-N/ kg N. [1]
        
    --------- Eq CO2 EFs --------- 
   
    N2O_CO2eq : TYPE, optional
        DESCRIPTION. The default is 273 kg CO2eq/kg CH4 [1].

    Returns
    -------
    Total Normalized N2O emissions from secondary treatment (kg CO2 eq./m3) [int].
    '''
    
    # source 1 (on-site)
    if influent_sc is None:
        influent_sc = [inf for inf in system.feeds if inf.phase == 'l']

    influent_sc_flow = np.array([inf.F_vol*24 for inf in influent_sc]) # in m3/day
    
    influent_N_sc = np.array([inf.TN for inf in influent_sc]) # in mg/L
    mass_influent_N = np.sum(influent_sc_flow*influent_N_sc/1000) # in kg/day
    
    N2O_emitted_sc = N2O_EF_sc*mass_influent_N
    
    N2O_CO2_eq_treatment = N2O_emitted_sc*N2O_CO2eq 
    
    normalized_CO2_eq_WRRF = N2O_CO2_eq_treatment/sum([24*s.F_vol for s in system.feeds])
    
    return normalized_CO2_eq_WRRF

def get_CH4_CO2_eq_discharge(system, effluent_sys =None, CH4_CO2eq=29.8, CH4_EF_discharge=0.009):
    
    '''
    Parameters
    ----------
    system : :class:`biosteam.System`
        The system for which normalized GHG emission is being determined.
        
    ----Discharge----
        
    effluent_sys : : iterable[:class:`WasteStream`], optional
        Effluent wastestreams from the system whose wastewater composition determine the potential for GHG emissions at discharge. The default is None.
    CH4_EF_discharge :  float, optional.
        The emission factor used to calculate methane emissions in discharge. The default is 0.009 kg CH4/ kg effluent COD. [1]

    --------- Eq CO2 EFs --------- 
   
    CH4_CO2eq : TYPE, optional
        DESCRIPTION. The default is 29.8 kg CO2eq/kg CH4 [1].

    Returns
    -------
    Total Normalized GHG emissions from onsite and offsite operations associated with WRRF (kg CO2 eq./m3) [int].

    '''
    
    # source 3 (off-site)
    effluent_flow = np.array([eff.F_vol*24 for eff in effluent_sys]) # in m3/day
    effluent_COD = np.array([eff.COD for eff in effluent_sys]) # in mg/L
    mass_effluent_COD = np.sum(effluent_flow*effluent_COD/1000) # in kg/day
    
    CH4_emitted_discharge = CH4_EF_discharge*mass_effluent_COD

    # source 3 (off-site)
    CH4_CO2_eq_discharge = CH4_emitted_discharge*CH4_CO2eq
    
    normalized_total_CO2_eq_WRRF = CH4_CO2_eq_discharge/sum([24*s.F_vol for s in system.feeds])
    
    return normalized_total_CO2_eq_WRRF

def get_N2O_CO2_eq_discharge(system, effluent_sys =None, N2O_CO2eq=273, F=0.5, N2O_EF_discharge=0.005):
    '''
    Parameters
    ----------
    system : :class:`biosteam.System`
        The system for which normalized GHG emission is being determined.
        
    ----Discharge----
        
    effluent_sys : : iterable[:class:`WasteStream`], optional
        Effluent wastestreams from the system whose wastewater composition determine the potential for GHG emissions at discharge. The default is None.
    N2O_EF_discharge : float, optional
        The emission factor used to calculate nitrous oxide emissions in discharge. The default is 0.005 kg N2O-N/ kg effluent N. [1]
    
    --------- Eq CO2 EFs --------- 

    N2O_CO2eq : TYPE, optional
        DESCRIPTION. The default is 273 kg CO2eq/kg CH4 [1].

    Returns
    -------
    Total Normalized GHG emissions from onsite and offsite operations associated with WRRF (kg CO2 eq./m3) [int].
    

    ''' 
    
    # source 3 (off-site)
    effluent_flow = np.array([eff.F_vol*24 for eff in effluent_sys]) # in m3/day
    effluent_N = np.array([eff.TN for eff in effluent_sys]) # in mg/L
    mass_effluent_N = np.sum(effluent_flow*effluent_N/1000) # in kg/day
    
    N2O_emitted_discharge = N2O_EF_discharge*mass_effluent_N
    N2O_CO2_eq_discharge = N2O_emitted_discharge*N2O_CO2eq 
    normalized_total_CO2_eq_WRRF = N2O_CO2_eq_discharge/sum([24*s.F_vol for s in system.feeds])
    
    return normalized_total_CO2_eq_WRRF

def get_CH4_emitted_during_pl(system, sludge=None, CH4_CO2eq=29.8, F=0.5, DOC_f = 0.45, MCF = 0.8, k = 0.06, pl=30):
    '''
    Parameters
    ----------
    system : :class:`biosteam.System`
        The system for which normalized GHG emission is being determined.
        
    ----Sludge disposal---
    
    sludge : : iterable[:class:`WasteStream`], optional
        Effluent sludge from the system for which GHG emissions are being calculated. The default is None.
    DOC_f : float, optional
        fraction of DOC that can decompose (fraction). The default value is 0.5.
    MCF : float, optional
        CH4 correction factor for aerobic decomposition in the year of deposition (fraction). The default is 0.8.
    k : TYPE, optional
        Methane generation rate (k). The default is 0.185. (1/year)
        
        The decomposition of carbon is assumed to follow first-order kinetics 
        (with rate constant k), and methane generation is dependent on the amount of remaining decomposable carbon 
        in the waste. For North America (boreal and temperate climate) the default values for wet and dry climate 
        are:
            k (dry climate) = 0.06
            k (wet climate) = 0.185
    F : float, optional
        Fraction of methane in generated landfill gas (volume fraction).  The default is 0.5.
    pl : float, optional
        The project lifetime over which methane emissions would be calculated. (years)
        The default is 30 years.
        
    --------- Eq CO2 EFs --------- 
   
    CH4_CO2eq : TYPE, optional
        DESCRIPTION. The default is 29.8 kg CO2eq/kg CH4 [1].

    Returns
    -------
    Total Normalized CH4 emissions from sludge disposal during project lifetime (kg CO2 eq./m3) [int].
    
    '''
    DOC_mass_flow = 0
    
    for slg in sludge:
        DOC_mass_flow += slg.composite("C", flow=True, exclude_gas=True, 
                      subgroup=None, particle_size=None,
                      degradability="b", organic=True, volatile=None,
                      specification=None, unit="kg/day")

    annual_DOC_mass = 365*DOC_mass_flow # in kg/year

    annual_DDOC = annual_DOC_mass*DOC_f*MCF
        
    t_vary = np.arange(pl)
    decomposed_DOC = annual_DDOC * (1 - np.exp(-1 * k * t_vary))
    total_decomposed_DOC = np.sum(decomposed_DOC)
    CH4_emitted_during_pl = total_decomposed_DOC*F*16/12
    
    days_in_year = 365

    CH4_CO2_eq_sludge_disposal_pl = (CH4_emitted_during_pl/(pl*days_in_year))*CH4_CO2eq
    
    normalized_total_CO2_eq_WRRF = CH4_CO2_eq_sludge_disposal_pl/sum([24*s.F_vol for s in system.feeds])
    
    return normalized_total_CO2_eq_WRRF

def get_CH4_emitted_after_pl(system, sludge=None, CH4_CO2eq=29.8, N2O_CO2eq=273, F=0.5,
                     # uncertain parameters 
                     DOC_f = 0.45, MCF = 0.8, k = 0.06, pl=30
                     ):
    
    '''

    Parameters
    ----------
    system : :class:`biosteam.System`
        The system for which normalized GHG emission is being determined.
    
    ----Sludge disposal---
    
    sludge : : iterable[:class:`WasteStream`], optional
        Effluent sludge from the system for which GHG emissions are being calculated. The default is None.
    DOC_f : float, optional
        fraction of DOC that can decompose (fraction). The default value is 0.5.
    MCF : float, optional
        CH4 correction factor for aerobic decomposition in the year of deposition (fraction). The default is 0.8.
    k : TYPE, optional
        Methane generation rate (k). The default is 0.185. (1/year)
        
        The decomposition of carbon is assumed to follow first-order kinetics 
        (with rate constant k), and methane generation is dependent on the amount of remaining decomposable carbon 
        in the waste. For North America (boreal and temperate climate) the default values for wet and dry climate 
        are:
            k (dry climate) = 0.06
            k (wet climate) = 0.185
    F : float, optional
        Fraction of methane in generated landfill gas (volume fraction).  The default is 0.5.
    pl : float, optional
        The project lifetime over which methane emissions would be calculated. (years)
        The default is 30 years.
        
        
    --------- Eq CO2 EFs --------- 
   
    CH4_CO2eq : TYPE, optional
        DESCRIPTION. The default is 29.8 kg CO2eq/kg CH4 [1].

    Returns
    -------
    Total Normalized CH4 emissions from acculumulated sludge disposal after project lifetime (kg CO2 eq./m3) [int].
    

    '''
    # source 4 (off-site)

    DOC_mass_flow = 0
    
    for slg in sludge:
        DOC_mass_flow += slg.composite("C", flow=True, exclude_gas=True, 
                      subgroup=None, particle_size=None,
                      degradability="b", organic=True, volatile=None,
                      specification=None, unit="kg/day")

    annual_DOC_mass = 365*DOC_mass_flow # in kg/year

    annual_DDOC = annual_DOC_mass*DOC_f*MCF
        
    accumulated_DOC_at_pl = annual_DDOC* (1 - np.exp(-1 * k * (pl-1))) / (1 - np.exp(-1 * k)) 
    CH4_emitted_after_pl = accumulated_DOC_at_pl*F*16/12
    
    days_in_year = 365

    CH4_CO2_eq_sludge_disposal_pl = (CH4_emitted_after_pl/(pl*days_in_year))*CH4_CO2eq
    
    normalized_total_CO2_eq_WRRF = CH4_CO2_eq_sludge_disposal_pl/sum([24*s.F_vol for s in system.feeds])
    
    return normalized_total_CO2_eq_WRRF

def get_CO2_eq_WRRF (system, GHG_treatment, GHG_discharge, GHG_electricity, 
                     GHG_sludge_disposal, CH4_CO2eq=29.8, N2O_CO2eq=273):
    '''

    Parameters
    ----------
    system : :class:`biosteam.System`
        The system for which normalized GHG emission is being determined.
    GHG_treatment : tuple[int], optional
        The amount of methane and nitrous oxide emitted during secondary treatment (kg/day).
    GHG_discharge : tuple[int], optional
        The amount of methane and nitrous oxide emitted during effluent discharge (kg/day).
    GHG_electricity : float
        The amount of eq. CO2 emitted due to electrity consumption (kg-CO2-Eq/day).
    GHG_sludge_disposal : int
        The average amount of methane emitted during sludge disposal (kg/day).
    CH4_CO2eq : TYPE, optional
        DESCRIPTION. The default is 29.8 kg CO2eq/kg CH4 [1].
    N2O_CO2eq : TYPE, optional
        DESCRIPTION. The default is 273 kg CO2eq/kg CH4 [1].

    Returns
    -------
    Normalized GHG emissions from onsite and offsite operations associated with WRRF (kg CO2 eq./m3). [numpy array] 
    
    References
    ----------
    [1] IPCC 2021 – 6th Assessment Report Values. 

    '''
    
    # source 1 (on-site)
    CH4_CO2_eq_treatment = GHG_treatment[0]*CH4_CO2eq
    N2O_CO2_eq_treatment = GHG_treatment[1]*N2O_CO2eq 
    
    
    # source 3 (off-site)
    CH4_CO2_eq_discharge = GHG_discharge[0]*CH4_CO2eq
    N2O_CO2_eq_discharge = GHG_discharge[1]*N2O_CO2eq 
    
    # source 4 (off-site)
    CH4_CO2_eq_sludge_disposal_pl = GHG_sludge_disposal[0]*CH4_CO2eq
    CH4_CO2_eq_sludge_disposal_after_pl = GHG_sludge_disposal[1]*CH4_CO2eq
    
    # source 5 (off-site)
    CO2_eq_electricity = GHG_electricity*1
    
    CO2_eq_WRRF = np.array([CH4_CO2_eq_treatment, N2O_CO2_eq_treatment, #1
                            CH4_CO2_eq_discharge, N2O_CO2_eq_discharge, #3
                            CH4_CO2_eq_sludge_disposal_pl,              #4
                            CH4_CO2_eq_sludge_disposal_after_pl,        #4
                            CO2_eq_electricity])                        #5
    
    normalized_CO2_eq_WRRF = CO2_eq_WRRF/sum([24*s.F_vol for s in system.feeds])
    
    return normalized_CO2_eq_WRRF 


def get_CO2_eq_electricity(system, q_air, active_unit_IDs=None, p_atm=101.325, K=0.283, T=20,
                     # uncertain parameters 
                     P_inlet_loss=1, P_diffuser_loss=7, h_submergance=5.18, efficiency=0.7, CO2_EF=0.675):
    '''

    Parameters
    ----------
    system : :class:`biosteam.System`
        The system for which normalized GHG emission is being determined.
    
    ----Electricity---
    
    --blower power-
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
        Blower efficiency. Default is 0.7. 
    K : float, optional
        Equal to (kappa - 1)/kappa, where kappa is the adiabatic exponent of air.
        Default is 0.283, equivalent to kappa = 1.3947.
        
    --pump power--
    active_unit_IDs : tuple[str], optional
        IDs of all units whose power needs to be accounted for. The default is None.
    

    Returns
    -------
    Total Normalized CO2 emissions from electricity consumption (kg CO2 eq./m3) [int].
    

    '''    
    # source 5 (off-site)
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
    
    blower_power = WP/1000
    
    pumping_power = 0
        
    for y in system.flowsheet.unit:
        if y.ID in active_unit_IDs:
            pumping_power += y.power_utility.power # in kW
        
    total_energy_consumed = (blower_power + pumping_power)*24 # in kWh/day
    
    GHG_electricity = total_energy_consumed*CO2_EF # in kg-CO2-Eq/day
    
    # source 5 (off-site)
    CO2_eq_electricity = GHG_electricity*1
    
    normalized_total_CO2_eq_WRRF = CO2_eq_electricity/sum([24*s.F_vol for s in system.feeds])
    
    return normalized_total_CO2_eq_WRRF

def get_total_CO2_eq(system, q_air, influent_sc =None, effluent_sc = None, effluent_sys =None, active_unit_IDs=None, sludge=None, 
                      p_atm=101.325, K=0.283, CH4_CO2eq=29.8, N2O_CO2eq=273, F=0.5, 
                      
                      CH4_EF_sc =0.0075, N2O_EF_sc =0.016, CH4_EF_discharge=0.009, N2O_EF_discharge=0.005,
                      T=20, 
                     
                     # uncertain parameters 
                     P_inlet_loss=1, P_diffuser_loss=7, h_submergance=5.18, efficiency=0.7,
                     
                     CO2_EF=0.675, DOC_f = 0.45, MCF = 0.8, k = 0.06, pl=30
                     ):
    
    '''

    Parameters
    ----------
    system : :class:`biosteam.System`
        The system for which normalized GHG emission is being determined.
        
    ----Secondary treatment----
     
     influent_sc : : iterable[:class:`WasteStream`], optional
         Influent wastestreams to secondary treatment whose wastewater composition determine the potential for GHG emissions. The default is None.
     effluent_sc  : : iterable[:class:`WasteStream`], optional
         The wastestreams which represents the effluent from the secondary treatment process. The default is None.
     CH4_EF_sc :  float, optional.
         The emission factor used to calculate methane emissions in secondary treatment. The default is 0.0075 kg CH4/ kg rCOD. [1]
     N2O_EF_sc : float, optional
         The emission factor used to calculate nitrous oxide emissions in secondary treatment. The default is 0.016 kg N2O-N/ kg N. [1]
        
    ----Discharge----
        
    effluent_sys : : iterable[:class:`WasteStream`], optional
        Effluent wastestreams from the system whose wastewater composition determine the potential for GHG emissions at discharge. The default is None.
    CH4_EF_discharge :  float, optional.
        The emission factor used to calculate methane emissions in discharge. The default is 0.009 kg CH4/ kg effluent COD. [1]
    N2O_EF_discharge : float, optional
        The emission factor used to calculate nitrous oxide emissions in discharge. The default is 0.005 kg N2O-N/ kg effluent N. [1]
        
    
    ----Electricity---
    
    --blower power-
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
        Blower efficiency. Default is 0.7. 
    K : float, optional
        Equal to (kappa - 1)/kappa, where kappa is the adiabatic exponent of air.
        Default is 0.283, equivalent to kappa = 1.3947.
        
    --pump power--
    active_unit_IDs : tuple[str], optional
        IDs of all units whose power needs to be accounted for. The default is None.
        
    
    ----Sludge disposal---
    
    sludge : : iterable[:class:`WasteStream`], optional
        Effluent sludge from the system for which GHG emissions are being calculated. The default is None.
    DOC_f : float, optional
        fraction of DOC that can decompose (fraction). The default value is 0.5.
    MCF : float, optional
        CH4 correction factor for aerobic decomposition in the year of deposition (fraction). The default is 0.8.
    k : TYPE, optional
        Methane generation rate (k). The default is 0.185. (1/year)
        
        The decomposition of carbon is assumed to follow first-order kinetics 
        (with rate constant k), and methane generation is dependent on the amount of remaining decomposable carbon 
        in the waste. For North America (boreal and temperate climate) the default values for wet and dry climate 
        are:
            k (dry climate) = 0.06
            k (wet climate) = 0.185
    F : float, optional
        Fraction of methane in generated landfill gas (volume fraction).  The default is 0.5.
    pl : float, optional
        The project lifetime over which methane emissions would be calculated. (years)
        The default is 30 years.
        
        
    --------- Eq CO2 EFs --------- 
   
    CH4_CO2eq : TYPE, optional
        DESCRIPTION. The default is 29.8 kg CO2eq/kg CH4 [1].
    N2O_CO2eq : TYPE, optional
        DESCRIPTION. The default is 273 kg CO2eq/kg CH4 [1].

    Returns
    -------
    Total Normalized GHG emissions from onsite and offsite operations associated with WRRF (kg CO2 eq./m3) [int].
    

    '''
    
    # source 1 (on-site)
    if influent_sc is None:
        influent_sc = [inf for inf in system.feeds if inf.phase == 'l']

    influent_sc_flow = np.array([inf.F_vol*24 for inf in influent_sc]) # in m3/day
    influent_sc_COD = np.array([inf.COD for inf in influent_sc]) # in mg/L
    mass_influent_COD_sc = np.sum(influent_sc_flow*influent_sc_COD/1000) # in kg/day
    
    effluent_sc_flow = np.array([eff.F_vol*24 for eff in effluent_sc])
    effluent_sc_COD =  np.array([eff.COD for eff in effluent_sc]) # in mg/L
    mass_effluent_COD_sc = np.sum(effluent_sc_flow*effluent_sc_COD/1000) # in kg/day
    
    mass_removed_COD_sc = mass_influent_COD_sc - mass_effluent_COD_sc
    CH4_emitted_sc = CH4_EF_sc*mass_removed_COD_sc
    
    influent_N_sc = np.array([inf.TN for inf in influent_sc]) # in mg/L
    mass_influent_N = np.sum(influent_sc_flow*influent_N_sc/1000) # in kg/day
    
    N2O_emitted_sc = N2O_EF_sc*mass_influent_N
    
    # source 1 (on-site)
    CH4_CO2_eq_treatment = CH4_emitted_sc*CH4_CO2eq
    N2O_CO2_eq_treatment = N2O_emitted_sc*N2O_CO2eq 
    
    # source 3 (off-site)
    effluent_flow = np.array([eff.F_vol*24 for eff in effluent_sys]) # in m3/day
    effluent_COD = np.array([eff.COD for eff in effluent_sys]) # in mg/L
    mass_effluent_COD = np.sum(effluent_flow*effluent_COD/1000) # in kg/day
    
    CH4_emitted_discharge = CH4_EF_discharge*mass_effluent_COD
    
    effluent_N = np.array([eff.TN for eff in effluent_sys]) # in mg/L
    mass_effluent_N = np.sum(effluent_flow*effluent_N/1000) # in kg/day
    
    N2O_emitted_discharge = N2O_EF_discharge*mass_effluent_N
    
    # source 3 (off-site)
    CH4_CO2_eq_discharge = CH4_emitted_discharge*CH4_CO2eq
    N2O_CO2_eq_discharge = N2O_emitted_discharge*N2O_CO2eq 
    
    
    # source 5 (off-site)
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
    
    blower_power = WP/1000
    
    pumping_power = 0
        
    for y in system.flowsheet.unit:
        if y.ID in active_unit_IDs:
            pumping_power += y.power_utility.power # in kW
        
    total_energy_consumed = (blower_power + pumping_power)*24 # in kWh/day
    
    GHG_electricity = total_energy_consumed*CO2_EF # in kg-CO2-Eq/day
    
    # source 5 (off-site)
    CO2_eq_electricity = GHG_electricity*1
    
    # source 4 (off-site)

    if sludge == None:
        CH4_CO2_eq_sludge_disposal_pl = 0
    else:
        DOC_mass_flow = 0
        
        for slg in sludge:
            DOC_mass_flow += slg.composite("C", flow=True, exclude_gas=True, 
                          subgroup=None, particle_size=None,
                          degradability="b", organic=True, volatile=None,
                          specification=None, unit="kg/day")
    
        annual_DOC_mass = 365*DOC_mass_flow # in kg/year
    
        annual_DDOC = annual_DOC_mass*DOC_f*MCF
            
        t_vary = np.arange(pl)
        decomposed_DOC = annual_DDOC * (1 - np.exp(-1 * k * t_vary))
        total_decomposed_DOC = np.sum(decomposed_DOC)
        CH4_emitted_during_pl = total_decomposed_DOC*F*16/12
        
        accumulated_DOC_at_pl = annual_DDOC* (1 - np.exp(-1 * k * (pl-1))) / (1 - np.exp(-1 * k)) 
        CH4_emitted_after_pl = accumulated_DOC_at_pl*F*16/12
        
        days_in_year = 365
    
        GHG_sludge_disposal = (CH4_emitted_during_pl + CH4_emitted_after_pl)/(pl*days_in_year)
    
        CH4_CO2_eq_sludge_disposal_pl = GHG_sludge_disposal*CH4_CO2eq
    
    CO2_eq_WRRF = np.array([CH4_CO2_eq_treatment, N2O_CO2_eq_treatment, #1
                            CH4_CO2_eq_discharge, N2O_CO2_eq_discharge, #3
                            CH4_CO2_eq_sludge_disposal_pl,              #4
                            CO2_eq_electricity])                        #5
    
    normalized_total_CO2_eq_WRRF = np.sum(CO2_eq_WRRF)/sum([24*s.F_vol for s in system.feeds])
    
    return normalized_total_CO2_eq_WRRF

def estimate_ww_treatment_energy_demand(daily_energy_demand, daily_flow = 37854, ww_pcpd = 0.175):
    '''
    Parameters
    ----------
    daily_energy_demand : float
        Energy consumption at a centralized WRRF. (kWh/day)
    daily_flow : TYPE
        Daily wastewater flow. (m3/day)
    ww_pcpd : TYPE, float
        Average wastewater generated per capita per day. The default is 0.175. [1]

    Returns
    -------
    None.
    
    [1] Jones, E. R., Van Vliet, M. T., Qadir, M., & Bierkens, M. F. (2021). Country-level and 
    gridded estimates of wastewater production, collection, treatment and reuse. Earth System Science Data, 
    13(2), 237-254.

    '''
    ed_pcpd = (daily_energy_demand/daily_flow)*ww_pcpd
    
    return ed_pcpd


def estimate_N_removal_energy_demand(daily_energy_demand, effluent_N_conc, daily_flow = 37854, influent_N_conc = 40, 
                           per_capita_protein_intake = 68.6, N_in_pro = 0.13, N_excreted = 1):
    '''
    Parameters
    ----------
    daily_energy_demand : float
        Energy consumption at a centralized WRRF. (kWh/day)
    daily_flow : TYPE
        Daily wastewater flow. (m3/day)
    influent_N_conc : float
        TN concentration in the influent. (mg/L)
    effluent_N_conc : float
        N concentration in the effluent. (mg/L)
    per_capita_protein_intake : float
        Per capita protein intake. (g/cap/day)
    N_in_pro : float
        % of N in protein. (%)
    N_excreted : float
        % of N intake that is excreted. (%)

    Returns
    -------
    None.

    '''
    daily_ed_N_removal = daily_energy_demand/(daily_flow*(influent_N_conc - effluent_N_conc))
    
    per_capita_per_day_N = per_capita_protein_intake*N_in_pro*N_excreted
    
    return daily_ed_N_removal*per_capita_per_day_N

def estimate_P_removal_energy_demand(daily_energy_demand, effluent_TP_conc, 
                                     daily_flow = 37854, influent_TP_conc = 7, 
                                     pc_ani_protein_intake = 12.39, 
                                     pc_plant_protein_intake = 40.29, 
                                     P_ani_pro = 0.011, 
                                     P_plant_pro = 0.022, 
                                     P_excreted = 1):
    '''
    Parameters
    ----------
    daily_energy_demand : float
        Energy consumption at a centralized WRRF. (kWh/day)
    effluent_TP_conc : float
        TN concentration in the effluent (mg/L). 
    daily_flow : TYPE, optional
        Daily wastewater flow (m3/day) The default is 37854.
    influent_TP_conc : float, optional
        TN concentration in the influent (mg/L). The default is 7.
    pc_ani_protein_intake : float, optional
        DESCRIPTION. The default is 12.39.
    pc_plant_protein_intake : TYPE, optional
        DESCRIPTION. The default is 40.29.
    P_ani_pro : TYPE, optional
        DESCRIPTION. The default is 0.011.
    P_plant_pro : TYPE, optional
        DESCRIPTION. The default is 0.022.
    P_excreted : TYPE, optional
        DESCRIPTION. The default is 1.

    Returns
    -------
    None.

    '''
    daily_ed_P_removal = daily_energy_demand/(daily_flow*(influent_TP_conc - effluent_TP_conc))
    per_capita_per_day_P = (pc_ani_protein_intake*P_ani_pro + pc_plant_protein_intake*P_plant_pro)*P_excreted
    
    return daily_ed_P_removal*per_capita_per_day_P
    