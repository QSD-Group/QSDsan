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
           'get_total_CO2_eq')


#%%
    
def get_SRT(system, biomass_IDs, wastage=None, active_unit_IDs=None):
    """
    Estimate sludge residence time (SRT) [day] of an activated sludge system.

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

def get_oxygen_heterotrophs(flow, influent_COD, eff_COD_soluble, 
                            f_d=0.1, b_H=0.4, SRT=10, Y_H=0.625):
    """
    Estimates oxygen requirement [kg-O2/day] for heterotrophic biological processes in 
    an activated sludge system given design and performance assumptions, 
    following equation 10.10 in [1].

    Parameters
    ----------
    flow : float
        Volumetric flow rate through the system [m3/d].
    influent_COD : float
        Influent COD concentration [mg/L].
    eff_COD_soluble : float
        Maximum effluent soluble COD concentration [mg/L]. 
    f_d : float, optional
        Fraction of biomass that remains as cell debris after decay. 
        Default value is 0.1 gCOD/gCOD, based on ASM2d 'f_XI_H'. 
    b_H : float, optional
        Decay rate constant of heterotrophs [d^(-1)]. The default is 0.4, based on ASM2d.
    SRT : float, optional
        Design sludge retention time of the system. Default is 10 days. 
    Y_H : float, optional
        Yield of heterotrophs [gCOD/gCOD]. The default is 0.625.
            
    References
    ----------
    [1] Grady Jr, C.P.L.; Daigger, G.T.; Love, N.G.; Filipe, C.D.M. Biological 
    Wastewater Treatment. 3rd ed. Boca Raton, FL: CRC Press 2011.
    """    
    mass_COD_treated = flow * (influent_COD - eff_COD_soluble) * 1e-3 # kg/day
    aeration_factor = 1 - (1 + f_d*b_H*SRT)*Y_H/(1 + b_H*SRT)
    
    return mass_COD_treated*aeration_factor

def get_oxygen_autotrophs(flow, influent_COD, eff_COD_soluble, influent_TKN,
                          f_d=0.1, b_H=0.4, b_AUT=0.15, SRT=10, 
                          Y_H=0.625, Y_AUT=0.24, K_NH=1, U_AUT=1, SF_DO=1.375):
    """
    Estimates oxygen requirement [kg-O2/day] for autotrophic biological processes in 
    an activated sludge system given design and performance assumptions, 
    following equations 11.16-11.19 in [1].
    
    Parameters
    ----------
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
        0.1 gCOD/gCOD, based on ASM2d 'f_XI_A'. 
    b_H : float
        Decay rate constant of heterotrophs [d^(-1)]. The default is 0.4, based on ASM2d.
    b_AUT : float
        Decay rate constant of autotrophs [d^(-1)]. The default is 0.15, based on ASM2d.
    SRT : float
        Design sludge retention time of the system. Default is 10 days. 
    Y_H : float
        Yield of heterotrophs [gCOD/gCOD]. The default is 0.625.
    Y_AUT : float
        Yield of autotrophs [g COD/g N]. The default is 0.24.
    K_NH: float
        Ammonium (nutrient) half saturation coefficient for autotrophs, in [g N/m^3].
        The default is 1.0.
    U_AUT: float
        Autotrophic maximum specific growth rate, in [d^(-1)]. The default is 1.0.
    SF_DO: float
        Safety factor for dissolved oxygen. The default is 1.375. 

    References
    ----------
    [1] Grady Jr, C.P.L.; Daigger, G.T.; Love, N.G.; Filipe, C.D.M. 
    Biological Wastewater Treatment. 3rd ed. Boca Raton, FL: CRC Press 2011.
    """
    
    NR = 0.087*(1 + f_d*b_H*SRT)*Y_H/(1 + b_H*SRT)
    TKN = influent_TKN
    S_N_a = TKN - NR*(influent_COD - eff_COD_soluble)
    S_NH = K_NH*(1/SRT  + b_AUT)/(U_AUT/SF_DO - (1 + b_AUT/SRT))
    aeration_factor = 4.57 - (1 + f_d*b_AUT*SRT)*Y_AUT/(1 + b_AUT*SRT)
    mass_N_removed = flow*(S_N_a - S_NH)/1000 # kg/day
    
    return mass_N_removed*aeration_factor

def get_airflow(oxygen_heterotrophs, oxygen_autotrophs, oxygen_transfer_efficiency=12):
    """
    Estimates diffused aeration air flow rate [m3/min] of an activated sludge 
    system based on oxygen requirements, following equation 11.2 in [1].
    
    Parameters
    ----------
    oxygen_heterotrophs : float
        Oxygen requirement for heterotrophic biological processes, in kg-O2/day.
    oxygen_autotrophs : float
        Oxygen requirement for autotrophic biological processes, in kg-O2/day.
    oxygen_transfer_efficiency : float
        Field oxygen transfer efficiency in percentage. The default is 12. 
    
    References
    ----------
    [1] Grady Jr, C.P.L.; Daigger, G.T.; Love, N.G.; Filipe, C.D.M. 
    Biological Wastewater Treatment. 3rd ed. Boca Raton, FL: CRC Press 2011.

    """
    
    required_oxygen = (oxygen_heterotrophs + oxygen_autotrophs)/24 # in kg/hr
    Q_air = 6*required_oxygen/oxygen_transfer_efficiency
    
    return Q_air

def get_P_blower(q_air, T=20, P_atm=101.325, P_inlet_loss=1, 
                 P_diffuser_loss=7, h_submergance=5.18, efficiency=0.7,
                 K=0.283):
    """
    Estimates blower power requirement [kW] for diffused aeration, following
    equation 5-77 in [1] and equation 4.27 in [2].
    
    Parameters
    ----------
    q_air : float
        Air volumetric flow rate [m3/min].
    T : float
        Air temperature [degree Celsius].
    P_atm : float
        Atmostpheric pressure [kPa]
    P_inlet_loss : float
        Head loss at inlet [kPa].
    P_diffuser_loss : float
        Head loss due to piping and diffuser [kPa].
    h_submergance : float
        Diffuser submergance depth in m. The default is 17 feet (5.18 m)
    efficiency : float
        Blower efficiency. The default is 0.7, usual range is 0.7 to 0.9. 
    K : float, optional
        Equal to (kappa - 1)/kappa, where kappa is the adiabatic exponent of air,
        i.e., the specific heat ratio. For single-stage centrifugal blower, 
        the default is 0.283, equivalent to kappa = 1.3947 for dry air.

    References
    ----------
    [1] Metcalf & Eddy, Wastewater Engineering: Treatment and 
        Resource Recovery. 5th ed. New York: McGraw-Hill Education. 2014.
    [2] Mueller, J.; William C.B.; and Popel, H.J. Aeration: 
        Principles and Practice, Volume 11. CRC press, 2002.
    """
    
    T += 273.15
    air_molar_vol = 22.4e-3 * (T/273.15)/(P_atm/101.325)   # m3/mol
    R = 8.31446261815324 # J/mol/K, universal gas constant
    
    p_in = P_atm - P_inlet_loss
    p_out = P_atm + 9.81*h_submergance + P_diffuser_loss
    
    Q_air = q_air/60 # m3/min to m3/s
    WP = (Q_air / air_molar_vol) * R * T / K * ((p_out/p_in)**K - 1) / efficiency  # in W
    return WP/1000

def get_power_utility(system, active_unit_IDs=None):
    '''
    Total power of the specified unit operations [kW].
    
    Parameters
    ----------
    system : :class:`biosteam.System`
        The system whose power will be calculated.
    active_unit_IDs : tuple[str], optional
        IDs of all units whose power needs to be accounted for. The default is None.
    '''

    return sum([y.power_utility.power for y in system.flowsheet.unit if y.ID in active_unit_IDs])        

def get_cost_sludge_disposal(sludge, unit_weight_disposal_cost=375):
    '''
    Returns the daily operating cost of sludge treatment and disposal [USD/day].
    Typical sludge disposal unit costs:
    
        Land application: 300 - 800 USD/US ton. [2]
        Landfill: 100 - 650 USD/US ton. [2]
        Incineration: 300 - 500 USD/US ton. [2]
    
    Parameters
    ----------
    sludge : iterable[:class:`WasteStream`]
        Effluent sludge from the system for which treatment and disposal costs are being calculated.
        The default is None.
    unit_weight_disposal_cost : float, optional
        The sludge treatment and disposal cost per unit weight [USD/ton].
        Feasible range for this value lies between 100-800 USD/ton [1]. 
        The default is 375 USD/US ton, which is the close to average of landfill.
    
    References
    ----------
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
    return np.sum(sludge_prod)*unit_weight_disposal_cost   #in USD/day
    

def get_normalized_energy(system, aeration_power, pumping_power, miscellaneous_power):
    '''
    Normalized operational energy consumption associated with WRRF [kWh/m3].
    
    Parameters
    ----------
    system : :class:`biosteam.System`
        The system for which normalized energy consumption is being determined.
    aeration_power : float
        Power of blower [kW].
    pumping_power : float
        Power rquired for sludge pumping and other equipments [kW].
    miscellaneous_power : float
        Any other power requirement in the system [kW].

    '''
    Q = sum([s.F_vol for s in system.feeds])
    return np.array([aeration_power, pumping_power, miscellaneous_power])/Q


def get_daily_operational_costs(aeration_power, pumping_power, miscellaneous_power, \
                                sludge_disposal_cost, unit_electricity_cost=0.161):
    '''
    Normalized daily operational costs associated with WRRF [USD/day], in the
    order of aeration, sludge pumping, sludge disposal, and miscellaneous. 
    
    Parameters
    ----------
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

    References
    ----------
    [1] https://www.bls.gov/regions/midwest/data/averageenergyprices_selectedareas_table.htm
    
    '''
    aeration_cost = aeration_power*24*unit_electricity_cost # in (kWh/day)*(USD/kWh) = USD/day
    pumping_cost = pumping_power*24*unit_electricity_cost # in (kWh/day)*(USD/kWh) = USD/day
    miscellaneous_cost = miscellaneous_power*24*unit_electricity_cost # in (kWh/day)*(USD/kWh) = USD/day
    
    return np.array([aeration_cost, pumping_cost, sludge_disposal_cost, miscellaneous_cost])                        #5
    
get_daily_operational_cost = get_daily_operational_costs

def get_total_operational_cost(q_air, # aeration (blower) power 
                               sludge,  # sludge disposal costs 
                               system, active_unit_IDs=None,  # sludge pumping power 
                               T=20, P_atm=101.325, P_inlet_loss=1, P_diffuser_loss=7, 
                               h_submergance=5.18, efficiency=0.7, K=0.283, # aeration (blower) power 
                               miscellaneous_power=0, 
                               unit_weight_disposal_cost=375, # sludge disposal costs 
                               unit_electricity_cost=0.161): 
    '''
    Normalized daily operational cost associated with WRRF [USD/day].
    
    Parameters
    ----------
    
    q_air : float
        Air volumetric flow rate for diffused aeration [m3/min].
    sludge : iterable[:class:`WasteStream`], optional
        Effluent sludge from the system for which treatment and disposal costs are being calculated.
        The default is None.
    system : :class:`biosteam.System`
        The system whose power will be calculated.
    active_unit_IDs : tuple[str], optional
        IDs of all units whose power needs to be accounted for. The default is None.
    T : float
        Air temperature [degree Celsius].
    P_atm : float
        Atmostpheric pressure [kPa]
    P_inlet_loss : float
        Head loss at aeration blower inlet [kPa]. The default is 1 kPa. 
    P_diffuser_loss : float
        Head loss due to piping and diffuser [kPa]. The default is 7 kPa. 
    h_submergance : float
        Diffuser submergance depth in m. The default is 17 feet (5.18 m)
    efficiency : float
        Blower efficiency. Default is 0.7. 
    K : float, optional
        Equal to (kappa - 1)/kappa, where kappa is the adiabatic exponent of air.
        Default is 0.283, equivalent to kappa = 1.3947.  
    unit_weight_disposal_cost : float
        The sludge treatment and disposal cost per unit weight [USD/ton].
        Feasible range for this value lies between 100-800 USD/ton. 
        The default is 375 USD/US ton, which is the close to average of landfill.
    miscellaneous_power : float, optional
        Any other power requirement in the system [kW].
    unit_electricity_cost : float
        Unit cost of electricity. Default value is taken as $0.161/kWh. 
    
    '''
    aeration_power = get_P_blower(q_air, T, P_atm, P_inlet_loss, P_diffuser_loss,
                                  h_submergance, efficiency, K)           
    pumping_power = get_power_utility(system, active_unit_IDs)
    sludge_disposal_cost = get_cost_sludge_disposal(sludge, unit_weight_disposal_cost)
    opex = get_daily_operational_costs(aeration_power, pumping_power, 
                                       miscellaneous_power,
                                       sludge_disposal_cost, 
                                       unit_electricity_cost)
    return sum(opex)

def get_GHG_emissions_sec_treatment(system=None, influent=None, effluent=None,
                                    CH4_EF=0.0075, N2O_EF=0.016):
    '''    
    Returns a 2-tuple of the fugitive emissions of CH4 and N2O [kg/day] 
    during secondary treatment.
    
    Parameters
    ----------
    system : :class:`biosteam.System`
        The system for which emissions during secondary treatment are being 
        calculated. The default is None.
    influent : iterable[:class:`WasteStream`], optional
        Influent wastewater to secondary treatment. The default is None.
    effluent  : iterable[:class:`WasteStream`], optional
        Effluent wastewater from the secondary treatment process. The default is None.
    CH4_EF : float, optional.
        The emission factor used to calculate methane emissions in secondary 
        treatment. The default is 0.0075 kg CH4/ kg COD removed. [1]
    N2O_EF : float, optional
        The emission factor used to calculate nitrous oxide emissions in
        secondary treatment. The default is 0.016 kg N2O-N/ kg N. [1]
        
    References
    ----------
    [1] IPCC (2019). Chapter 6 in 2019 Refinement to the 2006 IPCC Guidelines 
        for National Greenhouse Gas Inventories.

    '''
    
    if influent is None:
        influent = [inf for inf in system.feeds if inf.phase == 'l']
    mass_influent_COD = sum(inf.F_vol*24*inf.COD for inf in influent)/1000 # in kg/day
    mass_effluent_COD = sum(eff.F_vol*24*eff.COD for eff in effluent)/1000 # in kg/day

    CH4_emitted = CH4_EF*(mass_influent_COD - mass_effluent_COD)
    
    mass_influent_N = sum(inf.F_vol*24*inf.TN for inf in influent)/1000 # in kg/day    
    N2O_emitted = N2O_EF*mass_influent_N
    
    return CH4_emitted, N2O_emitted

def get_GHG_emissions_discharge(effluent=None, CH4_EF=0.009, N2O_EF=0.005):
    '''   
    Returns a 2-tuple of the fugitive emissions of CH4 and N2O [kg/day] 
    at discharge.
    
    Parameters
    ----------
    effluent : : iterable[:class:`WasteStream`], optional
        Effluent wastewater discharged from the system. The default is None.
    CH4_EF_discharge :  float, optional.
        The emission factor used to calculate methane emissions in discharge. 
        The default is 0.009 kg CH4/ kg effluent COD. [1]
    N2O_EF_discharge : float, optional
        The emission factor used to calculate nitrous oxide emissions in discharge. 
        The default is 0.005 kg N2O-N/ kg effluent N. [1]
        
    References
    ----------
    [1] IPCC, 2019. Chapter 6 in 2019 Refinement to the 2006 IPCC Guidelines 
        for National Greenhouse Gas Inventories.
    '''
    mass_effluent_COD = sum(eff.F_vol*24*eff.COD for eff in effluent)/1000 # in kg/day
    CH4_emitted = CH4_EF*mass_effluent_COD
    mass_effluent_N = sum(eff.F_vol*24*eff.TN for eff in effluent)/1000 # in kg/day    
    N2O_emitted = N2O_EF*mass_effluent_N
    
    return CH4_emitted, N2O_emitted
    
def get_GHG_emissions_electricity(system, power_blower, power_pump, CO2_EF=0.675):
    '''
    Returns GHG emission associated with operational electricity consumption [kg CO2-eq/day].
    
    Parameters
    ----------
    system : :class:`biosteam.System`
        The system for which tier-2 GHG emissions due to electricity consumption 
        are being calculated. 
    power_blower : float
        Power of blower [kW].
    power_pump : float
        Power required for pumping and other utilities at treatment facility [kW].
    CO2_EF : float
        The emission factor used to calculate scope-2 CO2 emissions due to electricity consumption. 
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

    '''
    
    total_energy_consumed = (power_blower + power_pump)*24 # in kWh/day
    CO2_emissions = total_energy_consumed*CO2_EF # in kg-CO2-Eq/day
    
    return CO2_emissions

def get_GHG_emissions_sludge_disposal(sludge=None, DOC_f=0.38, MCF=0.8, 
                                      k=0.06, F=0.5, pl=30):
    '''
    The average amount of methane emitted from sludge disposed in landfill [kg/day],
    returned in a 2-tuple representing emissions during and after project lifetime,
    respectively.
    
    Parameters
    ----------
    sludge : iterable[:class:`WasteStream`], optional
        Effluent sludge from the system for which GHG emissions are being 
        calculated. The default is None.
    DOC_f : float, optional
        fraction of DOC that can decompose. The default value is 0.5.
    MCF : float, optional
        CH4 correction factor for aerobic decomposition in the year of 
        deposition. The default is 0.8.
    k : float, optional
        Methane generation rate [yr^(-1)]. The default is 0.185. 
        The decomposition of carbon is assumed to follow 1st-order kinetics 
        (with rate constant k), and methane generation is dependent on the 
        amount of remaining decomposable carbon in the waste. 
        For North America (boreal and temperate climate) the default values are:
            k (dry climate) = 0.06
            k (wet climate) = 0.185
    F : float, optional
        Volume fraction of methane in generated landfill gas. The default is 0.5.
    pl : float, optional
        The project lifetime [yr] over which methane emissions would be calculated. 
        The default is 30 years.
        
    References
    ----------
    [1] IPCC, 2019. Chapter 3: Solid Waste Disposal, in 2019 Refinement to
        the 2006 IPCC Guidelines for National Greenhouse Gas Inventories.
    '''
    DOC_mass_flow = 0
    for slg in sludge:
        DOC_mass_flow += slg.composite("C", flow=True, exclude_gas=True, 
                      degradability="b", organic=True, unit="kg/day")

    annual_DOC_mass = 365*DOC_mass_flow # in kg/year
    annual_DDOC = annual_DOC_mass*DOC_f*MCF
        
    t_vary = np.arange(pl)
    decomposed_DOC = annual_DDOC * (1 - np.exp(-k * t_vary))
    CH4_emitted_during_pl = sum(decomposed_DOC)*F*16/12
    
    accumulated_DOC_at_pl = annual_DDOC* (1 - np.exp(-k * (pl-1))) / (1 - np.exp(-k)) 
    CH4_emitted_after_pl = accumulated_DOC_at_pl*F*16/12
    
    return CH4_emitted_during_pl/(pl*365), CH4_emitted_after_pl/(pl*365)

def get_CO2_eq_WRRF(system, GHG_treatment, GHG_discharge, GHG_electricity, 
                    GHG_sludge_disposal, CH4_CO2eq=29.8, N2O_CO2eq=273):
    '''
    Normalized GHG emissions from onsite and offsite operations associated 
    with WRRF [kg CO2 eq./m3].

    Parameters
    ----------
    system : :class:`biosteam.System`
        The system for which normalized GHG emission is being determined.
    GHG_treatment : tuple[float], optional
        The amount of methane and nitrous oxide emitted during secondary treatment (kg/day).
    GHG_discharge : tuple[float], optional
        The amount of methane and nitrous oxide emitted during effluent discharge (kg/day).
    GHG_electricity : float
        The amount of eq. CO2 emitted due to electrity consumption (kg-CO2-Eq/day).
    GHG_sludge_disposal : int
        The average amount of methane emitted during sludge disposal (kg/day).
    CH4_CO2eq : float, optional
        Conversion factor of CH4 to equivalent CO2. The default is 29.8 kg CO2eq/kg CH4 [1].
    N2O_CO2eq : float, optional
        Conversion factor of N2O to equivalent CO2. The default is 273 kg CO2eq/kg CH4 [1].
    
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

def get_total_CO2_eq(
        system, q_air, influent_sc=None, effluent_sc=None, effluent_sys=None, 
        active_unit_IDs=None, sludge=None, P_atm=101.325, K=0.283, 
        CH4_CO2eq=29.8, N2O_CO2eq=273, CH4_EF_sc=0.0075, N2O_EF_sc=0.016, 
        CH4_EF_discharge=0.009, N2O_EF_discharge=0.005, T=20, F=0.5, 
        P_inlet_loss=1, P_diffuser_loss=7, h_submergance=5.18, efficiency=0.7,
        CO2_EF=0.675, DOC_f=0.38, MCF=0.8, k=0.06, pl=30
        ):
    
    '''
    Returns the total normalized GHG emissions from onsite and offsite operations 
    associated with WRRF [kg CO2 eq./m3].
    
    Parameters
    ----------
    system : :class:`biosteam.System`
        The system for which normalized GHG emission is being determined.
        
    ----Secondary treatment----
     
    influent_sc : iterable[:class:`WasteStream`], optional
        Influent wastewater to secondary treatment. The default is None.
    effluent_sc : iterable[:class:`WasteStream`], optional
        Effluent wastewater from the secondary treatment process. The default is None.
    CH4_EF_sc :  float, optional.
        The emission factor used to calculate methane emissions in secondary 
        treatment. The default is 0.0075 kg CH4/ kg rCOD.
    N2O_EF_sc : float, optional
        The emission factor used to calculate nitrous oxide emissions in 
        secondary treatment. The default is 0.016 kg N2O-N/ kg N.
        
    ----Discharge----
        
    effluent_sys : iterable[:class:`WasteStream`], optional
        Effluent wastewater discharged from the system. The default is None.
    CH4_EF_discharge : float, optional.
        The emission factor used to calculate methane emissions in discharge. 
        The default is 0.009 kg CH4/ kg effluent COD.
    N2O_EF_discharge : float, optional
        The emission factor used to calculate nitrous oxide emissions in 
        discharge. The default is 0.005 kg N2O-N/ kg effluent N.
    
    ----Electricity---
    CO2_EF : float
        The emission factor used to calculate scope-2 CO2 emissions due to 
        electricity consumption. The default is 0.675 kg-CO2-Eq/kWh.
    
    --blower power-
    q_air : float
        Air volumetric flow rate for diffused aeration [m3/min].
    T : float
        Air temperature [degree Celsius].
    P_atm : float
        Atmostpheric pressure [kPa]
    P_inlet_loss : float
        Head loss at aeration blower inlet [kPa]. The default is 1 kPa. 
    P_diffuser_loss : float
        Head loss due to piping and diffuser [kPa]. The default is 7 kPa. 
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
    
    sludge : iterable[:class:`WasteStream`], optional
        Effluent sludge from the system for which GHG emissions are being 
        calculated. The default is None.
    DOC_f : float, optional
        fraction of DOC that can decompose. The default value is 0.5.
    MCF : float, optional
        CH4 correction factor for aerobic decomposition in the year of 
        deposition. The default is 0.8.
    k : float, optional
        Methane generation rate [yr^(-1)]. The default is 0.185. 
        The decomposition of carbon is assumed to follow 1st-order kinetics 
        (with rate constant k), and methane generation is dependent on the 
        amount of remaining decomposable carbon in the waste. 
        For North America (boreal and temperate climate) the default values are:
            k (dry climate) = 0.06
            k (wet climate) = 0.185
    F : float, optional
        Volume fraction of methane in generated landfill gas. The default is 0.5.
    pl : float, optional
        The project lifetime [yr] over which methane emissions would be calculated. 
        The default is 30 years.
        
    --------- Eq CO2 EFs --------- 
   
    CH4_CO2eq : float, optional
        Conversion factor of CH4 to equivalent CO2. The default is 29.8 kg CO2eq/kg CH4.
    N2O_CO2eq : float, optional
        Conversion factor of N2O to equivalent CO2. The default is 273 kg CO2eq/kg CH4.

    '''
    
    # source 1 (on-site)
    CH4_treatment, N2O_treatment = get_GHG_emissions_sec_treatment(
        system, influent_sc, effluent_sc, CH4_EF_sc, N2O_EF_sc
        )
    
    # source 3 (off-site)
    CH4_discharge, N2O_discharge = get_GHG_emissions_discharge(
        effluent_sys, CH4_EF_discharge, N2O_EF_discharge
        )
    
    # source 5 (off-site)
    blower_power = get_P_blower(q_air, T, P_atm, P_inlet_loss, P_diffuser_loss,
                                h_submergance, efficiency, K)           
    pumping_power = get_power_utility(system, active_unit_IDs)
    CO2_eq_electricity = (blower_power + pumping_power)*24*CO2_EF # in kg-CO2-Eq/day
    
    # source 4 (off-site)
    CH4_sludge_disposal = get_GHG_emissions_sludge_disposal(
        sludge, DOC_f, MCF, k, F, pl
        )
    
    CO2_eq_WRRF = np.sum([CH4_treatment*CH4_CO2eq, N2O_treatment*N2O_CO2eq, #1
                          CH4_discharge*CH4_CO2eq, N2O_discharge*N2O_CO2eq, #3
                          sum(CH4_sludge_disposal)*CH4_CO2eq,               #4
                          CO2_eq_electricity])                              #5
    
    return CO2_eq_WRRF/sum([24*s.F_vol for s in system.feeds])