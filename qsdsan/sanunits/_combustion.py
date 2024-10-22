#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    
    Yalin Li <mailto.yalin.li@gmail.com>
    
    Jianan Feng <jiananf2@illinois.edu>

Part of this module is based on the BioSTEAM package:
https://github.com/bioSTEAMDevelopmentGroup/biosteam

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''


# %%

import biosteam as bst
from warnings import warn
from flexsolve import IQ_interpolation
from biosteam import HeatUtility, Facility
from thermosteam.reaction import ParallelReaction
from .. import SanUnit, Construction
from ..utils import sum_system_utility

__all__ = ('BiogasCombustion', 'CombinedHeatPower',)


class BiogasCombustion(SanUnit):
    '''
    Use of biogas in combustion.

    Note that this unit does not include cost calculation and generation of electricity,
    use `CHP` instead for such purposes.

    Parameters
    ----------
    if_combustion : bool
        If include combustion reaction during simulation.
    biogas_loss : float
        Fraction of biogas loss (e.g., leaked).
    biogas_eff : float
        Combustion efficiency of biogas as a fraction of CH4.

    Examples
    --------
    `bwaise systems <https://github.com/QSD-Group/EXPOsan/blob/main/exposan/bwaise/systems.py>`_
    '''
    _N_ins = 2
    _N_outs = 3

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 if_combustion=False, biogas_loss=0.1, biogas_eff=0.55):

        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.if_combustion = if_combustion
        self._biogas_loss = biogas_loss
        self._biogas_eff = biogas_eff


    def _run(self):
        biogas, air = self.ins
        biogas.phase = air.phase = 'g'
        for i in self.outs:
            i.copy_like(biogas)
        used, lost, wasted = self.outs
        lost.mass *= self.biogas_loss
        used.mass -= lost.mass
        wasted.mass = used.mass * (1-self.biogas_eff)
        used.mass -= wasted.mass
        if self.if_combustion:
            rxn = self.chemicals.get_combustion_reactions()
            rxn.force_reaction(used.mol)
            air.imol['O2'] = -used.imol['O2']
            used.imol['O2'] = 0.
            air.imol['N2'] = 0.79/0.21 * air.imol['O2']
        else:
            air.empty()


    @property
    def biogas_loss(self):
        '''[float] Fraction of biogas loss (i.e., leaked).'''
        return self._biogas_loss
    @biogas_loss.setter
    def biogas_loss(self, i):
        self._biogas_loss = i

    @property
    def biogas_eff(self):
        '''[float] Combustion efficiency of biogas as a fraction of CH4.'''
        return self._biogas_eff
    @biogas_eff.setter
    def biogas_eff(self, i):
        self._biogas_eff = i


class CombinedHeatPower(SanUnit, Facility):
    '''
    Combustion of all feed streams with simple estimation of the capital cost
    of a combined heat and power (CHP) unit based on Shoener et al. [1]_
    and Havukainen et al. [2]_
    
    This unit is designed to always satisfy the system's heating demand,
    if there is not enough energy in the feed stream, natural gas will be purchased
    to fill the gap.
    
    The unit can also be set to satisfy the system power (i.e., electricity) needs
    if `supplement_power_utility` is set to True.

    Parameters
    ----------
    ins : Iterable(obj):
        Feed stream to be combusted, optional natural gas stream, air.
    outs : Iterable(obj):
        Gas emission, solid ash.
    unit_CAPEX : float
        Capital cost of the CHP per kW of power generated, $/kW.
    CHP_type : str
        Type of the CHP to adjust the heat-to-power efficiency.
        Can be "Fuel cell" (40.5%), "Microturbine" (27%),
        "Internal combustion" (36%), "Combustion gas" (31.5%).
        Note that when `CHP_type` is provided, the set `combined_efficiency`
        will be ignored.
    combustion_eff : float
        Efficiency of the boiler to convert the energy embedded in the feed
        to the heating steam.
    combined_eff : float
        Combined heat (the total embedded energy in feed streams)-to-power efficiency
        (i.e., combustion efficiency * power generation efficiency).
    system : obj
        The linked system whose heating/power utility needs will be supplied
        by this CHP unit.
    supplement_power_utility : bool
        Whether to purchase additional natural gas to supplement energy demand.
    
        .. note::
            
            If energy in the feed (without natural gas) is sufficient enough,
            the system could still be producing power even when
            `supplement_power_utility` is False.

    References
    ----------
    .. [1] Shoener et al., Design of Anaerobic Membrane Bioreactors for the
        Valorization of Dilute Organic Carbon Waste Streams.
        Energy Environ. Sci. 2016, 9 (3), 1102–1112.
        https://doi.org/10.1039/C5EE03715H.
    .. [2] Havukainen, J.; Nguyen, M. T.; Väisänen, S.; Horttanainen, M.
        Life Cycle Assessment of Small-Scale Combined Heat and Power Plant:
        Environmental Impacts of Different Forest Biofuels and Replacing
        District Heat Produced from Natural Gas. Journal of Cleaner
        Production 2018, 172, 837–846.
        https://doi.org/10.1016/j.jclepro.2017.10.241.
    '''
    _N_ins = 3
    _N_outs = 2
    network_priority = 0
    default_combined_eff = {
        'Fuel cell': 0.405,
        'Microturbine': 0.27,
        'Internal combustion': 0.36,
        'Combustion gas': 0.315,
        }
    _units = {'Steel': 'kg',
              'Furnace': 'kg',
              'Concrete': 'kg',
              'Reinforcing steel': 'kg'}

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 include_construction=True, unit_CAPEX=1225, CHP_type='Fuel cell',
                 combustion_eff=0.8, combined_eff=None,
                 system=None, supplement_power_utility=False, F_BM={'CHP': 1.}):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with, F_BM=F_BM, include_construction=include_construction)
        self.unit_CAPEX = unit_CAPEX
        self.CHP_type = CHP_type
        self.combustion_eff = combustion_eff
        self.combined_eff = combined_eff
        self.system = None
        self.supplement_power_utility = supplement_power_utility
        self._sys_heating_utilities = ()
        self._sys_steam_utilities = ()
        self._sys_power_utilities = ()
        
    def _init_lca(self):
        self.construction = [
            Construction('carbon_steel', linked_unit=self, item='Carbon_steel', quantity_unit='kg'),
            Construction('furnace', linked_unit=self, item='Furnace', quantity_unit='kg'),
            Construction('concrete', linked_unit=self, item='Concrete', quantity_unit='kg'),
            Construction('reinforcing_steel', linked_unit=self, item='Reinforcing_steel', quantity_unit='kg'),
            ]

    # Only simulate in design stage to ensure capturing all system-wise utility
    def _run(self): pass


    def _design(self):
        feed, natural_gas, air = self.ins
        emission, ash = self.outs
        for i in (natural_gas, air, ash):
            i.empty()
        feed.phase = natural_gas.phase = air.phase = emission.phase = 'g'
        ash.phase = 's'
        emission.P = ash.P = 101325
        emission.T = ash.T = 298.15
        self._refresh_sys()

        cmps = self.components
        rxns = []
        for cmp in cmps:
            if cmp.locked_state in ('l', 's') and (not cmp.organic or cmp.degradability=='Undegradable'):
                continue
            rxn = cmp.get_combustion_reaction()
            if rxn:
                rxns.append(rxn)
        combustion_rxns = self.combustion_reactions = ParallelReaction(rxns)

        def react(natural_gas_flow=0):
            emission.copy_flow(feed)
            emission.imol['CH4'] += natural_gas_flow
            natural_gas.imol['CH4'] = natural_gas_flow
            combustion_rxns.force_reaction(emission.mol)
            air.imol['O2'] = -emission.imol['O2']
            emission.imol['N2'] = air.imol['N2'] = air.imol['O2']/0.21*0.79
            emission.imol['O2'] = 0
            H_net_feed = feed.H + feed.HHV - emission.H # subtracting the energy in emission
            if natural_gas.imol['CH4'] != 0: # add natural gas H and HHV
                H_net_feed += natural_gas.H + natural_gas.HHV
            return H_net_feed

        # Calculate the amount of energy in the feed (no natural gas) and needs
        self.H_net_feeds_no_natural_gas = react(0)
        
        # Calculate all energy needs in kJ/hr as in H_net_feeds
        kwds = dict(system=self.system, operating_hours=self.system.operating_hours, exclude_units=(self,))
        pu = self.power_utility
        H_steam_needs = sum_system_utility(**kwds, utility='steam', result_unit='kJ/hr')/self.combustion_eff
        H_power_needs = sum_system_utility(**kwds, utility='power', result_unit='kJ/hr')/self.combined_eff
        
        # Calculate the amount of energy needs to be provided
        H_supp = H_steam_needs+H_power_needs if self.supplement_power_utility else H_steam_needs
                      
        # Objective function to calculate the heat deficit at a given natural gas flow rate
        def H_deficit_at_natural_gas_flow(flow):
            return H_supp-react(flow)
        # Initial lower and upper bounds for the solver
        lb = 0
        ub = react()/cmps.CH4.LHV*2
        if H_deficit_at_natural_gas_flow(0) > 0: # energy in the feeds is not enough
            while H_deficit_at_natural_gas_flow(ub) > 0: # increase bounds if not enough energy
                lb = ub
                ub *= 2
            natural_gas_flow = IQ_interpolation(
                H_deficit_at_natural_gas_flow,
                x0=lb, x1=ub, xtol=1e-3, ytol=1,
                checkbounds=False)
            H_net_feeds = react(natural_gas_flow)
        else: # enough energy in the feed, set natural_gas_flow to 0
            H_net_feeds = react(0)

        # Update heating utilities
        self.steam_utilities = HeatUtility.sum_by_agent(sum(self.sys_steam_utilities.values(), []))
        ngu = self.natural_gas_utilities = HeatUtility.sum_by_agent(sum(self.sys_natural_gas_utilities.values(), []))
        self.heat_utilities = HeatUtility.sum_by_agent(sum(self.sys_heating_utilities.values(), []))
        natural_gas.imol['CH4'] += sum(i.flow for i in ngu) # natural gas is added on separately
        for hu in self.heat_utilities: hu.reverse()
        
        # Power production if there is sufficient energy
        if H_net_feeds <= H_steam_needs:
            pu.production = 0
        else:
            pu.production = (H_net_feeds-H_steam_needs)/3600*self.combined_eff

        self.H_steam_needs = H_steam_needs
        self.H_power_needs = H_power_needs
        self.H_net_feeds = H_net_feeds

        ash_IDs = [i.ID for i in cmps if not i.formula]
        ash.copy_flow(emission, IDs=tuple(ash_IDs), remove=True)
        
        if self.include_construction:
            D = self.design_results
            constr = self.construction
            
            # material calculation based on [2], linearly scaled on power (kW)
            # in [2], a 580 kW CHP:
            # steel: 20098 kg
            # furnace: 12490 kg
            # reinforced concrete: 15000 kg (concrete + reinforcing steel)
            # 1 m3 reinforced concrete: 98 v/v% concrete with a density of 2500 kg/m3 (2450 kg)
            #                            2 v/v% reinforcing steel with a density of 7850 kg/m3 (157kg)
            factor = self.H_net_feeds/3600/580
            constr[0].quantity = D['Steel'] = factor*20098
            constr[1].quantity = D['Furnace'] = factor*12490
            constr[2].quantity = D['Concrete'] = factor*15000*2450/(2450 + 157)
            constr[3].quantity = D['Reinforcing steel'] = factor*15000*157/(2450 + 157)


    def _cost(self):
        unit_CAPEX = self.unit_CAPEX
        unit_CAPEX /= 3600 # convert to $ per kJ/hr
        self.baseline_purchase_costs['CHP'] = unit_CAPEX * self.H_net_feeds
        # Update biosteam utility costs
        uprices = bst.stream_utility_prices
        uprices['Fuel'] = uprices['Natural gas'] = self.ins[1].price
        uprices['Ash disposal'] = self.outs[1].price

    def _refresh_sys(self):
        sys = self._system
        ng_dct = self._sys_natural_gas_utilities = {}
        steam_dct = self._sys_steam_utilities = {}
        pu_dct = self._sys_power_utilities = {}
        if sys:
            units = [u for u in sys.units if u is not self]
            for u in units:
                pu = u.power_utility
                if pu: pu_dct[u.ID] = pu
                steam_utilities = []
                for hu in u.heat_utilities:
                    if hu.flow*hu.duty <= 0: continue # cooling utilities
                    if hu.ID=='natural_gas': ng_dct[u.ID] = [hu]                    
                    elif 'steam' in hu.ID: steam_utilities.append(hu)
                    else: raise ValueError(f'The heating utility {hu.ID} is not recognized by the CHP.')
                if steam_utilities: steam_dct[u.ID] = steam_utilities
        sys_hus = {k:ng_dct.get(k,[])+steam_dct.get(k, [])
                   for k in list(ng_dct.keys())+list(steam_dct.keys())}
        self._sys_heating_utilities = sys_hus

    @property
    def fuel_price(self):
        '''
        [Float] Price of fuel (natural gas), set to be the same as the price of ins[1]
        and `bst.stream_utility_prices['Natural gas']`.
        '''
        return self.ins[1].price

    natural_gas_price = fuel_price
    
    @property
    def ash_disposal_price(self):
        '''
        [Float] Price of ash disposal, set to be the same as the price of outs[1]
        and `bst.stream_utility_prices['Ash disposal']`.
        Negative means need to pay for ash disposal.
        '''
        
        """[Float] Price of ash disposal, same as `bst.stream_utility_prices['Ash disposal']`."""
        return self.outs[1].price

    @property
    def CHP_type(self):
        '''
        [str] Type of the combined heat and power (CHP) type
        to adjust the heat-to-power efficiency.
        Can be "Fuel cell" (40.5%), "Microturbine" (27%),
        "Internal combustion" (36%), "Combustion gas" (31.5%),
        or None and define an efficiency.
        '''
        return self._CHP_type
    @CHP_type.setter
    def CHP_type(self, i):
        eff = self.combined_eff if hasattr(self, '_combined_eff') else None
        if i is None and eff:
            self._CHP_type = i
            return
        if eff:
            warn(f'With the provided `CHP_type` ({i}), '
                 f'the current `combined_eff` {eff} will be ignored.')
        self._CHP_type = i

    @property
    def combustion_eff(self):
        '''
        [float] Combined heat (the total embedded energy in feed streams)-to-power efficiency
        (i.e., combustion efficiency * power generation efficiency).
        '''
        return self._combustion_eff
    @combustion_eff.setter
    def combustion_eff(self, i):
        if i is None:
            self._combustion_eff = i
            return
        if hasattr(self, '_combined_eff'):
            combined_eff = self.combined_eff
            if combined_eff:
                if i < combined_eff:
                    raise ValueError(f'`combustion_eff` ({i}) should be larger than '
                                     f'`combined_eff` ({combined_eff}).')
        self._combustion_eff = i

    @property
    def combined_eff(self):
        '''
        [float] Combined heat (the total embedded energy in feed streams)-to-power efficiency
        (i.e., combustion efficiency * power generation efficiency).
        Note that when `CHP_type` is provided, the set `combined_efficiency`
        will be ignored.
        '''
        if self.CHP_type:
            return self.default_combined_eff[self.CHP_type]
        return self._combined_eff
    @combined_eff.setter
    def combined_eff(self, i):
        if i is None:
            self._combined_eff = i
            return
        CHP_type = self.CHP_type
        if CHP_type:
            warn(f'With the provided `combined_eff` ({i}), '
                 f'the `CHP_type` ("{CHP_type}") is ignored.')
            self._combined_eff = i
            return
        if hasattr(self, '_combustion_eff'):
            combustion_eff = self.combustion_eff
            if combustion_eff:
                if i > combustion_eff:
                    raise ValueError(f'`combined_eff` ({i}) should be smaller than '
                                     f'`combustion_eff` ({combustion_eff}).')
        self._combined_eff = i

    @property
    def system(self):
        '''
        [obj] The linked system whose heating/power utility needs will be supplied
        by this CHP unit.
        '''
        return self._system
    @system.setter
    def system(self, i):
        self._system = i
        self._refresh_sys()

    @property
    def sys_heating_utilities(self):
        '''[dict] Heating utilities (steams and natural gas) of the given system (excluding this CHP unit).'''
        return self._sys_heating_utilities

    @property
    def sys_natural_gas_utilities(self):
        '''[dict] Steam utilities of the given system (excluding this CHP unit).'''
        return self._sys_natural_gas_utilities

    @property
    def sys_steam_utilities(self):
        '''[dict] Steam utilities of the given system (excluding this CHP unit).'''
        return self._sys_steam_utilities

    @property
    def sys_power_utilities(self):
        '''[dict] Power utilities of the given system (excluding this CHP unit).'''
        return self._sys_power_utilities

    @property
    def H_needs(self):
        '''
        Total energy needs of the system in kJ/hr,
        already divided by the combustion/combined efficiency.
        '''
        return self.H_heating_needs + self.H_power_needs