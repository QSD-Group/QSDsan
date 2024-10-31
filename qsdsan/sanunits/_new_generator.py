#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems
This module is developed by:
    Shion Watabe <shionwatabe@gmail.com>
    Yalin Li <mailto.yalin.li@gmail.com>
    Hannah Lohman <hlohman94@gmail.com>
    Lewis Rowles <stetsonsc@gmail.com>
    Tori Morgan <vlmorgan@illinois.edu>

This module contains unit operations used in the NEWgenerator system
as described in `Shyu et al. <https://doi.org/10.1016/j.jenvman.2021.112921>`_

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

import numpy as np
from warnings import warn
from .. import SanUnit, Construction
from ._non_reactive import Copier
from ..processes._decay import Decay
from ..utils import ospath, load_data, data_path, price_ratio

__all__ = (
    'NEWgeneratorAnMBR',
    'NEWgeneratorChlorination',
    'NEWgeneratorControls',
    'NEWgeneratorFoundation',
    'NEWgeneratorGridtied',
    'NEWgeneratorHousing',
    'NEWgeneratorIonExchange',
    'NEWgeneratorPhotovoltaic',
    'NEWgeneratorPretreatment',
    )

ng_data_path = ospath.join(data_path, 'sanunit_data/ng')


# %%

anmbr_path = ospath.join(ng_data_path, '_ng_anmbr.csv')

class NEWgeneratorAnMBR(SanUnit, Decay):
    '''
    Anaerobic membrane bioreactor of the NEWgenerator
    consisting of a anaerobic baffled reactor (ABR)
    and an ultrafiltration membrane, the ABR is based on Trimmer et al.


    Parameters
    ----------
    ins : Iterable(stream)
        waste: waste for treatment.
    outs : Iterable(stream)
        treated: treated waste.
        sludge: sludge (biomass) produced from anaerobic reactor.
        biogas: biogas produced from anaerobic reactor.
        CH4: fugitive methane emissions from biomass decay.
        N2O: fugitive nitrous oxide emissions from biomass decay.
        soluble CH4: fugitive methane emissions from CH4 dissolved in liquid effluent.
    sludge_moisture_content : float
        Moisture content of the sludge (mass of water/total mass).
    if_gridtied : bool
        If electricity source is grid-tied energy.
    user_scale_up : float
        Scaling up of consumables, electricity demand, capital parts,
        and replacement parts based on number of input users compared to the
        baseline design number of users

    References
    ----------
    [1] Trimmer et al., Navigating Multidimensional Social–Ecological System
    Trade-Offs across Sanitation Alternatives in an Urban Informal Settlement.
    Environ. Sci. Technol. 2020, 54 (19), 12641–12653.
    https://doi.org/10.1021/acs.est.0c03296.
    [2] Shyu et al., The NEWgeneratorTM Non-Sewered Sanitation System:
    Long-Term Field Testing at an Informal Settlement Community in
    EThekwini Municipality, South Africa.
    Journal of Environmental Management 2021, 296, 112921.
    https://doi.org/10.1016/j.jenvman.2021.112921.

    See Also
    --------
    :ref:`qsdsan.processes.Decay <processes_Decay>`
    '''
    _N_ins = 1
    _N_outs = 5

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 MCF_decay=0.1, N_max_decay=0.8, N2O_EF_decay=0.0005,
                 sludge_moisture_content=0.93,
                 if_gridtied=True, if_HFMC_membrane=True,
                 user_scale_up=1, **kwargs):           
        Decay.__init__(self, ID, ins, outs, thermo=thermo,
                       init_with=init_with, F_BM_default=1,
                       if_capture_biogas=True,
                       if_N2O_emission=True,)
        
        self.MCF_decay = MCF_decay
        self.N_max_decay = N_max_decay
        self.N2O_EF_decay = N2O_EF_decay
        self.sludge_moisture_content = sludge_moisture_content
        self.if_gridtied = if_gridtied
        self.price_ratio = 1
        self.user_scale_up = user_scale_up
        self.if_HFMC_membrane = if_HFMC_membrane
        
        data = load_data(path=anmbr_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data
        
        self._refres_lca()
        
        for attr, value in kwargs.items():
            setattr(self, attr, value)
            
    def _refres_lca(self):
        self.construction = [
            Construction(item='PE', linked_unit=self, quantity_unit='kg'),
            Construction(item='PVC', linked_unit=self, quantity_unit='kg'),
            Construction(item='Plastic', linked_unit=self, quantity_unit='kg'),
            Construction(item='Pump', linked_unit=self, quantity=1, quantity_unit='ea', lifetime=self.anmbr_permeate_pump_lifetime),
            Construction(item='Pump', linked_unit=self, quantity=2, quantity_unit='ea', lifetime=self.anmbr_membrane_pump_lifetime),
            Construction(item='Membrane', linked_unit=self, quantity_unit='kg', lifetime=self.anmbr_membrane_module_lifetime),
            ]
        if self.if_HFMC_membrane:
            self.construction.append(
                Construction(item='Membrane', linked_unit=self, quantity_unit='kg', lifetime=self.HFMC_lifetime),
                )


    def _run(self):
        waste = self.ins[0]
        treated, sludge, biogas, CH4, N2O = self.outs
        treated.copy_like(waste)
        sludge.empty()
        biogas.phase = CH4.phase = N2O.phase = 'g'

        ##### COD removal #####
        COD_removal = self.anmbr_COD_removal
        removed_COD = waste.COD/1e3 * waste.F_vol * COD_removal  # kg/hr
        sludge_prcd = self.anmbr_sludge_yield * removed_COD  # produced biomass

        sludge.copy_flow(treated, ('Mg', 'Ca', 'OtherSS', 'Tissue', 'WoodAsh'), remove=True)
        sludge.imass['OtherSS'] += sludge_prcd
        sludge.imass['H2O'] = sludge.imass['OtherSS']/(1-self.sludge_moisture_content) - sludge.imass['OtherSS']

        # Adjust water
        treated.imass['H2O'] -= sludge.imass['H2O']
        if treated.imass['H2O'] <= 0:  # not enough water
            sludge.imass['H2O'] = waste.imass['H2O']
            treated.imass['H2O'] = 0

        # Methane
        CH4_prcd = removed_COD * self.anmbr_methane_yield * self.methane_density  # kg CH4 produced/hr
        CH4_soluble = self.anmbr_soluble_methane_fraction * CH4_prcd
        if self.if_HFMC_membrane:
            treated.imass['SolubleCH4'] = CH4_soluble * (1 - self.dissolved_methane_HFMC_removed)
            recovered_soluble_CH4 = CH4_soluble * self.dissolved_methane_HFMC_removed
        else:
            treated.imass['SolubleCH4'] = CH4_soluble
            recovered_soluble_CH4 = 0
        CH4_gas = CH4_prcd - CH4_soluble + recovered_soluble_CH4
        if self.if_capture_biogas:
            biogas.imass['CH4'] = CH4_gas
            CH4.empty()
        else:
            CH4.imass['CH4'] = CH4_gas
            biogas.empty()
        treated.imass['SolubleCH4'] = CH4_soluble  # Soluble CH4 left in liquid stream

        ##### Nutrient removal #####
        # N2O
        N_removal = self.anmbr_TN_removal
        if self.if_N2O_emission:
            # Assume the removal part covers both degradation loss
            # and other unspecified removal as well
            N_loss = self.first_order_decay(k=self.decay_k_N,
                                            t=self.anmbr_tau/365,
                                            max_decay=self.N_max_decay)
            if N_loss > N_removal:
                warn(f'Nitrogen degradation loss ({N_loss:.2f}) larger than the given removal ({N_removal:.2f})), '
                      'the degradation loss is assumed to be the same as the removal.')
                N_loss = N_removal

            # N2O only from the degraded part
            N_loss_tot = N_loss*waste.TN/1e3*waste.F_vol
            N2O.imass['N2O'] = N_loss_tot*self.N2O_EF_decay*44/28
        else:
            N2O.empty()

        total_solubles = np.array([
            CH4_soluble,
            waste.imass['NH3']*(1-N_loss), # assume the same removal for NH3 and NonNH3
            waste.imass['NonNH3']*(1-N_loss),
            waste.imass['P'],
            waste.imass['K'],
            ])

        # Removed solubles in the sludge, assume minimal used for growth
        sludge_solubles = total_solubles * np.array([
            0, N_removal, N_removal, self.anmbr_TP_removal, 0,
            ])

        liquid_solubles = total_solubles - sludge_solubles

        # Allocate remaining nutrients between effluent and sludge based on the water ratio
        ratio = treated.imass['H2O']/waste.imass['H2O']
        treated.imass['SolubleCH4', 'NH3', 'NonNH3', 'P', 'K'] = liquid_solubles * ratio
        sludge.imass['SolubleCH4', 'NH3', 'NonNH3', 'P', 'K'] = liquid_solubles*(1-ratio) + sludge_solubles

        # COD as in the liquid part
        treated._COD = sludge._COD = waste.COD * (1-self.anmbr_COD_removal)

    def _design(self):
        design = self.design_results
        constr = self.construction
        design['PE'] = constr[0].quantity = self.anmbr_PE_weight
        design['PVC'] = constr[1].quantity = (self.anmbr_PVC_weight + self.anmbr_piping_PVC_weight)
        design['Plastic'] = constr[2].quantity = self.anmbr_plastic_weight
        design['Membrane'] = constr[-2].quantity = self.anmbr_membrane_weight
        design['Membrane'] = constr[-1].quantity = self.HFMC_weight
        self.add_construction(add_cost=False)

    def _cost(self):
        C = self.baseline_purchase_costs
        user_scale_up = self.user_scale_up
        C['Reactor'] = (self.anmbr_reactor_tank + self.anmbr_reactor_media + self.anmbr_reactor_fittings)
        C['MembraneModule'] = (self.anmbr_membrane_module + self.anmbr_membrane_fittings)
        C['MembranePumps'] = (self.anmbr_permeate_pump + self.anmbr_backwash_pump + self.anmbr_water_pump)*user_scale_up
        C['Biogas'] = (self.anmbr_biogas_gas_bag + self.anmbr_biogas_fittings + self.anmbr_biogas_ventilation)*user_scale_up
        C['Piping'] = self.anmbr_piping*user_scale_up
        if self.if_HFMC_membrane:
            C['HFMC'] = self.HFMC
        else:
            C['HFMC'] = 0
        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio

        self.add_OPEX = self._calc_replacement_cost() + self._calc_maintenance_labor_cost()

        if self.if_gridtied:
            power_demand = (
                self.anmbr_water_pump_energy_percycle +
                self.anmbr_permeate_pump_energy_percycle +
                self.anmbr_backwash_pump_energy_percycle
                ) * self.anmbr_batch_cycle_perday * self.user_scale_up
            power_demand /= 24  # convert from kWh/d to kW
        else:
            power_demand = 0

        self.power_utility(power_demand)  # kW

    def _calc_replacement_cost(self):
        if self.if_HFMC_membrane:
            anmbr_replacement_cost = ((self.anmbr_membrane_module / self.anmbr_membrane_module_lifetime)
                                      + (self.anmbr_permeate_pump*self.user_scale_up/ self.anmbr_permeate_pump_lifetime)
                                      + (self.anmbr_backwash_pump*self.user_scale_up/ self.anmbr_membrane_pump_lifetime)
                                      + (self.anmbr_water_pump*self.user_scale_up/ self.anmbr_membrane_pump_lifetime)
                                      + (self.HFMC/self.HFMC_lifetime))
        else:
            anmbr_replacement_cost = ((self.anmbr_membrane_module / self.anmbr_membrane_module_lifetime)
                                      + (self.anmbr_permeate_pump*self.user_scale_up/ self.anmbr_permeate_pump_lifetime)
                                      + (self.anmbr_backwash_pump*self.user_scale_up/ self.anmbr_membrane_pump_lifetime)
                                      + (self.anmbr_water_pump*self.user_scale_up/ self.anmbr_membrane_pump_lifetime))
        return anmbr_replacement_cost / (365 * 24) * self.price_ratio  # USD/hr

    def _calc_maintenance_labor_cost(self):
        anmbr_maintenance_labor_cost = ((self.anmbr_labor_maintenance_membrane_chem_cleaning * self.wages)
                                        + (self.anmbr_labor_maintenance_sludge_removal * self.wages)
                                        + (self.anmbr_labor_replacement_membrane_pump * self.wages*self.user_scale_up)
                                        + (self.anmbr_labor_replacement_membrane_pump * self.wages*self.user_scale_up)
                                        + (self.anmbr_labor_replacement_permeate_pump * self.wages*self.user_scale_up)
                                        + (self.anmbr_labor_replacement_membrane_module * self.wages)
                                        + (self.anmbr_labor_replacement_misc_repairs * self.wages))
        return anmbr_maintenance_labor_cost/(365 * 24)  # USD/hr
    
    @property
    def if_HFMC_membrane(self):
        return self._if_HFMC_membrane
    @if_HFMC_membrane.setter
    def if_HFMC_membrane(self, i):
        self._if_HFMC_membrane = i
        constr = self.construction
        if constr:
            for i in constr: i.registry.discard(i)
            self._refres_lca()


# %%

chlorination_path = ospath.join(ng_data_path, '_ng_chlorination.tsv')

@price_ratio()
class NEWgeneratorChlorination(SanUnit):
    '''
    Electrochlorination of the NEWgenerator consisting a chloroalkaline generator cell,
    where the brine (NaCl) solution produces chlorine gas.

    Parameters
    ----------
    ins : Iterable(stream)
        waste: waste for treatment.
        NaCl: brine solution from NaCl and water as chloroalkaline cell input.
    outs : Iterable(stream)
        treated: treated waste.
    if_gridtied : bool
        If electricity source is grid-tied energy.
    user_scale_up : float
        Scaling up of consumables, electricity demand, capital parts,
        and replacement parts based on number of input users compared to the
        baseline design number of users.

    References
    ----------
    [1] Shyu et al., The NEWgeneratorTM Non-Sewered Sanitation System:
    Long-Term Field Testing at an Informal Settlement Community in
    EThekwini Municipality, South Africa.
    Journal of Environmental Management 2021, 296, 112921.
    https://doi.org/10.1016/j.jenvman.2021.112921.
    '''
    _N_ins = 2
    _N_outs = 2

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 if_gridtied=True, user_scale_up=1,  **kwargs):

        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with,
                         F_BM_default=1)
        self.if_gridtied = if_gridtied
        self.user_scale_up = user_scale_up

        data = load_data(path=chlorination_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)
            
    def _init_lca(self):
        self.construction = [
            Construction(item='PE', linked_unit=self, quantity_unit='kg'),
            Construction(item='PVC', linked_unit=self, quantity_unit='kg'),
            Construction(item='ElectronicsActive', linked_unit=self, quantity_unit='kg'),
            ]

    def _run(self):
        waste, NaCl = self.ins
        NaCl.phase = 'l'
        treated, solubleCH4 = self.outs
        treated.copy_like(waste)
        solubleCH4.empty()

        # Moving soluble CH4 in effluent to separate output
        solubleCH4.copy_flow(treated, 'SolubleCH4', remove=True)

        # COD removal
        treated._COD *= (1-self.chlorination_COD_removal)

        # TN removal
        treated.imass['N'] = waste.imass['N']*(1-self.chlorination_TN_removal)

        # TP removal
        treated.imass['P'] = waste.imass['P']*(1-self.chlorination_TP_removal)

        treated.imass['OtherSS'] *= (1-self.chlorination_COD_removal)

        # NaCl
        self.NaCl_demand_time = (self.chlorination_NaCl_dosage * self.chlorination_throughput_flow)/1000000/24
        NaCl.imass['SodiumChloride'] = self.NaCl_demand_time*self.user_scale_up

    def _design(self):
        design = self.design_results
        constr = self.construction
        design['PE'] = constr[0].quantity = self.chlorination_PE_weight
        design['PVC'] = constr[1].quantity = (self.chlorination_PVC_weight + self.chlorination_piping_PVC_weight)
        design['ElectronicsActive'] = constr[2].quantity = self.chlorination_electronics_active_weight * (self.user_scale_up ** 0.6)
        self.add_construction(add_cost=False)

    def _cost(self):
        C = self.baseline_purchase_costs
        user_scale_up = self.user_scale_up
        C['RecyclePump'] = self.chlorination_recycle_pump * (user_scale_up ** 0.6)
        C['Electrochlorinator'] = self.chlorination_electochlorinator * (user_scale_up ** 0.6)
        C['PoweredValve'] = self.chlorination_powered_valve
        C['Tanks'] = self.chlorination_tanks
        C['Fittings'] = self.chlorination_fittings
        C['Piping'] = self.chlorination_piping
        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio

        self.add_OPEX = self._calc_replacement_cost() + self._calc_maintenance_labor_cost()

        if self.if_gridtied:
            power_demand = (
                (self.chlorination_chlorinator_energy_percycle * self.user_scale_up) +
                (self.chlorination_pump_energy_percycle * self.user_scale_up) +
                self.chlorination_ventilation_energy_percycle +
                self.chlorination_valve_energy_percycle +
                (self.chlorination_product_water_pump_energy_percycle * self.user_scale_up)
            ) * self.chlorination_batch_cycle_per_day
            power_demand /= 24  # convert from kWh/d to kW
        else:
            power_demand = 0

        self.power_utility(power_demand)  # kW

    def _calc_replacement_cost(self):
        chlorination_replacement_cost = (
            self.chlorination_recycle_pump*(self.user_scale_up**0.6) / self.chlorination_other_pump_lifetime +
            self.chlorination_electochlorinator*(self.user_scale_up**0.6) / self.chlorination_electrochlorinator_lifetime +
            self.chlorination_powered_valve / self.chlorination_solenoid_valve_lifetime
            )
        return chlorination_replacement_cost / (365 * 24) * self.price_ratio  # USD/hr

    def _calc_maintenance_labor_cost(self):
        chlorination_maintenance_labor = (
            (self.chlorination_labor_maintenance_chlorinator_recharging * self.user_scale_up) +
            (self.chlorination_labor_replacement_chlorinator_replacement * self.user_scale_up) +
            (self.chlorination_labor_replacement_other_pump * self.user_scale_up) +
            self.chlorination_labor_replacement_misc_repairs
            ) * self.wages
        return chlorination_maintenance_labor/(365 * 24) # USD/hr



# %%

controls_path = ospath.join(ng_data_path, '_ng_controls.tsv')

@price_ratio()
class NEWgeneratorControls(Copier):
    '''
    Controlling devices of the NEWgenerator, which is composed of
    a controls box, a microcontroller, a sensors box,
    an ORP (Oxidation-reduction potential) box,
    a pump, and additional supplies/sensors.

    This is a non-reactive unit (i.e., the effluent is copied from the influent).

    Parameters
    ----------
    if_gridtied : bool
        If electricity source is grid-tied energy.

    References
    ----------
    [1] Shyu et al., The NEWgeneratorTM Non-Sewered Sanitation System:
    Long-Term Field Testing at an Informal Settlement Community in
    EThekwini Municipality, South Africa.
    Journal of Environmental Management 2021, 296, 112921.
    https://doi.org/10.1016/j.jenvman.2021.112921.

    See Also
    --------
    :class:`~.sanunits.Copier`
    '''

    def __init__(self, ID='', ins=None, outs=(),  thermo=None, init_with='WasteStream',
                 if_gridtied=False, **kwargs):
        data = load_data(path=controls_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        Copier.__init__(self, ID, ins, outs, thermo, init_with)
        self.if_gridtied = if_gridtied
        
        for attr, value in kwargs.items():
            setattr(self, attr, value)

            
    def _init_lca(self):
        self.construction = [
            Construction(item='ControlUnits', linked_unit=self, quantity_unit='kg',
                         lifetime=self.control_system_PLC_lifetime),
            Construction(item='Polycarbonate', linked_unit=self, quantity_unit='kg'),
            Construction(item='Aluminum', linked_unit=self, quantity_unit='kg'),
            Construction(item='ElectricCables', linked_unit=self, quantity_unit='m'),
            Construction(item='ElectronicsPassive', linked_unit=self, quantity_unit='kg'),
            Construction(item='ElectronicsActive', linked_unit=self, quantity_unit='kg'),
            ]

    def _design(self):
        design = self.design_results
        constr = self.construction
        design['ControlUnits'] = constr[0].quantity = self.control_control_units_weight
        design['Polycarbonate'] = constr[1].quantity = self.control_polycarbonate_weight
        design['Aluminum'] = constr[2].quantity = self.control_aluminum_weight
        design['ElectricCables'] = constr[3].quantity = self.control_cable_length
        design['ElectronicsPassive'] = constr[4].quantity = self.control_electronics_passive_weight
        design['ElectronicsActive'] = constr[5].quantity = self.control_electronics_active_weight
        self.add_construction(add_cost=False)

    def _cost(self):
        C = self.baseline_purchase_costs
        C['MicrocontrollerBox'] = (
            self.control_box_enclosure_panel +
            self.control_box_arduino_microcontroller +
            self.control_box_microcontroller_din_rails +
            self.control_box_wifi +
            self.control_box_solid_state_relays +
            self.control_box_other_relays +
            self.control_box_conductors +
            self.control_box_din_rails
            )
        C['SensorsBox'] = (
            self.control_sensors_enclosure_panel +
            self.control_sensors_signal_conditioners +
            self.control_sensors_din_rails +
            self.control_sensors_voltage_converter
            )
        C['ORPBox'] = (
            self.control_orp_enclosure_panel +
            self.control_orp_voltage_converter +
            self.control_orp_controller +
            self.control_orp_probe +
            self.control_orp_cables
            )
        C['PumpBox'] = (
            self.control_pump_enclosure_panel +
            self.control_pump_voltage_converter +
            self.control_pump_voltage_rectifier +
            self.control_pump_box_relay
            )
        C['Sensors'] = (
            self.control_sensors_level_sensor +
            self.control_sensors_level_switch +
            self.control_sensors_flow_switches +
            self.control_sensors_proximity_sensor
            )
        C['AdditionalSupplies'] = (
            self.control_stranded_cables +
            self.control_sensor_cables +
            self.control_cable_glands +
            self.control_gauge_ferrules +
            self.control_zip_ties +
            self.control_terminal_block_jumper
            )
        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio

        self.add_OPEX =  self._calc_replacement_cost() + self._calc_maintenance_labor_cost()

        if self.if_gridtied:
            power_demand = (
                self.control_system_energy_percycle * self.control_batch_cycle_perday +
                self.control_system_ORP_energy_percycle * self.control_batch_cycle_perday +
                self.control_background_runtime_energy_perday
                )
            power_demand = power_demand / 24  # convert from kWh/d to kW
        else:
            power_demand = 0

        self.power_utility(power_demand) # kW

    def _calc_replacement_cost(self):
        control_system_replacement_cost = (
            self.control_box_arduino_microcontroller / self.control_system_PLC_lifetime +
            self.control_orp_controller / self.control_system_PLC_lifetime +
            self.control_sensors_voltage_converter / self.control_voltage_converter_lifetime +
            self.control_orp_voltage_converter / self.control_voltage_converter_lifetime +
            self.control_pump_voltage_converter / self.control_voltage_converter_lifetime +
            self.control_sensors_level_sensor / self.control_level_sensor_lifetime
            )
        return control_system_replacement_cost / (365 * 24) * self.price_ratio # USD/hr

    def _calc_maintenance_labor_cost(self):
        control_system_maintenance_labor = self.control_labor_replacement_misc_repairs * self.wages
        return control_system_maintenance_labor / (365 * 24) # USD/hr


# %%

foundation_path = ospath.join(ng_data_path, '_ng_foundation.csv')

@price_ratio()
class NEWgeneratorFoundation(Copier):
    '''
    Concrete slab foundation for the NEWgenerator.

    This is a non-reactive unit (i.e., the effluent is copied from the influent).

    Parameters
    ----------
    if_larger_NCS : bool
        If baseline design zeolite capacity is increased 3 fold and no zeolite
        is regenerated, zeolite replacement only. Accounts for change in cost
        and emissions from larger concrete foundation to support system.
    user_scale_up : float
        Scaling up of consumables, electricity demand, capital parts,
        and replacement parts based on number of input users compared to the
        baseline design number of users.

    References
    ----------
    [1] Shyu et al., The NEWgeneratorTM Non-Sewered Sanitation System:
    Long-Term Field Testing at an Informal Settlement Community in
    EThekwini Municipality, South Africa.
    Journal of Environmental Management 2021, 296, 112921.
    https://doi.org/10.1016/j.jenvman.2021.112921.

    See Also
    --------
    :class:`~.sanunits.Copier`
    '''

    def __init__(self, ID='', ins=None, outs=(),  thermo=None, init_with='WasteStream',
                 if_larger_NCS=False, user_scale_up=1, **kwargs):
        Copier.__init__(self, ID, ins, outs, thermo, init_with)
        self.if_larger_NCS = if_larger_NCS
        self.user_scale_up = user_scale_up

        data = load_data(path=foundation_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)

    def _init_lca(self):
        self.construction = [
            Construction(item='Concrete', linked_unit=self, quantity_unit='m3'),
            ]

    def _design(self):
        if self.if_larger_NCS:
            q = self.foundation_concrete_thickness*self.foundation_largerncs_concrete_area
        else:
            q = self.foundation_concrete_thickness*self.foundation_concrete_area
        self.design_results['Concrete'] = self.construction[0].quantity = q
        self.add_construction(add_cost=False)

    def _cost(self):
        C = self.baseline_purchase_costs
        C['Concrete'] = self.design_results['Concrete'] * self.foundation_concrete_unit_cost
        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio

        self.add_OPEX = self._calc_replacement_cost() + self._calc_maintenance_labor_cost()

    def _calc_replacement_cost(self):
        return 0

    def _calc_maintenance_labor_cost(self):
        return 0


# %%

grid_path = ospath.join(ng_data_path, '_ng_gridtied.tsv')

@price_ratio()
class NEWgeneratorGridtied(Copier):
    '''
    The grid-tied configuration of the NEWgenerator,
    which is composed of AC conversion system, DC distribution box,
    and miscellaneous parts.

    This is a non-reactive unit (i.e., the effluent is copied from the influent).

    Parameters
    ----------
    user_scale_up : float
        Scaling up of consumables, electricity demand, capital parts,
        and replacement parts based on number of input users compared to the
        baseline design number of users.

    References
    ----------
    [1] Shyu et al., The NEWgeneratorTM Non-Sewered Sanitation System:
    Long-Term Field Testing at an Informal Settlement Community in
    EThekwini Municipality, South Africa.
    Journal of Environmental Management 2021, 296, 112921.
    https://doi.org/10.1016/j.jenvman.2021.112921.

    See Also
    --------
    :class:`~.sanunits.Copier`
    '''

    def __init__(self, ID='', ins=None, outs=(),  thermo=None, init_with='WasteStream',
                 user_scale_up=1, **kwargs):
        Copier.__init__(self, ID, ins, outs, thermo, init_with)
        self.user_scale_up = user_scale_up

        data = load_data(path=grid_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)
            
    def _init_lca(self):
        self.construction = [
            Construction(item='Polycarbonate', linked_unit=self, quantity_unit='kg'),
            Construction(item='Aluminum', linked_unit=self, quantity_unit='kg'),
            Construction(item='ElectronicsPassive', linked_unit=self, quantity_unit='kg'),
            Construction(item='ElectronicsActive', linked_unit=self, quantity_unit='kg'),
            ]

    def _design(self):
        design = self.design_results
        constr = self.construction
        design['Polycarbonate'] = constr[0].quantity = self.grid_polycarbonate_weight
        design['Aluminum'] = constr[1].quantity = self.grid_aluminium_weight
        design['ElectronicsPassive'] = constr[2].quantity = self.grid_electronics_passive_weight*(self.user_scale_up*0.6)
        design['ElectronicsActive'] = constr[3].quantity = self.grid_electronics_active_weight*(self.user_scale_up*0.6)
        self.add_construction(add_cost=False)

    def _cost(self):
        C = self.baseline_purchase_costs
        C['ACDCBlock'] = (self.grid_power_ac_conversion)*(self.user_scale_up*0.6)
        C['DCdistribution'] = (self.grid_power_dc_distribution)*(self.user_scale_up*0.6)
        C['Misc'] = self.grid_power_misc
        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio

        self.add_OPEX = self._calc_replacement_cost() + self._calc_maintenance_labor_cost()

    def _calc_replacement_cost(self):
        return 0

    def _calc_maintenance_labor_cost(self):
        gridtied_maintenance_labor = self.grid_labor_replacement_misc_repairs * self.wages
        return gridtied_maintenance_labor/ (365 * 24)  # USD/hr


# %%

housing_path = ospath.join(ng_data_path, '_ng_housing.csv')

@price_ratio()
class NEWgeneratorHousing(Copier):
    '''
    Framework and fittings for the NEWgenerator.

    This is a non-reactive unit (i.e., the effluent is copied from the influent).

    Parameters
    ----------
    if_larger_NCS_or_greater_150users : bool
        If baseline design zeolite capacity is increased 3 fold and no zeolite
        is regenerated, zeolite replacement only. Accounts for change in cost
        and emissions from larger housing framework (lower cost).
        Or users is greater or equal to 150 users, accounting for larger
        housing framework (lower cost).
    if_cheaper_housing : bool
        If the housing framework cost is reduced by 75% of baseline design
        BOM cost of housing framework.
    user_scale_up : float
        Scaling up of consumables, electricity demand, capital parts,
        and replacement parts based on number of input users compared to the
        baseline design number of users.

    References
    ----------
    [1] Shyu et al., The NEWgeneratorTM Non-Sewered Sanitation System:
    Long-Term Field Testing at an Informal Settlement Community in
    EThekwini Municipality, South Africa.
    Journal of Environmental Management 2021, 296, 112921.
    https://doi.org/10.1016/j.jenvman.2021.112921.

    See Also
    --------
    :class:`~.sanunits.Copier`
    '''

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 if_larger_NCS_or_greater_150users=False,
                 if_cheaper_housing=False, user_scale_up=1, **kwargs):
        Copier.__init__(self, ID, ins, outs, thermo, init_with)
        self.if_larger_NCS_or_greater_150users = if_larger_NCS_or_greater_150users
        self.if_cheaper_housing = if_cheaper_housing
        self.user_scale_up = user_scale_up

        data = load_data(path=housing_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)

        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)
            
    def _init_lca(self):
        self.construction = [
            Construction(item='Steel', linked_unit=self, quantity_unit='kg'),
            Construction(item='ZincCoat', linked_unit=self, quantity_unit='m2'),
            Construction(item='Aluminum', linked_unit=self, quantity_unit='kg'),
            Construction(item='Plastic', linked_unit=self, quantity_unit='kg'),
            ]

    def _design(self):
        design = self.design_results
        constr = self.construction
        design['Steel'] = constr[0].quantity = (self.housing_steel_weight + self.housing_zinc_coated_steel_weight)
        design['ZincCoat'] = constr[1].quantity = self.housing_zinc_coat
        design['Aluminum'] = constr[2].quantity = self.housing_aluminum_weight
        design['Plastic'] = constr[3].quantity = self.housing_PE_weight
        self.add_construction(add_cost=False)

    def _cost(self):
        C = self.baseline_purchase_costs
        if self.if_cheaper_housing:
            C['Framework'] = self.cheaper_housing_framework
        elif self.if_larger_NCS_or_greater_150users:
            C['Framework'] = self.housing_largerncs_framework
        else:
            C['Framework'] = self.housing_framework

        C['Fittings'] = (self.housing_fittings)
        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio

        self.add_OPEX =  self._calc_replacement_cost() + self._calc_maintenance_labor_cost()


    def _calc_replacement_cost(self):
        return 0


    def _calc_maintenance_labor_cost(self):
        return 0


# %%

ix_path = ospath.join(ng_data_path, '_ng_ion_exchange.tsv')

@price_ratio()
class NEWgeneratorIonExchange(SanUnit):
    '''
    Ion exchange in the NEWgenerator for N recovery from liquid stream.
    Concentrated NH3 is recovered.

    Parameters
    ----------
    ins : Iterable(stream)
        waste: waste for treatment.
        Zeolite_in: used for NH3 ion exchange adsorption, residual COD and color [kg/h].
        GAC_out: used for additional COD removal [kg/h].
        NaCl: regenerant solution for zeolite beds [kg/h].
        NaOH: elevate regenerant solution to pH 10-11 [kg/h].
    outs : Iterable(stream)
        treated: treated waste.
        Zeolite_out: Used zeolite as waste [kg/h].
        GAC_out: Used GAC as waste [kg/h].
        concentrated NH3: Concentration of NH3 in regenerant solution and desorped from zeolite [kg/h].
    if_gridtied : bool
        If electricity source is grid-tied energy.
    if_larger_NCS : bool
        If baseline design zeolite capacity is increased 3 fold and no zeolite
        is regenerated, zeolite replacement only. Accounts for change in cost
        and emissions from housing, tanks, consumable zeolite, labor covered.
    if_zeolite_replacement_only : bool
        If the baseline design zeolite is replaced only, and the zeolite is NOT
        regenerated. Accounts for change in cost and emissions from
        consumable zeolite and labor covered.
    user_scale_up : float
        Scaling up of consumables, electricity demand, capital parts,
        and replacement parts based on number of input users compared to the
        baseline design number of users

    References
    ----------
    [1] Shyu et al., The NEWgeneratorTM Non-Sewered Sanitation System:
    Long-Term Field Testing at an Informal Settlement Community in
    EThekwini Municipality, South Africa.
    Journal of Environmental Management 2021, 296, 112921.
    https://doi.org/10.1016/j.jenvman.2021.112921.
    [2] Castro et al., Performance and Onsite Regeneration of Natural Zeolite
    for Ammonium Removal in a Field-Scale Non-Sewered Sanitation System.
    Science of The Total Environment 2021, 776, 145938.
    https://doi.org/10.1016/j.scitotenv.2021.145938.
    [3] Lohman et al., Advancing Sustainable Sanitation and Agriculture
    through Investments in Human-Derived Nutrient Systems.
    Environ. Sci. Technol. 2020, 54, (15), 9217-9227.
    https://dx.doi.org/10.1021/acs.est.0c03764
    [4] Tarpeh et al., Evaluating ion exchange for nitrogen recovery from
    source-separated urine in Nairobi, Kenya. Development Engineering. 2018,
    3, 188–195.
    https://doi.org/10.1016/j.deveng.2018.07.002
    '''
    _N_ins = 5
    _N_outs = 4

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 if_gridtied=True, if_larger_NCS=False,
                 if_zeolite_replacement_only=False, user_scale_up=1, **kwargs):

        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with,
                         F_BM_default=1)
        self.if_gridtied = if_gridtied
        self.if_larger_NCS = if_larger_NCS
        self.if_zeolite_replacement_only = if_zeolite_replacement_only
        self.user_scale_up = user_scale_up

        data = load_data(path=ix_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)
    
    def _init_lca(self):
        self.construction = [
            Construction(item='Plastic', linked_unit=self, quantity_unit='kg'),
            Construction(item='PE', linked_unit=self, quantity_unit='kg'),
            Construction(item='PVC', linked_unit=self, quantity_unit='kg'),
            Construction(item='Polycarbonate', linked_unit=self, quantity_unit='kg'),
            Construction(item='Steel', linked_unit=self, quantity_unit='kg'),
            ]

    def _run(self):
        waste, zeolite_in, gac_in, NaCl, NaOH = self.ins
        treated, zeolite_out, gac_out, conc_NH3 = self.outs
        treated.copy_like(waste)
        for stream in (zeolite_in, zeolite_out, gac_in, gac_out, NaCl): stream.phase = 's'

        user_scale_up = self.user_scale_up
        hr_per_yr = 365*24
        ion_gac_lifetime_hr = self.ion_gac_lifetime * hr_per_yr
        NaCl_demand_time = NaOH_demand_time = 0
        if self.if_larger_NCS:
            zeolite_demand_time = self.ion_largerncs_scenario_zeolite_replacement_mass*user_scale_up/hr_per_yr
            gac_demand_time = self.ion_gac_weight*user_scale_up/ion_gac_lifetime_hr
        elif self.if_zeolite_replacement_only:
            zeolite_demand_time = self.ion_replacement_scenario_zeolite_massflow*user_scale_up/hr_per_yr
            gac_demand_time = self.ion_gac_weight*user_scale_up/ion_gac_lifetime_hr
        else:
            zeolite_demand_time = self.ion_zeolite_weight*user_scale_up / (self.ion_zeolite_lifetime*hr_per_yr)
            NaCl_demand_time = self.ion_NaCl_weight*user_scale_up * (self.ion_regen_freq_per_yr/hr_per_yr)
            NaOH_demand_time = self.ion_NaOH_weight*user_scale_up * (self.ion_regen_freq_per_yr/hr_per_yr)
            gac_demand_time = self.ion_gac_weight*user_scale_up / (self.ion_gac_lifetime*hr_per_yr)

        zeolite_in.imass['Zeolite'] = zeolite_out.imass['Zeolite'] = zeolite_demand_time
        NaCl.imass['SodiumChloride'] = NaCl_demand_time
        NaOH.imass['SodiumHydroxide'] = NaOH_demand_time
        gac_in.imass['GAC'] = gac_out.imass['GAC'] = gac_demand_time

        N_removed = waste.imass['NH3'] * self.ion_TN_removal
        N_recovered = N_removed * self.ion_desorption_recovery_efficiency
        treated.imass['NH3'] = waste.imass['NH3'] - N_removed
        conc_NH3.imass['NH3'] = N_recovered # kg N / hr

        # Not sure why NaCl and NaOH was added to the concentrated NH3 stream - Hannah Lohman 6/6/2022
        # conc_NH3.imass['SodiumChloride'] = self.NaCl_demand_time
        # conc_NH3.imass['SodiumHydroxide'] = self.NaOH_demand_time
        zeolite_out.imass['NH3'] = N_removed - N_recovered

    def _design(self):
        design = self.design_results
        constr = self.construction
        design['Plastic'] = constr[0].quantity = self.ion_plastic_weight
        design['PE'] = constr[1].quantity = self.ion_PE_weight
        design['PVC'] = constr[2].quantity = (self.ion_PVC_weight + self.ion_piping_PVC_weight)
        design['Polycarbonate'] = constr[3].quantity = self.ion_polycarbonate_weight
        design['Steel'] = constr[4].quantity = self.ion_steel_weight
        self.add_construction(add_cost=False)

    def _cost(self):
        C = self.baseline_purchase_costs

        if self.if_larger_NCS:
            C['EqualizationTank2'] = (
                self.ion_eq2_tank +
                self.ion_eq2_pump*(self.user_scale_up ** 0.6) +
                self.ion_eq2_powered_valve +
                self.ion_eq2_fittings
                )
            C['ZeoliteTanks'] = self.ion_largerncs_scenario_zeolite_tanks
            C['Fittings'] = self.ion_zeolite_fittings_largerncs
            C['ZeoliteGACTreatment'] = (
                self.ion_filtration_media +
                self.ion_gac +
                self.ion_zeolite
                )
            C['EqualizationTank3'] = (
                self.ion_eq3_tank +
                self.ion_eq3_pump*(self.user_scale_up*0.6) +
                self.ion_eq3_filter +
                self.ion_eq3_fittings
                )
            C['Piping'] = self.ion_piping
        else:
            C['EqualizationTank2'] = (
                self.ion_eq2_tank +
                self.ion_eq2_pump*(self.user_scale_up*0.6) +
                self.ion_eq2_powered_valve +
                self.ion_eq2_fittings
                )
            C['ZeoliteTanks'] = self.ion_zeolite_tanks
            C['Fittings'] = self.ion_zeolite_fittings
            C['ZeoliteGACTreatment'] = (
                self.ion_filtration_media +
                self.ion_gac +
                self.ion_zeolite
                )
            C['EqualizationTank3'] = (
                self.ion_eq3_tank +
                self.ion_eq3_pump * (self.user_scale_up*0.6) +
                self.ion_eq3_filter + self.ion_eq3_fittings
                )
            C['Piping'] = self.ion_piping

        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio

        self.add_OPEX = self._calc_replacement_cost() + self._calc_maintenance_labor_cost()

        if self.if_gridtied:
            pump_power = self.ion_backwash_pump_energy_percycle*self.ion_batch_cycle_perday
            power_demand = 2 * pump_power * self.user_scale_up # 2 pumps
            power_demand = power_demand / 24  # convert from kWh/d to kW
        else:
            power_demand = 0

        self.power_utility(power_demand) # kW


    def _calc_replacement_cost(self):
        ion_exchange_replacement_cost = (
            self.ion_eq2_pump*(self.user_scale_up*0.6)/self.ion_other_pump_lifetime +
            self.ion_eq2_powered_valve/self.ion_solenoid_valve_lifetime +
            self.ion_eq3_pump*(self.user_scale_up*0.6)/self.ion_other_pump_lifetime
            )
        return ion_exchange_replacement_cost/ (365 * 24) * self.price_ratio # USD/hr


    def _calc_maintenance_labor_cost(self):
        labor = (
            self.ion_labor_maintenance_GAC_replacement * self.user_scale_up +
            self.ion_labor_replacement_misc_repairs +
            self.ion_labor_replacement_other_pump * 2 * self.user_scale_up  # 2 pumps (i.e., eq2 and eq3)
            )
        if self.if_larger_NCS: labor += (self.ion_largerncs_scenario_replacement * self.user_scale_up)
        elif self.if_zeolite_replacement_only: # replacement only and no regeneration
            labor += (self.ion_replacement_scenario_replacement * self.user_scale_up)
        else: # replacement and regeneration
            labor += ((self.ion_replacement_scenario_replacement+self.ion_labor_maintenance_zeolite_regeneration) * self.user_scale_up)

        return labor*self.wages/(365 * 24) # USD/hr


# %%

pv_path = ospath.join(ng_data_path, '_ng_photovoltaic.tsv')

@price_ratio()
class NEWgeneratorPhotovoltaic(Copier):
    '''
    The photovoltaic configuration of the NEWgenerator,
    which is composed of solar panels, batteries, DC distribution box,
    and other parts.

    This is a non-reactive unit (i.e., the effluent is copied from the influent).

    Parameters
    ----------
    if_lithium_battery : bool
        If alternative lithium battery is used.
    user_scale_up : float
        Scaling up of consumables, electricity demand, capital parts,
        and replacement parts based on number of input users compared to the
        baseline design number of users

    References
    ----------
    [1] Shyu et al., The NEWgeneratorTM Non-Sewered Sanitation System:
    Long-Term Field Testing at an Informal Settlement Community in
    EThekwini Municipality, South Africa.
    Journal of Environmental Management 2021, 296, 112921.
    https://doi.org/10.1016/j.jenvman.2021.112921.

    See Also
    --------
    :class:`~.sanunits.Copier`
    '''

    def __init__(self, ID='', ins=None, outs=(),  thermo=None, init_with='WasteStream',
                 if_lithium_battery=False, user_scale_up=1, **kwargs):
        Copier.__init__(self, ID, ins, outs, thermo, init_with)
        self.if_lithium_battery = if_lithium_battery
        self.user_scale_up = user_scale_up

        data = load_data(path=pv_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        self._refres_lca()

        for attr, value in kwargs.items():
            setattr(self, attr, value)


    def _refres_lca(self):
        if self.if_lithium_battery:
            battery_construction = Construction(
                item='LithiumBattery', linked_unit=self, quantity_unit='kg',
                lifetime=self.pv_lithium_photovoltaic_battery_lifetime)
        else:
            battery_construction = Construction(
                item='Battery', linked_unit=self, quantity_unit='kg',
                lifetime=self.pv_battery_lifetime)
        
        self.construction = [
            Construction(item='PhotovoltaicPanel', linked_unit=self, quantity_unit='m2'),
            battery_construction,
            Construction(item='ElectricConnectors', linked_unit=self, quantity_unit='kg'),
            Construction(item='Switches', linked_unit=self, quantity_unit='kg'),
            Construction(item='Polycarbonate', linked_unit=self, quantity_unit='kg'),
            Construction(item='ElectronicsPassive', linked_unit=self, quantity_unit='kg'),
            Construction(item='ElectricCables', linked_unit=self, quantity_unit='m'),
            Construction(item='ControlUnits', linked_unit=self, quantity_unit='kg'),
            Construction(item='Aluminum', linked_unit=self, quantity_unit='kg'),
            ]

    def _design(self):
        design = self.design_results
        constr = self.construction
        user_scale_up = self.user_scale_up
        design['PhotovoltaicPanel'] = constr[0].quantity = self.pv_photovoltaic_panel_area * (user_scale_up**0.6)
        design['Battery'] = constr[1].quantity = self.pv_battery_quant * (user_scale_up**0.6)
        design['ElectricConnectors'] = constr[2].quantity = self.pv_connectors_weight
        design['Switches'] = constr[3].quantity = self.pv_switch_weight
        design['Polycarbonate'] = constr[4].quantity = self.pv_polycarbonate_weight
        design['ElectronicsPassive'] = constr[5].quantity = self.pv_electronics_passive_weight
        design['ElectricCables'] = constr[6].quantity = self.pv_cable_length
        design['ControlUnits'] = constr[7].quantity = self.pv_control_units_weight
        design['Aluminum'] = constr[8].quantity = self.pv_aluminum_weight
        self.add_construction(add_cost=False)

    def _cost(self):
        C = self.baseline_purchase_costs
        C['Battery'] = (
            self.pv_battery_connectors +
            self.pv_battery_switch +
            self.pv_battery_terminal
            )
        if self.if_lithium_battery: C['Battery'] += self.pv_lithium_battery * (self.user_scale_up**0.6)
        else: C['Battery'] += self.pv_battery * (self.user_scale_up**0.6)

        C['PhotovoltaicPanel'] = (
            self.pv_panel * (self.user_scale_up**0.6) +
            self.pv_racking_system +
            self.pv_charge_controllers +
            self.pv_connectors +
            self.pv_combiner_box +
            self.pv_circuit_breaker
            )
        C['DCdistribution'] = (
            self.pv_dc_enclosure_panel +
            self.pv_dc_din_rail +
            self.pv_dc_din_circuit_breaker +
            self.pv_dc_distribution_block +
            self.pv_dc_voltage_converter
            )
        C['Misc'] = self.pv_misc
        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio

        self.add_OPEX = self._calc_replacement_cost() + self._calc_maintenance_labor_cost()

    def _calc_replacement_cost(self):
        if self.if_lithium_battery:
            battery_replacememt_parts_annual_cost = self.pv_lithium_battery * (self.user_scale_up**0.6) / self.pv_lithium_photovoltaic_battery_lifetime
        else:
            battery_replacememt_parts_annual_cost = self.pv_battery * (self.user_scale_up**0.6) / self.pv_battery_lifetime

        photovoltaic_replacement_cost = (
            battery_replacememt_parts_annual_cost +
            self.pv_dc_voltage_converter / self.pv_voltage_converter_lifetime
            )
        return photovoltaic_replacement_cost/ (365 * 24) * self.price_ratio # USD/hr

    def _calc_maintenance_labor_cost(self):
        photovoltaic_maintenance_labor = (
            self.pv_labor_replacement_batteries * self.user_scale_up +
            self.pv_labor_replacement_misc_repairs
            ) * self.wages
        return photovoltaic_maintenance_labor/ (365 * 24) # USD/hr


    @property
    def if_lithium_battery(self):
        return self._if_lithium_battery
    @if_lithium_battery.setter
    def if_lithium_battery(self, i):
        self._if_lithium_battery = i
        constr = self.construction
        if constr:
            for i in constr: i.registry.discard(i)
            self._init_lca()
            

# %%

pretreatment_path = ospath.join(ng_data_path, '_ng_pretreatment.csv')

@price_ratio()
class NEWgeneratorPretreatment(Copier):
    '''
    Pretreatment unit for the NEWgenerator,
    which is composed of a concrete equalization tank and a feed system pump.

    This is a non-reactive unit (i.e., the effluent is copied from the influent).

    Parameters
    ----------
    if_gridtied : bool
        If electricity source is grid-tied energy.
    user_scale_up : float
        Scaling up of consumables, electricity demand, capital parts,
        and replacement parts based on number of input users compared to the
        baseline design number of users.

    References
    ----------
    [1] Shyu et al., The NEWgeneratorTM Non-Sewered Sanitation System:
    Long-Term Field Testing at an Informal Settlement Community in
    EThekwini Municipality, South Africa.
    Journal of Environmental Management 2021, 296, 112921.
    https://doi.org/10.1016/j.jenvman.2021.112921.

    See Also
    --------
    :class:`~.sanunits.Copier`
    '''

    def __init__(self, ID='', ins=None, outs=(),  thermo=None, init_with='WasteStream',
                 if_gridtied=True, user_scale_up=1, **kwargs):
        Copier.__init__(self, ID, ins, outs, thermo, init_with)
        self.if_gridtied = if_gridtied
        self.user_scale_up = user_scale_up

        data = load_data(path=pretreatment_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)

    def _init_lca(self):
        self.construction = [
            Construction(item='Plastic', linked_unit=self, quantity_unit='kg'),
            Construction(item='Polycarbonate', linked_unit=self, quantity_unit='kg'),
            Construction(item='Pump', linked_unit=self, quantity=1., quantity_unit='ea'),
            Construction(item='Steel', linked_unit=self, quantity_unit='kg')
            ]


    def _design(self):
        design = self.design_results
        constr = self.construction
        design['Plastic'] = constr[0].quantity = self.pretreatment_tank_PE_weight
        design['Polycarbonate'] = constr[1].quantity = self.pretreatment_piping_PVC_weight
        design['Steel'] = constr[-1].quantity = self.pretreatment_barscreen_weight
        self.add_construction(add_cost=False)

    def _cost(self):
        C = self.baseline_purchase_costs
        user_scale_up = self.user_scale_up
        C['Tank'] = self.pretreatment_eq_tank
        C['Piping'] = self.pretreatment_piping
        C['Barscreen'] = self.pretreatment_barscreen
        C['Pump'] = self.pretreatment_feed_pump * (user_scale_up ** 0.6)
        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio

        self.add_OPEX = self._calc_replacement_cost() + self._calc_maintenance_labor_cost()

        if self.if_gridtied:
            power_demand = self.pretreatment_feed_pump_energy_percycle*self.pretreatment_batch_cycle_perday*self.user_scale_up
            power_demand = power_demand / 24  # convert from kWh/d to kW
        else:
            power_demand = 0

        self.power_utility(power_demand)  # kW

    def _calc_replacement_cost(self):
        pretreatment_replacement_cost = self.pretreatment_feed_pump*(self.user_scale_up ** 0.6) / self.pretreatment_other_pump_lifetime
        return pretreatment_replacement_cost/(365 * 24)* self.price_ratio # USD/hr

    def _calc_maintenance_labor_cost(self):
        pretreatment_maintenance_labor_cost = (
            self.pretreatment_labor_replacement_other_pump * self.user_scale_up +
            self.pretreatment_labor_screen_pit_cleaning
            ) * self.wages
        return pretreatment_maintenance_labor_cost/(365 * 24) # USD/hr