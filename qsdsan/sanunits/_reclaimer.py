#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:

    Tori Morgan <vlmorgan@illinois.edu>

    Yalin Li <mailto.yalin.li@gmail.com>

    Hannah Lohman <hlohman94@gmail.com>

    Lewis Rowles <stetsonsc@gmail.com>

This module contains unit operations used in the Reclaimer system
developed by Duke University as described in
https://washaid.pratt.duke.edu/work/water-sanitation/reinvent-toilet-challenge
and `Trotochaud et al. <https://doi.org/10.1021/acs.est.0c02755>`_

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

from math import ceil
from qsdsan import SanUnit, Construction
from ..utils import ospath, load_data, data_path, price_ratio

__all__ = (
    'ReclaimerECR',
    'ReclaimerHousing',
    'ReclaimerIonExchange',
    'ReclaimerSolar',
    'ReclaimerSystem',
    'ReclaimerUltrafiltration',
    )

re_su_data_path = ospath.join(data_path, 'sanunit_data/re')


# %%

electrochemical_path = ospath.join(re_su_data_path, '_re_ecr.csv')

class ReclaimerECR(SanUnit):
    '''
    Electrochemical reactor (ECR) in the Reclaimer system with chlorine dosing.

    This is a non-reactive unit (i.e., the effluent is copied from the influent).

    The following impact items should be pre-constructed for life cycle assessment:
    Titanium.

    Parameters
    ----------
    ppl: int
        Total number of users for scaling of costs.
    if_gridtied: bool
        If using grid electricity instead of photovoltaic electricity.

    References
    ----------
    [1] Trotochaud et al., Laboratory Demonstration and Preliminary Techno-Economic Analysis of an Onsite
    Wastewater Treatment System Environ. Sci. Technol. 2020, 54, (24), 16147–16155.
    https://dx.doi.org/10.1021/acs.est.0c02755

    [2] Duke Center for WaSH-AID Reclaimer design team data and guidance
    https://washaid.pratt.duke.edu/work/water-sanitation/reinvent-toilet-challenge
    '''

    # Constants
    baseline_ppl = 30  # baseline population served by Reclaimer
    exponent_scale = 0.6  # exponential scaling constant

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 if_gridtied=True, ppl=1, **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with, F_BM_default=1)
        self.ppl = ppl
        self.if_gridtied = if_gridtied

        data = load_data(path=electrochemical_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)

    _N_ins = 1
    _N_outs = 1


    def _design(self):
        design = self.design_results
        design['Titanium'] = electrode_quant = self.Titanium_weight * (self.ppl / self.baseline_ppl)  # linear scale
        self.construction = Construction(item='Titanium', quantity=electrode_quant, quantity_unit='kg')
        self.add_construction(add_cost=False)

    def _cost(self):
        C = self.baseline_purchase_costs
        C['EC_brush'] = self.EC_brush
        C['EC_cell'] = self.EC_cell

        # Exponentially scale capital cost with number of users
        scale = (self.ppl / self.baseline_ppl) ** self.exponent_scale
        for equipment, cost in C.items():
            C[equipment] = cost * scale

        self.add_OPEX = self._calc_replacement_cost()

        if self.if_gridtied:  # scale linearly with number of units
            power_demand = (self.power_demand_ecr / 1000) * self.N_reclaimers  # [W/day][1 kW/1000 W] = [kW/d]
        else:
            power_demand = 0
        self.power_utility(power_demand)  # kW

    def _calc_replacement_cost(self):
        scale = (self.ppl / self.baseline_ppl) ** self.exponent_scale
        ec_brush_replacement_cost = scale * (self.EC_brush / self.EC_brush_lifetime)
        ec_cell_replacement_cost = scale * (self.EC_cell / self.EC_cell_lifetime)
        ec_replacement_cost = ec_brush_replacement_cost + ec_cell_replacement_cost
        ec_replacement_cost = ec_replacement_cost / (365 * 24)  # convert from USD/year to USD/hour
        return ec_replacement_cost

    @property
    def N_reclaimers(self):
        '''[int] Number of the reclaimer units needed, calculated by `ppl`/`baseline_ppl`.'''
        return ceil(self.ppl / self.baseline_ppl)


# %%

housing_path = ospath.join(re_su_data_path,'_re_housing.csv')

@price_ratio()
class ReclaimerHousing(SanUnit):
    '''
    Structural housing for the Reclaimer system.

    This is a non-reactive unit (i.e., the effluent is copied from the influent).

    The following impact items should be pre-constructed for life cycle assessment:
    Steel.

    Parameters
    ----------
    ppl: int
        Total number of users for scaling of costs.

    References
    ----------
    [1] Trotochaud et al., Laboratory Demonstration and Preliminary Techno-Economic Analysis of an Onsite
    Wastewater Treatment System Environ. Sci. Technol. 2020, 54, (24), 16147–16155.
    https://dx.doi.org/10.1021/acs.est.0c02755

    [2] Duke Center for WaSH-AID Reclaimer design team data and guidance
    https://washaid.pratt.duke.edu/work/water-sanitation/reinvent-toilet-challenge
    
    [3] Eco-san water recycling toilet Reinvented Toilet design team bill of materials
    https://sanitation.ansi.org/EcoSanToilet
    '''

    baseline_ppl = 30  # baseline population served by Reclaimer
    ppl_per_MURT = 25 # assume 25 people per MURT

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream', ppl=1, **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with, F_BM_default=1)
        self.ppl = ppl

        data = load_data(path=housing_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)


    def _design(self):
        design = self.design_results
        design['Steel'] = steel_quant = (
            self.steel_weight +
            self.framework_weight/4 +
            self.fittings_weight
            ) * (self.ppl / self.baseline_ppl)  # linear scale
        design['Metal'] = metal_quant = self.aluminum_weight * (self.ppl / self.baseline_ppl)  # linear scale
        self.construction = (
            Construction(item='Steel', quantity=steel_quant, quantity_unit='kg'),
            Construction(item='Metal', quantity=metal_quant, quantity_unit='kg')
            )
        self.add_construction(add_cost=False)

    def _cost(self):
        C = self.baseline_purchase_costs
        C['Housing'] = (
            self.frame +
            self.extrusion +
            self.angle_frame +
            self.angle +
            self.door_sheet +
            self.plate_valve +
            self.powder +
            self.container
            ) * (1 + 0.1 * (self.N_reclaimers-1))

        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio


    @property
    def N_reclaimers(self):
        '''[int] Number of the reclaimer units needed, calculated by `ppl`/`baseline_ppl`.'''
        return ceil(self.ppl / self.baseline_ppl)

    @property
    def N_toilets(self):
        '''[int] Number of the MURT units, calculated by `ppl`/`ppl_per_MURT`.'''
        return ceil(self.ppl / self.ppl_per_MURT)


# %%

ion_exchange_path = ospath.join(re_su_data_path, '_re_ion_exchange.csv')

class ReclaimerIonExchange(SanUnit):
    '''
    Ion exchange in the Reclaimer system is used for N recovery from liquid stream.
    Concentrated NH3 is recovered.

    The following impact items should be pre-constructed for life cycle assessment:
    Plastic, PVC, Steel.

    Parameters
    ----------
    ins : Iterable(stream)
        waste: liquid waste stream to be treated by ion exchange unit.
        zeolite_in: zeolite input.
        gac_in: GAC input.
        KCl: KCl input.
    outs : Iterable(stream)
        treated: treated liquid leaving ion exchange unit.
        zeolite_out: spent zeolite.
        gac_out: spent GAC.
        conc_NH3: concentrated NH3.
    ppl: int
        Total number of users for scaling of costs.

    References
    ----------
    [1] Lohman et al., Advancing Sustainable Sanitation and Agriculture
    through Investments in Human-Derived Nutrient Systems.
    Environ. Sci. Technol. 2020, 54, (15), 9217-9227.
    https://dx.doi.org/10.1021/acs.est.0c03764

    [2] Tarpeh et al., Evaluating ion exchange for nitrogen recovery from
    source-separated urine in Nairobi, Kenya. Development Engineering. 2018,
    3, 188–195.
    https://doi.org/10.1016/j.deveng.2018.07.002
    
    [3] Trotochaud et al., Laboratory Demonstration and Preliminary Techno-Economic Analysis of an Onsite
    Wastewater Treatment System Environ. Sci. Technol. 2020, 54, (24), 16147–16155.
    https://dx.doi.org/10.1021/acs.est.0c02755
    
    [4] Duke Center for WaSH-AID Reclaimer design team data and guidance
    https://washaid.pratt.duke.edu/work/water-sanitation/reinvent-toilet-challenge
    '''

    # Constants
    baseline_ppl = 30  # baseline population served by Reclaimer
    exponent_scale = 0.6  # exponential scaling constant

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream', ppl=1, **kwargs):

        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with, F_BM_default=1)
        self.ppl = ppl

        data = load_data(path=ion_exchange_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)

    _N_ins = 4
    _N_outs = 4

    def _run(self):
        waste, zeolite_in, gac_in, KCl = self.ins
        treated, zeolite_out, gac_out, conc_NH3 = self.outs
        treated.copy_like(self.ins[0])
        for s in (zeolite_in, zeolite_out, gac_in, gac_out, KCl): s.phase = 's'

        N = self.N_reclaimers
        # Zeolite
        zeolite_demand_time = (self.zeolite_weight / (self.zeolite_lifetime * 365 * 24)) * N  # kg zeolite/hr
        zeolite_in.imass['Zeolite'] = zeolite_out.imass['Zeolite'] = zeolite_demand_time

        # GAC
        gac_demand_time = (self.gac_weight / (self.gac_lifetime * 365 * 24)) * N  # kg GAC/hr
        gac_in.imass['GAC'] = gac_out.imass['GAC'] = gac_demand_time

        # KCl
        KCl_demand_time = (self.KCl_weight / (self.KCl_regeneration_freq * 365 * 24)) * N  # kg KCl/hr
        KCl.imass['PotassiumChloride'] = KCl_demand_time

        N_removed = waste.imass['NH3'] * self.TN_removal
        N_recovered = N_removed * self.desorption_recovery_efficiency  # kg N / hr
        treated.imass['NH3'] = waste.imass['NH3'] - N_removed  # kg N / hr
        conc_NH3.imass['NH3'] = N_recovered  # kg N / hr

        # Not sure why KCl was added to the concentrated NH3 stream - Hannah Lohman 6/6/2022
        # conc_NH3.imass['PotassiumChloride'] = self.KCl_demand_time
        zeolite_out.imass['NH3'] = N_removed - N_recovered

    def _design(self):
        design = self.design_results
        factor = self.ppl / self.baseline_ppl # linear scale
        design['Plastic'] = P_quant = self.Plastic_weight * factor
        design['PVC'] = PVC_quant = self.PVC_weight * factor
        design['Steel'] = S_quant = self.Steel_weight * factor

        self.construction = (
            Construction(item='Plastic', quantity=P_quant, quantity_unit='kg'),
            Construction(item='PVC', quantity=PVC_quant, quantity_unit='kg'),
            Construction(item='Steel', quantity=S_quant, quantity_unit='kg'),
            )
        self.add_construction(add_cost=False)

    def _cost(self):
        C = self.baseline_purchase_costs
        C['Pipes'] = self.four_in_pipe_SCH40 + self.four_in_pipe_SCH80
        C['Fittings'] = (
            self.four_in_pipe_SCH80_endcap +
            self.NRV +
            self.connector +
            self.ball_valve +
            self.three_eight_elbow +
            self.ten_ten_mm_tee +
            self.OD_tube +
            self.four_in_pipe_clamp
            )
        C['GAC_Zeolite'] = self.GAC_zeolite_mesh + (self.GAC_cost * self.gac_weight) + (self.Zeolite_cost * self.zeolite_weight)
        C['Regeneration Solution'] = self.KCl_cost * self.KCl_weight

        # Exponentially scale capital cost with number of users
        scale = (self.ppl / self.baseline_ppl) ** self.exponent_scale
        for equipment, cost in C.items():
            C[equipment] = cost * scale

        self.add_OPEX = self._calc_replacement_cost() + self._calc_maintenance_labor_cost()

    def _calc_replacement_cost(self):
        scale = (self.ppl / self.baseline_ppl) ** self.exponent_scale
        zeolite_replacement_cost = self.Zeolite_cost * self.zeolite_weight / self.zeolite_lifetime  # USD/year
        gac_replacement_cost = self.GAC_cost * self.gac_weight / self.gac_lifetime  # USD/year
        kcl_replacement_cost = self.KCl_cost * self.KCl_weight / self.KCl_regeneration_freq  # USD/year

        # Replacement parts for everything else, assume 2-6% of capital as annual cost
        C = self.baseline_purchase_costs
        other_replacement_cost = (C['Pipes']+C['Fittings']+self.GAC_zeolite_mesh)*self.om_capital_ratio
        ion_exchange_replacement_cost = (
            zeolite_replacement_cost +
            gac_replacement_cost +
            kcl_replacement_cost
            + other_replacement_cost) * scale / (365 * 24)  # USD/hr
        return ion_exchange_replacement_cost

    def _calc_maintenance_labor_cost(self):
        scale = (self.ppl / self.baseline_ppl) ** self.exponent_scale
        labor_cost = (self.wages * self.labor_maintenance_zeolite_regeneration) * scale  # USD/year
        return labor_cost / (365 * 24)  # USD/hr


    @property
    def N_reclaimers(self):
        '''[int] Number of the reclaimer units needed, calculated by `ppl`/`baseline_ppl`.'''
        return ceil(self.ppl / self.baseline_ppl)


# %%

solar_path = ospath.join(re_su_data_path, '_re_solar.csv')

@price_ratio()
class ReclaimerSolar(SanUnit):
    '''
    Photovoltaic system for solar power generation in the Reclaimer system.

    This is a non-reactive unit (i.e., the effluent is copied from the influent).

    The following impact items should be pre-constructed for life cycle assessment:
    Battery, Solar.

    Parameters
    ----------
    ppl: int
        Total number of users for scaling of costs.

    References
    ----------
    [1] Trotochaud et al., Laboratory Demonstration and Preliminary Techno-Economic Analysis of an Onsite
    Wastewater Treatment System Environ. Sci. Technol. 2020, 54, (24), 16147–16155.
    https://dx.doi.org/10.1021/acs.est.0c02755

    [2] Duke Center for WaSH-AID Reclaimer design team data and guidance
    https://washaid.pratt.duke.edu/work/water-sanitation/reinvent-toilet-challenge
    
    [3] Eco-san water recycling toilet Reinvented Toilet design team bill of materials
    https://sanitation.ansi.org/EcoSanToilet

    '''

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream', **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with, F_BM_default=1)

        data = load_data(path=solar_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)


    def _design(self):
        design = self.design_results
        design['Solar'] = solar_quant = self.solar_capacity
        design['Battery'] = battery_quant = self.battery_kg

        self.construction = Construction(item='Solar', quantity=solar_quant, quantity_unit='m2')
        self.construction = Construction(item='Battery', quantity=battery_quant, quantity_unit='kg')

        self.add_construction(add_cost=False)

    def _cost(self):
        C = self.baseline_purchase_costs
        C['Battery System'] = self.battery_storage_cost + self.battery_holder_cost
        C['Solar Cost'] = (
            self.solar_cost * self.power_demand_30users +
            self.solar_module_system +
            self.inverter_cost
            )

        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio

        self.add_OPEX = self._calc_replacement_cost() + self._calc_maintenance_labor_cost()

    def _calc_maintenance_labor_cost(self):
        labor_cost = (self.wages * self.pannel_cleaning)  # USD/year
        return labor_cost / (365 * 24)  # USD/hr

    def _calc_replacement_cost(self):
        solar_panel_replacement_cost = (self.solar_cost * self.power_demand_30users) / self.solar_lifetime  # USD/year
        battery_replacement_cost = self.battery_storage_cost / self.battery_lifetime  # USD/year
        replacement_cost = (solar_panel_replacement_cost + battery_replacement_cost) / (365 * 24) * self.price_ratio  # USD/hr
        return replacement_cost


# %%

system_path = ospath.join(re_su_data_path, '_re_system.csv')

@price_ratio()
class ReclaimerSystem(SanUnit):
    '''
    System connection components for the Reclaimer system.

    This is a non-reactive unit (i.e., the effluent is copied from the influent).

    The following impact items should be pre-constructed for life cycle assessment:
    Steel.

    Parameters
    ----------
    ppl: int
        Total number of users for scaling of costs.
    if_gridtied: bool
        If using grid electricity instead of photovoltaic electricity.

    References
    ----------
    [1] Trotochaud et al., Laboratory Demonstration and Preliminary Techno-Economic Analysis of an Onsite
    Wastewater Treatment System Environ. Sci. Technol. 2020, 54, (24), 16147–16155.
    https://dx.doi.org/10.1021/acs.est.0c02755

    [2] Duke Center for WaSH-AID Reclaimer design team data and guidance
    https://washaid.pratt.duke.edu/work/water-sanitation/reinvent-toilet-challenge
    '''

    baseline_ppl = 30 # baseline population served by the Reclaimer system

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 if_gridtied=True, ppl=1, **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with, F_BM_default=1)
        self.ppl = ppl
        self.if_gridtied = if_gridtied

        data = load_data(path=system_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)


    def _design(self):
        design = self.design_results
        design['Steel'] = steel_quant = self.steel_weight * (self.ppl / self.baseline_ppl)  # linear scale
        self.construction = Construction(item='Steel', quantity=steel_quant, quantity_unit='kg')
        self.add_construction(add_cost=False)

    def _cost(self):
        C = self.baseline_purchase_costs
        N = self.N_reclaimers
        C['System'] = (
            self.T_nut +
            self.die_cast_hinge +
            self.SLS_locks +
            self.DC_round_key +
            self.handle_rod +
            self.eight_mm_bolt +
            self.button_headed_nut +
            self.twelve_mm_bolt +
            self.ten_mm_CSK +
            self.sixteen_mm_bolt +
            self.coupling_brass +
            self.socket +
            self.onehalf_tank_nipple +
            self.onehalf_in_coupling_brass +
            self.onehalf_in_fitting +
            self.plate +
            self.pump +
            self.three_way_valve +
            self.lofted_tank) * (1 + 0.05 * (N - 1))

        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio

        self.add_OPEX = self._calc_replacement_cost()

        if self.if_gridtied:
            power_demand = (self.power_demand_system / 1000) * N  # [W/day][1 kW/1000 W] = [kWh/d]
        else:
            power_demand = 0

        self.power_utility(power_demand)  # kW

    def _calc_replacement_cost(self):
        system_replacement_cost = self.baseline_purchase_costs['System'] * self.om_capital_ratio
        system_replacement_cost = system_replacement_cost / (365 * 24)  # convert from USD/year to USD/hour
        return system_replacement_cost

    @property
    def N_reclaimers(self):
        '''[int] Number of the reclaimer units needed, calculated by `ppl`/`baseline_ppl`.'''
        return ceil(self.ppl / self.baseline_ppl)


# %%

ultrafiltration_path = ospath.join(re_su_data_path, '_re_ultrafiltration.csv')

class ReclaimerUltrafiltration(SanUnit):
    '''
    Ultrafiltration in the Reclaimer system is used for removing suspended solids
    with automated backwash.

    The following impact items should be pre-constructed for life cycle assessment:
    Plastic, Steel.

    Parameters
    ----------
    ins : Iterable(stream)
        waste: liquid waste stream to be treated by ultrafiltration unit.
    outs : Iterable(stream)
        treated: treated liquid leaving ultrafiltration unit.
        retentate: concentrated retentate leaving ultrafiltration unit.
    ppl: int
        Total number of users for scaling of costs.
    if_gridtied: bool
        If using grid electricity instead of photovoltaic electricity.

    References
    ----------
    [1] Trotochaud et al., Laboratory Demonstration and Preliminary Techno-Economic Analysis of an Onsite
    Wastewater Treatment System Environ. Sci. Technol. 2020, 54, (24), 16147–16155.
    https://dx.doi.org/10.1021/acs.est.0c02755
    
    [2] Duke Center for WaSH-AID Reclaimer design team data and guidance
    https://washaid.pratt.duke.edu/work/water-sanitation/reinvent-toilet-challenge
    '''

    # Constants
    baseline_ppl = 30  # baseline population served by Reclaimer
    exponent_scale = 0.6  # exponential scaling constant

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 if_gridtied=True, ppl=1, **kwargs):

        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with, F_BM_default=1)
        self.if_gridtied = if_gridtied
        self.ppl = ppl

        data = load_data(path=ultrafiltration_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)

    _N_ins = 1
    _N_outs = 2

    def _run(self):
        waste = self.ins[0]
        treated, retentate = self.outs
        treated.copy_like(self.ins[0])

        self.retentate_prcd = (waste.imass['OtherSS'] * self.TSS_removal)  # mg/L
        retentate.imass['OtherSS'] = self.retentate_prcd

    def _design(self):
        design = self.design_results
        factor = self.ppl / self.baseline_ppl # linear scale
        design['Plastic'] = plastic_quant = self.Plastic_weight * factor
        design['Steel'] = steel_quant = self.Steel_weight * factor

        self.construction = (Construction(item='Plastic', quantity=plastic_quant, quantity_unit='kg'),
                             Construction(item='Steel', quantity=steel_quant, quantity_unit='kg'))

        self.add_construction(add_cost=False)

    def _cost(self):
        C = self.baseline_purchase_costs
        C['Pipes'] = self.one_in_pipe_SCH40 + self.onehalf_in_pipe_SCH40 + self.three_in_pipe_SCH80
        C['Fittings'] = (
            self.one_in_elbow_SCH80 +
            self.one_in_tee_SCH80 +
            self.one_in_SCH80 +
            self.one_onehalf_in_SCH80 +
            self.onehalf_in_SCH80 +
            self.three_in_SCH80_endcap +
            self.one_one_NB_MTA +
            self.one_onehalf_NB_MTA +
            self.foot_valve +
            self.one_onehalf_in_SCH80_threadedtee +
            self.three_in_pipe_clamp +
            self.one_in_pipe_clamp +
            self.onehalf_in_pipe_clamp +
            self.two_way_valve +
            self.UF_brush
            )
        C['UF_unit'] = self.UF_unit

        # Exponentially scale capital cost with number of users
        scale = (self.ppl / self.baseline_ppl) ** self.exponent_scale
        for equipment, cost in C.items():
            C[equipment] = cost * scale

        self.add_OPEX = self._calc_replacement_cost()

        N = self.N_reclaimers
        if self.if_gridtied:
            if N <=3:
                power_demand = getattr(self, f'power_demand_{N}') / 1000  # [W/day][1 kW/1000 W] = [kW/d]
            else: power_demand = self.power_demand_4 / 1000  # [W/day][1 kW/1000 W] = [kW/d]
        else:
            power_demand = 0

        self.power_utility(power_demand)

    def _calc_replacement_cost(self):
        scale = (self.ppl / self.baseline_ppl) ** self.exponent_scale

        pipe_replacement_cost = (
            self.one_in_pipe_SCH40 / self.one_in_pipe_SCH40_lifetime +
            self.onehalf_in_pipe_SCH40 / self.onehalf_in_pipe_SCH40_lifetime +
            self.three_in_pipe_SCH80 / self.three_in_pipe_SCH80_lifetime
            )

        fittings_replacement_cost = (
            self.one_in_elbow_SCH80 / self.one_in_elbow_SCH80_lifetime +
            self.one_in_tee_SCH80 / self.one_in_tee_SCH80_lifetime +
            self.one_in_SCH80 / self.one_in_SCH80_lifetime +
            self.one_onehalf_in_SCH80 / self.one_onehalf_in_SCH80_lifetime +
            self.onehalf_in_SCH80 / self.onehalf_in_SCH80_lifetime +
            self.three_in_SCH80_endcap / self.three_in_SCH80_endcap_lifetime +
            self.one_one_NB_MTA / self.one_one_NB_MTA_lifetime +
            self.one_onehalf_NB_MTA / self.one_onehalf_NB_MTA_lifetime +
            self.foot_valve / self.foot_valve_lifetime +
            self.one_onehalf_in_SCH80_threadedtee / self.one_onehalf_in_SCH80_threadedtee_lifetime +
            self.three_in_pipe_clamp / self.three_in_pipe_clamp_lifetime +
            self.one_in_pipe_clamp / self.one_in_pipe_clamp_lifetime +
            self.onehalf_in_pipe_clamp / self.onehalf_in_pipe_clamp_lifetime +
            self.two_way_valve / self.two_way_valve_lifetime +
            self.UF_brush / self.UF_brush_lifetime
            )

        uf_replacement_cost = self.UF_unit / self.UF_unit_lifetime

        total_replacement_cost = scale * (pipe_replacement_cost + fittings_replacement_cost + uf_replacement_cost)  # USD/year
        return total_replacement_cost / (365 * 24)  # USD/hr

    @property
    def N_reclaimers(self):
        '''[int] Number of the reclaimer units needed, calculated by `ppl`/`baseline_ppl`.'''
        return ceil(self.ppl / self.baseline_ppl)