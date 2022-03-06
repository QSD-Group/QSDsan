#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>
    Lewis Rowles <stetsonsc@gmail.com>
    Joy Zhang <joycheung1994@gmail.com>
    Lane To <lane20@illinois.edu>

Part of this module is based on the biosteam package:
https://github.com/BioSTEAMDevelopmentGroup/biosteam

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

from warnings import warn
from biosteam.units import HXprocess, HXutility
from .. import SanUnit, Construction
from ..utils import ospath, load_data, data_path, price_ratio

__all__ = (
    'HXprocess', 'HXutility',
    'HydronicHeatExchanger', 'HHXdryer',
    'OilHeatExchanger',
    )


# %%

# =============================================================================
# BioSTEAM-inherited ones
# =============================================================================

class HXprocess(SanUnit, HXprocess):
    '''
    Similar to :class:`biosteam.units.HXprocess`,
    but can be initilized with :class:`qsdsan.SanStream` and :class:`qsdsan.WasteStream`.

    See Also
    --------
    `biosteam.units.HXprocess <https://biosteam.readthedocs.io/en/latest/units/heat_exchange.html>`_
    '''

    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='Stream', F_BM_default=None,
                 *, U=None, dT=5., T_lim0=None, T_lim1=None,
                 material="Carbon steel/carbon steel",
                 heat_exchanger_type="Floating head",
                 N_shells=2, ft=None,
                 phase0=None,
                 phase1=None,
                 H_lim0=None,
                 H_lim1=None):
        SanUnit.__init__(self, ID, ins, outs, thermo,
                         init_with=init_with, F_BM_default=F_BM_default)

        #: [float] Enforced overall heat transfer coefficent (kW/m^2/K)
        self.U = U

        #: [float] Total heat transfered.
        self.Q = None

        #: Number of shells for LMTD correction factor method.
        self.N_shells = N_shells

        #: User imposed correction factor.
        self.ft = ft

        #: [float] Pinch temperature difference.
        self.dT = dT

        #: [float] Temperature limit of outlet stream at index 0.
        self.T_lim0 = T_lim0

        #: [float] Temperature limit of outlet stream at index 1.
        self.T_lim1 = T_lim1

        #: [float] Temperature limit of outlet stream at index 0.
        self.H_lim0 = H_lim0

        #: [float] Temperature limit of outlet stream at index 1.
        self.H_lim1 = H_lim1

        #: Enforced phase of outlet at index 0
        self.phase0 = phase0

        #: Enforced phase of outlet at index 1
        self.phase1 = phase1

        self.material = material
        self.heat_exchanger_type = heat_exchanger_type


class HXutility(SanUnit, HXutility):
    '''
    Similar to :class:`biosteam.units.HXutility`,
    but can be initilized with :class:`qsdsan.SanStream` and :class:`qsdsan.WasteStream`.

    See Also
    --------
    `biosteam.units.HXutility <https://biosteam.readthedocs.io/en/latest/units/heat_exchange.html>`_
    '''

    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='Stream', F_BM_default=None,
                 *, T=None, V=None, rigorous=False, U=None, H=None,
                 heat_exchanger_type="Floating head",
                 material="Carbon steel/carbon steel",
                 N_shells=2, ft=None, heat_only=None, cool_only=None):
        SanUnit.__init__(self, ID, ins, outs, thermo,
                         init_with=init_with, F_BM_default=F_BM_default)
        self.T = T #: [float] Temperature of outlet stream (K).
        self.V = V #: [float] Vapor fraction of outlet stream.
        self.H = H #: [float] Enthalpy of outlet stream.

        #: [bool] If true, calculate vapor liquid equilibrium
        self.rigorous = rigorous

        #: [float] Enforced overall heat transfer coefficent (kW/m^2/K)
        self.U = U

        #: Number of shells for LMTD correction factor method.
        self.N_shells = N_shells

        #: User imposed correction factor.
        self.ft = ft

        #: [bool] If True, heat exchanger can only heat.
        self.heat_only = heat_only

        #: [bool] If True, heat exchanger can only cool.
        self.cool_only = cool_only

        self.material = material
        self.heat_exchanger_type = heat_exchanger_type


# %%

hhx_path = ospath.join(data_path, 'sanunit_data/_hydronic_heat_exchanger.tsv')

@price_ratio(default_price_ratio=1)
class HydronicHeatExchanger(SanUnit):
    '''
    Hydronic heat exchanger is used for applications that require drying
    of the feedstock before the refinery is capable of processing the
    material. The heat is exchanged between the exhaust gas and water,
    which is then pumped into radiators connected to a dryer.
    The refinery monitors the temperature of the water to ensure that
    the feedstock is being sufficiently dried before entering the refinery.

    This class should be used together with :class:`~.sanunits.DryerFromHHX`.

    The following components should be included in system thermo object for simulation:
    N2O.

    The following impact items should be pre-constructed for life cycle assessment:
    OilHeatExchanger, Pump.

    Parameters
    ----------
    ins : WasteStream
        Hot gas.
    outs : WasteStream
        Hot gas.

    References
    ----------
    #!!! Reference?

    See Also
    --------
    :class:`~.sanunits.DryerFromHHX`

    '''

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with, F_BM_default=1)

        self.construction = (
            Construction('stainless_steel', linked_unit=self, item='StainlessSteel', quantity_unit='kg'),
            Construction('steel', linked_unit=self, item='Steel', quantity_unit='kg'),
            Construction('hydronic_heat_exchanger', linked_unit=self, item='HydronicHeatExchanger', quantity_unit='ea'),
            Construction('pump', linked_unit=self, item='Pump', quantity_unit='ea'),
            )

        data = load_data(path=hhx_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)


    _N_ins = 1
    _N_outs = 1

    def _run(self):
        hot_gas_in = self.ins[0]
        hot_gas_out = self.outs[0]
        # #!!! The following should probably b
        # hot_gas_out.copy_like(hot_gas_in)
        # hot_gas_out.phase = hot_gas_in.phase = 'g'
        hot_gas_out.phase = 'g'

        # Set temperature
        hot_gas_out.T = self.hhx_temp
        # Calculate the heat that was delivered to the HHX
        self.heat_output_water = \
            self.water_flowrate * \
            (self.water_out_temp-self.water_in_temp) * \
            self.water_density_kg_m_3 * \
            self.water_heat_capacity_k_j_kg_k * (60/1000)  # MJ/hr
        # Calculate losses through water pipe
        self.heat_loss_water_pipe = \
            (self.water_out_temp-self.ambient_temp) / \
            (self.water_r_pipe_k_k_w+self.water_r_insulation_k_k_w) * 3.6 # MJ/hr


    def _design(self):
        design = self.design_results
        constr = self.construction
        design['StainlessSteel'] = constr[0].quantity = self.heat_exchanger_hydronic_stainless
        design['Steel'] = constr[1].quantity  = self.heat_exchanger_hydronic_steel
        design['HydronicHeatExchanger'] = constr[2].quantity = 1
        design['Pump'] = constr[3].quantity = 17.2/2.72
        self.add_construction()


    def _cost(self):
        C = self.baseline_purchase_costs
        D = self.design_results
        C['Stainless steel'] = self.stainless_steel_cost * D['StainlessSteel']
        C['Steel'] = self.steel_cost * D['Steel']
        C['Misc. parts'] = \
            self.hhx_stack + \
            self.hhx_stack_thermocouple + \
            self.hhx_oxygen_sensor + \
            self.hhx_inducer_fan + \
            self.hhx_flow_meter + \
            self.hhx_pump + \
            self.hhx_water_in_thermistor + \
            self.hhx_water_out_thermistor + \
            self.hhx_load_tank + \
            self.hhx_expansion_tank + \
            self.hhx_heat_exchanger + \
            self.hhx_values + \
            self.hhx_thermal_well + \
            self.hhx_hot_water_tank + \
            self.hhx_overflow_tank

        ratio = self.price_ratio
        for equipment, cost in C.items():
           C[equipment] = cost * ratio

        # O&M cost converted to annual basis, labor included,
        # USD/yr only accounts for time running
        num = 1 / self.frequency_corrective_maintenance
        annual_maintenance = \
            self.service_team_adjustdoor_hhx*12 + \
            num * (self.service_team_replacewaterpump_hhx+self.service_team_purgewaterloop_hhx)

        self.add_OPEX =  annual_maintenance * self.service_team_wages / 60 / (365 * 24) # USD/hr (all items are per hour)

        self.power_utility(self.water_pump_power+self.hhx_inducer_fan_power) # kWh/hr


hhx_dryer_path = ospath.join(data_path, 'sanunit_data/_hhx_dryer.tsv')

@price_ratio(default_price_ratio=1)
class HHXdryer(SanUnit):
    '''
    Dryer is used in combination with :class:`~.sanunits.HydronicHeatExchanger`.

    The following components should be included in system thermo object for simulation:
    H2O, N, CH4, N2O.

    Parameters
    ----------
    ins : WasteStream
        Dewatered solids, heat.
    outs : WasteStream
        Dried solids, fugitive N2O, fugitive CH4.
    #!!! Get rid of these?
            set both as WasteStream
            others could be:
            Gases consist of SO2_emissions, NOx_emissions, CO_emissions,
            Hg_emissions, Cd_emissions, As_emissions, Dioxin_Furans_emissions.
    moisture_content_out : float
        Desired moisture content of the effluent.

    References
    ----------
    #!!! Reference?

    See Also
    --------
    :class:`~.sanunits.HydronicHeatExchanger`

    '''

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 moisture_content_out=0.35, **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with, F_BM_default=1)
        self.moisture_content_out = moisture_content_out

        data = load_data(path=hhx_dryer_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)


    _N_ins = 2
    _N_outs = 3

    def _run(self):
        waste_in, heat_in = self.ins
        waste_out, N2O, CH4  = self.outs
        waste_out.copy_like(self.ins[0])
        heat_in.phase = N2O.phase = CH4.phase = 'g'

        # Calculate heat needed to dry to the desired moisture content
        mc_in = waste_in.imass['H2O'] / waste_in.F_mass # fraction
        mc_out = self.moisture_content_out
        if mc_in < mc_out:
            warn(f'Moisture content of the influent stream ({mc_in:.2f}) '
                f'is smaller than the desired moisture content ({mc_out:.2f}).')
        TS_in = waste_in.F_mass - waste_in.imass['H2O'] # kg TS dry/hr
        waste_out.imass['H2O'] = TS_in/(1-mc_out)*mc_out
        water_to_dry = waste_in.imass['H2O'] - waste_out.imass['H2O'] # kg water/hr
        heat_needed_to_dry_35 = water_to_dry * self.energy_required_to_dry_sludge # MJ/hr

        #!!! heat_supplied isn't calculated from heat_in and heat loss?
        heat_supplied = self.dryer_heat_transfer_coeff * self.area_surface * (self.water_out_temp-self.feedstock_temp)

        if heat_needed_to_dry_35 > heat_supplied:
            warn('Heat required exceeds heat supplied.')

        #!!! Still needs this?
        # set to use COD or C for carbon based on influent composition
        # if waste_in.COD == 0:
        #     C_in = waste_in.imass['C']

        # else:
        #     C_in = self.carbon_COD_ratio * waste_in.COD

        # Emissions
        drying_CO2_to_air = (self.drying_CO2_emissions * self.carbon_COD_ratio
                             * waste_in.COD * waste_in.F_vol / 1000) # kg CO2 /hr

        #!!! Add or not?
        # add conversion factor for COD to TC?
        CH4.imass['CH4'] = drying_CH4_to_air = \
            self.drying_CH4_emissions * self.carbon_COD_ratio * \
            waste_in.COD * waste_in.F_vol / 1000 # kg CH4 /hr

        drying_NH3_to_air = self.drying_NH3_emissions * waste_in.TN * waste_in.F_vol / 1000 # kg NH3 /hr
        N2O.imass['N2O'] = drying_NH3_to_air * self.NH3_to_N2O # kg N2O /hr

        #!!! This is not correct, need to convert from CH4/CO2 to C
        # Reduce COD and TN in waste_out based on emissions
        waste_out._COD = (waste_in.COD * (waste_in.F_vol/waste_out.F_vol)) - ((drying_CO2_to_air + drying_CH4_to_air) / self.carbon_COD_ratio)
        # # Probably should be like this
        # waste_out._COD = (waste_in.COD*waste_in.F_vol - (drying_CO2_to_air/44*12+drying_CH4_to_air/16*12) / self.carbon_COD_ratio)

        #!!! Why this is commented out?
        # waste_out.imass['C'] -= ((drying_CO2_to_air + drying_CH4_to_air))

        #!!! Should do the conversion from NH3 to N
        waste_out.imass['N'] -= drying_NH3_to_air

        #!!! Are these used?
        # Jacket loss data and calculate losses due to convection and radiation
        self.jacket_heat_loss_conv = self.heat_transfer_coeff * (self.water_air_hx_temp1 - self.ambient_temp) * self.water_air_hx_area1 / 1000 / 0.2778 # MJ/hr
        self.jacket_heat_loss_radiation = self.radiative_emissivity * 5.67e-8 * ((self.water_air_hx_temp1 + 273)**4 - (self.ambient_temp + 273)**4) / 1000 / 0.2778 # MJ/hr
        self.jacket_heat_loss_sum = self.jacket_heat_loss_conv + self.jacket_heat_loss_radiation # MJ/hr
        #these can be calculated for pyrolysis and catalytic converter also


# %%

oilhx_path = ospath.join(data_path,'sanunit_data/_oil_heat_exchanger.tsv')

@price_ratio(default_price_ratio=1)
class OilHeatExchanger(SanUnit):
    '''
    Oil heat exchanger utilizes an organic Rankin cycle. This CHP system is
    used to generate additional electricity that the refinery and/or
    facility can use to decrease the units electrical demand on the
    electrical grid. This type of system is required for ISO 31800
    certification as the treatment unity needs to be energy independent
    when processing faecal sludge.

    The following components should be included in system thermo object for simulation:
    N2O.

    The following impact items should be pre-constructed for life cycle assessment:
    OilHeatExchanger, Pump.

    Parameters
    ----------
    ins : stream obj
        Hot gas.
    outs : stream obj
        Hot gas.

    References
    ----------
    #!!! Reference?

    '''
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with, F_BM_default=1)

        self.construction = (
            Construction('oil_heat_exchanger', linked_unit=self, item='OilHeatExchanger', quantity_unit='ea'),
            Construction('pump', linked_unit=self, item='Pump', quantity_unit='ea'),
            )

        data = load_data(path=oilhx_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)

    _N_ins = 1
    _N_outs = 1

    def _run(self):
        heat_in = self.ins[0]
        heat_out = self.outs[0]
        heat_out.copy_like(heat_in)
        # #!!! The following should probably b
        # hot_gas_out.phase = hot_gas_in.phase = 'g'
        heat_out.phase = 'g'

        # Set temperature
        heat_out.T = self.ohx_temp

        #!!! Below aren't used anywhere
        # Calculate the power that was delivered to the ORC
        self.power_delivery_orc = \
            self.oil_flowrate * \
            (self.oil_temp_out-self.oil_temp_in) * \
            self.oil_density * self.oil_specific_heat * (60/1000)  # MJ/hr
        # Calculate losses through pipe
        self.pipe_heat_loss = \
            (self.oil_temp_out-self.amb_temp) / \
            (self.oil_r_pipe_k_k_w+self.oil_r_insulation_k_k_w) * 3.6  # MJ/hr


    def _design(self):
        design = self.design_results
        constr = self.construction
        design['OilHeatExchanger'] = constr[0].quantity = 4/200
        design['Pump'] = constr[1].quantity = 2.834/2.27
        self.add_construction()


    def _cost(self):
        self.baseline_purchase_costs['Oil Heat Exchanger'] = self.orc_cost * self.price_ratio
        self.power_utility(self.oil_pump_power-self.oil_electrical_energy_generated) # kWh/hr