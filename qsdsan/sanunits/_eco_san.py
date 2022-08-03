#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:

    Tori Morgan <vlmorgan@illinois.edu>

    Yalin Li <mailto.yalin.li@gmail.com>

    Hannah Lohman <hlohman94@gmail.com>

This module contains unit operations used in the Eco-San system as described in
Li, M.; Xiaokang, Z. Technical Report of Eco-san Water Recycling System. Version 20190504-V3. 2019.
Note that the report is not publicly available, but information about the Eco-San system can be found at:
https://sanitation.ansi.org/EcoSanToilet
http://www.eco-san.cn/e_main.html

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

from qsdsan import SanUnit, Construction
from ._decay import Decay
from ._septic_tank import SepticTank
from ..utils import ospath, load_data, data_path, price_ratio

__all__ = (
    'EcoSanAerobic',
    'EcoSanAnaerobic',
    'EcoSanAnoxic',
    'EcoSanBioCost',
    'EcoSanECR',
    'EcoSanMBR',
    'EcoSanPrimary',
    'EcoSanSolar',
    'EcoSanSystem',
    )

es_su_data_path = ospath.join(data_path, 'sanunit_data/es')


# %%

aerobic_path = ospath.join(es_su_data_path, '_es_aerobic.tsv')

class EcoSanAerobic(SanUnit, Decay):
    '''
    Aerobic treatment unit in the Eco-San system.

    Note that costs and environmental impacts associated with unit is
    considered in :class:`EcoSanBioCost`.

    Parameters
    ----------
    ins : Iterable(stream)
        waste: waste stream to be treated.
    outs : Iterable(stream)
        treated: treated liquid leaving septic tank.
        CH4: fugitive CH4 emissions.
        N2O: fugitive N2O emissions.

    References
    ----------
    [1] Ma Li and Zhou Xiaokang. Technical Report of Eco-san Water Recycling System. Version 20190504-V3. 2019.

    See Also
    --------
    :ref:`qsdsan.sanunits.Decay <sanunits_Decay>`

    :class:`qsdsan.sanunits.EcoSanBioCost`
    '''
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 degraded_components=('OtherSS',), **kwargs):
        Decay.__init__(self, ID, ins, outs, thermo=thermo,
                       init_with=init_with, F_BM_default=1,
                       degraded_components=degraded_components)

        data = load_data(path=aerobic_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)

    _N_ins = 1
    _N_outs = 3

    def _run(self):
        Decay._first_order_run(self, if_capture_biogas=False, if_N2O_emission=True)


# %%

anaerobic_path = ospath.join(es_su_data_path, '_es_anaerobic.tsv')

class EcoSanAnaerobic(SanUnit, Decay):
    '''
    Anaerobic treatment unit in the Eco-San system.

    Note that costs and environmental impacts associated with unit is
    considered in :class:`EcoSanBioCost`.

    Parameters
    ----------
    ins : Iterable(stream)
        waste: waste stream to be treated.
    outs : Iterable(stream)
        treated: treated liquid leaving septic tank.
        CH4: fugitive CH4 emissions.
        N2O: fugitive N2O emissions.

    References
    ----------
    [1] Ma Li and Zhou Xiaokang. Technical Report of Eco-san Water Recycling System. Version 20190504-V3. 2019.

    See Also
    --------
    :ref:`qsdsan.sanunits.Decay <sanunits_Decay>`

    :class:`qsdsan.sanunits.EcoSanBioCost`
    '''
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 degraded_components=('OtherSS',), **kwargs):
        Decay.__init__(self, ID, ins, outs, thermo=thermo,
                       init_with=init_with, F_BM_default=1,
                       degraded_components=degraded_components)

        data = load_data(path=anaerobic_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)

    _N_ins = 1
    _N_outs = 3

    def _run(self):
        Decay._first_order_run(self, if_capture_biogas=False, if_N2O_emission=True)


# %%

anoxic_path = ospath.join(es_su_data_path, '_es_anoxic.tsv')

class EcoSanAnoxic(SanUnit, Decay):
    '''
    Anoxic treatment unit in the Eco-San system.

    Note that costs and environmental impacts associated with unit is
    considered in :class:`EcoSanBioCost`.

    Parameters
    ----------
    ins : Iterable(stream)
        waste: waste stream to be treated.
    outs : Iterable(stream)
        treated: treated liquid leaving septic tank.
        CH4: fugitive CH4 emissions.
        N2O: fugitive N2O emissions.

    References
    ----------
    [1] Ma Li and Zhou Xiaokang. Technical Report of Eco-san Water Recycling System. Version 20190504-V3. 2019.

    See Also
    --------
    :ref:`qsdsan.sanunits.Decay <sanunits_Decay>`

    :class:`qsdsan.sanunits.EcoSanBioCost`
    '''
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 degraded_components=('OtherSS',), **kwargs):
        Decay.__init__(self, ID, ins, outs, thermo=thermo,
                       init_with=init_with, F_BM_default=1,
                       degraded_components=degraded_components)

        data = load_data(path=anoxic_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)

    _N_ins = 1
    _N_outs = 3

    def _run(self):
        Decay._first_order_run(self, if_capture_biogas=False, if_N2O_emission=True)


# %%

bio_cost_path = ospath.join(es_su_data_path, '_es_bio_cost.tsv')

class EcoSanBioCost(SanUnit):
    '''
    A non-reactive unit to account for all costs and environmental impacts
    associated with the biological treatment units within the Eco-San system
    (other than EcoSanPrimary and EcoSanMBR).

    The following impact items should be pre-constructed for life cycle assessment:
    FRP, Pump.

    References
    ----------
    [1] Ma Li and Zhou Xiaokang. Technical Report of Eco-san Water Recycling System. Version 20190504-V3. 2019.
    '''
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream', **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with, F_BM_default=1)

        data = load_data(path=bio_cost_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)


    def _design(self):
        design = self.design_results
        design['FRP'] = FRP_quant = self.FRP
        design['Pump'] = pump_quant = self.pump_ammount
        self.construction = (
            Construction(item='Pump', quantity = pump_quant, quantity_unit = 'ea'),
            Construction(item='FRP', quantity = FRP_quant, quantity_unit = 'kg')
            )
        self.add_construction(add_cost=False)

    def _cost(self):
        C = self.baseline_purchase_costs
        C['Tanks'] = self.tank_cost
        C['Misc. parts'] = (
            self.middle_tank_pump_cost +
            self.fan_cost +
            self.level_guage_cost +
            self.pipe_cost +
            self.regulating_tank_pump_cost
            )

        bio_replacement_parts_annual_cost = (
            self.guage_life * self.level_guage_cost +
            self.pump_life * self.pump_replacement_cost +
            self.bio_tube_diffuser_replacement_cost * self.diffuser_life
            )

        self.add_OPEX =  bio_replacement_parts_annual_cost / (365 * 24) # USD/hr


# %%

ecr_path = ospath.join(es_su_data_path, '_es_ecr.xlsx')

@price_ratio()
class EcoSanECR(SanUnit, Decay):
    '''
    Electrochemical reactor (ECR) in the Eco-San system with chlorine dosing.

    The following impact items should be pre-constructed for life cycle assessment:
    Metal, Pump.

    Parameters
    ----------
    ins : Iterable(stream)
        waste: waste stream to be treated.
        salt: NaCl to be added for treatment.
        HCl: HCl to be added for treatment.
    outs : Iterable(stream)
        treated: treated liquid leaving septic tank.
        CH4: fugitive CH4 emissions.
        N2O: fugitive N2O emissions.
    if_after_MBR: bool
        If this unit is used after a membrane bioreactor (MBR),
        ECR after an MBR will have lower costs compared to the scenario without an MBR.

    References
    ----------
    [1] Ma Li and Zhou Xiaokang. Technical Report of Eco-san Water Recycling System. Version 20190504-V3. 2019.

    See Also
    --------
    :ref:`qsdsan.sanunits.Decay <sanunits_Decay>`
    '''

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 degraded_components=('OtherSS',),
                 if_after_MBR=False, **kwargs):
        Decay.__init__(self, ID, ins, outs, thermo=thermo,
                       init_with=init_with, F_BM_default=1,
                       degraded_components=degraded_components)

        self._if_after_MBR = None # initialize the attr
        self.if_after_MBR = if_after_MBR

        for attr, value in kwargs.items():
            setattr(self, attr, value)

    def _refresh_data(self):
        sheet_name = 'standalone' if not self.if_after_MBR else 'after_MBR'
        data = load_data(path=ecr_path, sheet_name=sheet_name)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

    _N_ins = 3
    _N_outs = 3

    def _run(self):
        self.ins[1].imass['NaCl'] = self.salt_dosing / 24 # NaCl
        HCl_density = 1.2 # g/ml
        self.ins[2].imass['HCl'] =  self.HCl_life/52/24/7 * HCl_density * 1000 # kg/h

        Decay._first_order_run(self, degraded_components=('OtherSS',),
                               if_capture_biogas=False, if_N2O_emission=True)


    def _design(self):
        design = self.design_results
        design['Metal'] = electrode_quant = self.electrode
        design['Pump'] = pump_quant = self.pump
        self.construction = (
            Construction(item='Pump', quantity = pump_quant, quantity_unit = 'ea'),
            Construction(item='Metal', quantity = electrode_quant, quantity_unit = 'kg')
        )
        self.add_construction(add_cost=False)


    def _cost(self):
        C = self.baseline_purchase_costs
        C['level_guage_cost'] = self.level_guage_cost
        C['pump_cost'] = self.pump_cost
        C['fan_cost'] = self.fan_cost
        C['salt_dosing_device_cost'] = self.salt_dosing_device_cost
        C['UPVC_pipe_cost'] = self.UPVC_pipe_cost
        C['UPVC_electric_ball_cost'] = self.UPVC_electric_ball_cost
        C['GAC_cost'] = self.GAC_cost
        C['electrode_cost'] = self.electrode_cost
        C['ECR_reactor_cost'] = self.ECR_reactor_cost
        C['power_supply_cost'] = self.power_supply_cost
        C['plastic_spraying_cabinent_cost'] = self.plastic_spraying_cabinent_cost

        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio

        ECR_replacement_parts_annual_cost = (
            self.HCl_replacement_cost * self.HCl_life +
            self.salt_replacement_cost * 52
        )
        if not self.if_after_MBR:
            ECR_replacement_parts_annual_cost += (
                self.electrode_replacement_cost * self.electrode_life +
                self.pump_cost * self.pump_life +
                self.level_guage_replacement_cost * self.level_guage_life +
                self.GAC_cost * self.GAC_filter_life
            )
        else:
            ECR_replacement_parts_annual_cost += self.electrode_replacement_cost

        self.add_OPEX =  ECR_replacement_parts_annual_cost*self.price_ratio / (365 * 24) # USD/hr

        self.power_utility(self.power_demand * self.working_time) # kWh/hr


    @property
    def if_after_MBR(self):
        '''
        [bool] If this unit is used after a membrane bioreactor (MBR),
        ECR after an MBR will have lower costs compared to the scenario without an MBR.
        '''
        return self._if_after_MBR
    @if_after_MBR.setter
    def if_after_MBR(self, i):
        if i != self._if_after_MBR:
            self._if_after_MBR = i
            self._refresh_data()
        return self._if_after_MBR


# %%

mbr_path = ospath.join(es_su_data_path, '_es_mbr.tsv')

@price_ratio()
class EcoSanMBR(SanUnit, Decay):
    '''
    Membrane bioreactor for the Eco-San system.

    The following impact items should be pre-constructed for life cycle assessment:
    FRP.

    Parameters
    ----------
    ins : Iterable(stream)
        waste: waste stream to be treated.
    outs : Iterable(stream)
        treated: treated liquid leaving septic tank.
        CH4: fugitive CH4 emissions.
        N2O: fugitive N2O emissions.

    References
    ----------
    [1] Ma Li and Zhou Xiaokang. Technical Report of Eco-san Water Recycling System. Version 20190504-V3. 2019.

    See Also
    --------
    :ref:`qsdsan.sanunits.Decay <sanunits_Decay>`
    '''

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 degraded_components=('OtherSS',), **kwargs):
        Decay.__init__(self, ID, ins, outs, thermo=thermo,
                       init_with=init_with, F_BM_default=1,
                       degraded_components=degraded_components)

        data = load_data(path=mbr_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)

    _N_ins = 1
    _N_outs = 3

    def _run(self):
        Decay._first_order_run(self, if_capture_biogas=False, if_N2O_emission=True)


    def _design(self):
        self.design_results['FRP'] = FRP_quant = self.FRP_per_tank
        self.construction = (Construction(item='FRP', quantity = FRP_quant, quantity_unit = 'kg'),)
        self.add_construction(add_cost=False)

    def _cost(self):
        C = self.baseline_purchase_costs
        C['MBR tank'] = self.MBR_cost
        C['Bio tank'] = self.tank_cost
        C['Primary tank'] = self.primary_tank_cost
        C['Regulating pump'] = self.regulating_pump_cost
        C['Fan'] = self.fan_cost
        C['Gauge'] = self.gauge_cost
        C['Pipe'] = self.pipe_cost
        C['Mid tank pump'] = self.mid_tank_pump

        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio

        MBR_replacement_parts_annual_cost = (
            self.MBR_cost * self.MBR_replacement +
            self.gauge_cost * self.gauge_life
            )

        self.add_OPEX =  MBR_replacement_parts_annual_cost / (365 * 24) * self.price_ratio # USD/hr
        self.power_utility(self.power_demand * self.working_time)


# %%

primary_path = ospath.join(es_su_data_path, '_es_primary.xlsx')

class EcoSanPrimary(SepticTank):
    '''
    The primary treatment of the Eco-San system uses anaerobic digestion
    to treat wastes (similar to a septic tank).

    It can be used in conjunction with a membrane bioreactor (MBR)
    to recovery N and P as struvite.

    The following impact items should be pre-constructed for life cycle assessment:
    FRP.

    Parameters
    ----------
    ins : Iterable(stream)
        Waste for treatment, Mg(OH2) used if struvite is produced.
    outs : Iterable(stream)
        Treated waste, fugitive CH4, fugitive N2O, generated sludge, and struvite.

    References
    ----------
    [1] Ma Li and Zhou Xiaokang. Technical Report of Eco-san Water Recycling System. Version 20190504-V3. 2019.

    See Also
    --------
    :ref:`qsdsan.sanunits.Decay <sanunits_Decay>`

    :class:`qsdsan.sanunits.SepticTank`
    '''

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 if_with_MBR=False, sludge_moisture_content=0.95, **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with, F_BM_default=1)
        self.sludge_moisture_content = sludge_moisture_content
        self._if_with_MBR = None # initialize the attr
        self.if_with_MBR = if_with_MBR

        for attr, value in kwargs.items():
            setattr(self, attr, value)

    def _refresh_data(self):
        sheet_name = 'with_MBR' if not self.if_with_MBR else 'without_MBR'
        data = load_data(path=primary_path, sheet_name=sheet_name)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data


    def _design(self):
        self.design_results['FRP'] = FRP_quant = self.FRP_per_tank
        self.construction = (Construction(item='FRP', quantity=FRP_quant, quantity_unit='kg'),)
        self.add_construction(add_cost=False)


    def _cost(self):
        # No replacement parts for the tank and cleaning performed
        # is considered in TEA
        ratio = self.price_ratio
        self.baseline_purchase_costs['Tanks'] = self.FRP_tank_cost * ratio
        self.add_OPEX =  (self.Mg_dose * self.MgOH2_cost / 365 / 24) * ratio
        self.power_utility(self.power_demand)


    @property
    def user_scale_up(self):
        raise AttributeError('No scaling is considered for `EcoSanPrimary`.')


    @property
    def if_generate_struvite(self):
        '''[bool] If generating struvite (True when with MBR).'''
        if self.if_with_MBR: return True
        return False

    @property
    def if_struvite_in_sludge(self):
        '''
        [bool] If the generated struvite is in sludge, always False
        (when MBR is included, struvite is generated as a separate stream).
        '''
        return False

    @property
    def if_with_MBR(self):
        '''
        [bool] If this unit has a membrane bioreactor (MBR),
        primary treatment with an MBR is able to recovery N and P as struvite.
        '''
        return self._if_with_MBR
    @if_with_MBR.setter
    def if_with_MBR(self, i):
        if i != self._if_with_MBR:
            self._if_with_MBR = i
            self._refresh_data()
        return self._if_with_MBR


# %%

solar_path = ospath.join(es_su_data_path, '_es_solar.tsv')

@price_ratio()
class EcoSanSolar(SanUnit):
    '''
    A non-reactive unit to account for all costs and environmental impacts
    associated with the solar energy system.

    The following impact items should be pre-constructed for life cycle assessment:
    FRP, Pump.

    References
    ----------
    [1] Ma Li and Zhou Xiaokang. Technical Report of Eco-san Water Recycling System. Version 20190504-V3. 2019.
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
        design['Solar system'] = solar_quant = self.solar_capacity
        design['Battery'] = battery_quant = self.battery_kg
        self.construction = (
            Construction(item='Solar', quantity = solar_quant, quantity_unit = 'm2'),
            Construction(item='Battery', quantity = battery_quant, quantity_unit = 'kg')
            )
        self.add_construction(add_cost=False)


    def _cost(self):
        C = self.baseline_purchase_costs
        C['Battery System'] = self.battery_storage_cost + self.battery_holder_cost
        C['Solar Cost'] = self.solar_cost + self.solar_module_system + self.inverter_cost

        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio

        # The following was originally commented out, but OPEX should be included - Yalin Li 2022-08-03
        self.add_OPEX =  self.solar_replacement * self.solar_cost / (365 * 24) # USD/hr


# %%

system_path = ospath.join(es_su_data_path, '_es_system.tsv')

@price_ratio()
class EcoSanSystem(SanUnit):
    '''
    A non-reactive unit to account for all costs and environmental impacts
    associated with controlling, recycling, and other miscellaneous parts.

    The following impact items should be pre-constructed for life cycle assessment:
    Pump.

    References
    ----------
    [1] Ma Li and Zhou Xiaokang. Technical Report of Eco-san Water Recycling System. Version 20190504-V3. 2019.
    '''
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream', **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with, F_BM_default=1)

        data = load_data(path=system_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)


    def _design(self):
        self.design_results['Pump'] = pump_quant = self.pump_amount
        self.construction = (Construction(item='Pump', quantity = pump_quant, quantity_unit = 'ea'),)
        self.add_construction(add_cost=False)


    def _cost(self):
        C = self.baseline_purchase_costs
        C['Tanks'] = self.high_level_tank_cost
        C['Misc. parts'] = (
            self.booster_pump_cost +
            self.level_guage_cost +
            self.UPVC_pipe_cost +
            self.filter_cost +
            self.control_system_cost +
            self.shell_cost +
            self.container_cost +
            self.packaging +
            self.door_cost +
            self.heat_preservation
            )
        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio

        # The following was originally commented out, but OPEX should be included - Yalin Li 2022-08-03
        self.add_OPEX =  self.filter_replacement * self.filter_cost / (365 * 24) # USD/hr