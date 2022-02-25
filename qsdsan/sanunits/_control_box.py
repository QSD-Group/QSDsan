#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Lewis Rowles <stetsonsc@gmail.com>
    Yalin Li <zoe.yalin.li@gmail.com>
    Lane To <lane20@illinois.edu>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

# %%

from .. import SanUnit, Construction
from ..utils import ospath, load_data, data_path, price_ratio

__all__ = ('ControlBoxOP',)

op_path = ospath.join(data_path,'sanunit_data/_control_box_op.tsv')

@price_ratio(default_price_ratio=1)
class ControlBoxOP(SanUnit):
    '''
    Control box (industrial control panel) for the omni processor.
    No process algorithm is included, only design (including) cost algorithms are included.

    The following impact items should be pre-constructed for life cycle assessment:
    Electronics, ElectricConnectors, ElectricCables.

    Parameters
    ----------
    ins : stream obj
        Influent stream.
    outs : stream obj
        Effluent stream, is copied from the influent stream.

    References
    ----------
    #!!! Reference?
    '''

    _N_ins = _N_outs = 1

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 **kwargs):
        SanUnit.__init__(self, ID, ins, outs, F_BM_default=1)

        self.construction = (
            Construction('electronics', linked_unit=self, item='Electronics', quantity_unit='kg'),
            Construction('electric_connectors', linked_unit=self, item='ElectricConnectors', quantity_unit='kg'),
            Construction('electric_cables', linked_unit=self, item='ElectricCables', quantity_unit='m'),
            )

        data = load_data(path=op_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)

    def _design(self):
        design = self.design_results
        constr = self.construction
        design['Electronics'] = constr[0].quantity = 2
        design['ElectricConnectors'] = constr[1].quantity = 0.5
        design['ElectricCables'] = constr[2].quantity = 3
        self.add_construction()

    # Cost based on amount of steel and stainless plus individual components
    def _cost(self):
        self.baseline_purchase_costs['Misc. parts'] = (
            self.icp_controller_board +
            self.icp_variable_frequence_drives +
            self.icp_power_meter +
            self.icp_line_filter +
            self.icp_transformer +
            self.icp_power_meter_transformer +
            self.icp_AC_to_DC +
            self.icp_DC_to_AC +
            self.icp_touch_screen) * self.price_ratio

        cb_annual_maintenance = \
            (self.electrician_replacecables_icp+self.electrician_replacewires_icp)*self.certified_electrician_wages + \
            self.service_team_replacetouchscreen_icp*self.service_team_wages + \
            self.facility_manager_configurevariable_icp*self.facility_manager_wages + \
            (self.biomass_controls_replaceboard_icp+self.biomass_controls_codemalfunctioning_icp)*self.biomass_controls_wages

        self.add_OPEX =  cb_annual_maintenance/60/self.frequency_corrective_maintenance/(365*24) # USD/hr
        # kWh/hr
        self.power_utility(self.icp_controller_board_power+self.icp_variable_frequence_drives_power)