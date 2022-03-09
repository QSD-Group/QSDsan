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
from ..utils import ospath, load_data, data_path

__all__ = ('Grinder',)

grinder_path = ospath.join(data_path, 'sanunit_data/_grinder.tsv')

#!!! Need `price_ratio`?
class Grinder(SanUnit):
    '''
    Grinder is used to break up solids.

    .. note:

        Moisture content of the effluent is adjusted to be 65%, although the grinder
        itself can't change the moisture. This assumption was made based on pilot experiments.

    The following components should be included in system thermo object for simulation:
    H2O, OtherSS.

    The following impact items should be pre-constructed for life cycle assessment:
    Steel.

    Parameters
    ----------
    ins : stream obj
        Influent stream.
    outs : stream obj
        Effluent stream.
    moisture_content_out : float
        Moisture content of the effleunt solids stream.

    References
    ----------
    #!!! Reference?

    '''

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 moisture_content_out=0.65, **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with,
                         F_BM_default=1)
        self.moisture_content_out = moisture_content_out
        self.construction = (
            Construction('steel', linked_unit=self, item='Steel', quantity_unit='kg'),
            )

        data = load_data(path=grinder_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)

    _N_ins = 1
    _N_outs = 1

    def _run(self):
        waste_in = self.ins[0]
        waste_out = self.outs[0]
        waste_out.copy_like(waste_in)

        mc_in = waste_in.imass['H2O'] / waste_in.F_mass # fraction
        mc_out = self.moisture_content_out
        if mc_in < mc_out:
            raise RuntimeError(f'Moisture content of the influent stream ({mc_in:.2f}) '
                               f'is smaller than the desired moisture content ({mc_out:.2f}).')
        TS_in = waste_in.F_mass - waste_in.imass['H2O'] # kg TS dry/hr
        waste_out.imass['H2O'] = TS_in/(1-mc_out)*mc_out
        waste_out._COD = waste_in.COD * waste_in.F_vol / waste_out.F_vol

    def _design(self):
        self.design_results['Steel'] = self.construction[0].quantity = self.grinder_steel
        self.add_construction(add_cost=False)


    def _cost(self):
        self.baseline_purchase_costs['Grinder'] = self.grinder
        self.power_utility(self.grinder_electricity * self.ins[0].imass['OtherSS']) # kWh/hr