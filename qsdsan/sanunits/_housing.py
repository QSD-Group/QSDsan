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
from ..utils import load_data, data_path, ospath, price_ratio

__all__ = ('HousingBiogenicRefinery',)

br_path = ospath.join(data_path, 'sanunit_data/_housing_biogenic_refinery.tsv')

@price_ratio(default_price_ratio=1)
class HousingBiogenicRefinery(SanUnit):
    '''
    Housing for the biogenic refinery which is composed of the casing around the system,
    containers, and the concrete slab.
    No process algorithm is included, only design (including) cost algorithms are included.

    The following impact items should be pre-constructed for life cycle assessment:
    Steel, StainlessSteelSheet, Concrete.

    Parameters
    ----------
    ins : stream obj
        Influent stream.
    outs : stream obj
        Effluent stream, is copied from the influent stream.
    '''

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 const_wage=15, const_person_days=100, **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with, F_BM_default=1)
        self.const_wage = const_wage
        self.const_person_days = const_person_days

        self.construction = (
            Construction('steel', linked_unit=self, item='Steel', quantity_unit='kg'),
            Construction('stainlesssteelsheet', linked_unit=self, quantity_unit='kg'),
            Construction('concrete', item='Concrete', linked_unit=self, quantity_unit='m3'),
            )

        data = load_data(path=br_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)

    
    def _design(self):
        design = self.design_results
        constr = self.construction
        design['Steel'] = constr[0].quantity = 2000 + 4000
        design['StainlessSteelSheet'] = constr[1].quantity = 4.88 * 2 * 3 * 4.5
        design['Concrete'] = constr[2].quantity = self.concrete_thickness * self.footprint
        self.add_construction()

    
    def _cost(self):
        D = self.design_results
        C = self.baseline_purchase_costs
        C['Containers'] = self.container20ft_cost + self.container40ft_cost
        C['Equip Housing'] = D['StainlessSteelSheet'] / 4.88 * self.stainless_steel_housing
        C['Concrete'] = D['Concrete'] * self.concrete_cost

        #!!! Why is labor included in the capital cost dict? Why it is multiplied by the `price_ratio`?
        C['Labor'] = self.const_wage * self.const_person_days
        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio