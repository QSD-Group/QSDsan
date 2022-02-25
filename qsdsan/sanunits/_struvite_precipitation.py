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

from math import ceil
from .. import SanUnit, Construction
from ..utils.loading import ospath, load_data, data_path

__all__ = ('StruvitePrecipitation',)

struvite_path = ospath.join(data_path, '/sanunit_data/_struvite_precipitation.tsv')

#!!! No `price_ratio`?
class StruvitePrecipitation(SanUnit):
    '''
    Stuvite precipitation for P recovery from liquid stream as solid stuvite.

    The following components should be included in system thermo object for simulation:
    P, NH3, K, MagnesiumHydroxide, MagnesiumCarbonate, Struvite, FilterBag.

    The following impact items should be pre-constructed for life cycle assessment:
    StainlessSteel, PVC.


    Parameters
    ----------
    ins : Iterable (stream obj)
        Liquid waste, Mg(OH)2, MgCO3, filter bag.
    outs : Iterable (stream obj)
        Treated waste, struvite.
    Mg_molar_split : Iterable(float)
        The molar split between Mg(OH)2 and MgCO3.
        (1, 0) means all Mg is added as Mg(OH)2 and (0,1) means all MgCO3.

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

    '''

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 Mg_molar_split=(1,0), **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with,
                         F_BM_default=1)
        self.Mg_molar_split = Mg_molar_split

        self.construction = (
            Construction('stainless_steel', linked_unit=self, item='StainlessSteel', quantity_unit='kg'),
            Construction('pvc', linked_unit=self, item='PVC', quantity_unit='kg'),
            )

        data = load_data(path=data_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)


    _N_ins = 4
    _N_outs = 2

# in _run: define influent and effluent streams and treatment processes
    def _run(self):
        waste, magnesium_hydroxide, magnesium_carbonate, bag_filter = self.ins
        treated, struvite = self.outs
        treated.copy_like(waste)
        for ws in self.ins[1:]+self.outs[1:]:
            ws.phase = 's'

        Mg_molar_split = self.Mg_molar_split
        Mg_molar_split /= sum(Mg_molar_split) # normalize to a sum of 1
        Mg_dose, P_rec_1 = self.Mg_dose, self.P_rec_1
        # Total P recovered [kmol P/hr]
        P_recovered = waste.imol['P'] * P_rec_1
        treated.imol['P'] =  waste.imol['P'] - P_recovered

        # Total N recovered, 1:1 mol/mol N:P in struvite ((NH4)MgPO4·6(H2O))
        treated.imol['NH3'] =  waste.imol['NH3'] - P_recovered

        # Total K recovered [kg K/hr]
        treated.imol['K'] =  waste.imol['K'] * (1-self.K_rec_1)

        magnesium_hydroxide.empty()
        magnesium_carbonate.empty()
        magnesium_hydroxide.imol['MagnesiumHydroxide'] = waste.imol['P']*Mg_dose*Mg_molar_split[0]
        magnesium_carbonate.imol['MagnesiumCarbonate'] = waste.imol['P']*Mg_dose*Mg_molar_split[1]
        struvite.imol['Struvite'] = P_recovered
        bag_filter.imass['FilterBag'] = self.N_tank * self.cycles_per_day / self.filter_reuse / 24 # bags/hr


    def _design(self):
        design = self.design_results
        constr = self.construction
        N_tank = self.N_tank
        design['StainlessSteel'] = constr[0].quantity = N_tank * self.reactor_weight # kg SS
        design['PVC'] = constr[1].quantity = N_tank * self.material_P_pipe * self.pvc_mass # kg PVC
        self.add_construction()

    def _cost(self):
        C = self.baseline_purchase_costs
        N_tank = self.N_tank
        C['Reactor'] = N_tank * self.cost_P_reactor
        C['Stirrer'] = N_tank * self.cost_P_stirrer
        C['PVC'] = N_tank * self.material_P_pipe * self.cost_P_pipe
        self.add_OPEX =  0.35/4 * self.baseline_purchase_cost / (365 * 24) # USD/hr (all items are per hour)

        #!!! Still need this?
        # costs associated with full time operators can be added in the TEA as staff


    @property
    def N_tank(self):
        '''[int] Number of reactor tanks.'''
        return ceil(self.ins[0].F_vol*1000*24/self.cycles_per_day/self.reactor_volume)