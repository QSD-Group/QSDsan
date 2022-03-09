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
from ..utils import ospath, load_data, data_path

__all__ = ('IonExchangeNH3',)

ix_NH3_path = ospath.join(data_path, 'sanunit_data/_ion_exchange_NH3.tsv')

class IonExchangeNH3(SanUnit):
    '''
    Ion Exchange for N recovery from liquid stream. Concentrated NH3 is recovered.

    The following components should be included in system thermo object for simulation:
    NH3, Polystyrene, H2SO4.

    The following impact items should be pre-constructed for life cycle assessment:
    PVC, PE.

    Parameters
    ----------
    ins : Iterable (stream obj)
        Liquid waste, fresh resin, H2SO4.
    outs : Iterable (stream obj)
        Treated waste, spent resin, concentrated NH3.

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
                 **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with,
                         F_BM_default=1)

        self.construction = (
            Construction('pvc', linked_unit=self, item='PVC', quantity_unit='kg'),
            Construction('pe', linked_unit=self, item='PE', quantity_unit='kg'),
            )

        data = load_data(path=ix_NH3_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)


    _N_ins = 3
    _N_outs = 3

    def _run(self):
        waste, resin_in, H2SO4 = self.ins
        treated, resin_out, conc_NH3 = self.outs
        treated.copy_like(waste)

        #!!! During storage most N as urea goes to NH3, should that
        # conversion be added or just use total N here?
        N_recovered = waste.imass['NH3'] * self.N_rec # kg N / hr
        treated.imass['NH3'] =  waste.imass['NH3'] - N_recovered # kg N / hr
        conc_NH3.imass['NH3'] = N_recovered # kg N / hr

        # following Tarpeh et al. 2018 estimates for regenerating resin
        # assume volume of eluent = volume of urine
        # !!! cost associated with this influent stream neg compared to cost of resin and acid?
        #conc_NH3.F_vol = waste.F_vol
        # 0.1 M H2SO4 (0.65%) used to regenerate resin

        # !!! need to add SO4 as a component?
        # conc_NH3.imass['SO4'] = 0.1 * 98 * conc_NH3.F_vol # kg SO4 / hr

        resin_demand_influent = waste.TN / self.resin_lifetime / self.ad_density / 14 # kg resin / m3 treated
        resin_demand_time = resin_demand_influent * waste.F_vol # kg resin / hr
        resin_in.imass['Polystyrene'] = resin_out.imass['Polystyrene'] = resin_demand_time

        acid_demand_influent = waste.TN * self.vol_H2SO4 / self.ad_density / 14 # L acid / L treated
        acid_demand_time = acid_demand_influent * waste.F_vol * 1000 * 1.83 # kg acid / hr
        H2SO4.imass['H2SO4'] = acid_demand_time


    def _design(self):
        design = self.design_results
        constr = self.construction
        N_column = self.N_column

        design['PVC'] = constr[0].quantity = \
            N_column * self.column_length * self.pvc_mass # kg PVC
        tubing_quant = N_column * self.tubing_length * self.tubing_mass # kg PE
        design['PE'] = constr[1].quantity = tubing_quant + self.N_tank*self.tank_mass

        self.add_construction()


    def _cost(self):
        C = self.baseline_purchase_costs
        N_column = self.N_column
        C['PVC'] = self.cost_PVC_column * self.column_length * N_column
        C['Tubing'] = self.cost_tubing * self.tubing_length*N_column
        C['Tank'] = N_column/3 * self.tank_cost # one tank has three columns

        #!!! Still need this?
        # costs associated with full time operators can be added in the TEA as staff

    @property
    def N_column(self):
        '''[int] Number of resin columns.'''
        return ceil(self.ins[0].F_vol*1000*24/self.column_daily_loading_rate)

    @property
    def N_tank(self):
        '''[int] Number of tanks for cost estimation, might be float (instead of int).'''
        return self.N_column/3