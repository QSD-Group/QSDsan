#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Lewis Rowles <stetsonsc@gmail.com>
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''


# %%

from . import Toilet
from ... import WasteStream, Construction
from ...utils import ospath, load_data, data_path

__all__ = ('MURTToilet',)


murt_path = ospath.join(data_path, 'sanunit_data/_murt_toilet.tsv')


# %%

class MURTToilet(Toilet):
    '''
    Multi-Unit Reinvented Toilet.

    To enable life cycle assessment, the following impact items should be pre-constructed:
    `Ceramic`, `Fan`.

    Parameters
    ----------
    ins : :class:`~.WasteStream`
        Excreta.
    outs : Iterable(:class:`~.WasteStream`)
        Recyclable mixed excreta, stream leached to soil, fugitive CH4, and fugitive N2O.
    lifetime : int
        Lifetime of this UDDT, [yr].
    N_squatting_pans : int
        Number of squatting pans.
    N_urinal : int
        Number of urinals.
    N_fan : int
        Number of fans.

    References
    ----------
    #!!! What should be used for the reference?

    See Also
    --------
    :ref:`qsdsan.sanunits.Toilet <sanunits_Toilet>`
    '''

    #change to be 1 squatting pan and 2 urinals1
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 N_user=1, N_toilet=1, lifetime=25,
                 if_toilet_paper=True, if_flushing=True, if_cleansing=False,
                 if_desiccant=False, if_air_emission=True, if_ideal_emptying=True,
                 CAPEX=0, OPEX_over_CAPEX=0.075,
                 N_squatting_pan=1, N_urinal=1, N_fan=1,
                 **kwargs):

        Toilet.__init__(self, ID, ins, outs, thermo, init_with, N_user, N_toilet,
                        if_toilet_paper, if_flushing, if_cleansing, if_desiccant,
                        if_air_emission, if_ideal_emptying, CAPEX, OPEX_over_CAPEX)
        self.lifetime = lifetime
        #!!! Legacy note: change to be 1 squatting pan and 2 urinals
        self.N_squatting_pan = N_squatting_pan
        self.N_urinal = N_urinal
        self.N_fan = N_fan
        self.construction = (
            Construction('ceramic', linked_unit=self, item='Cement', quantity_unit='kg'),
            Construction('fan', linked_unit=self, item='Sand', quantity_unit='ea'),
            )

        data = load_data(path=data_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)

        # Initialize a mixed stream
        self._mixed = WasteStream(f'{self.ID}_mixed')

    _N_outs = 3

    def _run(self):
        Toilet._run(self)
        waste, CH4, N2O = self.outs
        CH4.phase = N2O.phase = 'g'

        mixed = self._mixed
        mixed.mix_from(self.ins)
        tot_COD_kg = sum(float(getattr(i, 'COD'))*i.F_vol for i in self.ins)/1e3

        # Air emission
        if self.if_air_emission:
            # N loss due to ammonia volatilization
            NH3_rmd, NonNH3_rmd = \
                self.allocate_N_removal(mixed.TN/1e3*mixed.F_vol*self.N_volatilization,
                                        mixed.imass['NH3'])
            mixed.imass ['NH3'] -= NH3_rmd
            mixed.imass['NonNH3'] -= NonNH3_rmd
            # Energy/N loss due to degradation
            COD_loss = self.first_order_decay(k=self.decay_k_COD,
                                              t=self.collection_period/365,
                                              max_decay=self.COD_max_decay)
            COD_loss_kg = tot_COD_kg * COD_loss
            CH4.imass['CH4'] = COD_loss_kg * self.max_CH4_emission * self.MCF_decay
            mixed.imass['OtherSS'] *= 1 - COD_loss

            N_loss = self.first_order_decay(k=self.decay_k_N,
                                            t=self.collection_period/365,
                                            max_decay=self.N_max_decay)
            N_loss_tot = mixed.TN/1e3*mixed.F_vol*N_loss
            NH3_rmd, NonNH3_rmd = \
                self.allocate_N_removal(N_loss_tot,
                                        mixed.imass['NH3'])
            mixed.imass ['NH3'] -= NH3_rmd
            mixed.imass['NonNH3'] -= NonNH3_rmd
            N2O.imass['N2O'] = N_loss_tot * self.N2O_EF_decay * 44/28
            mixed._COD = (tot_COD_kg-COD_loss_kg)*1e3/mixed.F_vol
        else:
            CH4.empty()
            N2O.empty()

        # Non-ideal emptying
        if not self.if_ideal_emptying:
            mixed, CH4, N2O = self.get_emptying_emission(
                waste=mixed, CH4=CH4, N2O=N2O,
                empty_ratio=self.empty_ratio,
                CH4_factor=self.COD_max_decay*self.MCF_aq*self.max_CH4_emission,
                N2O_factor=self.N2O_EF_decay*44/28)
        waste.copy_like(mixed)

        # Scale up the effluent based on the number of user per toilet and
        # toilet number
        tot_user = self.N_user * self.N_toilet
        for i in self.outs:
            if not i.F_mass == 0:
                i.F_mass *= tot_user
    _units = {
        'Collection period': 'd',
        }

    def _design(self):
        design = self.design_results
        design['Number of users per toilet'] = self.N_user
        design['Parallel toilets'] = N = self.N_toilet
        design['Collection period'] = self.collection_period
        design['Ceramic (single MURT)'] = Ceramic_quant = \
            self.squatting_pan_weight*self.N_squatting_pan + self.urinal_weight*self.N_urinal
        design['Fan (single MURT)'] = N_fan = self.N_fan

        constr = self.construction
        constr[0].quantity = Ceramic_quant * N
        constr[1].quantity = N_fan * N
        self.add_construction()

    def _cost(self):
        C = self.baseline_purchase_costs
        N_toilet, N_squatting_pan, N_urinal = self.N_toilet
        C['Ceramic Toilets'] = \
            (self.squatting_pan_cost*N_squatting_pan+self.urinal_cost*N_urinal)*N_toilet
        C['Fan'] = self.fan_cost * N_toilet
        C['Misc. parts'] = \
            (self.led_cost+self.anticor_floor_cost+self.circuit_change_cost+self.pipe_cost)*N_toilet

    @property
    def collection_period(self):
        '''[float] Time interval between storage tank collection, [d].'''
        return self._collection_period
    @collection_period.setter
    def collection_period(self, i):
        self._collection_period = float(i)