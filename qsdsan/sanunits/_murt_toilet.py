#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems
This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>
This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

# %%


import numpy as np
from warnings import warn
from qsdsan import SanUnit, Construction, WasteStream
from ._decay import Decay
from ..utils import ospath, load_data, data_path, price_ratio
from qsdsan.sanunits._toilets import Toilet


__all__ = ('MURTToilet',)

murt_toilet_path = ospath.join(data_path, 'sanunit_data/_murt_toilet.tsv')

# %%


@price_ratio(default_price_ratio=1)
class MURTToilet(Toilet):
    '''
    Single toilet based on data from Eco-San, a subclass of :class:`~.Toilet`.

    Parameters
    ----------


    Returns
    -------
    waste : WasteStream
        Recyclable mixed excreta.
    CH4 : WasteStream
        Fugitive CH4.
    N2O : WasteStream
        Fugitive N2O.

    References
    ----------
    ..

    See Also
    --------
    :ref:`qsdsan.sanunits.Toilet <sanunits_Toilet>`

    '''
    
    # Legacy code to add checkers
    # _P_leaching = Frac_D(name='P_leaching')

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 N_user=1, N_tot_user=1, N_toilet=None, lifetime=8, if_include_front_end=True,
                 if_toilet_paper=True, if_flushing=True, if_cleansing=False,
                 if_desiccant=True, if_air_emission=True, if_ideal_emptying=True,
                 CAPEX=0, OPEX_over_CAPEX=0,
                 T=273.15 + 24, safety_factor=1, **kwargs):

        Toilet.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with, N_user=N_user, N_tot_user=N_tot_user, N_toilet=N_toilet,
                        if_toilet_paper=if_toilet_paper, if_flushing=if_flushing, if_cleansing=if_cleansing,
                        if_desiccant=if_desiccant, if_air_emission=if_air_emission, if_ideal_emptying=if_ideal_emptying,
                        CAPEX=CAPEX, OPEX_over_CAPEX=OPEX_over_CAPEX)
        self.lifetime = lifetime
        self.T = T
        self._safety_factor = safety_factor
        self.if_include_front_end = if_include_front_end

        data = load_data(path=murt_toilet_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        self._tank_V = 60 / 1e3  # m3
        for attr, value in kwargs.items():
            setattr(self, attr, value)

    _N_outs = 3

    def _run(self):
        Toilet._run(self)
        waste, CH4, N2O = self.outs
        CH4.phase = N2O.phase = 'g'

        mixed = WasteStream(f'{self.ID}_mixed')
        mixed.mix_from(self.ins)
        tot_COD_kg = sum(float(getattr(i, 'COD')) * i.F_vol for i in self.ins) / 1e3
        
        # All composite variables in mg/L
        # Leaching
        # Here COD change due to leaching not considered
        
        
        if self.if_air_emission:
            # N loss due to ammonia volatilization
            NH3_rmd, NonNH3_rmd = \
                self.allocate_N_removal(mixed.TN / 1e3 * mixed.F_vol * self.N_volatilization,
                                        mixed.imass['NH3'])
            mixed.imass['NH3'] -= NH3_rmd
            mixed.imass['NonNH3'] -= NonNH3_rmd
            # Energy/N loss due to degradation
            COD_loss = self.first_order_decay(k=self.decay_k_COD,
                                              t=self.collection_period / 365,
                                              max_decay=self.COD_max_decay)
            COD_loss_kg = tot_COD_kg * COD_loss
            CH4.imass['CH4'] = COD_loss_kg * self.max_CH4_emission * self.MCF_decay
            mixed.imass['OtherSS'] *= 1 - COD_loss

            N_loss = self.first_order_decay(k=self.decay_k_N,
                                            t=self.collection_period / 365,
                                            max_decay=self.N_max_decay)
            N_loss_tot = mixed.TN / 1e3 * mixed.F_vol * N_loss
            NH3_rmd, NonNH3_rmd = \
                self.allocate_N_removal(N_loss_tot,
                                        mixed.imass['NH3'])
            mixed.imass['NH3'] -= NH3_rmd
            mixed.imass['NonNH3'] -= NonNH3_rmd
            N2O.imass['N2O'] = N_loss_tot * self.N2O_EF_decay * 44 / 28
            mixed._COD = (tot_COD_kg - COD_loss_kg) * 1e3 / mixed.F_vol
        else:
            CH4.empty()
            N2O.empty()

        # Non-ideal emptying
        if not self.if_ideal_emptying:
            mixed, CH4, N2O = self.get_emptying_emission(
                waste=mixed, CH4=CH4, N2O=N2O,
                empty_ratio=self.empty_ratio,
                CH4_factor=self.COD_max_decay * self.MCF_aq * self.max_CH4_emission,
                N2O_factor=self.N2O_EF_decay * 44 / 28)

        waste.copy_like(mixed)
        
        # Scale up the effluent based on the number of user per toilet and
        # toilet number
        tot_user = self.N_tot_user
        for i in self.outs:
            if not i.F_mass == 0:
                i.F_mass *= tot_user

    _units = {
        'Collection period': 'd',
        }

    def _design(self):
        if self.if_include_front_end:
            design = self.design_results
            design['Number of users per toilet'] = self.N_user
            design['Parallel toilets'] = N = self.N_toilet
            design['Collection period'] = self.collection_period
            design['Ceramic'] = Ceramic_quant = (self.squatting_pan_weight +
                                                 self.urinal_weight)
            design['Fan'] = Fan_quant = 1  # assume fan is 1 kg

            self.construction = (
                Construction(item='Ceramic', quantity=Ceramic_quant * N, quantity_unit='kg'),
                Construction(item='Fan', quantity=Fan_quant * N, quantity_unit='kg'),
            )

            self.add_construction()
        else:
            self.design_results.clear()
            self.construction = ()

    def _cost(self):
        if self.if_include_front_end:
            C = self.baseline_purchase_costs
            C['Ceramic Toilets'] = (self.squatting_pan_cost + self.urinal_cost) * self.N_toilet
            C['Fan'] = self.fan_cost * self.N_toilet
            C['Misc. parts'] = (self.led_cost + self.anticor_floor_cost + self.circuit_change_cost +self.pipe_cost) * self.N_toilet
            ratio = self.price_ratio
            for equipment, cost in C.items():
                C[equipment] = cost * ratio
        else:
            self.baseline_purchase_costs.clear()

    @property
    def collection_period(self):
        '''[float] Time interval between storage tank collection, [d].'''
        return self._collection_period

    @collection_period.setter
    def collection_period(self, i):
        self._collection_period = float(i)
