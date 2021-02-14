#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/master/LICENSE.txt
for license details.
'''


# %%

from .. import WasteStream, Construction
from ._toilet import Toilet
from ..utils.loading import load_data, data_path

# Note that here two different methods are used to check the user-input values
from ..utils.checkers import Fraction as Frac_C
from ..utils.descriptors import Fraction as Frac_D

__all__ = ('PitLatrine',)

data_path += 'sanunit_data/_pit_latrine.csv'


# %%

class PitLatrine(Toilet):
    '''
    Single pit latrine based on Trimmer et al. [1]_, a subclass of :class:`~.Toilet`.
    
    Parameters
    ----------
    if_leaching : bool
        If infiltration to soil occurs
        (i.e., if the pit walls and floors are permeable).
    if_pit_above_water_table : bool
        If the pit is above local water table.
        
    Returns
    -------
    waste : WasteStream
        Recyclable mixed excreta.
    leachate : WasteStream
        Leached to soil.
    CH4 : WasteStream
        Fugitive CH4.
    N2O : WasteStream
        Fugitive N2O.
    
    References
    ----------
    .. [1] Trimmer et al., Navigating Multidimensional Social–Ecological System
        Trade-Offs across Sanitation Alternatives in an Urban Informal Settlement.
        Environ. Sci. Technol. 2020, 54 (19), 12641–12653.
        https://doi.org/10.1021/acs.est.0c03296.
    
    See Also
    --------
    :ref:`qsdsan.sanunits.Toilet <sanunits_Toilet>`
    
    '''

    _P_leaching = Frac_D(name='P_leaching')

    def __init__(self, ID='', ins=None, outs=(), N_user=1, N_toilet=1, lifetime=8,
                 if_toilet_paper=True, if_flushing=True, if_cleansing=False,
                 if_desiccant=False, if_air_emission=True, if_ideal_emptying=True, 
                 CAPEX=449, OPEX_over_CAPEX=0.05,
                 if_leaching=True, if_shared=True,
                 if_pit_above_water_table=True, **kwargs):

        Toilet.__init__(self, ID, ins, outs, N_user, N_toilet,
                        if_toilet_paper, if_flushing, if_cleansing, if_desiccant,
                        if_air_emission, if_ideal_emptying, CAPEX, OPEX_over_CAPEX)
        self.lifetime = lifetime
        self.if_leaching = if_leaching
        self.if_pit_above_water_table = if_pit_above_water_table
        self.if_shared = if_shared
        data = load_data(path=data_path)
        for para in data.index:
            if para in ('MCF_decay', 'N2O_EF_decay'):
                value = eval(data.loc[para]['expected'])
            else:
                value = float(data.loc[para]['expected'])
            setattr(self, '_'+para, value)
        del data
        self._pit_depth = 4.57 # m
        self._pit_area = 0.8 # m2
        self._liq_leaching = None
        for attr, value in kwargs.items():
            setattr(self, attr, value)


    _N_outs = 4
        
    def _run(self):
        Toilet._run(self)
        waste, leachate, CH4, N2O = self.outs
        CH4.phase = N2O.phase = 'g'

        mixed = WasteStream()
        mixed.mix_from(self.ins)
        tot_COD_kg = sum(float(getattr(i, 'COD'))*i.F_vol for i in self.ins)/1e3
        
        # All composite variables in mg/L
        # Leaching
        # Here COD change due to leaching not considered
        if self.if_leaching:
            # Additional assumption not in ref [1]
            leachate.imass['H2O'] = mixed.imass['H2O'] * self.liq_leaching
            leachate.imass['NH3'], leachate.imass['NonNH3'] = \
                self.allocate_N_removal(mixed.TN/1e3*mixed.F_vol*self.N_leaching,
                                        mixed.imass['NH3'])
            leachate.imass['P'] = mixed.imass['P'] * self.P_leaching
            leachate.imass['K'] = mixed.imass['K'] * self.K_leaching
            mixed.mass -= leachate.mass
        
        # Air emission
        #!!! Based on the logic, COD won't degrade without air emission?
        if self.if_air_emission:
            # N loss due to ammonia volatilization
            NH3_rmd, NonNH3_rmd = \
                self.allocate_N_removal(mixed.TN/1e3*mixed.F_vol*self.N_volatilization,
                                        mixed.imass['NH3'])
            mixed.imass ['NH3'] -= NH3_rmd
            mixed.imass['NonNH3'] -= NonNH3_rmd
            # Energy/N loss due to degradation
            COD_loss = self.first_order_decay(k=self.decay_k_COD,
                                              t=self.emptying_period,
                                              max_decay=self.COD_max_decay)
            COD_loss_kg = tot_COD_kg * COD_loss
            CH4.imass['CH4'] = COD_loss_kg * self.max_CH4_emission * self.MCF_decay
            mixed.imass['OtherSS'] *= 1 - COD_loss

            N_loss = self.first_order_decay(k=self.decay_k_N,
                                            t=self.emptying_period,
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

        # Aquatic emission when not ideally emptied        
        if not self.if_ideal_emptying:
            mixed, CH4, N2O = self.get_emptying_emission(
                waste=mixed, CH4=CH4, N2O=N2O,
                empty_ratio=self.empty_ratio,
                CH4_factor=self.COD_max_decay*self.MCF_aq*self.max_CH4_emission,
                N2O_factor=self.N2O_EF_decay*44/28)

        # Drain extra water, assume density of water (1 kg/L)
        sludge = self.sludge_accum_rate/(365*24)
        diff = mixed.F_mass - sludge
        if diff > 0:
            mixed_COD = mixed._COD * mixed.F_vol
            mixed.imass['H2O'] -= diff
            mixed.imass['H2O'] = max(0, mixed.imass['H2O'])
            mixed._COD = mixed_COD / mixed.F_vol
        
        waste.copy_like(mixed)
        
        # Scale up the effluent based on the number of user per toilet and
        # toilet number
        tot_user = self.N_user * self.N_toilet
        for i in self.outs:
            if not i.F_mass == 0:
                i.F_mass *= tot_user

    _units = {
        'Emptying period': 'yr',
        'Single pit volume': 'm3',
        'Single pit area': 'm2',
        'Singple pit depth': 'm'
        }

    def _design(self):
        design = self.design_results
        design['Number of users per toilet'] = self.N_user
        design['Parallel toilets'] = N = self.N_toilet
        design['Emptying period'] = self.emptying_period
        design['Single pit volume'] = self.pit_V
        design['Single pit area'] = self.pit_area
        design['Single pit depth'] = self.pit_depth
        
        density = self.density_dct
        self.construction = (
            Construction(item='Cement', quantity=700*N, unit='kg'),
            Construction(item='Sand', quantity=2.2*density['Sand']*N, unit='kg'),
            Construction(item='Gravel', quantity=0.8*density['Gravel']*N, unit='kg'),
            Construction(item='Brick', quantity=54*0.0024*density['Brick']*N, unit='kg'),
            Construction(item='Plastic', quantity=16*density['Plastic']*N, unit='kg'),
            Construction(item='Steel', quantity=0.00425*density['Steel']*N, unit='kg'),
            Construction(item='Wood', quantity=0.19*N, unit='m3'),
            Construction(item='Excavation', quantity=self.pit_V*N, unit='m3'),
            )

        self.add_construction(add_cost=False)


    @property
    def pit_depth(self):
        '''[float] Depth of the pit, [m].'''
        return self._pit_depth
    @pit_depth.setter
    def pit_depth(self, i):
        self._pit_depth = float(i)
        
    @property
    def pit_area(self):
        '''[float] Area of the pit, [m2].'''
        return self._pit_area
    @pit_area.setter
    def pit_area(self, i):
        self._pit_area = float(i)
    
    #!!! Should add some contraints, maybe in _run, to make sure the pit is big
    # enough for the amount of excreta
    @property
    def pit_V(self):
        '''[float] Volume of the pit, [m3].'''
        return self.pit_area*self.pit_depth     

    @property
    def emptying_period(self):
        '''[float] Time interval between pit emptying, [yr].'''
        return self._emptying_period
    @emptying_period.setter
    def emptying_period(self, i):
        self._emptying_period = float(i)

    @property
    def sludge_accum_rate(self):
        '''[float] Sludge accumulation rate, [L/cap/yr].'''
        return self._sludge_accum_rate
    @sludge_accum_rate.setter
    def sludge_accum_rate(self, i):
        self._sludge_accum_rate = float(i)

    @property
    def liq_leaching(self):
        '''
        [float] Fraction of input water that leaches to the soil
        (if if_leaching is True). If not set, then return the maximum of
        fraction of N, P, K leaching
        '''
        return self._liq_leaching or \
            max(self.N_leaching, self.P_leaching, self.K_leaching)
    @liq_leaching.setter
    def liq_leaching(self, i):
        self._liq_leaching = float(i)

    @property
    def N_leaching(self):
        '''
        [float] Fraction of input N that leaches to the soil
        (if if_leaching is True).
        '''
        return self._N_leaching
    @N_leaching.setter
    def N_leaching(self, i):
        @Frac_C(self)
        def N_leaching(): return i

    @property
    def P_leaching(self):
        '''
        [float] Fraction of input P that leaches to the soil
        (if if_leaching is True).
        '''
        return self._P_leaching
    @P_leaching.setter
    def P_leaching(self, i):
        self._P_leaching = float(i)

    @property
    def K_leaching(self):
        '''
        [float] Fraction of input K that leaches to the soil
        (if if_leaching is True).
        '''
        return self._K_leaching
    @K_leaching.setter
    # This is faster than using descriptors or decorators, but (I think) less elegant
    def K_leaching(self, i):
        if i < 0:
            raise ValueError('Value for K_leaching cannot be negative')
        self._K_leaching = float(i)

    def _return_MCF_EF(self):
        # self._MCF and self._N2O_EF are dict for
        # single_above_water, communal_above_water, below_water
        if self.if_pit_above_water_table:
            if not self.if_shared:
                return 'single_above_water'
            else:
                return 'communal_above_water'
        else:
            return 'below_water'

    @property
    def MCF_decay(self):
        '''[float] Methane correction factor for COD degraded during storage.'''
        return float(self._MCF_decay[self._return_MCF_EF()])
    @MCF_decay.setter
    def MCF_decay(self, i):
        self._MCF_decay[self._return_MCF_EF()] = float(i)

    @property
    def N2O_EF_decay(self):
        '''[float] Fraction of N emitted as N2O during storage.'''
        return float(self._N2O_EF_decay[self._return_MCF_EF()])
    @N2O_EF_decay.setter
    def N2O_EF_decay(self, i):
        self._N2O_EF_decay[self._return_MCF_EF()] = float(i)





