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

import numpy as np
from .. import Construction
from ._toilet import Toilet
from ..utils.loading import load_data, data_path

__all__ = ('UDDT',)

data_path += 'sanunit_data/_uddt.csv'


# %%

class UDDT(Toilet):
    '''
    Urine-diverting dry toilet with liquid storage tank and dehydration vault 
    for urine and feces storage, respectively, based on Trimmer et al. [1]_,
    a subclass of qsdsan.sanunits.Toilet.
    
    Parameters
    ----------
    T : float
        Temperature, [K].
    safety_factor : float
        Safety factor for pathogen removal during onsite treatment,
        must be larger than 1.            
    if_treatment : bool
        If has onsite treatment.
    if_pit_above_water_table : bool
        If the pit is above local water table.
    
    Returns
    -------
    liq : WasteStream
        Recyclable liquid urine.
    sol : WasteStream
        Recyclable solid feces.
    struvite : WasteStream
        Struvite scaling (irrecoverable).
    HAP : WasteStream
        Hydroxyapatite scaling (irrecoverable).
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
    
    def __init__(self, ID='', ins=None, outs=(), N_user=1, N_toilet=1, lifetime=8,
                 if_toilet_paper=True, if_flushing=True, if_cleansing=False,
                 if_desiccant=True, if_air_emission=True, if_ideal_emptying=True,
                 CAPEX=553, OPEX_over_CAPEX=0.1,
                 T=273.15+24, safety_factor=1, if_prep_loss=True, if_treatment=False,
                 **kwargs):

        Toilet.__init__(self, ID, ins, outs, N_user, N_toilet,
                        if_toilet_paper, if_flushing, if_cleansing, if_desiccant,
                        if_air_emission, if_ideal_emptying, CAPEX, OPEX_over_CAPEX)
        self.lifetime = lifetime
        self.T = T
        self._safety_factor = safety_factor
        self.if_prep_loss = if_prep_loss
        self.if_treatment = if_treatment

        data = load_data(path=data_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, '_'+para, value)
        del data

        self._tank_V = 60/1e3 # m3
        for attr, value in kwargs.items():
            setattr(self, attr, value)
    
    _N_outs = 6

    def _run(self):
        Toilet._run(self)
        liq, sol, struvite, HAP, CH4, N2O = self.outs
        liq.copy_like(self.ins[0])
        # Assume all of additives (toilet paper and desiccant) and water 
        # (cleansing and flushing) go to feces
        sol_COD = self.ins[1]._COD*self.ins[1].F_vol/1e3
        sol.mix_from(self.ins[1:])
        sol._COD = sol_COD*1e3/sol.F_vol
        struvite.phase = HAP.phase = 's'
        CH4.phase = N2O.phase = 'g'
        
        #!!! Modified from ref [1], assume this only happens when air emission occurs
        if self.if_air_emission:
            # N loss due to ammonia volatilization
            NH3_rmd, NonNH3_rmd = \
                self.allocate_N_removal(liq.TN/1e3*liq.F_vol*self.N_volatilization,
                                        liq.imass['NH3'])
            liq.imass ['NH3'] -= NH3_rmd
            liq.imass['NonNH3'] -= NonNH3_rmd
            # Energy/N loss due to degradation
            COD_loss = self.first_order_decay(k=self.decay_k_COD,
                                              t=self.collection_period/365,
                                              max_decay=self.COD_max_decay)
            CH4.imass['CH4'] = sol.COD/1e3*sol.F_vol*COD_loss * \
                self.max_CH4_emission*self.MCF_decay # COD in mg/L (g/m3)
            sol._COD *= 1 - COD_loss
            sol.imass['OtherSS'] *= 1 - COD_loss

            N_loss = self.first_order_decay(k=self.decay_k_N,
                                            t=self.collection_period/365,
                                            max_decay=self.N_max_decay)
            N_loss_tot = sol.TN/1e3*sol.F_vol*N_loss
            NH3_rmd, NonNH3_rmd = \
                self.allocate_N_removal(N_loss_tot, sol.imass['NH3'])
            sol.imass ['NH3'] -= NH3_rmd
            sol.imass['NonNH3'] -= NonNH3_rmd
            N2O.imass['N2O'] = N_loss_tot * self.N2O_EF_decay * 44/28
        else:
            CH4.empty()
            N2O.empty()
            
        # N and P losses due to struvite and hydroxyapatite (HAp)
        if self.if_prep_loss:
            # Struvite
            NH3_mol = liq.imol['NH3']
            P_mol = liq.imol['P']
            Mg_mol = liq.imol['Mg']
            Ksp = 10**(-self.struvite_pKsp)
            # Ksp = (initial N - struvite)(initial P - struvite)(initial Mg - struvite)
            coeff = [1, -(NH3_mol+P_mol+Mg_mol),
                     (NH3_mol*P_mol + P_mol*Mg_mol + Mg_mol*NH3_mol),
                     (Ksp - NH3_mol*P_mol*Mg_mol)]
            struvite_mol = 0
            for i in np.roots(coeff):
                if 0 < i < min(NH3_mol, P_mol, Mg_mol):
                    struvite_mol = i
            struvite.imol['Struvite'] = \
                max(0, min(NH3_mol, P_mol, Mg_mol, struvite_mol))
            liq.imol['NH3'] -= struvite_mol
            liq.imol['P'] -= struvite_mol
            liq.imol['Mg'] -= struvite_mol
            # HAP
            left_P = liq.imol['P'] - 3*(liq.imol['Ca']/5)
            # Remaining P enough to precipitate all Ca as HAP
            if left_P > 0:
                HAP.imol['HAP'] = liq.imol['Ca']/5
                liq.imol['P'] = left_P
                liq.imol['Ca'] = 0
            else:
                HAP.imol['HAP'] = liq.imol['P']/3
                liq.imol['Ca'] -= 5*(liq.imol['P']/3)
                liq.imol['P'] = 0
        else:
            struvite.empty()
            HAP.empty()
        
        # Onsite treatment
        if self.if_treatment:
            NH3_mmol = liq.imol['NH3'] * 1e3
            ur_DM = 1 - liq.imass['H2O']/liq.F_mass
            pKa = 0.09018 + (2729.92/self.T)
            f_NH3_Emerson = 1 / (10**(pKa-self.ur_pH)+1)
            alpha = 0.82 - 0.011*np.sqrt(NH3_mmol+1700*ur_DM)
            beta = 1.17 + 0.02 * np.sqrt(NH3_mmol+1100*ur_DM)
            f_NH3_Pitzer = f_NH3_Emerson * \
                (alpha + ((1-alpha)*(f_NH3_Emerson**beta)))
            NH3_conc = NH3_mmol * f_NH3_Pitzer

            #!!! Shouldn't the collection period and design be affected by this as well?            
            # Time (in days) to reach desired inactivation level
            self.treatment_tau = ((3.2 + self.log_removal) \
                             / (10**(-3.7+0.062*(self.T-273.15)) * (NH3_conc**0.7))) \
                        * 1.14*self.safety_factor
            # Total volume in m3
            self.treatment_V = self.treatment_tau * liq.F_vol*24
        else:
            self.treatment_tau = self.treatment_V = 0

        # Feces water loss if desiccant is added
        if self.if_desiccant:
            sol_COD = sol._COD*sol.F_vol/1e3
            MC_min = self.fec_moi_min
            r = self.fec_moi_red_rate
            t = self.collection_period
            mixed = sol.copy()
            mixed.mix_from(self.ins[1:])
            fec_moi_int = mixed.imass['H2O']/mixed.F_mass
            fec_moi = MC_min + (fec_moi_int-MC_min)/(r*t)*(1-np.exp(-r*t))
            sol.imass['H2O'] = (sol.F_mass-sol.imass['H2O']) / (1-fec_moi)
            sol._COD = sol_COD*1e3/sol.F_vol
        
        #!!! Maybe don't need this, only add this to the design_results dict,
        # remove setter if shouldn't be set
        self.vault_V = sol.F_vol*self.collection_period*24 # in day

        # Non-ideal emptying
        if not self.if_ideal_emptying:
            liq, CH4, N2O = self.get_emptying_emission(
                waste=liq, CH4=CH4, N2O=N2O,
                empty_ratio=self.empty_ratio,
                CH4_factor=self.COD_max_decay*self.MCF_aq*self.max_CH4_emission,
                N2O_factor=self.N2O_EF_decay*44/28)
            sol, CH4, N2O = self.get_emptying_emission(
                waste=sol, CH4=CH4, N2O=N2O,
                empty_ratio=self.empty_ratio,
                CH4_factor=self.COD_max_decay*self.MCF_aq*self.max_CH4_emission,
                N2O_factor=self.N2O_EF_decay*44/28)

        # Scale up the effluent based on the number of user per toilet and
        # toilet number
        tot_user = self.N_user * self.N_toilet
        for i in self.outs:
            if not i.F_mass == 0:
                i.F_mass *= tot_user

    _units = {
        'Collection period': 'd',
        'Single tank volume': 'm3',
        'Single vault volume': 'm3',
        'Treatment time': 'd',
        'Treatment volume': 'm3'
        }

    def _design(self):
        design = self.design_results
        design['Number of users per toilet'] = self.N_user
        design['Parallel toilets'] = N = self.N_toilet
        design['Collection period'] = self.collection_period
        design['Single tank volume'] = self.tank_V
        design['Single vault volume'] = self.vault_V
        design['Treatment time'] = self.treatment_tau
        design['Treatment volume'] = self.treatment_V
        
        density = self.density_dct
        self.construction = (
            Construction(item='Cement', quantity=200*N, unit='kg'),
            Construction(item='Sand', quantity=0.6*density['Sand']*N, unit='kg'),
            Construction(item='Gravel', quantity=0.2*density['Gravel']*N, unit='kg'),
            Construction(item='Brick', quantity=682*0.0024*density['Brick']*N, unit='kg'),
            Construction(item='Plastic', quantity=4*density['Plastic']*N, unit='kg'),
            Construction(item='Steel', quantity=0.00351*density['Steel']*N, unit='kg'),
            Construction(item='StainlessSteelSheet',
                         quantity=28.05*density['StainlessSteelSheet']*N, unit='kg'),
            Construction(item='Wood', quantity=0.222*N, unit='m3'),
            )

        self.add_construction(add_cost=False)


    @property
    def safety_factor(self):
        return self._safety_factor
    @safety_factor.setter
    def safety_factor(self, i):
        if i < 1:
            raise ValueError(f'safety_factor must be larger than 1, not {i}.')
        self._safety_factor = float(i)

    @property
    def collection_period(self):
        '''[float] Time interval between storage tank collection, [d].'''
        return self._collection_period
    @collection_period.setter
    def collection_period(self, i):
        self._collection_period = float(i)

    @property
    def treatment_tau(self):
        '''[float] Time for onsite treatment (if treating), [d].'''
        return self._treatment_tau
    @treatment_tau.setter
    def treatment_tau(self, i):
        self._treatment_tau = float(i)

    @property
    def treatment_V(self):
        '''[float] Volume needed to achieve treatment target (if treating), [d].'''
        return self._treatment_V
    @treatment_V.setter
    def treatment_V(self, i):
        self._treatment_V = float(i)

    @property
    def tank_V(self):
        '''[float] Volume of the urine storage tank, [m3].'''
        return self._tank_V
    @tank_V.setter
    def tank_V(self, i):
        self._tank_V = float(i)
        
    @property
    def vault_V(self):
        '''[float] Volume of the feces dehydration vault, [m3].'''
        return self._vault_V
    @vault_V.setter
    def vault_V(self, i):
        self._vault_V = float(i)

    @property
    def struvite_pKsp(self):
        '''[float] Precipitation constant of struvite.'''
        return self._struvite_pKsp
    @struvite_pKsp.setter
    def struvite_pKsp(self, i):
        self._struvite_pKsp = float(i)

    @property
    def prep_sludge(self):
        '''
        [float] Fraction of total precipitate appearing as sludge that can
        settle and be removed.
        '''
        return self._prep_sludge
    @prep_sludge.setter
    def prep_sludge(self, i):
        self._prep_sludge = float(i)

    @property
    def log_removal(self):
        '''Desired level of pathogen inactivation.'''
        return self._log_removal
    @log_removal.setter
    def log_removal(self, i):
        self._log_removal = float(i)

    @property
    def ur_pH(self):
        '''Urine pH.'''
        return self._ur_pH
    @ur_pH.setter
    def ur_pH(self, i):
        self._ur_pH = float(i)

    @property
    def fec_moi_min(self):
        '''[float] Minimum moisture content of feces.'''
        return self._fec_moi_min
    @fec_moi_min.setter
    def fec_moi_min(self, i):
        self._fec_moi_min = float(i)

    @property
    def fec_moi_red_rate(self):
        '''[float] Exponential reduction rate of feces moisture.'''
        return self._fec_moi_red_rate
    @fec_moi_red_rate.setter
    def fec_moi_red_rate(self, i):
        self._fec_moi_red_rate = float(i)




























