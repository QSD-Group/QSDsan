#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    
    Yalin Li <mailto.yalin.li@gmail.com>
    
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

# %%

from .. import SanUnit
from ..utils import ospath, load_data, data_path
from warnings import warn
# from scipy.linalg import solve as la_solve
import numpy as np

__all__ = ('Excretion', 'ExcretionmASM2d')

excretion_path = ospath.join(data_path, 'sanunit_data/_excretion.tsv')


# %%

class Excretion(SanUnit):
    '''
    Estimation of N, P, K, and COD in urine and feces based on dietary intake
    for one person based on `Trimmer et al. <https://doi.org/10.1021/acs.est.0c03296>`_

    Parameters
    ----------
    waste_ratio : float
        A ratio in [0, 1] to indicate the amount of intake calories and nutrients
        (N, P, K) that is wasted.

    Examples
    --------
    `bwaise systems <https://github.com/QSD-Group/EXPOsan/blob/main/exposan/bwaise/systems.py>`_

    References
    ----------
    [1] Trimmer et al., Navigating Multidimensional Social–Ecological System
    Trade-Offs across Sanitation Alternatives in an Urban Informal Settlement.
    Environ. Sci. Technol. 2020, 54 (19), 12641–12653.
    https://doi.org/10.1021/acs.est.0c03296
    '''

    _N_ins = 0
    _N_outs = 2

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 waste_ratio=0, **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.waste_ratio = waste_ratio

        data = load_data(path=excretion_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            # value = float(data.loc[para]['low'])
            # value = float(data.loc[para]['high'])
            setattr(self, '_'+para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)

    def _run(self):
        ur, fec = self.outs
        ur.empty()
        fec.empty()

        not_wasted = 1 - self.waste_ratio
        factor = 24 * 1e3 # from g per person per day to kg per hour

        ur_N = (self.p_veg+self.p_anim)/factor*self.N_prot \
           * self.N_exc*self.N_ur*not_wasted
        ur.imass['NH3'] = ur_N * self.N_ur_NH3
        ur.imass['NonNH3'] = ur_N - ur.imass['NH3']

        ur.imass['P'] = (self.p_veg*self.P_prot_v+self.p_anim*self.P_prot_a)/factor \
            * self.P_exc*self.P_ur*not_wasted

        e_cal = self.e_cal / 24 * not_wasted
        ur.imass['K'] = e_cal/1e3 * self.K_cal/1e3 * self.K_exc*self.K_ur
        ur.imass['Mg'] = self.Mg_ur / factor
        ur.imass['Ca'] = self.Ca_ur / factor

        ur_exc = self.ur_exc / factor
        ur.imass['H2O'] = self.ur_moi * ur_exc
        ur.imass['OtherSS'] = ur_exc - ur.F_mass

        fec_exc = self.fec_exc / factor
        fec_N = (1-self.N_ur)/self.N_ur * ur_N
        fec.imass['NH3'] = fec_N * self.N_fec_NH3
        fec.imass['NonNH3'] = fec_N - fec.imass['NH3']
        fec.imass['P'] = (1-self.P_ur)/self.P_ur * ur.imass['P']
        fec.imass['K'] = (1-self.K_ur)/self.K_ur * ur.imass['K']
        fec.imass['Mg'] = self.Mg_fec / factor
        fec.imass['Ca'] = self.Ca_fec / factor
        fec.imass['H2O'] = self.fec_moi * fec_exc
        fec.imass['OtherSS'] = fec_exc - fec.F_mass

        # 14 kJ/g COD, the average lower heating value of excreta
        tot_COD = e_cal*self.e_exc*4.184/14/1e3 # in kg COD/hr
        ur._COD = tot_COD*(1-self.e_fec) / (ur.F_vol/1e3) # in mg/L
        fec._COD = tot_COD*self.e_fec / (fec.F_vol/1e3) # in mg/L

    @property
    def e_cal(self):
        '''[float] Caloric intake, [kcal/cap/d].'''
        return self._e_cal
    @e_cal.setter
    def e_cal(self, i):
        self._e_cal = i

    @property
    def p_veg(self):
        '''[float] Vegetal protein intake, [g/cap/d].'''
        return self._p_veg
    @p_veg.setter
    def p_veg(self, i):
        self._p_veg = i

    @property
    def p_anim(self):
        '''[float] Animal protein intake, [g/cap/d].'''
        return self._p_anim
    @p_anim.setter
    def p_anim(self, i):
        self._p_anim = i

    @property
    def N_prot(self):
        '''[float] Nitrogen content in protein, [wt%].'''
        return self._N_prot
    @N_prot.setter
    def N_prot(self, i):
        self._N_prot = i

    @property
    def P_prot_v(self):
        '''[float] Phosphorus content in vegetal protein, [wt%].'''
        return self._P_prot_v
    @P_prot_v.setter
    def P_prot_v(self, i):
        self._P_prot_v = i

    @property
    def P_prot_a(self):
        '''[float] Phosphorus content in animal protein, [wt%].'''
        return self._P_prot_a
    @P_prot_a.setter
    def P_prot_a(self, i):
        self._P_prot_a = i

    @property
    def K_cal(self):
        '''[float] Potassium intake relative to caloric intake, [g K/1000 kcal].'''
        return self._K_cal
    @K_cal.setter
    def K_cal(self, i):
        self._K_cal = i

    @property
    def N_exc(self):
        '''[float] Nitrogen excretion factor, [% of intake].'''
        return self._N_exc
    @N_exc.setter
    def N_exc(self, i):
        self._N_exc = i

    @property
    def P_exc(self):
        '''[float] Phosphorus excretion factor, [% of intake].'''
        return self._P_exc
    @P_exc.setter
    def P_exc(self, i):
        self._P_exc = i

    @property
    def K_exc(self):
        '''[float] Potassium excretion factor, [% of intake].'''
        return self._K_exc
    @K_exc.setter
    def K_exc(self, i):
        self._K_exc = i

    @property
    def e_exc(self):
        '''[float] Energy excretion factor, [% of intake].'''
        return self._e_exc
    @e_exc.setter
    def e_exc(self, i):
        self._e_exc = i

    @property
    def N_ur(self):
        '''[float] Nitrogen recovered in urine, [wt%].'''
        return self._N_ur
    @N_ur.setter
    def N_ur(self, i):
        self._N_ur = i

    @property
    def P_ur(self):
        '''[float] Phosphorus recovered in urine, [wt%].'''
        return self._P_ur
    @P_ur.setter
    def P_ur(self, i):
        self._P_ur = i

    @property
    def K_ur(self):
        '''[float] Potassium recovered in urine, [wt%].'''
        return self._K_ur
    @K_ur.setter
    def K_ur(self, i):
        self._K_ur = i

    @property
    def e_fec(self):
        '''[float] Percent of excreted energy in feces, [%].'''
        return self._e_fec
    @e_fec.setter
    def e_fec(self, i):
        self._e_fec = i

    @property
    def N_ur_NH3(self):
        '''[float] Reduced inorganic nitrogen in urine, modeled as NH3, [% of total urine N].'''
        return self._N_ur_NH3
    @N_ur_NH3.setter
    def N_ur_NH3(self, i):
        self._N_ur_NH3 = i

    @property
    def N_fec_NH3(self):
        '''[float] Reduced inorganic nitrogen in feces, modeled as NH3, [% of total feces N].'''
        return self._N_fec_NH3
    @N_fec_NH3.setter
    def N_fec_NH3(self, i):
        self._N_fec_NH3 = i

    @property
    def ur_exc(self):
        '''[float] Urine generated per day, [g/cap/d].'''
        return self._ur_exc
    @ur_exc.setter
    def ur_exc(self, i):
        self._ur_exc = i

    @property
    def fec_exc(self):
        '''[float] Feces generated per day, [g/cap/d].'''
        return self._fec_exc
    @fec_exc.setter
    def fec_exc(self, i):
        self._fec_exc = i

    @property
    def ur_moi(self):
        '''[float] Moisture (water) content of urine, [wt%].'''
        return self._ur_moi
    @ur_moi.setter
    def ur_moi(self, i):
        self._ur_moi = i

    @property
    def fec_moi(self):
        '''[float] Moisture (water) content of feces, [wt%].'''
        return self._fec_moi
    @fec_moi.setter
    def fec_moi(self, i):
        self._fec_moi = i

    @property
    def Mg_ur(self):
        '''[float] Magnesium excreted in urine, [g Mg/cap/d].'''
        return self._Mg_ur
    @Mg_ur.setter
    def Mg_ur(self, i):
        self._Mg_ur = i

    @property
    def Mg_fec(self):
        '''[float] Magnesium excreted in feces, [g Mg/cap/d].'''
        return self._Mg_fec
    @Mg_fec.setter
    def Mg_fec(self, i):
        self._Mg_fec = i

    @property
    def Ca_ur(self):
        '''[float] Calcium excreted in urine, [g Ca/cap/d].'''
        return self._Ca_ur
    @Ca_ur.setter
    def Ca_ur(self, i):
        self._Ca_ur = i

    @property
    def Ca_fec(self):
        '''[float] Calcium excreted in feces, [g Ca/cap/d].'''
        return self._Ca_fec
    @Ca_fec.setter
    def Ca_fec(self, i):
        self._Ca_fec = i

    @property
    def waste_ratio(self):
        '''
        [float] The amount of intake calories and nutrients
        (N, P, K) that is wasted.

        .. note::
            Not considered for Mg and Ca.
        '''
        return self._waste_ratio
    @waste_ratio.setter
    def waste_ratio(self, i):
        self._waste_ratio = i


#%%

class ExcretionmASM2d(Excretion):
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 waste_ratio=0, **kwargs):
        super().__init__(ID, ins, outs, thermo, init_with, waste_ratio, **kwargs)
        isdyn = kwargs.pop('isdynamic', False)
        if isdyn: self._init_dynamic()
    
    def _run(self):
        ur, fec = self.outs
        ur.empty()
        fec.empty()
        cmps = ur.components
        sf_iN = cmps.S_F.i_N
        xs_iN = cmps.X_S.i_N
        xb_iN = cmps.X_H.i_N
        sxi_iN = cmps.S_I.i_N
        i_mass = cmps.i_mass
        i_P = cmps.i_P
        hco3_imass = cmps.S_IC.i_mass

        not_wasted = 1 - self.waste_ratio
        factor = 24 * 1e3 # from g/cap/d to kg/hr(/cap)
        e_cal = self.e_cal / 24 * not_wasted # kcal/cap/d --> kcal/cap/hr
        ur_exc = self.ur_exc / factor
        fec_exc = self.fec_exc / factor

        # 14 kJ/g COD, the average lower heating value of excreta
        tot_COD = e_cal*self.e_exc*4.184/14/1e3 # in kg COD/hr
        fec_COD = tot_COD*self.e_fec
        ur_COD = tot_COD - fec_COD
        
        tot_N = (self.p_veg+self.p_anim)*self.N_prot/factor \
            * self.N_exc*not_wasted
        ur_N = tot_N*self.N_ur
        fec_N = tot_N - ur_N
        
        tot_P = (self.p_veg*self.P_prot_v+self.p_anim*self.P_prot_a)/factor \
            * self.P_exc*not_wasted
        ur_P = tot_P*self.P_ur
        fec_P = tot_P - ur_P
        
        # breakpoint()
        ur.imass['S_NH4'] = ur_nh4 = ur_N * self.N_ur_NH3
        req_sf_cod = (ur_N - ur_nh4) / sf_iN
        if req_sf_cod <= ur_COD:
            ur.imass['S_F'] = sf = req_sf_cod
            ur.imass['S_A'] = ur_COD - sf  # contains no N or P
        else:
            req_si_cod = (ur_N - ur_nh4) / sxi_iN
            if req_si_cod <= ur_COD:
                ur.imass['S_F'] = sf = (sxi_iN * ur_COD - (ur_N - ur_nh4))/(sxi_iN - sf_iN)
                ur.imass['S_I'] = ur_COD - sf
            else:
                ur.imass['S_F'] = sf = ur_COD
                ur_other_n = ur_N - ur_nh4 - sf * sf_iN
                warn(f"Excess non-NH3 nitrogen cannot be accounted for by organics "
                     f"in urine: {ur_other_n} kg/hr. Added to NH3-N.")
                ur.imass['S_NH4'] += ur_other_n # debatable, has negative COD # raise warning/error
        
        ur.imass['S_PO4'] = ur_P - sum(ur.mass * i_P)
        ur.imass['S_K'] = e_cal/1e3 * self.K_cal/1e3 * self.K_exc*self.K_ur
        ur.imass['S_Mg'] = self.Mg_ur / factor
        ur.imass['S_Ca'] = self.Ca_ur / factor

        ur.imass['H2O'] = self.ur_moi * ur_exc
        ur_others = ur_exc - sum(ur.mass * i_mass)
        ur.imass['S_IC'] = ur_others * 0.34 / hco3_imass
        ur.imass['S_Na'] = ur_others * 0.35
        ur.imass['S_Cl'] = ur_others * 0.31

        fec.imass['S_NH4'] = fec_nh4 = fec_N * self.N_fec_NH3
        req_xs_cod = (fec_N - fec_nh4) / xs_iN
        if req_xs_cod <= fec_COD:
            fec.imass['X_S'] = xs = req_xs_cod
            fec.imass['S_A'] = fec_COD - xs
        else:
            req_xi_cod = (fec_N - fec_nh4) / sxi_iN
            if req_xi_cod <= fec_COD:
                fec.imass['X_S'] = xs = (sxi_iN * fec_COD - (fec_N - fec_nh4))/(sxi_iN - xs_iN)
                fec.imass['X_I'] = fec_COD - xs
            else:
                req_xb_cod = (fec_N - fec_nh4) / xb_iN
                if req_xb_cod <= fec_COD:
                    fec.imass['X_S'] = xs = (xb_iN * fec_COD - (fec_N - fec_nh4))/(xb_iN - xs_iN)
                    fec.imass['X_H'] = fec_COD - xs
                else:
                    fec.imass['X_S'] = xs = fec_COD
                    fec_other_n = fec_N - fec_nh4 - xs * xs_iN
                    warn(f"Excess non-NH3 nitrogen cannot be accounted for by organics "
                         f"in feces: {fec_other_n} kg/hr. Added to NH3-N.")
                    fec.imass['S_NH4'] += fec_other_n # debatable, has negative COD
        
        fec.imass['S_PO4'] = fec_P - sum(fec.mass * i_P)
        fec.imass['S_K'] = (1-self.K_ur)/self.K_ur * ur.imass['S_K']
        fec.imass['S_Mg'] = self.Mg_fec / factor
        fec.imass['S_Ca'] = self.Ca_fec / factor
        fec.imass['H2O'] = self.fec_moi * fec_exc
        
        fec_others = fec_exc - sum(fec.mass * i_mass)
        fec.imass['S_IC'] = fec_others * 0.34 / hco3_imass
        fec.imass['S_Na'] = fec_others * 0.35
        fec.imass['S_Cl'] = fec_others * 0.31
        
    
    @property
    def AE(self):
        if self._AE is None:
            self._compile_AE()
        return self._AE
    
    def _compile_AE(self):
        def yt(t, QC_ins, dQC_ins):
            pass
        self._AE = yt
        
    def _init_state(self):
        ur, fec = self.outs
        self._state = np.append(ur.mass, fec.mass)
        for ws in self.outs:
            ws.state = np.append(ws.conc, ws.F_vol * 24)
            ws.dstate = np.zeros_like(ws.state)

    def _update_state(self):
        pass

    def _update_dstate(self):
        pass