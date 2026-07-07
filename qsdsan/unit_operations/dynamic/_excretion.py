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

Note: `unit_operations/static/_excretion.py` is a separate module holding the
base `Excretion` class that `ExcretionmASM2d` here subclasses.
'''

# %%

from warnings import warn
import numpy as np
from ..static._excretion import Excretion

__all__ = ('ExcretionmASM2d',)


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
