# -*- coding: utf-8 -*-
'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
        
    Joy Zhang <joycheung1994@gmail.com>

Part of this module is based on the biosteam package:
https://github.com/BioSTEAMDevelopmentGroup/biosteam

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

import numpy as np
from warnings import warn
from math import exp
from .. import SanUnit
import flexsolve as flx


__all__ = ('MetalDosage')

    
def _coagulation_mass_balance(S_out, S_in, Pme, Fmax, Fmin, ka):
    '''
    Coagulation mass balance at equilibrium

    Parameters
    ----------
    S_out : float
        Soluble component concentration at equilibrium.
    S_in : float
        Influent soluble component concentration.
    Pme : float
        Metal dosage.
    Fmax : float
        Maximum required metal dosage for unit removal.
    Fmin : float
        Minimum required metal dosage for unit removal.
    ka : float
        Affinity factor.

    '''
    return Fmax/ka * (exp(-ka*S_out) - exp(-ka*S_in)) + Fmin*(S_in - S_out) - Pme

def _precipitation_mass_balance(SP, Me_in, SP_in, Ksp_mass, x, y, i, j, alpha):
    '''
    Precipitation mass balance at equilibrium.

    Parameters
    ----------
    SP : float
        Ortho-phosphate concentration at equilibrium.
    Me_in : float
        Influent metal concentration.
    SP_in : float
        Influent ortho-P concentration.
    Ksp_mass : float
        Metal phosphate solubility product, based on mass concentration of measured substances.
    x : int
        Molar stoichiometric coefficient for metal to precipitate one unit of metal phosphate.
    y : int
        Molar stoichiometric coefficient for ortho-P.
    i : float
        Mass based stoichiometric coefficient for metal.
    j : float
        Mass based stoichiometric coefficient for ortho-P.
    alpha : float
        Dissociation factor of PO4(3-).

    '''
    Me = Me_in + i/j*(SP-SP_in)                 # mass-based stoichiometry
    return  Me**x * (SP*alpha)**y - Ksp_mass    # solubility product

class MetalDosage(SanUnit):
    
    _N_ins = 1
    _N_outs = 1
    _ins_size_is_fixed = False
    _outs_size_is_fixed = True

    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', F_BM_default=None, isdynamic=False,
                 P_precipitation_stoichiometry=None, metal_dosage=0.,
                 metal_ID='', orthoP_ID='', metalP_ID='', 
                 metalP_pKsp=13.75, alpha=1.8e-6,
                 soluble_substrate_ID='S_I', particulate_substrate_ID='X_I',
                 soluble_inert_ID='S_F', particulate_inert_ID='X_S', 
                 Fmin_substrate=4, Fmax_substrate=20, ka_substrate=0.5,
                 Fmin_inert=4, Fmax_inert=20, ka_inert=0.5,        
                 ):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with,
                         F_BM_default=F_BM_default, isdynamic=isdynamic)
        self.metal_dosage = metal_dosage
        self.metal_ID = metal_ID
        self.orthoP_ID = orthoP_ID
        self.metalP_ID = metalP_ID
        self.P_precipitation_stoichiometry = P_precipitation_stoichiometry
        self.metalP_pKsp = metalP_pKsp
        self.alpha = alpha
        self.soluble_substrate_ID = soluble_substrate_ID
        self.particulate_substrate_ID = particulate_substrate_ID
        self._check_composition(soluble_substrate_ID, particulate_substrate_ID)
        self.soluble_inert_ID = soluble_inert_ID
        self.particulate_inert_ID = particulate_inert_ID
        self._check_composition(soluble_inert_ID, particulate_inert_ID)
        self.Fmin_substrate = Fmin_substrate
        self.Fmax_substrate = Fmax_substrate
        self.ka_substrate = ka_substrate
        self.Fmin_inert = Fmin_inert
        self.Fmax_inert = Fmax_inert
        self.ka_inert = ka_inert
        self._cmps_idx = self.chemicals.indices([metal_ID, orthoP_ID, metalP_ID])
        self._org_idx = self.chemicals.indices([
            soluble_substrate_ID, particulate_substrate_ID, 
            soluble_inert_ID, particulate_inert_ID])

    @property
    def metal_dosage(self):
        '''[float] Metal dosage, in g metal component / m3'''
        return self._Pme
    @metal_dosage.setter
    def metal_dosage(self, P):
        self._Pme = P  

    @property
    def P_precipitation_stoichiometry(self):
        '''[dict] Metal phsophate precipitation process stoichiometry (molar).'''
        return self._stoichio
    @P_precipitation_stoichiometry.setter
    def P_precipitation_stoichiometry(self, stoichio):
        metalP = stoichio[self.metalP_ID]
        if metalP != 1:
            stoichio = {k: v/metalP for k,v in stoichio.items()}
        self._stoichio = stoichio
        self._Me_stoi = -stoichio[self.metal_ID]
        self._SP_stoi = -stoichio[self.orthoP_ID]
        cmps = self.chemicals
        unit_conv = cmps.chem_MW / cmps.i_mass  # convert from mol to mass
        stoichio_arr = cmps.kwarray(stoichio) * unit_conv # mass-based stoichiometry
        self._mstoichio = stoichio_arr[self._cmps_idx]

    @property
    def metalP_pKsp(self):
        '''[float] Apparent solubility product (molar) of metal phosphate.'''
        return self._pKsp
    @metalP_pKsp.setter
    def metalP_pKsp(self, pKsp):
        self._pKsp = pKsp
        cmps = self.chemicals
        unit_conv = cmps.chem_MW / cmps.i_mass * 1e3  # convert from M to mg (measured_as) / L
        i_Me, i_SP, i_MeP = unit_conv[self._cmps_idx]
        self._Ksp_mass = 10**(-pKsp) * (i_Me ** self._Me_stoi) * (i_SP ** self._SP_stoi)
    
    @property
    def alpha(self):
        '''[float] Dissociation factor for estimating the precipitating phosphoric
        species concentrationn from total ortho-phosphate concentration.'''
        return self._alpha
    @alpha.setter
    def alpha(self, a):
        self._alpha = a

    @property
    def Fmin_substrate(self):
        '''[float] Minimum (i.e., at high soluble substrate concentration) 
        required metal dose for unit soluble substrate coagulation, 
        in g metal / g soluble substrate component. '''
        return self._Fmin_ss
    @Fmin_substrate.setter
    def Fmin_substrate(self, F):
        self._Fmin_ss = F
        
    @property
    def Fmax_substrate(self):
        '''[float] Maximum (i.e., at low soluble substrate concentration) 
        required metal dose for unit soluble substrate coagulation, 
        in g metal / g soluble substrate component. '''
        return self._Fmax_ss
    @Fmax_substrate.setter
    def Fmax_substrate(self, F):
        self._Fmax_ss = F

    @property
    def ka_substrate(self):
        '''Soluble substrate affinity factor for coagulation, in m3/g.'''
        return self._ka_ss
    @ka_substrate.setter
    def ka_substrate(self, ka):
        self._ka_ss = ka

    @property
    def Fmin_inert(self):
        '''[float] Minimum (i.e., at high soluble inert concentration) 
        required metal dose for unit soluble inert coagulation, 
        in g metal / g soluble inert component. '''
        return self._Fmin_si
    @Fmin_inert.setter
    def Fmin_inert(self, F):
        self._Fmin_si = F
        
    @property
    def Fmax_inert(self):
        '''[float] Maximum (i.e., at low soluble inert concentration) 
        required metal dose for unit soluble inert coagulation, 
        in g metal / g soluble inert component. '''
        return self._Fmax_si
    @Fmax_inert.setter
    def Fmax_inert(self, F):
        self._Fmax_si = F

    @property
    def ka_inert(self):
        '''Soluble inert affinity factor for coagulation, in m3/g.'''
        return self._ka_si
    @ka_inert.setter
    def ka_inert(self, ka):
        self._ka_si = ka
        
    def _check_composition(self, soluble, particulate):
        get = getattr
        cmps = self.chemicals
        s = get(cmps, soluble)
        x = get(cmps, particulate)
        for attr in ('i_COD', 'i_C', 'i_N', 'i_P', 'i_mass'):
            if get(s, attr) != get(x, attr):
                warn(f'ignored unequal {attr} between {soluble} and {particulate} in coagulation')

    
    def _run(self):
        out, = self.outs
        out.mix_from(self.ins)
        Q = out.F_vol
        out.imass[self.metal_ID] += self._Pme * Q * 1e-3
        Me_in, SP_in, MeP_in = out.conc[self._cmps_idx] # mg/L
        i,j,k = self._mstoichio
        SP = flx.IQ_interpolation(_precipitation_mass_balance, 0, SP_in, args=(
            Me_in, SP_in, self._Ksp_mass, self._Me_stoi, self._SP_stoi, i, j, self.alpha
            ))
        Me = Me_in + i/j*(SP-SP_in)
        MeP = MeP_in + k/j*(SP-SP_in)
        SS_in, XS_in, SI_in, XI_in = out.conc[self._org_idx]
        SS = flx.IQ_interpolation(_coagulation_mass_balance, 0, SS_in, args=(
            SS_in, self._Pme, self._Fmax_ss, self._Fmin_ss, self._ka_ss
            ))
        XS = XS_in + (SS_in - SS)
        SI = flx.IQ_interpolation(_coagulation_mass_balance, 0, SI_in, args=(
            SI_in, self._Pme, self._Fmax_si, self._Fmin_si, self._ka_si            
            ))
        XI = XI_in + (SI_in - SI)
        out.mass[[*self._cmps_idx, *self._org_idx]] = \
            [c*Q*1e-3 for c in [Me, SP, MeP, SS, XS, SI, XI]]

    def _init_state(self):
        out, = self.outs
        self._state = np.append(out.conc, out.F_vol*24)
        self._dstate = self._state * 0.

    def _update_state(self):
        self._outs[0].state = self._state

    def _update_dstate(self):
        self._outs[0].dstate = self._dstate

    @property
    def AE(self):
        if self._AE is None:
            self._compile_AE()
        return self._AE

    def _compile_AE(self):
        _state = self._state
        _update_state = self._update_state
        chem_idx = self._cmps_idx
        org_idx = self._org_idx
        Ksp_mass = self._Ksp_mass
        x = self._Me_stoi
        y = self._SP_stoi
        i,j,k = self._mstoichio
        alpha = self.alpha
        Pme = self._Pme
        Fmax_ss = self._Fmax_ss
        Fmin_ss = self._Fmin_ss
        ka_ss = self._ka_ss
        Fmax_si = self._Fmax_si
        Fmin_si = self._Fmin_si
        ka_si = self._ka_si

        def solve_sp(Me_in, SP_in):
            sp = flx.IQ_interpolation(
                _precipitation_mass_balance, 0, SP_in, args=(
                Me_in, SP_in, Ksp_mass, x, y, i, j, alpha
                ))
            return sp
        
        def solve_ss(SS_in):
            ss = flx.IQ_interpolation(
                _coagulation_mass_balance, 0, SS_in, args=(
                SS_in, Pme, Fmax_ss, Fmin_ss, ka_ss
                ))
            return ss
        
        def solve_si(SI_in):
            si = flx.IQ_interpolation(
                _coagulation_mass_balance, 0, SI_in, args=(
                SI_in, Pme, Fmax_si, Fmin_si, ka_si
                ))
            return si
            
        def yt(t, QC_ins, dQC_ins):
            Q_ins = QC_ins[:, -1]
            C_ins = QC_ins[:, :-1]
            Q = Q_ins.sum()
            C = Q_ins @ C_ins / Q
            _state[-1] = Q
            _state[:-1] = C
            Me_in, SP_in, MeP_in = C[chem_idx]
            Me_in += Pme
            SP = solve_sp(Me_in, SP_in)
            Me = Me_in + i/j*(SP-SP_in)
            MeP = MeP_in + k/j*(SP-SP_in)
            SS_in, XS_in, SI_in, XI_in = C[org_idx]
            SS = solve_ss(SS_in)
            XS = XS_in + (SS_in - SS)
            SI = solve_si(SI_in)
            XI = XI_in + (SI_in - SI)
            _state[[*chem_idx, *org_idx]] = [Me, SP, MeP, SS, XS, SI, XI]
            _update_state()

        self._AE = yt