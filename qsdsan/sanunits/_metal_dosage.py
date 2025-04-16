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
from ..processes import mASM2d, ion_speciation
import flexsolve as flx
from chemicals.elements import molecular_weight as mw


__all__ = ('MetalDosage',)

    
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
    '''
    In-line metal dosage for chemical phosphorus removal and soluble organic
    removal through coagulation.

    Parameters
    ----------
    metal_dosage : float
        Concentration-based metal dosage, in g (metal component) / m3.
    metal_ID : str
        Component ID for the dosed metal.
    orthoP_ID : str
        Component ID for ortho-phosphate.
    metalP_ID : str, optional
        Component ID for metal phosphate.
    soluble_substrate_ID : str, optional
        Component ID for soluble organic substate. The default is 'S_F'.
    particulate_substrate_ID : str, optional
        Component ID for particulate organic substrate. The default is 'X_S'.
    soluble_inert_ID : str, optional
        Component ID for soluble inert organic matter. The default is 'S_I'.
    particulate_inert_ID : str, optional
        Component ID for particulate inert organic matter. The default is 'X_I'.
    metal_price : float, optional
        Price of the dosed metal, in USD/kg metal component.    
    '''    
    _N_ins = 1
    _N_outs = 1
    _ins_size_is_fixed = False
    _outs_size_is_fixed = True
    
    _units = {
        'Coagulant': 'kg/hr',
              }
    
    _default_flocculant_prices = {
        'aluminum sulfate': (16, 0.32, 0.081),  # Al2(SO4)3·18H2O, purity in %, price in USD/kg, Al content -> 24.69 USD/kg Al
        'PAC': (28, 0.5, 0.26),                 # Al2(OH)nCl(6-n) or "polyaluminum chloride", 0.2174-0.3093 gAl/g depending on n -> 5.77~8.21 USD/kg Al
        'sodium aluminate': (48, 1.2, 0.329),   # NaAlO2 -> 7.60 USD/kg Al
        'ferric chloride': (21, 0.4, 0.3443),   # FeCl3 -> 5.53 USD/kg Fe
        'ferric sulfate': (20, 0.5, 0.1987),    # Fe2(SO4)3·9H2O -> 12.58 USD/kg Fe
        }

    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', F_BM_default=None, isdynamic=False,
                 P_precipitation_stoichiometry=None, metal_dosage=0.,
                 metal_ID='', orthoP_ID='', metalP_ID='', metal_price=None,
                 metalP_pKsp=13.75, alpha=1.8e-6,
                 soluble_substrate_ID='S_F', particulate_substrate_ID='X_S',
                 soluble_inert_ID='S_I', particulate_inert_ID='X_I', 
                 Fmin_substrate=4, Fmax_substrate=20, ka_substrate=0.5,
                 Fmin_inert=4, Fmax_inert=20, ka_inert=0.5,        
                 ):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with,
                         F_BM_default=F_BM_default, isdynamic=isdynamic)
        self.metal_dosage = metal_dosage
        self._cmps_idx = self.chemicals.indices([metal_ID, orthoP_ID, metalP_ID])
        self._org_idx = self.chemicals.indices([
            soluble_substrate_ID, particulate_substrate_ID, 
            soluble_inert_ID, particulate_inert_ID])
        self.metal_ID = metal_ID
        self.orthoP_ID = orthoP_ID
        self.metalP_ID = metalP_ID
        self.soluble_substrate_ID = soluble_substrate_ID
        self.particulate_substrate_ID = particulate_substrate_ID
        self._check_composition(soluble_substrate_ID, particulate_substrate_ID)
        self.soluble_inert_ID = soluble_inert_ID
        self.particulate_inert_ID = particulate_inert_ID
        self._check_composition(soluble_inert_ID, particulate_inert_ID)
        self.P_precipitation_stoichiometry = P_precipitation_stoichiometry
        self.metalP_pKsp = metalP_pKsp
        self.alpha = alpha
        self.Fmin_substrate = Fmin_substrate
        self.Fmax_substrate = Fmax_substrate
        self.ka_substrate = ka_substrate
        self.Fmin_inert = Fmin_inert
        self.Fmax_inert = Fmax_inert
        self.ka_inert = ka_inert
        self.metal_price = metal_price

    @classmethod
    def from_mASM2d(cls, ID, ins=None, outs=(), model=None, 
                    metal_ID='X_FeOH', metal_dosage=0., pH=7, 
                    Fmin_substrate=4, Fmax_substrate=20, ka_substrate=0.5,
                    Fmin_inert=4, Fmax_inert=20, ka_inert=0.5, 
                    thermo=None, init_with='WasteStream', 
                    F_BM_default=None, isdynamic=True, **kwargs):
        """
        Creating a metal dosage unit using data from a `mASM2d` process model.

        Parameters
        ----------
        model : :class:`mASM2d`
            mASM2d process model.
        pH : float, optional
            Affects orpho-phosphate speciation. The default is 7.

        Examples
        --------
        >>> from qsdsan import processes as pc, sanunits as su
        >>> cmps = pc.create_masm2d_cmps()
        >>> inf = pc.create_masm2d_inf('inf', 10)
        >>> asm = pc.mASM2d()
        >>> MD = su.MetalDosage.from_mASM2d('MD', inf, 'eff', metal_dosage=50, 
        ...                                 model=asm, isdynamic=False)
        >>> inf.show()
        WasteStream: inf to <MetalDosage: MD>
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (g/hr): S_N2   7.5
                     S_NH4  10.4
                     S_PO4  3.33
                     S_F    35.8
                     S_I    8.96
                     S_IC   35
                     S_K    11.7
                     S_Mg   20.8
                     X_I    23.3
                     X_S    111
                     S_Ca   58.3
                     S_Na   36.2
                     S_Cl   177
                     H2O    4.15e+05
         WasteStream-specific properties:
          pH         : 7.0
          Alkalinity : 7.0 mmol/L
          COD        : 430.0 mg/L
          BOD        : 216.3 mg/L
          TC         : 224.3 mg/L
          TOC        : 140.3 mg/L
          TN         : 41.5 mg/L
          TP         : 10.5 mg/L
          TK         : 28.0 mg/L
          TSS        : 241.9 mg/L
         Component concentrations (mg/L):
          S_N2     18.0
          S_NH4    25.0
          S_PO4    8.0
          S_F      86.0
          S_I      21.5
          S_IC     84.0
          S_K      28.0
          S_Mg     50.0
          X_I      55.9
          X_S      266.6
          S_Ca     140.0
          S_Na     87.0
          S_Cl     425.0
          H2O      995387.5

        >>> MD.simulate()
        >>> eff, = MD.outs
        >>> eff.show()
        WasteStream: eff from <MetalDosage: MD>
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (g/hr): S_N2     7.5
                     S_NH4    10.4
                     S_PO4    0.000407
                     S_F      30.6
                     S_I      3.79
                     S_IC     35
                     S_K      11.7
                     S_Mg     20.8
                     X_I      28.5
                     X_S      116
                     S_Ca     58.3
                     X_FeOH   9.33
                     X_FePO4  16.2
                     S_Na     36.2
                     S_Cl     177
                     ...      4.15e+05
         WasteStream-specific properties:
          pH         : 7.0
          Alkalinity : 7.0 mmol/L
          COD        : 430.0 mg/L
          BOD        : 214.6 mg/L
          TC         : 224.3 mg/L
          TOC        : 140.3 mg/L
          TN         : 41.5 mg/L
          TP         : 10.5 mg/L
          TK         : 28.0 mg/L
          TSS        : 321.9 mg/L
         Component concentrations (mg/L):
          S_N2     18.0
          S_NH4    25.0
          S_PO4    0.0
          S_F      73.5
          S_I      9.1
          S_IC     84.0
          S_K      28.0
          S_Mg     50.0
          X_I      68.3
          X_S      279.1
          S_Ca     140.0
          X_FeOH   22.4
          X_FePO4  38.9
          S_Na     87.0
          S_Cl     425.0
          ...
        """
        if not isinstance(model, mASM2d):
            raise TypeError(f'model must be an instance of `mASM2d` class, not {type(model)}')
        if metal_ID == 'X_FeOH': 
            metalP_ID = 'X_FePO4'
            metalP_pKsp = kwargs.pop('metalP_pKsp', 26.4)
        elif metal_ID == 'X_AlOH': 
            metalP_ID = 'X_AlPO4'
            metalP_pKsp = kwargs.pop('metalP_pKsp', 18.2)
        else: 
            metalP_ID = kwargs.pop('metalP_ID')
            metalP_pKsp = kwargs.pop('metalP_pKsp', 13.75)
        stoichio = model.mmp_stoichio[metalP_ID]
        Kas = model.rate_function.params['Ka'][4:7]
        alpha = kwargs.pop('alpha', ion_speciation(10**(-pH), *Kas)[-1])
        self = cls(ID, ins, outs, thermo, init_with, F_BM_default, isdynamic,
                   P_precipitation_stoichiometry=stoichio, 
                   metal_dosage=metal_dosage,
                   metal_ID=metal_ID, orthoP_ID='S_PO4', metalP_ID=metalP_ID, 
                   metalP_pKsp=metalP_pKsp, alpha=alpha,
                   soluble_substrate_ID='S_F', particulate_substrate_ID='X_S',
                   soluble_inert_ID='S_I', particulate_inert_ID='X_I', 
                   Fmin_substrate=Fmin_substrate, 
                   Fmax_substrate=Fmax_substrate, ka_substrate=ka_substrate,
                   Fmin_inert=Fmin_inert, Fmax_inert=Fmax_inert, ka_inert=ka_inert)
        return self

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
        metalP = stoichio.pop(self.metalP_ID, 1)
        if metalP != 1:
            stoichio = {k: v/metalP for k,v in stoichio.items()}
        stoichio[self.metalP_ID] = metalP
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
    
    @property
    def metal_price(self):
        '''[float] Price of the dosed metal, in USD/kg metal component.'''
        return self._price
    @metal_price.setter
    def metal_price(self, p):
        if p is None:
            cmps = self.chemicals
            mtid = self.metal_ID
            mt = cmps[mtid]
            if 'Al' in mtid:
                i_metal = mw(dict(Al=mt.atoms['Al'])) / mt.chem_MW * mt.i_mass
                purity, price, content = self._default_flocculant_prices['PAC']
            elif 'Fe' in mtid:
                i_metal = mw(dict(Fe=mt.atoms['Fe']))/mt.chem_MW * mt.i_mass
                purity, price, content = self._default_flocculant_prices['ferric chloride']
            else:
                raise ValueError(f'unrecognized metal element in {mtid},'
                                 'must provide metal price.')
            p = price / (purity/100) / content * i_metal
        self._price = p
        
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
        SP = flx.bisection(_precipitation_mass_balance, 0, SP_in, args=(
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
            try: 
                sp = flx.bisection(
                    _precipitation_mass_balance, 0, SP_in, args=(
                    Me_in, SP_in, Ksp_mass, x, y, i, j, alpha
                    ))
            except:
                # sp = flx.aitken_secant(
                #     _precipitation_mass_balance, SP_in, args=(
                #     Me_in, SP_in, Ksp_mass, x, y, i, j, alpha
                #     ))
                # if sp < 0 or sp > SP_in:
                #     warn(f'sp = {sp}; sp_in = {SP_in}')
                #     sp = max(0, min(sp, SP_in))
                sp = SP_in
            return sp
        
        def solve_ss(SS_in):
            try:
                ss = flx.IQ_interpolation(
                    _coagulation_mass_balance, 0, SS_in, args=(
                    SS_in, Pme, Fmax_ss, Fmin_ss, ka_ss
                    ))
            except:
                ss = SS_in
            return ss
        
        def solve_si(SI_in):
            try:
                si = flx.IQ_interpolation(
                    _coagulation_mass_balance, 0, SI_in, args=(
                    SI_in, Pme, Fmax_si, Fmin_si, ka_si
                    ))
            except:
                si = SI_in
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
        
    def _design(self):
        D = self.design_results
        out, = self.outs
        D['Coagulant'] = self.metal_dosage / 1000 * out.F_vol # kg/m3 * m3/hr
        
    def _cost(self):
        D = self.design_results
        opex = self.add_OPEX = {}
        opex['Coagulant'] = D['Coagulant'] * self.metal_price # kg/hr * USD/kg