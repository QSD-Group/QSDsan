# -*- coding: utf-8 -*-
'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    
    Joy Zhang <joycheung1994@gmail.com>
    
    Yalin Li <mailto.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

import numpy as np
from warnings import warn
from math import isclose
from biosteam.units import Junction as BSTjunction
from .. import SanUnit, processes as pc

__all__ = (
    'Junction',
    'ADMjunction', 'ADMtoASM', 'ASMtoADM',
    )

class Junction(SanUnit):
    '''
    A non-reactive class that serves to convert a stream with one set of components
    and into another.
    
    Thermal conditions of the downstream (T, P) will be copied from those of the upstream.

    Parameters
    ----------
    upstream : stream or str
        Influent stream.
    downstream : stream or str
        Effluent stream.
    reactions : iterable(dict) | callable
        Iterable of dict that has the conversion of upstream components to
        downstream components,
        or a function that will return the concentration of the effluent
        with influent concentration as the input.
        
        If given as an iterable of dict, keys of each of the dict should be 
        the ID or alias of components,
        values should be the conversion/yield,
        which should be negative for reactants and positive for products.
    '''
    _graphics = BSTjunction._graphics

    def __init__(self, ID='', upstream=None, downstream=(), thermo=None,
                 init_with='WasteStream', F_BM_default=None, isdynamic=False,
                 reactions=None, **kwargs):
        thermo = downstream.thermo if downstream else thermo
        SanUnit.__init__(self, ID, ins=upstream, outs=downstream, thermo=thermo,
                         init_with=init_with,
                         F_BM_default=F_BM_default, isdynamic=isdynamic,
                         skip_property_package_check=True)
        if reactions: self.reactions = reactions
        else: self.reactions = None
        
        for key, val in kwargs.items():
            setattr(self, key, val)

    def _no_parse_reactions(self, rxns):
        if rxns is None: return
        raise RuntimeError('Reactions are automatically compiled.')


    def _parse_reactions(self, rxns):
        cmps_in = self.ins[0].components
        cmps_outs = self.outs[0].components

        num_rxns = len(rxns)
        num_upcmps = len(cmps_in)
        num_downcmps = len(cmps_outs)
        RX = np.zeros(shape=(num_rxns, num_upcmps)) # each row is a reaction
        RY = np.zeros(shape=(num_rxns, num_downcmps))

        for n, dct in enumerate(rxns):
            RXn, RYn = RX[n], RY[n]
            for cmp, val in dct.items():
                if val < 0: RXn[cmps_in.index(cmp)] = val
                elif val > 0: RYn[cmps_outs.index(cmp)] = val

        # Transfer overlapped components
        overlapped = set(cmps_in.IDs).intersection(set(cmps_outs.IDs))
        RXsum = RX.sum(axis=0)
        RXnew = []
        RYnew = []
        for ID in overlapped:
            idx_up = cmps_in.index(ID)
            idx_down= cmps_outs.index(ID)
            if RXsum[idx_up] != -1:
                newRXn = np.zeros(num_upcmps)
                newRYn = np.zeros(num_downcmps)
                newRXn[idx_up] = -1 - RXsum[idx_up]
                newRYn[idx_down] = -newRXn[idx_up]
                RXnew.append(newRXn)
                RYnew.append(newRYn)
        RX = np.concatenate((RX, np.array(RXnew)))
        RY = np.concatenate((RY, np.array(RYnew)))

        # Check if all upstream components are converted
        RXsum = RX.sum(axis=0)
        if np.where(RXsum!=-1)[0].any():
            index = np.where(RXsum!=-1)[0][0]
            raise ValueError(f'Conversion for Component "{cmps_in.IDs[index]}" '
                             f'is {abs(RXsum[index])}, not 100%.')
        self._RX = RX
        self._RY = RY
        
        
    def _compile_reactions(self):
        def reactions(X):
            X.reshape(1, len(X))
            Yarr = -(self._RX*X).T @ self._RY # _RX: (num_rxns, num_upcmps); _RY: (num_rxns, num_downcmps)
            Y = Yarr.sum(axis=0) # Yarr: (num_upcmps, num_downcmps)
            return Y.reshape(Y.shape[1],)
        self._reactions = reactions
        

    def _run(self):
        ins = self.ins[0]
        rxns = self.reactions
        X = ins.conc.value
        Y = rxns(X)
        outs = self.outs[0]
        outs.thermal_condition.copy_like(ins.thermal_condition)
        concs = dict(zip(outs.components.IDs, Y))
        concs.pop('H2O', None)
        outs.set_flow_by_concentration(
            flow_tot=ins.F_vol,
            concentrations=concs,
            units=('m3/hr', 'mg/L'))


    # Below are dynamic simulation-related properties
    @property
    def state(self):
        '''The state of the Junction, including component concentrations [mg/L] and flow rate [m^3/d].'''
        if self._state is None: return None
        else:
            return dict(zip(list(self.components.IDs) + ['Q'], self._state))

    def _init_dynamic(self):
        super()._init_dynamic()
        # Need to use ins' components, otherwise _ins_QC will follow the shape of
        # the unit's (i.e., downstream) components
        self._ins_QC = np.zeros((len(self._ins), len(self.ins[0].components)+1))
        self._ins_dQC = self._ins_QC.copy()

    def _init_state(self):
        '''
        Initialize state by specifying or calculating component concentrations
        based on influents. Total flow rate is always initialized as the sum of
        influent wastestream flows.
        '''
        self._state = np.append(self.outs[0].conc, self.outs[0].F_vol*24)
        self._dstate = self._state * 0.
        self._cached_state = self._state.copy()
        self._cached_t = 0

    def _update_state(self):
        '''
        Updates conditions of output stream based on conditions of the Junction.
        '''
        self._outs[0].state = self._state

    def _update_dstate(self):
        '''
        Updates rates of change of output stream from rates of change of the Junction.
        '''
        self._outs[0].dstate = self._dstate

    # The unit's state should be the same as the effluent state
    # react the state arr and dstate arr
    def _compile_AE(self):
        _state = self._state
        _dstate = self._dstate
        _update_state = self._update_state
        _update_dstate = self._update_dstate
        rxns = self.reactions
        def yt(t, QC_ins, dQC_ins):
            for i, j in zip((QC_ins, dQC_ins), (_state, _dstate)):
                X = i[0][:-1] # shape = (1, num_upcmps)
                Y = rxns(X)
                j[:-1] = Y
                j[-1] = i[0][-1] # volumetric flow of outs should equal that of ins
            _update_state()
            _update_dstate()
        self._AE = yt

    
    @property
    def upstream(self):
        '''[qsdsan.WasteStream] Influent.'''
        return self.ins[0]
    @upstream.setter
    def upstream(self, upstream):
        self.ins[0] = upstream

    @property
    def downstream(self):
        '''[qsdsan.WasteStream] Effluent.'''
        return self.outs[0]
    @downstream.setter
    def downstream(self, downstream):
        self.outs[0] = downstream
    
    @property
    def AE(self):
        if self._AE is None:
            self._compile_AE()
        return self._AE
    
    @property
    def reactions(self):
        '''
        [callable] Function that takes the concentration array of the influent
        and convert to the concentration array of the effluent.
        '''
        if not self._reactions: self._compile_reactions()
        return self._reactions
    @reactions.setter
    def reactions(self, i):
        if callable(i): self._reactions = i
        else:
            self._parse_reactions(i)
            self._compile_reactions()
            
            
# %%

#TODO: add a `rtol` kwargs for error checking
class ADMjunction(Junction):
    '''
    An abstract superclass holding common properties of ADM interface classes.
    Users should use its subclasses (e.g., ``ASMtoADM``, ``ADMtoASM``) instead.
    
    See Also
    --------
    :class:`qsdsan.sanunits.Junction`
    
    :class:`qsdsan.sanunits.ADMtoASM`
    
    :class:`qsdsan.sanunits.ASMtoADM`
    '''
    _parse_reactions = Junction._no_parse_reactions
    rtol = 1e-2
    atol = 1e-6
    
    def __init__(self, ID='', upstream=None, downstream=(), thermo=None,
                 init_with='WasteStream', F_BM_default=None, isdynamic=False,
                 adm1_model=None):
        self.adm1_model = adm1_model # otherwise there won't be adm1_model when `_compile_reactions` is called
        if thermo is None:
            warn('No `thermo` object is provided and is prone to raise error. '
                 'If you are not sure how to get the `thermo` object, '
                 'use `thermo = qsdsan.set_thermo` after setting thermo with the `Components` object.')
        super().__init__(ID=ID, upstream=upstream, downstream=downstream,
                         thermo=thermo, init_with=init_with, 
                         F_BM_default=F_BM_default, isdynamic=isdynamic)
        
   
    @property
    def T(self):
        '''[float] Temperature of the upstream/downstream [K].'''
        return self.ins[0].T
    @T.setter
    def T(self, T):
        self.ins[0].T = self.outs[0].T = T
    
    @property
    def pH(self):
        '''[float] pH of the upstream/downstream.'''
        return self.ins[0].pH
    
    @property
    def adm1_model(self):
        '''[qsdsan.Process] ADM process model.'''
        return self._adm1_model
    @adm1_model.setter
    def adm1_model(self, model):
        if not isinstance(model, pc.ADM1):
            raise ValueError('`adm1_model` must be an `AMD1` object, '
                             f'the given object is {type(model).__name__}.')
        self._adm1_model = model
        
    @property
    def T_base(self):
        '''[float] Base temperature in the ADM1 model.'''
        return self.adm1_model.rate_function.params['T_base']
    
    @property
    def pKa_base(self):
        '''[float] pKa of the acid-base pairs at the base temperature in the ADM1 model.'''
        Ka_base = self.adm1_model.rate_function.params['Ka_base']
        return -np.log10(Ka_base)
    
    @property
    def Ka_dH(self):
        '''[float] Heat of reaction for Ka.'''
        return self.adm1_model.rate_function.params['Ka_dH']
    
    @property
    def pKa(self):
        '''
        [numpy.array] pKa array of the following acid-base pairs:
        ('H+', 'OH-'), ('NH4+', 'NH3'), ('CO2', 'HCO3-'),
        ('HAc', 'Ac-'), ('HPr', 'Pr-'), ('HBu', 'Bu-'), ('HVa', 'Va-')
        '''
        return self.pKa_base-np.log10(np.exp(pc.T_correction_factor(self.T_base, self.T, self.Ka_dH)))
    
    @property
    def alpha_IC(self):
        '''[float] Charge per g of C.'''
        pH = self.pH
        pKa_IC = self.pKa[2]
        return -1/(1+10**(pKa_IC-pH))/12

    @property
    def alpha_IN(self):
        '''[float] Charge per g of N.'''
        pH = self.pH
        pKa_IN = self.pKa[1]
        return 10**(pKa_IN-pH)/(1+10**(pKa_IN-pH))/14
    
        
    def _compile_AE(self):
        _state = self._state
        _dstate = self._dstate
        _cached_state = self._cached_state
        _update_state = self._update_state
        _update_dstate = self._update_dstate
        rxn = self.reactions
               
        def yt(t, QC_ins, dQC_ins):
            before_vals = QC_ins[0,:-1]
            _state[:-1] = rxn(before_vals)
            _state[-1] = QC_ins[0,-1]
            if t > self._cached_t:
                _dstate[:] = (_state - _cached_state)/(t-self._cached_t)
            _cached_state[:] = _state
            self._cached_t = t
            _update_state()
            _update_dstate()
        
        self._AE = yt
    
# %%

class ADMtoASM(ADMjunction):
    '''
    Interface unit to convert anaerobic digestion model (ADM) components
    to activated sludge model (ASM) components.
    
    Parameters
    ----------
    upstream : stream or str
        Influent stream with ADM components.
    downstream : stream or str
        Effluent stream with ASM components.
    adm1_model : obj
        The anaerobic digestion process model (:class:`qsdsan.processes.ADM1`).
    bio_to_xs : float
        Split of the total biomass COD to slowly biodegradable substrate (X_S),
        the rest is assumed to be mapped into X_P.
    rtol : float
        Relative tolerance for COD and TKN balance.
    atol : float
        Absolute tolerance for COD and TKN balance.
    
    References
    ----------
    [1] Nopens, I.; Batstone, D. J.; Copp, J. B.; Jeppsson, U.; Volcke, E.; 
    Alex, J.; Vanrolleghem, P. A. An ASM/ADM Model Interface for Dynamic 
    Plant-Wide Simulation. Water Res. 2009, 43, 1913–1923.
    
    See Also
    --------
    :class:`qsdsan.sanunits.ADMjunction`
    
    :class:`qsdsan.sanunits.ASMtoADM`  
    
    `math.isclose <https://docs.python.org/3.8/library/math.html#math.isclose>`
    '''
    # User defined values
    bio_to_xs = 0.7
    
    # Should be constants
    cod_vfa = np.array([64, 112, 160, 208])
    
    def isbalanced(self, lhs, rhs_vals, rhs_i):
        rhs = sum(rhs_vals*rhs_i)
        error = rhs - lhs
        tol = max(self.rtol*lhs, self.rtol*rhs, self.atol)
        return abs(error) <= tol, error, tol, rhs
    
    def balance_cod_tkn(self, adm_vals, asm_vals):
        cmps_adm = self.ins[0].components
        cmps_asm = self.outs[0].components
        asm_i_COD = cmps_asm.i_COD
        adm_i_COD = cmps_adm.i_COD
        gas_idx = cmps_adm.indices(('S_h2', 'S_ch4'))
        asm_i_N = cmps_asm.i_N
        adm_i_N = cmps_adm.i_N
        adm_cod = sum(adm_vals*adm_i_COD) - sum(adm_vals[gas_idx])
        adm_tkn = sum(adm_vals*adm_i_N)
        cod_bl, cod_err, cod_tol, asm_cod = self.isbalanced(adm_cod, asm_vals, asm_i_COD)
        tkn_bl, tkn_err, tkn_tol, asm_tkn = self.isbalanced(adm_tkn, asm_vals, asm_i_N)
        if cod_bl:
            if tkn_bl: return asm_vals
            else:
                if tkn_err > 0: dtkn = -(tkn_err - tkn_tol)/asm_tkn
                else: dtkn = -(tkn_err + tkn_tol)/asm_tkn
                _asm_vals = asm_vals * (1 + (asm_i_N>0)*dtkn)
                _cod_bl, _cod_err, _cod_tol, _asm_cod = self.isbalanced(adm_cod, _asm_vals, asm_i_COD)
                if _cod_bl: return _asm_vals
                else: 
                    warn('cannot balance COD and TKN at the same '
                        f'time with rtol={self.rtol} and atol={self.atol}. '
                        f'influent (ADM) COD is {adm_cod}, '
                        f'effluent (ASM) COD is {asm_cod} or {_asm_cod}. '
                        f'influent TKN is {adm_tkn}, ' 
                        f'effluent TKN is {asm_tkn} or {asm_tkn*(1+dtkn)}. ')
                    return asm_vals
        else:
            if cod_err > 0: dcod = -(cod_err - cod_tol)/asm_cod
            else: dcod = -(cod_err + cod_tol)/asm_cod
            _asm_vals = asm_vals * (1 + (asm_i_COD>0)*dcod)
            _tkn_bl, _tkn_err, _tkn_tol, _asm_tkn = self.isbalanced(adm_tkn, _asm_vals, asm_i_N)
            if _tkn_bl: return _asm_vals
            else:
                warn('cannot balance COD and TKN at the same '
                    f'time with rtol={self.rtol} and atol={self.atol}. '
                    f'influent (ADM) COD is {adm_cod}, '
                    f'effluent (ASM) COD is {asm_cod} or {asm_cod*(1+dcod)}. '
                    f'influent TKN is {adm_tkn}, ' 
                    f'effluent TKN is {asm_tkn} or {_asm_tkn}. ')
                return asm_vals
    
    def _compile_reactions(self):
        # Retrieve constants
        ins = self.ins[0]
        outs = self.outs[0]
        rtol = self.rtol
        atol = self.atol
        
        cmps_adm = ins.components
        X_c_i_N = cmps_adm.X_c.i_N
        X_pr_i_N = cmps_adm.X_pr.i_N
        S_aa_i_N = cmps_adm.S_aa.i_N
        adm_X_I_i_N = cmps_adm.X_I.i_N
        adm_S_I_i_N = cmps_adm.S_I.i_N
        adm_i_N = cmps_adm.i_N
        adm_bio_N_indices = cmps_adm.indices(('X_su', 'X_aa', 'X_fa', 
                                              'X_c4', 'X_pro', 'X_ac', 'X_h2'))

        cmps_asm = outs.components
        X_P_i_N = cmps_asm.X_P.i_N
        X_S_i_N = cmps_asm.X_S.i_N
        S_S_i_N = cmps_asm.S_S.i_N
        asm_X_I_i_N = cmps_asm.X_I.i_N
        asm_S_I_i_N = cmps_asm.S_I.i_N
        asm_X_P_i_N = cmps_asm.X_P.i_N
        asm_ions_idx = cmps_asm.indices(('S_NH', 'S_ALK'))
        
        alpha_IN = self.alpha_IN
        alpha_IC = self.alpha_IC
        alpha_vfa = self.alpha_vfa
        f_corr = self.balance_cod_tkn

        def adm2asm(adm_vals):    
            S_su, S_aa, S_fa, S_va, S_bu, S_pro, S_ac, S_h2, S_ch4, S_IC, S_IN, S_I, \
                X_c, X_ch, X_pr, X_li, X_su, X_aa, X_fa, X_c4, X_pro, X_ac, X_h2, X_I, \
                S_cat, S_an, H2O = adm_vals
                       
            # Step 0: snapshot of charged components
            _ions = np.array([S_IN, S_IC, S_ac, S_pro, S_bu, S_va])
            
            # Step 1a: convert biomass into X_S+X_ND and X_P
            bio_cod = X_su + X_aa + X_fa + X_c4 + X_pro + X_ac + X_h2
            bio_n = sum((adm_vals*adm_i_N)[adm_bio_N_indices])
            xp_cod = bio_cod * (1-self.bio_to_xs)
            xp_ndm = xp_cod*X_P_i_N
            if xp_ndm > bio_n:
                warn('Not enough biomass N to map the specified proportion of '
                     'biomass COD into X_P. Mapped as much COD as possible, the rest '
                     'goes to X_S.')
                X_P = bio_n/asm_X_P_i_N
                bio_n = 0
            else:
                X_P = xp_cod
                bio_n -= xp_ndm
            X_S = bio_cod - X_P
            xs_ndm = X_S*X_S_i_N
            if xs_ndm <= bio_n:
                X_ND = bio_n - xs_ndm
                bio_n = 0
            elif xs_ndm <= bio_n + S_IN:
                X_ND = 0
                S_IN -= (xs_ndm - bio_n)
                bio_n = 0
            else:
                if isclose(xs_ndm,  bio_n + S_IN, rel_tol=rtol, abs_tol=atol):
                    X_ND = S_IN = bio_n = 0
                else:
                    raise RuntimeError('Not enough nitrogen (S_IN + biomass) to map '
                                       'all biomass COD into X_P and X_S')
            
            # Step 1b: convert particulate substrates into X_S + X_ND
            xsub_cod = X_c + X_ch + X_pr + X_li
            xsub_n = X_c*X_c_i_N + X_pr*X_pr_i_N
            X_S += xsub_cod
            X_ND += xsub_n - xsub_cod*X_S_i_N  # X_S.i_N should technically be zero
            if X_ND < 0:
                if isclose(X_ND, 0, rel_tol=rtol, abs_tol=atol): X_ND = 0
                else:
                    raise RuntimeError('Not enough nitrogen (substrate + excess X_ND) '
                                       'to map all particulate substrate COD into X_S')
            
            # Step 2: map all X_I from ADM to ASM           
            excess_XIn = X_I * (adm_X_I_i_N - asm_X_I_i_N)
            S_IN += excess_XIn
            if S_IN < 0:
                if isclose(S_IN, 0, rel_tol=rtol, abs_tol=atol): S_IN = 0
                raise RuntimeError('Not enough nitrogen (X_I + S_IN) to map '
                                   'all ADM X_I into ASM X_I')
            
            # Step 3: map ADM S_I into ASM S_I and S_NH
            excess_SIn = S_I * (adm_S_I_i_N - asm_S_I_i_N)
            if excess_SIn > 0:
                S_NH = excess_SIn
            else:
                S_NH = 0
                S_IN += excess_SIn
                if S_IN < 0:
                    if isclose(S_IN, 0, rel_tol=rtol, abs_tol=atol): S_IN = 0
                    raise RuntimeError('Not enough nitrogen (S_I + S_IN) to map '
                                       'all ADM S_I into ASM S_I')
            S_NH += S_IN
                
            # Step 4: map all soluble substrates into S_S and S_ND        
            ssub_cod = S_su + S_aa + S_fa + S_va + S_bu + S_pro + S_ac
            ssub_n = S_aa * S_aa_i_N
            if ssub_cod*S_S_i_N <= ssub_n:
                S_S = ssub_cod
                S_ND = ssub_n
                if S_S_i_N != 0: S_ND -= S_S/S_S_i_N # S_S.i_N should technically be zero
            else:
                if isclose(ssub_cod*S_S_i_N, ssub_n, rel_tol=rtol, abs_tol=atol): 
                    S_S = ssub_cod
                    S_ND = 0
                else:
                    raise RuntimeError('Not enough nitrogen to map all soluble '
                                       'substrates into ASM S_S')
                        
            # Step 6: check COD and TKN balance
            asm_vals = np.array(([
                S_I, S_S, X_I, X_S, 
                0, 0, # X_BH, X_BA, 
                X_P, 
                0, 0, # S_O, S_NO, 
                S_NH, S_ND, X_ND,  
                0, 0, # temporary S_ALK, S_N2, 
                H2O]))
            
            if S_h2 > 0 or S_ch4 > 0:
                warn('Ignored dissolved H2 or CH4.')
            
            asm_vals = f_corr(adm_vals, asm_vals)
            
            # Step 5: charge balance for alkalinity
            S_NH = asm_vals[asm_ions_idx[0]]
            S_ALK = (sum(_ions * np.append([alpha_IN, alpha_IC], alpha_vfa)) - S_NH/14)*(-12)
            asm_vals[asm_ions_idx[1]] = S_ALK
            
            return asm_vals
        
        self._reactions = adm2asm
    
    @property
    def alpha_vfa(self):
        return 1.0/self.cod_vfa*(-1.0/(1.0 + 10**(self.pKa[3:]-self.pH)))
        
        
# %%

class ASMtoADM(ADMjunction):
    '''
    Interface unit to convert activated sludge model (ASM) components
    to anaerobic digestion model (ADM) components.
    
    Parameters
    ----------
    upstream : stream or str
        Influent stream with ASM components.
    downstream : stream or str
        Effluent stream with ADM components.
    adm1_model : obj
        The anaerobic digestion process model (:class:`qsdsan.processes.ADM1`).
    xs_to_li : float
        Split of slowly biodegradable substrate COD to lipid, 
        after all N is mapped into protein.
    bio_to_li : float
        Split of biomass COD to lipid, after all biomass N is
        mapped into protein.
    frac_deg : float
        Biodegradable fraction of biomass COD.
    rtol : float
        Relative tolerance for COD and TKN balance.
    atol : float
        Absolute tolerance for COD and TKN balance.
    
    References
    ----------
    [1] Nopens, I.; Batstone, D. J.; Copp, J. B.; Jeppsson, U.; Volcke, E.; 
    Alex, J.; Vanrolleghem, P. A. An ASM/ADM Model Interface for Dynamic 
    Plant-Wide Simulation. Water Res. 2009, 43, 1913–1923.
    
    See Also
    --------
    :class:`qsdsan.sanunits.ADMjunction`
    
    :class:`qsdsan.sanunits.ADMtoASM` 
    
    `math.isclose <https://docs.python.org/3.8/library/math.html#math.isclose>`
    '''    
    # User defined values
    xs_to_li = 0.7
    bio_to_li = 0.4
    frac_deg = 0.68
    
    
    def isbalanced(self, lhs, rhs_vals, rhs_i):
        rhs = sum(rhs_vals*rhs_i)
        error = rhs - lhs
        tol = max(self.rtol*lhs, self.rtol*rhs, self.atol)
        return abs(error) <= tol, error, tol, rhs
    
    def balance_cod_tkn(self, asm_vals, adm_vals):
        cmps_asm = self.ins[0].components
        cmps_adm = self.outs[0].components
        asm_i_COD = cmps_asm.i_COD
        adm_i_COD = cmps_adm.i_COD
        non_tkn_idx = cmps_asm.indices(('S_NO', 'S_N2'))
        asm_i_N = cmps_asm.i_N
        adm_i_N = cmps_adm.i_N
        asm_cod = sum(asm_vals*asm_i_COD)
        asm_tkn = sum(asm_vals*asm_i_N) - sum(asm_vals[non_tkn_idx])
        cod_bl, cod_err, cod_tol, adm_cod = self.isbalanced(asm_cod, adm_vals, adm_i_COD)
        tkn_bl, tkn_err, tkn_tol, adm_tkn = self.isbalanced(asm_tkn, adm_vals, adm_i_N)
        if cod_bl:
            if tkn_bl: return adm_vals
            else:
                if tkn_err > 0: dtkn = -(tkn_err - tkn_tol)/adm_tkn
                else: dtkn = -(tkn_err + tkn_tol)/adm_tkn
                _adm_vals = adm_vals * (1 + (adm_i_N>0)*dtkn)
                _cod_bl, _cod_err, _cod_tol, _adm_cod = self.isbalanced(asm_cod, _adm_vals, adm_i_COD)
                if _cod_bl: return _adm_vals
                else: 
                    warn('cannot balance COD and TKN at the same '
                        f'time with rtol={self.rtol} and atol={self.atol}.\n '
                        f'influent (ASM) COD is {asm_cod}\n '
                        f'effluent (ADM) COD is {adm_cod} or {_adm_cod}\n '
                        f'influent TKN is {asm_tkn}\n ' 
                        f'effluent TKN is {adm_tkn} or {adm_tkn*(1+dtkn)}. ')
                    return adm_vals
        else:
            if cod_err > 0: dcod = -(cod_err - cod_tol)/adm_cod
            else: dcod = -(cod_err + cod_tol)/adm_cod
            _adm_vals = adm_vals * (1 + (adm_i_COD>0)*dcod)
            _tkn_bl, _tkn_err, _tkn_tol, _adm_tkn = self.isbalanced(asm_tkn, _adm_vals, adm_i_N)
            if _tkn_bl: return _adm_vals
            else:
                warn('cannot balance COD and TKN at the same '
                    f'time with rtol={self.rtol} and atol={self.atol}.\n '
                    f'influent (ASM) COD is {asm_cod}\n '
                    f'effluent (ADM) COD is {adm_cod} or {adm_cod*(1+dcod)}\n '
                    f'influent TKN is {asm_tkn}\n ' 
                    f'effluent TKN is {adm_tkn} or {_adm_tkn}. ')
                return adm_vals
                
    def _compile_reactions(self):
        # Retrieve constants
        ins = self.ins[0]
        outs = self.outs[0]
        rtol = self.rtol
        atol = self.atol

        cmps_asm = ins.components
        S_NO_i_COD = cmps_asm.S_NO.i_COD
        X_BH_i_N = cmps_asm.X_BH.i_N
        X_BA_i_N = cmps_asm.X_BA.i_N
        asm_X_I_i_N = cmps_asm.X_I.i_N
        X_P_i_N = cmps_asm.X_P.i_N
        if cmps_asm.X_S.i_N > 0: 
            warn(f'X_S in ASM has positive nitrogen content: {cmps_asm.X_S.i_N} gN/gCOD. '
                 'These nitrogen will be ignored by the interface model '
                 'and could lead to imbalance of TKN after conversion.')
        if cmps_asm.S_S.i_N > 0: 
            warn(f'S_S in ASM has positive nitrogen content: {cmps_asm.S_S.i_N} gN/gCOD. '
                 'These nitrogen will be ignored by the interface model '
                 'and could lead to imbalance of TKN after conversion.')
        
        cmps_adm = outs.components
        S_aa_i_N = cmps_adm.S_aa.i_N
        X_pr_i_N = cmps_adm.X_pr.i_N
        S_I_i_N = cmps_adm.S_I.i_N
        adm_X_I_i_N = cmps_adm.X_I.i_N
        adm_ions_idx = cmps_adm.indices(['S_IN', 'S_IC', 'S_cat', 'S_an'])
        
        frac_deg = self.frac_deg
        alpha_IN = self.alpha_IN
        alpha_IC = self.alpha_IC
        proton_charge = 10**(-self.pKa[0]+self.pH) - 10**(-self.pH) # self.pKa[0] is pKw
        f_corr = self.balance_cod_tkn

        def asm2adm(asm_vals):
            S_I, S_S, X_I, X_S, X_BH, X_BA, X_P, S_O, S_NO, S_NH, S_ND, X_ND, S_ALK, S_N2, H2O = asm_vals

            # Step 0: charged component snapshot
            _sno = S_NO
            _snh = S_NH
            _salk = S_ALK
            
            # Step 1: remove any remaining COD demand
            O_coddm = S_O
            NO_coddm = -S_NO*S_NO_i_COD
            cod_spl = S_S + X_S + X_BH + X_BA
            bioN = X_BH*X_BH_i_N + X_BA*X_BA_i_N
            
            if cod_spl <= O_coddm:
                S_O = O_coddm - cod_spl
                S_S = X_S = X_BH = X_BA = 0
            elif cod_spl <= O_coddm + NO_coddm:
                S_O = 0
                S_NO = -(O_coddm + NO_coddm - cod_spl)/S_NO_i_COD
                S_S = X_S = X_BH = X_BA = 0
            else:
                S_S -= O_coddm + NO_coddm
                if S_S < 0:
                    X_S += S_S
                    S_S = 0
                    if X_S < 0:
                        X_BH += X_S                        
                        X_S = 0    
                        if X_BH < 0:
                            X_BA += X_BH
                            X_BH = 0
                S_O = S_NO = 0
            
            # Step 2: convert any readily biodegradable 
            # COD and TKN into amino acids and sugars
            # Assumed S_S, X_S has no nitrogen
            req_scod = S_ND / S_aa_i_N
            if S_S < req_scod:
                S_aa = S_S
                S_su = 0
                S_ND -= S_aa * S_aa_i_N
            else:
                S_aa = req_scod
                S_su = S_S - S_aa
                S_ND = 0

            # Step 3: convert slowly biodegradable COD and TKN
            # into proteins, lipids, and carbohydrates
            req_xcod = X_ND / X_pr_i_N
            if X_S < req_xcod:
                X_pr = X_S
                X_li = X_ch = 0
                X_ND -= X_pr * X_pr_i_N
            else:
                X_pr = req_xcod
                X_li = self.xs_to_li * (X_S - X_pr)
                X_ch = (X_S - X_pr) - X_li
                X_ND = 0
            
            # Step 4: convert active biomass into protein, lipids, 
            # carbohydrates and potentially particulate TKN
            available_bioN = bioN - (X_BH+X_BA) * (1-frac_deg) * adm_X_I_i_N
            if available_bioN < 0:
                raise RuntimeError('Not enough N in X_BA and X_BH to fully convert '
                                   'the non-biodegradable portion into X_I in ADM1.')
            req_bioN = (X_BH+X_BA) * frac_deg * X_pr_i_N
            if available_bioN + X_ND >= req_bioN:
                X_pr += (X_BH+X_BA) * frac_deg
                X_ND += available_bioN - req_bioN
            else:
                bio2pr = (available_bioN + X_ND)/X_pr_i_N
                X_pr += bio2pr
                bio_to_split = (X_BH+X_BA) * frac_deg - bio2pr
                bio_split_to_li = bio_to_split * self.bio_to_li
                X_li += bio_split_to_li
                X_ch += (bio_to_split - bio_split_to_li)
                X_ND = 0
            
            # Step 5: map particulate inerts
            xi_nsp = X_P_i_N * X_P + asm_X_I_i_N * X_I
            xi_ndm = (X_P+X_I) * adm_X_I_i_N
            if xi_nsp + X_ND >= xi_ndm:
                deficit = xi_ndm - xi_nsp
                X_I += X_P + (X_BH+X_BA) * (1-frac_deg)
                X_ND -= deficit
            elif isclose(xi_nsp+X_ND, xi_ndm, rel_tol=rtol, abs_tol=atol):
                X_I += X_P + (X_BH+X_BA) * (1-frac_deg)
                X_ND = 0
            else:
                raise RuntimeError('Not enough N in X_I, X_P, X_ND to fully '
                                   'convert X_I and X_P into X_I in ADM1.')

            req_sn = S_I * S_I_i_N
            if req_sn <= S_ND:
                S_ND -= req_sn
            elif req_sn <= S_ND + X_ND:
                X_ND -= (req_sn - S_ND)
                S_ND = 0
            elif req_sn <= S_ND + X_ND + S_NH:
                S_NH -= (req_sn - S_ND - X_ND)
                S_ND = X_ND = 0
            else:
                warn('Additional soluble inert COD is mapped to S_su.')
                SI_cod = (S_ND + X_ND + S_NH)/S_I_i_N
                S_su += S_I - SI_cod
                S_I = SI_cod
                S_ND = X_ND = S_NH = 0
                
            # Step 6: map any remaining TKN
            S_IN = S_ND + X_ND + S_NH            
            
            # Step 8: check COD and TKN balance
            # has TKN: S_aa, S_IN, S_I, X_pr, X_I
            S_IC = S_cat = S_an = 0
            adm_vals = np.array([
                S_su, S_aa, 
                0, 0, 0, 0, 0, # S_fa, S_va, S_bu, S_pro, S_ac, 
                0, 0, # S_h2, S_ch4,
                S_IC, S_IN, S_I, 
                0, # X_c, 
                X_ch, X_pr, X_li, 
                0, 0, 0, 0, 0, 0, 0, # X_su, X_aa, X_fa, X_c4, X_pro, X_ac, X_h2,
                X_I, S_cat, S_an, H2O])
            
            adm_vals = f_corr(asm_vals, adm_vals)
            
            # Step 7: charge balance
            asm_charge_tot = _snh/14 - _sno/14 - _salk/12
            #!!! charge balance should technically include VFAs, 
            # but VFAs concentrations are assumed zero per previous steps??
            S_IN = adm_vals[adm_ions_idx[0]]
            S_IC = (asm_charge_tot -  S_IN*alpha_IN)/alpha_IC
            net_Scat = asm_charge_tot + proton_charge
            if net_Scat > 0:  
                S_cat = net_Scat
                S_an = 0
            else:
                S_cat = 0
                S_an = -net_Scat
            
            adm_vals[adm_ions_idx[1:]] = [S_IC, S_cat, S_an]
            
            return adm_vals
        
        self._reactions = asm2adm    