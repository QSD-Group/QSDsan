# -*- coding: utf-8 -*-
'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    
    Joy Zhang <joycheung1994@gmail.com>
    
    Yalin Li <mailto.yalin.li@gmail.com>
    
    Saumitra Rai <mailto.raisaumitra9@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

import numpy as np
from warnings import warn
from math import isclose
from biosteam.units import Junction as BSTjunction
from .. import Stream, SanUnit, processes as pc

__all__ = (
    'Junction',
    'ADMjunction', 
    'mADMjunction',
    'ADMtoASM', 
    'ASMtoADM', 
    'ASM2dtoADM1', 
    'ADM1toASM2d',
    'ASM2dtomADM1',
    'mADM1toASM2d',
    'A1junction',
    'ADM1ptomASM2d',
    'mASM2dtoADM1p'
          )

#%% Junction
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
        thermo = downstream.thermo if isinstance(downstream, Stream) else thermo
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
            
            
#%% ADMjunction

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
    # Should be constants
    cod_vfa = np.array([64, 112, 160, 208])
    
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
        return self.pKa_base-np.log10(pc.T_correction_factor(self.T_base, self.T, self.Ka_dH))
    
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
    
    @property
    def alpha_vfa(self):
        '''[float] charge per g of VFA-COD.'''
        return 1.0/self.cod_vfa*(-1.0/(1.0 + 10**(self.pKa[3:]-self.pH)))
    
    def isbalanced(self, lhs, rhs_vals, rhs_i):
        rhs = sum(rhs_vals*rhs_i)
        error = rhs - lhs
        tol = max(self.rtol*lhs, self.rtol*rhs, self.atol)
        return abs(error) <= tol, error, tol, rhs
    
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

#%% mADMjunction
class mADMjunction(ADMjunction):
    '''
    An abstract superclass holding common properties of modified ADM interface classes.
    Users should use its subclasses (e.g., ``ASM2dtomADM1``, ``mADM1toASM2d``) instead.
    
    See Also
    --------
    :class:`qsdsan.sanunits.ADMJunction`
    
    :class:`qsdsan.sanunits.mADM1toASM2d`
    
    :class:`qsdsan.sanunits.ASM2dtomADM1`
    '''
    _parse_reactions = Junction._no_parse_reactions
    rtol = 1e-2
    atol = 1e-6
    cod_vfa = np.array([64, 112, 160, 208])

    def __init__(self, ID='', upstream=None, downstream=(), thermo=None,
                 init_with='WasteStream', F_BM_default=None, isdynamic=False,
                 adm1_model=None, asm2d_model=None):
        self.asm2d_model = asm2d_model
        super().__init__(ID=ID, upstream=upstream, downstream=downstream,
                         thermo=thermo, init_with=init_with, 
                         F_BM_default=F_BM_default, isdynamic=isdynamic,
                         adm1_model=adm1_model)
        
    # @property
    # def T(self):
    #     '''[float] Temperature of the upstream/downstream [K].'''
    #     return self.ins[0].T
    # @T.setter
    # def T(self, T):
    #     self.ins[0].T = self.outs[0].T = T
    
    # @property
    # def pH(self):
    #     '''[float] pH of the upstream/downstream.'''
    #     return self.ins[0].pH
    
    @property
    def asm2d_model(self):
        '''[qsdsan.CompiledProcesses] ASM2d process model.'''
        return self._asm2d_model
    @asm2d_model.setter
    def asm2d_model(self, model):
        if not isinstance(model, (pc.ASM2d, pc.mASM2d)):
            raise ValueError('`asm2d_model` must be an `ASM2d` object, '
                              f'the given object is {type(model).__name__}.')
        self._asm2d_model = model

    @property
    def adm1_model(self):
        '''[qsdsan.CompiledProcesses] mADM1 process model.'''
        return self._adm1_model
    @adm1_model.setter
    def adm1_model(self, model):
        if not isinstance(model, (pc.ADM1_p_extension, pc.ADM1p)):
            raise ValueError('`adm1_model` must be an `ADM1_p_extension` object, '   #!!! update error message
                              f'the given object is {type(model).__name__}.')
        self._adm1_model = model
        
    # @property
    # def T_base(self):
    #     '''[float] Base temperature in the ADM1 model.'''
    #     return self.adm1_model.rate_function.params['T_base']
    
    # @property
    # def pKa_base(self):
    #     '''[float] pKa of the acid-base pairs at the base temperature in the ADM1 model.'''
    #     Ka_base = self.adm1_model.rate_function.params['Ka_base']
    #     return -np.log10(Ka_base)
    
    # @property
    # def Ka_dH(self):
    #     '''[float] Heat of reaction for Ka.'''
    #     return self.adm1_model.rate_function.params['Ka_dH']
    
    @property
    def pKa(self):
        '''
        [numpy.array] pKa array of the following acid-base pairs:
        ('H+', 'OH-'), ('NH4+', 'NH3'), ('H2PO4-', 'HPO4-2'), ('CO2', 'HCO3-'),
        ('HAc', 'Ac-'), ('HPr', 'Pr-'), ('HBu', 'Bu-'), ('HVa', 'Va-')
        '''
        return self.pKa_base-np.log10(pc.T_correction_factor(self.T_base, self.T, self.Ka_dH))

    @property
    def alpha_IN(self):
        '''[float] Charge per g of N.'''
        pH = self.pH
        pKa_IN = self.pKa[1]
        return 10**(pKa_IN-pH)/(1+10**(pKa_IN-pH))/14
    
    @property
    def alpha_IP(self):
        '''[float] Charge per g of P.'''
        pH = self.pH
        pKa_IP = self.pKa[2]
        return (-1/(1+10**(pKa_IP-pH))-1)/31 #!!! alpha IP should be negative
    
    @property
    def alpha_IC(self):
        '''[float] Charge per g of C.'''
        pH = self.pH
        pKa_IC = self.pKa[3]
        return -1/(1+10**(pKa_IC-pH))/12
    
    @property
    def alpha_vfa(self):
        return 1.0/self.cod_vfa*(-1.0/(1.0 + 10**(self.pKa[4:]-self.pH)))
    
    def check_component_properties(self, cmps_asm, cmps_adm):
        get = getattr
        setv = setattr
        for name in ('X_PHA', 'X_PP', 'X_PAO'):
            casm = get(cmps_asm, name)
            cadm = get(cmps_adm, name)
            for attr in ('measured_as', 'i_COD', 'i_C', 'i_N', 'i_P'):
                vasm = get(casm, attr)
                if get(cadm, attr) != vasm:
                    setv(cadm, attr, vasm)
                    warn(f"ADM component {name}'s {attr} is changed to match "
                         "the corresponding ASM component")
        
        for name in ('S_I', 'X_I'):
            casm = get(cmps_asm, name)
            cadm = get(cmps_adm, name)
            for attr in ('measured_as', 'i_C', 'i_N', 'i_P'):
                vadm = get(cadm, attr)
                if get(casm, attr) != vadm:
                    setv(casm, attr, vadm)
                    warn(f"ASM component {name}'s {attr} is changed to match "
                         "the corresponding ADM component")        
        
        for attr in ('i_N', 'i_P'):
            vadm = get(cmps_adm.S_ac, attr)
            if get(cmps_asm.S_A, attr) != vadm:
                cmps_asm.S_A.i_N = vadm
                warn(f"ASM component S_A's {attr} is changed to match "
                     "the ADM component S_ac.") 
        
        if cmps_asm.S_ALK.measured_as != cmps_adm.S_IC.measured_as:
            raise RuntimeError('S_ALK in ASM and S_IC in ADM must both be measured as "C".')
        if cmps_asm.S_NH4.measured_as != cmps_adm.S_IN.measured_as:
            raise RuntimeError('S_NH4 in ASM and S_IN in ADM must both be measured as "N".')
        if cmps_asm.S_PO4.measured_as != cmps_adm.S_IP.measured_as:
            raise RuntimeError('S_PO4 in ASM and S_IP in ADM must both be measured as "P".')
        cmps_asm.refresh_constants()
        cmps_adm.refresh_constants()


#%% ADMtoASM
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
    
    # whether to conserve the nitrogen split between soluble and particulate components
    conserve_particulate_N = False          
   
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
        
        f_corr = self.balance_cod_tkn
        conserve_particulate_N = self.conserve_particulate_N

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
            if conserve_particulate_N:
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
            else:
                xs_ndm = X_S*X_c_i_N  # requires X_S.i_N = 0
                if xs_ndm <= bio_n + S_IN:
                    X_ND = xs_ndm
                    S_IN += bio_n - X_ND
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
            alpha_IN = self.alpha_IN
            alpha_IC = self.alpha_IC
            alpha_vfa = self.alpha_vfa
            S_NH = asm_vals[asm_ions_idx[0]]
            S_ALK = (sum(_ions * np.append([alpha_IN, alpha_IC], alpha_vfa)) - S_NH/14)*(-12)
            asm_vals[asm_ions_idx[1]] = S_ALK
            
            return asm_vals
        
        self._reactions = adm2asm
        
#%% ASMtoADM
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
    
    def __init__(self, ID='', upstream=None, downstream=(), thermo=None,
                 init_with='WasteStream', F_BM_default=None, isdynamic=False,
                 adm1_model=None, T=298.15, pH=7):
        self._T = T
        self._pH = pH
        super().__init__(ID=ID, upstream=upstream, downstream=downstream,
                         thermo=thermo, init_with=init_with, 
                         F_BM_default=F_BM_default, isdynamic=isdynamic,
                         adm1_model=adm1_model)
    
    @property
    def T(self):
        '''[float] Temperature of the downstream [K].'''
        try: return self.outs[0].sink.T
        except: return self._T
    @T.setter
    def T(self, T):
        self._T = self.outs[0].T = T
    
    @property
    def pH(self):
        '''[float] downstream pH.'''
        if self._pH: return self._pH
        else:
            try: return self.outs[0].sink.outs[1].pH
            except: return 7.
    @pH.setter
    def pH(self, ph):
        self._pH = self.outs[0].pH = ph
    
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
            alpha_IN = self.alpha_IN
            alpha_IC = self.alpha_IC
            proton_charge = 10**(-self.pKa[0]+self.pH) - 10**(-self.pH) # self.pKa[0] is pKw
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
        
#%% ASM2dtoADM1

class ASM2dtoADM1(ADMjunction):
    '''
    Interface unit to convert activated sludge model No. (ASM2d) components
    to original anaerobic digestion model (ADM1) components.
    
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
    asm_X_I_i_N = 0.06
    
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
        non_tkn_idx = cmps_asm.indices(('S_N2', 'S_NO3'))
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
        
        S_NO3_i_COD = cmps_asm.S_NO3.i_COD
        X_H_i_N = cmps_asm.X_H.i_N
        X_AUT_i_N = cmps_asm.X_AUT.i_N
        X_PAO_i_N = cmps_asm.X_PAO.i_N
        S_F_i_N = cmps_asm.S_F.i_N
        X_S_i_N = cmps_asm.X_S.i_N
        asm_S_I_i_N = cmps_asm.S_I.i_N
        
        if self.asm_X_I_i_N == None:
            asm_X_I_i_N = cmps_asm.X_I.i_N
        else:
            asm_X_I_i_N = self.asm_X_I_i_N
            
        if cmps_asm.S_A.i_N > 0: 
            warn(f'S_A in ASM2d has positive nitrogen content: {cmps_asm.S_S.i_N} gN/gCOD. '
                 'These nitrogen will be ignored by the interface model '
                 'and could lead to imbalance of TKN after conversion.')
        
        cmps_adm = outs.components
        S_aa_i_N = cmps_adm.S_aa.i_N
        X_pr_i_N = cmps_adm.X_pr.i_N
        adm_S_I_i_N = cmps_adm.S_I.i_N
        adm_X_I_i_N = cmps_adm.X_I.i_N
        
        
        adm_ions_idx = cmps_adm.indices(['S_IN', 'S_IC', 'S_cat', 'S_an'])
        
        frac_deg = self.frac_deg
        alpha_IN = self.alpha_IN
        alpha_IC = self.alpha_IC
        proton_charge = 10**(-self.pKa[0]+self.pH) - 10**(-self.pH) # self.pKa[0] is pKw
        f_corr = self.balance_cod_tkn

        def asm2adm(asm_vals):
            
            S_O2, S_N2, S_NH4, S_NO3, S_PO4, S_F, S_A, S_I, S_ALK, X_I, X_S, X_H, \
                X_PAO, X_PP, X_PHA, X_AUT, X_MeOH, X_MeP, H2O = asm_vals

            # Step 0: charged component snapshot
            _sa = S_A
            _snh4 = S_NH4
            _sno3 = S_NO3
            _spo4 = S_PO4
            _salk = S_ALK
            _xpp = X_PP
            
            # Step 1: remove any remaining COD demand
            O2_coddm = S_O2
            NO3_coddm = -S_NO3*S_NO3_i_COD
            
            cod_spl = (S_A + S_F) + (X_S + X_PHA) + (X_H + X_AUT + X_PAO)
            
            
            bioN = X_H*X_H_i_N + X_AUT*X_AUT_i_N + X_PAO*X_PAO_i_N
            # To be used in Step 2
            S_F_N = S_F*S_F_i_N   #S_ND (in asm1) equals the N content in S_F
            # To be used in Step 3
            X_S_N = X_S*X_S_i_N   #X_ND (in asm1) equals the N content in X_S
            
            
            if cod_spl <= O2_coddm:
                S_O2 = O2_coddm - cod_spl
                S_F = S_A =  X_S = X_H = X_AUT = 0
            elif cod_spl <= O2_coddm + NO3_coddm:
                S_O2 = 0
                S_NO3 = -(O2_coddm + NO3_coddm - cod_spl)/S_NO3_i_COD
                S_A = S_F = X_S = X_H = X_AUT = 0
            else:
                S_A -= O2_coddm + NO3_coddm
                if S_A < 0:
                    S_F += S_A
                    S_A = 0
                    if S_F < 0:
                        X_S += S_F
                        S_F = 0
                        if X_S < 0:
                            X_PHA += X_S
                            X_S = 0
                            if X_PHA < 0:
                                X_H += X_PHA
                                X_PHA = 0
                                if X_H < 0:
                                    X_AUT += X_H
                                    X_H = 0
                                    if X_AUT < 0:
                                        X_PAO += X_AUT
                                        X_AUT = 0
                            
                S_O2 = S_NO3 = 0

            # Step 2: convert any readily biodegradable 
            # COD and TKN into amino acids and sugars
            
            # S_S (in asm1) equals to the sum of S_F and S_A (pg. 82 IWA ASM models handbook)
            S_S_asm1 = S_F + S_A 
            
            # First we calculate the amount of amino acid required in ADM1
            # if all available soluble organic N can be mapped to amino acid
            req_scod = S_F_N / S_aa_i_N
            
            # if available S_S is not enough to fulfill that amino acid requirement 
            if S_S_asm1 < req_scod: 
                # then all available S_S is mapped to amino acids 
                S_aa = S_S_asm1
                # and no S_S would be available for conversion to sugars
                S_su = 0
                # This needs to be followed by a corresponding loss in soluble organic N 
                S_F_N -= S_aa * S_aa_i_N
            # if available S_S is more than enough to fulfill that amino acid requirement 
            else:
                # All soluble organic N will be mapped to amino acid
                S_aa = req_scod
                # The line above implies that a certain portion of S_S would also be consumed to form amino acid
                # The S_S which is left would form sugar 
                # In simpler terms; S_S = S_S - S_aa; S_su = S_S 
                S_su = S_S_asm1 - S_aa
                # All soluble organic N would thus be consumed in amino acid formation
                S_F_N = 0
                
                
            # Step 3: convert slowly biodegradable COD and TKN
            # into proteins, lipids, and carbohydrates
            
            # First we calculate the amount of protein required in ADM1
            # if all available particulate organic N can be mapped to protein
            req_xcod = X_S_N / X_pr_i_N
            
            # if available X_S is not enough to fulfill that protein requirement
            if X_S < req_xcod:
                # then all available X_S is mapped to amino acids
                X_pr = X_S
                # and no X_S would be available for conversion to lipid or carbohydrates 
                
                X_li = self.xs_to_li * X_PHA
                X_ch = (1 - self.xs_to_li)*X_PHA 
               
                # This needs to be followed by a corresponding loss in particulate organic N 
                X_S_N -= X_pr * X_pr_i_N
                
            # if available X_S is more than enough to fulfill that protein requirement
            else:
                # All particulate organic N will be mapped to amino acid
                X_pr = req_xcod
                # The line above implies that a certain portion of X_S would also be consumed to form protein
                # The X_S which is left would form lipid and carbohydrates in a percentage define by the user  
                X_li = self.xs_to_li * (X_S + X_PHA - X_pr)
                X_ch = (X_S + X_PHA - X_pr) - X_li
                # All particulate organic N would thus be consumed in amino acid formation
                X_S_N = 0
            
            # Step 4: convert active biomass into protein, lipids, 
            # carbohydrates and potentially particulate TKN
            
            # First the amount of biomass N available for protein, lipid etc is determined
            # For this calculation, from total biomass N available the amount 
            # of particulate inert N expected in ADM1 is subtracted 
            
            available_bioN = bioN - (X_H + X_AUT + X_PAO) * (1-frac_deg) * adm_X_I_i_N
            
            if available_bioN < 0:
                raise RuntimeError('Not enough N in X_H, X_AUT and X_PAO to fully convert '
                                   'the non-biodegradable portion into X_I in ADM1.')
                
            # Then the amount of biomass N required for biomass conversion to protein is determined
            req_bioN = (X_H + X_AUT + X_PAO) * frac_deg * X_pr_i_N
            # req_bioP = (X_H + X_AUT) * frac_deg * X_pr_i_P
            
            # If available biomass N and particulate organic N is greater than 
            # required biomass N for conversion to protein
            if available_bioN + X_S_N >= req_bioN:
                # then all biodegradable biomass N (corrsponding to protein demand) is converted to protein
                X_pr += (X_H + X_AUT + X_PAO) * frac_deg
                # the remaining biomass N is transfered as organic N
                X_S_N += available_bioN - req_bioN 
            else:
                # all available N and particulate organic N is converted to protein
                bio2pr = (available_bioN + X_S_N)/X_pr_i_N
                X_pr += bio2pr
                # Biodegradable biomass available after conversion to protein is calculated 
                bio_to_split = (X_H + X_AUT + X_PAO) * frac_deg - bio2pr
                # Part of the remaining biomass is mapped to lipid based on user defined value 
                bio_split_to_li = bio_to_split * self.bio_to_li
                X_li += bio_split_to_li
                # The other portion of the remanining biomass is mapped to carbohydrates 
                X_ch += (bio_to_split - bio_split_to_li)
                # Since all organic N has been mapped to protein, none is left
                X_S_N = 0
            
            # Step 5: map particulate inerts
            
            # 5 (a)
            # First determine the amount of particulate inert N available from ASM2d
            xi_nsp_asm2d = X_I * asm_X_I_i_N
            
            # Then determine the amount of particulate inert N that could be produced 
            # in ADM1 given the ASM1 X_I
            xi_ndm = X_I * adm_X_I_i_N

            # if particulate inert N available in ASM1 is greater than ADM1 demand
            if xi_nsp_asm2d + X_S_N >= xi_ndm:
                deficit = xi_ndm - xi_nsp_asm2d
                # COD balance 
                X_I += (X_H + X_AUT + X_PAO) * (1-frac_deg)
                # N balance 
                X_S_N -= deficit
            elif isclose(xi_nsp_asm2d+X_S_N, xi_ndm, rel_tol=rtol, abs_tol=atol):
                # COD balance 
                X_I += (X_H + X_AUT + X_PAO) * (1-frac_deg)
                # N balance 
                X_S_N = 0
            else:
                raise RuntimeError('Not enough N in X_I, X_S to fully '
                                    'convert X_I in ASM2d into X_I in ADM1.')
                
            # 5(b)
            
            # Then determine the amount of soluble inert N that could be produced 
            # in ADM1 given the ASM1 X_I
            req_sn = S_I * adm_S_I_i_N
            supply_inert_n_asm2d = S_I * asm_S_I_i_N
            
            # N balance 
            if req_sn <= S_F_N + supply_inert_n_asm2d:
                S_F_N -= (req_sn - supply_inert_n_asm2d)
                supply_inert_n_asm2d = 0 
            # N balance
            elif req_sn <= S_F_N + X_S_N + supply_inert_n_asm2d:
                X_S_N -= (req_sn - S_F_N - supply_inert_n_asm2d)
                S_F_N = supply_inert_n_asm2d = 0
            # N balance
            elif req_sn <= S_F_N + X_S_N + S_NH4 + supply_inert_n_asm2d:
                S_NH4 -= (req_sn - S_F_N - X_S_N - supply_inert_n_asm2d)
                S_F_N = X_S_N = supply_inert_n_asm2d = 0
            else:
                warn('Additional soluble inert COD is mapped to S_su.')
                SI_cod = (S_F_N + X_S_N + S_NH4 + supply_inert_n_asm2d)/adm_S_I_i_N
                S_su += S_I - SI_cod
                S_I = SI_cod
                S_F_N = X_S_N = S_NH4 = supply_inert_n_asm2d = 0
            
            # Step 6: Step map any remaining TKN/P
            S_IN = S_F_N + X_S_N + S_NH4 + supply_inert_n_asm2d
            
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
            asm_charge_tot = - _sa/64 + _snh4/14 - _sno3/14 - 1.5*_spo4/31 - _salk - _xpp/31 #Based on page 84 of IWA ASM handbook
            S_IN = adm_vals[adm_ions_idx[0]]            
            S_IC = (asm_charge_tot -S_IN*alpha_IN)/alpha_IC
            
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
        
#%% ADM1toASM2d

class ADM1toASM2d(ADMjunction):
    '''
    Interface unit to convert anaerobic digestion model no. 1 (ADM1) components
    to activated sludge model no. 2 (ASM2d) components.
    
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
    bio_to_xs = 0.9
    
    # Should be constants
    cod_vfa = np.array([64, 112, 160, 208])
    
    # whether to conserve the nitrogen split between soluble and particulate components
    conserve_particulate_N = False
    
    
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
        
        X_S_i_N = cmps_asm.X_S.i_N
        S_F_i_N = cmps_asm.S_F.i_N
        
        asm_X_I_i_N = cmps_asm.X_I.i_N
        asm_S_I_i_N = cmps_asm.S_I.i_N
       
        asm_ions_idx = cmps_asm.indices(('S_A', 'S_NH4', 'S_NO3', 'S_PO4','S_ALK', 'X_PP'))
        
        alpha_IN = self.alpha_IN
        alpha_IC = self.alpha_IC
        alpha_vfa = self.alpha_vfa
        f_corr = self.balance_cod_tkn
        conserve_particulate_N = self.conserve_particulate_N

        def adm2asm(adm_vals):    
                
            S_su, S_aa, S_fa, S_va, S_bu, S_pro, S_ac, S_h2, S_ch4, S_IC, S_IN, S_I, X_c, \
            X_ch, X_pr, X_li, X_su, X_aa, X_fa, X_c4, X_pro, X_ac, X_h2, X_I, S_cat, S_an, H2O = adm_vals
                       
            # Step 0: snapshot of charged components
            _ions = np.array([S_IN, S_IC, S_ac, S_pro, S_bu, S_va])
            
            # Step 1a: convert biomass into X_S and X_I
            
            bio_cod = X_su + X_aa + X_fa + X_c4 + X_pro + X_ac + X_h2
            bio_n = sum((adm_vals*adm_i_N)[adm_bio_N_indices])
            
            #!!! In default ASM2d stoichiometry, biomass decay (cell lysis)
            #!!! yields 90% particulate substrate + 10% X_I
            #!!! so: convert both biomass and X_I in adm to X_S and X_I in asm
            xi_n = X_I*adm_X_I_i_N
            xs_cod = bio_cod * self.bio_to_xs
            xs_ndm = xs_cod * X_S_i_N
            
            xi_cod = bio_cod * (1-self.bio_to_xs) + X_I
            xi_ndm = xi_cod * asm_X_I_i_N
            
            if xs_ndm > bio_n:
                warn('Not enough biomass N to map the specified proportion of '
                     'biomass COD into X_S. Rest of the biomass COD goes to S_A')
                X_S = bio_n / X_S_i_N
                xs_cod -= X_S
                bio_n = 0
            else:
                X_S = xs_cod
                xs_cod = 0
                bio_n -= xs_ndm
                
            if xi_ndm > bio_n + xi_n + S_IN:
                warn('Not enough N in biomass and X_I to map the specified proportion of '
                     'biomass COD into X_I. Rest of the biomass COD goes to S_A')
                X_I = (bio_n + xi_n + S_IN) / asm_X_I_i_N
                xi_cod -= X_I
                bio_n = xi_n = S_IN = 0
            else:
                X_I = xi_cod
                xi_cod = 0
                xi_n -= xi_ndm
                if xi_n < 0:
                    bio_n += xi_n
                    xi_n = 0
                    if bio_n < 0:
                        S_IN += bio_n
                        bio_n = 0
        
            xsub_cod = X_c + X_ch + X_pr + X_li 
            xsub_n = X_c*X_c_i_N + X_pr*X_pr_i_N
            
            xs_ndm = xsub_cod * X_S_i_N
            
            if xs_ndm > xsub_n + bio_n:
                X_S_temp = (xsub_n + bio_n)/X_S_i_N
                X_S += X_S_temp
                xsub_cod -= X_S_temp
                xsub_n = bio_n = 0
            else:
                X_S += xsub_cod
                xsub_cod = 0
                xsub_n -= xs_ndm
                if xsub_n < 0: 
                    bio_n += xsub_n
                    xsub_n = 0
                    
            ssub_cod = S_su + S_aa + S_fa
            ssub_n = S_aa * S_aa_i_N
            sf_ndm = ssub_cod * S_F_i_N
            
            if sf_ndm > ssub_n + xsub_n + bio_n:
                S_F = (ssub_n + xsub_n + bio_n) / S_F_i_N
                ssub_cod -= S_F
                ssub_n = xsub_n = bio_n = 0
            else:
                S_F = ssub_cod
                ssub_cod = 0
                ssub_n -= sf_ndm
                if ssub_n < 0:
                    xsub_n += ssub_n
                    ssub_n = 0
                    if xsub_n < 0:
                        bio_n += xsub_n
                        xsub_n = 0
            
            S_A = S_ac + S_pro + S_bu + S_va
            
            si_cod = S_I    
            si_n = S_I * adm_S_I_i_N
            si_ndm = si_cod * asm_S_I_i_N
            if si_ndm > si_n + xi_n + S_IN:
                warn('Not enough N in S_I and X_I to map all S_I from ADM1 to ASM2d. '
                     'Rest of the S_I COD goes to S_A')
                S_I = (si_n + xi_n + S_IN) / asm_S_I_i_N
                si_cod -= S_I
                si_n = xi_n = S_IN = 0
            else:
                S_I = si_cod
                si_cod = 0
                si_n -= si_ndm
                if si_n < 0:
                    xi_n += si_n
                    si_n = 0
                    if xi_n < 0:
                        S_IN += xi_n
                        xi_n = 0
            
           
            S_NH4 = S_IN + si_n + ssub_n + xsub_n + xi_n + bio_n
            S_A += si_cod + ssub_cod + xsub_cod + xi_cod + xs_cod
            S_ALK = S_IC
                
            # Step 6: check COD and TKN balance
            asm_vals = np.array(([
                0, 0, # S_O2, S_N2,
                S_NH4, 
                0, 
                0, 
                S_F, S_A, S_I, S_ALK,
                X_I, X_S, 
                0,  # X_H,
                0, 0, 0, 
                0, # X_AUT,
                0, 0, H2O]))
            
            if S_h2 > 0 or S_ch4 > 0:
                warn('Ignored dissolved H2 or CH4.')
            
            asm_vals = f_corr(adm_vals, asm_vals)
            
            # Step 5: charge balance for alkalinity
            
            # asm_ions_idx = cmps_asm.indices(('S_A', 'S_NH4', 'S_NO3', 'S_PO4','S_ALK', 'X_PP'))
            
            _sa = asm_vals[asm_ions_idx[0]]
            _snh4 = asm_vals[asm_ions_idx[1]]
            _sno3 = asm_vals[asm_ions_idx[2]]
            _spo4 = asm_vals[asm_ions_idx[3]]
            _xpp = asm_vals[asm_ions_idx[5]]
            
            S_ALK = (sum(_ions * np.append([alpha_IN, alpha_IC], alpha_vfa)) - \
                     (- _sa/64 + _snh4/14 - _sno3/14 - 1.5*_spo4/31 - _xpp/31))*(-12)
            
            asm_vals[asm_ions_idx[4]] = S_ALK
            
            return asm_vals
        
        self._reactions = adm2asm
    
    @property
    def alpha_vfa(self):
        return 1.0/self.cod_vfa*(-1.0/(1.0 + 10**(self.pKa[3:]-self.pH)))
    
#%% mADM1toASM2d

class mADM1toASM2d(mADMjunction):
    '''
    Interface unit to convert modified anaerobic digestion model no. 1 (ADM1) components
    to activated sludge model no. 2d (ASM2d) components.
    
    Parameters
    ----------
    upstream : stream or str
        Influent stream with ADM components.
    downstream : stream or str
        Effluent stream with ASM components.
    adm1_model : obj
        The anaerobic digestion process model (:class:`qsdsan.processes.ADM1_p_extension`).
    rtol : float
        Relative tolerance for COD and TKN balance.
    atol : float
        Absolute tolerance for COD and TKN balance.
    
    References
    ----------
    [1] Nopens, I.; Batstone, D. J.; Copp, J. B.; Jeppsson, U.; Volcke, E.; 
    Alex, J.; Vanrolleghem, P. A. An ASM/ADM Model Interface for Dynamic 
    Plant-Wide Simulation. Water Res. 2009, 43, 1913–1923.
    
    [2] Flores-Alsina, X., Solon, K., Kazadi Mbamba, C., Tait, S., Gernaey, K. V., 
    Jeppsson, U., & Batstone, D. J. (2016). Modelling phosphorus (P), sulfur (S) 
    and iron (FE) interactions for dynamic simulations of anaerobic digestion processes. 
    Water Research, 95, 370–382.
    
    See Also
    --------
    :class:`qsdsan.sanunits.mADMjunction`
    
    :class:`qsdsan.sanunits.ASM2dtomADM1`
    '''
    
    # User defined values
    bio_to_xs = 0.9
    
    # Since we are matching PAOs directly from ASM2d to mADM1, it is important 
    # for PAOs to have identical N/P content across models
    # adm_X_PAO_i_N = 0.07 
    # adm_X_PAO_i_P = 0.02
    
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
    
    def balance_cod_tkn(self, adm_vals, asm_vals):
        cmps_adm = self.ins[0].components
        cmps_asm = self.outs[0].components
        asm_i_COD = cmps_asm.i_COD
        adm_i_COD = cmps_adm.i_COD
        gas_idx = cmps_adm.indices(('S_h2', 'S_ch4'))
        
        asm_i_N = cmps_asm.i_N
        adm_i_N = cmps_adm.i_N    
        asm_i_P = cmps_asm.i_P
        adm_i_P = cmps_adm.i_P
        adm_cod = sum(adm_vals*adm_i_COD) - sum(adm_vals[gas_idx])
        adm_tkn = sum(adm_vals*adm_i_N)
        adm_tp = sum(adm_vals*adm_i_P)
        
        cod_bl, cod_err, cod_tol, asm_cod = self.isbalanced(adm_cod, asm_vals, asm_i_COD)
        tkn_bl, tkn_err, tkn_tol, asm_tkn = self.isbalanced(adm_tkn, asm_vals, asm_i_N)
        tp_bl, tp_err, tp_tol, asm_tp = self.isbalanced(adm_tp, asm_vals, asm_i_P)
        
        if tkn_bl and tp_bl:
            if cod_bl:
                return asm_vals
            else:
                print('COD not balanced')
                breakpoint()
                if cod_err > 0: dcod = -(cod_err - cod_tol)/asm_cod
                else: dcod = -(cod_err + cod_tol)/asm_cod
                _asm_vals = asm_vals * (1 + (asm_i_COD>0)*dcod)
                _tkn_bl, _tkn_err, _tkn_tol, _asm_tkn = self.isbalanced(adm_tkn, _asm_vals, asm_i_N)
                _tp_bl, _tp_err, _tp_tol, _asm_tp = self.isbalanced(adm_tp, _asm_vals, asm_i_P)
                if _tkn_bl and _tp_bl: return _asm_vals
                else: 
                    warn('cannot balance COD, TKN, and TP at the same \n'
                        f'time with rtol={self.rtol} and atol={self.atol}.\n '
                        f'influent (ADM) TKN is {adm_tkn}\n '
                        f'effluent (ASM) TKN is {asm_tkn} or {_asm_tkn}\n '
                        f'influent (ADM) TP is {adm_tp}\n ' 
                        f'effluent (ASM) TP is {asm_tp} or {_asm_tp}. '
                        f'influent (ADM) COD is {adm_cod}\n ' 
                        f'effluent (ASM) COD is {asm_cod} or {asm_cod*(1+dcod)}. ')
                    return asm_vals
        elif cod_bl and tp_bl:
            if tkn_bl:
                return asm_vals
            else:
                print('TKN not balanced')
                breakpoint()
                if tkn_err > 0: dtkn = -(tkn_err - tkn_tol)/asm_tkn
                else: dtkn = -(tkn_err + tkn_tol)/asm_tkn
                _asm_vals = asm_vals * (1 + (asm_i_N>0)*dtkn)
                _cod_bl, _cod_err, _cod_tol, _asm_cod = self.isbalanced(adm_cod, _asm_vals, asm_i_COD)
                _tp_bl, _tp_err, _tp_tol, _asm_tp = self.isbalanced(adm_tp, _asm_vals, asm_i_P)
                if _cod_bl and _tp_bl: return _asm_vals
                else: 
                    warn('cannot balance COD, TKN, and TP at the same time'
                        f'time with rtol={self.rtol} and atol={self.atol}.\n '
                        f'influent (ADM) COD is {adm_cod}\n '
                        f'effluent (ASM) COD is {asm_cod} or {_asm_cod}\n '
                        f'influent (ADM) TP is {adm_tp}\n ' 
                        f'effluent (ASM) TP is {asm_tp} or {_asm_tp}. '
                        f'influent (ADM) TKN is {adm_tkn}\n ' 
                        f'effluent (ASM) TKN is {asm_tkn} or {asm_tkn*(1+dtkn)}. ')
                    return asm_vals
        elif cod_bl and tkn_bl:
            if tp_bl:
                return asm_vals
            else:
                print('TP not balanced')
                breakpoint()
                if tp_err > 0: dtp = -(tp_err - tp_tol)/asm_tp
                else: dtp = -(tp_err + tp_tol)/asm_tp
                _asm_vals = asm_vals * (1 + (asm_i_P>0)*dtp)
                _cod_bl, _cod_err, _cod_tol, _asm_cod = self.isbalanced(adm_cod, _asm_vals, asm_i_COD)
                _tkn_bl, _tkn_err, _tkn_tol, _asm_tkn = self.isbalanced(adm_tkn, _asm_vals, asm_i_N)
                if _cod_bl and _tkn_bl: return _asm_vals
                else: 
                    warn('cannot balance COD, TKN, and TP at the same time'
                        f'time with rtol={self.rtol} and atol={self.atol}.\n '
                        f'influent (ADM) COD is {adm_cod}\n '
                        f'effluent (ASM) COD is {asm_cod} or {_asm_cod}\n '
                        f'influent (ADM) TKN is {adm_tkn}\n ' 
                        f'effluent (ASM) TKN is {asm_tkn} or {_asm_tkn}. '
                        f'influent (ADM) TP is {adm_tp}\n ' 
                        f'effluent (ASM) TP is {asm_tp} or {asm_tp*(1+dtp)}. ')
                    return asm_vals
        else:
            print('At least two of COD, TKN, and TP not balanced')
            breakpoint()
            warn('cannot balance COD, TKN and TP at the same time. \n'
                 'Atleast two of the three COD, TKN, and TP are not balanced \n'
                f'time with rtol={self.rtol} and atol={self.atol}.\n '
                f'influent (ADM) COD is {adm_cod}\n '
                f'effluent (ASM) COD is {asm_cod}\n '
                f'influent (ADM) TP is {adm_tp}\n ' 
                f'effluent (ASM) TP is {asm_tp}'
                f'influent (ADM) TKN is {adm_tkn}\n ' 
                f'effluent (ASM) TKN is {asm_tkn}. ')
            return asm_vals
    
    def _compile_reactions(self):
        # Retrieve constants
        # rtol = self.rtol
        # atol = self.atol
        
        cmps_adm = self.ins[0].components
        cmps_asm = self.outs[0].components
        self.check_component_properties(cmps_asm, cmps_adm)
        
        # N balance 
        X_pr_i_N = cmps_adm.X_pr.i_N
        S_aa_i_N = cmps_adm.S_aa.i_N
        # adm_X_I_i_N = cmps_adm.X_I.i_N
        # adm_S_I_i_N = cmps_adm.S_I.i_N
        adm_i_N = cmps_adm.i_N
        adm_bio_indices = cmps_adm.indices(('X_su', 'X_aa', 'X_fa', 
                                            'X_c4', 'X_pro', 'X_ac', 'X_h2'))
        
        # P balance
        X_pr_i_P = cmps_adm.X_pr.i_P
        # adm_X_I_i_P = cmps_adm.X_I.i_P
        # adm_S_I_i_P = cmps_adm.S_I.i_P
        S_aa_i_P = cmps_adm.S_aa.i_P
        adm_i_P = cmps_adm.i_P
        
        # N balance 
        X_S_i_N = cmps_asm.X_S.i_N
        S_F_i_N = cmps_asm.S_F.i_N
        # S_A_i_N = cmps_asm.S_A.i_N
            
        asm_X_I_i_N = cmps_asm.X_I.i_N
        # asm_S_I_i_N = cmps_asm.S_I.i_N
        asm_ions_idx = cmps_asm.indices(('S_NH4', 'S_A', 'S_NO3', 'S_PO4', 'S_ALK'))
        
        # P balance 
        X_S_i_P = cmps_asm.X_S.i_P
        S_F_i_P = cmps_asm.S_F.i_P
        # S_A_i_P = cmps_asm.S_A.i_P
        asm_X_I_i_P = cmps_asm.X_I.i_P 
        # asm_S_I_i_P = cmps_asm.S_I.i_P
        
        #!!! All checks are now done in `check_component_property`
        # Checks for direct mapping of X_PAO, X_PP, X_PHA
        
        # # Check for X_PAO (measured as COD so i_COD = 1 in both ASM2d and ADM1)
        
        # asm_X_PAO_i_N = cmps_asm.X_PAO.i_N
        # adm_X_PAO_i_N = cmps_adm.X_PAO.i_N
        
        # if self.adm_X_PAO_i_N == None:
        #     adm_X_PAO_i_N = cmps_adm.X_PAO.i_N
        # else:
        #     adm_X_PAO_i_N = self.adm_X_PAO_i_N
        
        # if asm_X_PAO_i_N != adm_X_PAO_i_N:
        #     raise RuntimeError('X_PAO cannot be directly mapped as N content'
        #                        f'in asm2d_X_PAO_i_N = {asm_X_PAO_i_N} is not equal to'
        #                        f'adm_X_PAO_i_N = {adm_X_PAO_i_N}')
            
        # asm_X_PAO_i_P = cmps_asm.X_PAO.i_P
        # adm_X_PAO_i_P = cmps_adm.X_PAO.i_P
        
        # if self.adm_X_PAO_i_P == None:
        #     adm_X_PAO_i_P = cmps_adm.X_PAO.i_P
        # else:
        #     adm_X_PAO_i_P = self.adm_X_PAO_i_P
        
        # if asm_X_PAO_i_P != adm_X_PAO_i_P:
        #     raise RuntimeError('X_PAO cannot be directly mapped as P content'
        #                        f'in asm2d_X_PAO_i_P = {asm_X_PAO_i_P} is not equal to'
        #                        f'adm_X_PAO_i_P = {adm_X_PAO_i_P}')
        
        # # Checks not required for X_PP as measured as P in both, with i_COD = i_N = 0
        # # Checks not required for X_PHA as measured as COD in both, with i_N = i_P = 0
        
        #!!! alpha must be called within the function.
        # alpha_IN = self.alpha_IN
        # alpha_IC = self.alpha_IC
        # alpha_IP = self.alpha_IP
        # alpha_vfa = self.alpha_vfa
        f_corr = self.balance_cod_tkn

        # To convert components from mADM1 to ASM2d (madm1-2-asm2d)
        def madm12asm2d(adm_vals):    
            
            S_su, S_aa, S_fa, S_va, S_bu, S_pro, S_ac, S_h2, S_ch4, S_IC, S_IN, S_IP, S_I, \
                X_ch, X_pr, X_li, X_su, X_aa, X_fa, X_c4, X_pro, X_ac, X_h2, X_I, \
                X_PHA, X_PP, X_PAO, S_K, S_Mg, X_MeOH, X_MeP, S_cat, S_an, H2O = adm_vals
                       
            # Step 0: snapshot of charged components
            # Not sure about charge on X_PP, S_Mg, S_K (PHA and PAO would have zero charge)
            # X_PP in ADM1 is charge neutral
            _ions = np.array([S_IN, S_IC, S_IP, S_Mg, S_K, S_ac, S_pro, S_bu, S_va])
            
            # Step 1a: convert biomass and inert particulates into X_S and X_I 
            
            # What is available
            bio_cod = X_su + X_aa + X_fa + X_c4 + X_pro + X_ac + X_h2
            bio_n = sum((adm_vals*adm_i_N)[adm_bio_indices])
            bio_p = sum((adm_vals*adm_i_P)[adm_bio_indices])
            
            
            #!!! In default ASM2d stoichiometry, biomass decay (cell lysis)
            #!!! yields 90% particulate substrate + 10% X_I
            #!!! so: convert both biomass and X_I in adm to X_S and X_I in asm
            
            # ----------------------------Rai version---------------------------------------------
            # What is available
            # xi_n = X_I*adm_X_I_i_N
            # xi_p = X_I*adm_X_I_i_P
            # ----------------------------Rai version---------------------------------------------
        
            # What would be formed by X_S
            xs_cod = bio_cod * self.bio_to_xs
            xs_ndm = xs_cod * X_S_i_N
            xs_pdm = xs_cod * X_S_i_P
            
            # What would be formed by X_I (ASM2d)
            # xi_cod = bio_cod * (1 - self.bio_to_xs) + X_I
            X_I += bio_cod * (1 - self.bio_to_xs)
            
            # ----------------------------Rai version---------------------------------------------
            # xi_ndm = xi_cod * asm_X_I_i_N
            # xi_pdm = xi_cod * asm_X_I_i_P
            # ----------------------------Rai version---------------------------------------------
            
            # MAPPING OF X_S
            
            # Case I: Both bio_N and bio_P are sufficient 
            if xs_ndm <= bio_n and xs_pdm <= bio_p:
                X_S = xs_cod
                xs_cod = 0
                bio_n -= xs_ndm
                bio_p -= xs_pdm
            else:
            # Case II, III, and, IV: At least one of the two biological N/P is not sufficient
                if bio_p / X_S_i_P > bio_n / X_S_i_N:
                    warn('Not enough biomass N to map the specified proportion of '
                         'biomass COD into X_S. Rest of the biomass COD goes to S_A in last step')
                    X_S = bio_n / X_S_i_N
                    xs_cod -= X_S
                    bio_n = 0
                    bio_p -= X_S*X_S_i_P #mathematically, bio_p can become negative at this point
                    if bio_p < 0:
                        S_IP += bio_p
                        bio_p = 0
                else:
                    warn('Not enough biomass P to map the specified proportion of '
                         'biomass COD into X_S. Rest of the biomass COD goes to S_A in last step')
                    X_S = bio_p / X_S_i_P
                    xs_cod -= X_S
                    bio_p = 0
                    bio_n -= X_S*X_S_i_N #mathematically, bio_n can become negative at this point
                    if bio_n < 0:
                        S_IN += bio_n
                        bio_n = 0
               
            # Step 2: MAPPING OF X_I
            
            # Flores Alsina
            
            #!!! equality is enforced in `check_component_properties`
            # if asm_X_I_i_N == adm_X_I_i_N and asm_X_I_i_P == adm_X_I_i_P:
            #     X_I = xi_cod
            excess_N = (bio_cod * (1 - self.bio_to_xs))*asm_X_I_i_N
            excess_P = (bio_cod * (1 - self.bio_to_xs))*asm_X_I_i_P
            # else:
            #     raise RuntimeError('N and P content in X_I should be same')
            
            # ----------------------------Rai version---------------------------------------------
            # if xi_ndm < bio_n + xi_n + S_IN and xi_pdm < bio_p + xi_p + S_IP:
                
            #     X_I = xi_cod
            #     xi_cod = 0
                
            #     xi_n -= xi_ndm
            #     if xi_n < 0:
            #         bio_n += xi_n
            #         xi_n = 0
            #         if bio_n < 0:
            #             S_IN += bio_n
            #             bio_n = 0
                
            #     xi_p -= xi_pdm
            #     if xi_p < 0:
            #         bio_p += xi_p
            #         xi_p = 0
            #         if bio_p < 0:
            #             S_IP += bio_p
            #             bio_p = 0
                        
            # else:
            #     if (bio_p + xi_p + S_IP) / asm_X_I_i_P  >  (bio_n + xi_n + S_IN) / asm_X_I_i_N:
                    
            #         warn('Not enough N in biomass and X_I to map the specified proportion of'
            #              'biomass COD into X_I. Rest of the biomass COD goes to S_A')
            #         X_I =  (bio_n + xi_n + S_IN) / asm_X_I_i_N
            #         xi_cod -= X_I
                    
            #         bio_n = xi_n = S_IN = 0
                    
            #         xi_p -= X_I*asm_X_I_i_P
            #         if xi_p < 0:
            #             bio_p += xi_p
            #             xi_p = 0
            #             if bio_p < 0:
            #                 S_IP += bio_p
            #                 bio_p = 0
                
            #     else:
                    
            #         warn('Not enough P in biomass and X_I to map the specified proportion of'
            #              'biomass COD into X_I. Rest of the biomass COD goes to S_A')
            #         X_I =  (bio_p + xi_p + S_IP) / asm_X_I_i_P
            #         xi_cod -= X_I
                    
            #         bio_p = xi_p = S_IP = 0
                    
            #         xi_n -= X_I*asm_X_I_i_N
            #         if xi_n < 0:
            #             bio_n += xi_n
            #             xi_n = 0
            #             if bio_n < 0:
            #                 S_IN += bio_n
            #                 bio_n = 0
            #------------------------------Rai version---------------------------------------------
            
            # Step 1(b)
                
            xsub_cod = X_ch + X_pr + X_li 
            xsub_n = X_pr*X_pr_i_N
            xsub_p = X_pr*X_pr_i_P
            
            xs_ndm = xsub_cod * X_S_i_N
            xs_pdm = xsub_cod * X_S_i_P
            
            if xs_ndm <= xsub_n + bio_n and xs_pdm <= xsub_p + bio_p:
                X_S += xsub_cod
                xsub_cod = 0
                xsub_n -= xs_ndm
                xsub_p -= xs_pdm
            else:
                if (xsub_n + bio_n)/X_S_i_N < (xsub_p + bio_p)/X_S_i_P:
                    
                    X_S_temp = (xsub_n + bio_n)/X_S_i_N 
                    X_S += X_S_temp
                    xsub_cod -= X_S_temp
                    xsub_n = bio_n = 0
                    
                    xsub_p -= X_S_temp*X_S_i_P
                    if xsub_p < 0: 
                        bio_p += xsub_p
                        xsub_p = 0
                    
                else: 
                    X_S_temp = (xsub_p + bio_p)/X_S_i_P
                    X_S += X_S_temp
                    xsub_cod -= X_S_temp
                    xsub_p = bio_p = 0
                    
                    xsub_n -= X_S_temp*X_S_i_N
                    if xsub_n < 0: 
                        bio_n += xsub_n
                        xsub_n = 0
                        
            # Step 3(A)
            
            # P balance not required as S_su, S_aa, S_fa do not have P
            ssub_cod = S_su + S_aa + S_fa
            ssub_n = S_aa * S_aa_i_N
            ssub_p = S_aa * S_aa_i_P # which would be 0
            
            sf_ndm = ssub_cod * S_F_i_N
            sf_pdm = ssub_cod * S_F_i_P
            
            if sf_ndm <= ssub_n + xsub_n + bio_n and sf_pdm <= ssub_p + xsub_p + bio_p:
                
                S_F = ssub_cod
                ssub_cod = 0
                
                ssub_n -= sf_ndm
                if ssub_n < 0:
                    xsub_n += ssub_n
                    ssub_n = 0
                    if xsub_n < 0:
                        bio_n += xsub_n
                        xsub_n = 0
                        
                ssub_p -= sf_pdm
                if ssub_p < 0:
                    xsub_p += ssub_p
                    ssub_p = 0
                    if xsub_p < 0:
                        bio_p += xsub_p
                        xsub_p = 0
                        
            else:
                if (ssub_n + xsub_n + bio_n) / S_F_i_N < (ssub_p + xsub_p + bio_p) / S_F_i_P:
                    
                    S_F = (ssub_n + xsub_n + bio_n) / S_F_i_N
                    ssub_cod -= S_F
                    ssub_n = xsub_n = bio_n = 0
                    
                    ssub_p -= S_F*S_F_i_P
                    
                    if ssub_p < 0:
                        xsub_p += ssub_p
                        ssub_p = 0
                        if xsub_p < 0:
                            bio_p += xsub_p
                            xsub_p = 0
                            
                else:
                    
                    S_F = (ssub_p + xsub_p + bio_p) / S_F_i_P
                    ssub_cod -= S_F
                    ssub_p = xsub_p = bio_p = 0
                    
                    ssub_n -= S_F*S_F_i_N
                    
                    if ssub_n < 0:
                        xsub_n += ssub_n
                        ssub_n = 0
                        if xsub_n < 0:
                            bio_n += xsub_n
                            xsub_n = 0
                
            # N and P balance not required as S_ac, S_pro, S_bu, and S_va do not have N and P
            S_A = S_ac + S_pro + S_bu + S_va
            
            # Step 4
            #!!! equality if enforced in `check_component_properties`
            # if asm_S_I_i_N == adm_S_I_i_N and asm_S_I_i_P == adm_S_I_i_P:
            #     S_I = S_I
            # else:
            #     raise RuntimeError('N and P content in S_I should be same')
            
            #------------------------------Rai version---------------------------------------------
            # si_n = S_I * adm_S_I_i_N
            # si_p = S_I * adm_S_I_i_P
            
            # si_ndm = si_cod * asm_S_I_i_N
            # si_pdm = si_cod * asm_S_I_i_P
            
            # if si_ndm <= si_n + xi_n + S_IN and si_pdm <= si_p + xi_p + S_IP:
            #     S_I = si_cod
            #     si_cod = 0
            #     si_n -= si_ndm
            #     if si_n < 0:
            #         xi_n += si_n
            #         si_n = 0
            #         if xi_n < 0:
            #             S_IN += xi_n
            #             xi_n = 0
            #     si_p -= si_pdm
            #     if si_p < 0:
            #         xi_p += si_p
            #         si_p = 0
            #         if xi_p < 0:
            #             S_IP += xi_p
            #             xi_p = 0
            # else:
            #     if (si_n + xi_n + S_IN) / asm_S_I_i_N < (si_p + xi_p + S_IP) / asm_S_I_i_P:
            #         S_I = (si_n + xi_n + S_IN) / asm_S_I_i_N
            #         si_cod -= S_I
            #         si_n = xi_n = S_IN = 0
            #         si_p -= S_I * asm_S_I_i_P
            #         if si_p < 0:
            #             xi_p += si_p
            #             si_p = 0
            #             if xi_p < 0:
            #                 S_IP += xi_p
            #                 xi_p = 0
            #     else:
            #         S_I = (si_p + xi_p + S_IP) / asm_S_I_i_P
            #         si_cod -= S_I
            #         si_p = xi_p = S_IP = 0
            #         si_n -= S_I * asm_S_I_i_N
            #         if si_n < 0:
            #             xi_n += si_n
            #             si_n = 0
            #             if xi_n < 0:
            #                 S_IN += xi_n
            #                 xi_n = 0
            #------------------------------Rai version---------------------------------------------
            
            #------------------------------Rai version---------------------------------------------
            # S_NH4 = S_IN + si_n + ssub_n + xsub_n + xi_n + bio_n
            # S_PO4 = S_IP + si_p + ssub_p + xsub_p + xi_p + bio_p
            #------------------------------Rai version---------------------------------------------
            
            S_NH4 = S_IN + ssub_n + xsub_n + bio_n - excess_N 
            S_PO4 = S_IP + ssub_p + xsub_p + bio_p - excess_P
            
            #------------------------------Rai version---------------------------------------------
            # S_A += si_cod + ssub_cod + xsub_cod + xi_cod + xs_cod
            #------------------------------Rai version---------------------------------------------
            
            S_A += ssub_cod + xsub_cod + xs_cod
            
            # Step 6: check COD and TKN balance
            asm_vals = np.array(([
                0, 0, # S_O2, S_N2,
                S_NH4, 
                0, # S_NO3
                S_PO4, S_F, S_A, S_I, 
                0,  # S_ALK(for now)
                X_I, X_S, 
                0,  # X_H,
                X_PAO, X_PP, X_PHA, # directly mapped 
                0, # X_AUT,
                X_MeOH, X_MeP, H2O])) # directly mapped
            
            if S_h2 > 0 or S_ch4 > 0:
                warn('Ignored dissolved H2 or CH4.')
            
            asm_vals = f_corr(adm_vals, asm_vals)
            
            # Step 5: charge balance for alkalinity
            
            # asm_ions_idx = cmps_asm.indices(('S_NH4', 'S_A', 'S_NO3', 'S_PO4', 'S_ALK'))
            
            S_NH4, S_A, S_NO3, S_PO4 = asm_vals[asm_ions_idx[:4]]
            
            # Need to include S_K, S_Mg in the charge balance
            
            # _ions = np.array([S_IN, S_IC, S_IP, S_Mg, S_K, S_ac, S_pro, S_bu, S_va])
            
            adm_alphas = np.array([self.alpha_IN, self.alpha_IC, self.alpha_IP, 
                                   2/24, 1/39, *self.alpha_vfa]) #!!! should be in unit of charge per g
            adm_charge = np.dot(_ions, adm_alphas)
            #!!! X_PP in ASM2d has negative charge, to compensate for the absent variables S_K & S_Mg
            S_ALK = (adm_charge - (S_NH4/14 - S_A/64 - S_NO3/14 - 1.5*S_PO4/31 - X_PP/31))*(-12)
    
            asm_vals[asm_ions_idx[-1]] = S_ALK
            
            return asm_vals
        
        self._reactions = madm12asm2d


#%% ASM2dtomADM1
class ASM2dtomADM1(mADMjunction):
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
        The anaerobic digestion process model (:class:`qsdsan.processes.ADM1_p_extension`).
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
    
    [2] Flores-Alsina, X., Solon, K., Kazadi Mbamba, C., Tait, S., Gernaey, K. V., 
    Jeppsson, U., & Batstone, D. J. (2016). Modelling phosphorus (P), sulfur (S) 
    and iron (FE) interactions for dynamic simulations of anaerobic digestion processes. 
    Water Research, 95, 370–382. 
    
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
    
    def __init__(self, ID='', upstream=None, downstream=(), thermo=None,
                 init_with='WasteStream', F_BM_default=None, isdynamic=False,
                 adm1_model=None, asm2d_model=None, T=298.15, pH=7):
        self._T = T
        self._pH = pH
        super().__init__(ID=ID, upstream=upstream, downstream=downstream,
                         thermo=thermo, init_with=init_with, 
                         F_BM_default=F_BM_default, isdynamic=isdynamic,
                         adm1_model=adm1_model, asm2d_model=asm2d_model)
        
    @property
    def T(self):
        '''[float] Temperature of the downstream [K].'''
        try: return self.outs[0].sink.T
        except: return self._T
    @T.setter
    def T(self, T):
        self._T = self.outs[0].T = T
    
    @property
    def pH(self):
        '''[float] downstream pH.'''
        if self._pH: return self._pH
        else:
            try: return self.outs[0].sink.outs[1].pH
            except: return 7.
    @pH.setter
    def pH(self, ph):
        self._pH = self.outs[0].pH = ph
    
    # def isbalanced(self, lhs, rhs_vals, rhs_i):
    #     rhs = sum(rhs_vals*rhs_i)
    #     error = rhs - lhs
    #     tol = max(self.rtol*lhs, self.rtol*rhs, self.atol)
    #     return abs(error) <= tol, error, tol, rhs
    
    def balance_cod_tkn_tp(self, asm_vals, adm_vals):
        cmps_asm = self.ins[0].components
        cmps_adm = self.outs[0].components
        asm_i_COD = cmps_asm.i_COD
        adm_i_COD = cmps_adm.i_COD
        non_tkn_idx = cmps_asm.indices(('S_N2', 'S_NO3'))
        asm_i_N = cmps_asm.i_N
        
        adm_i_N = cmps_adm.i_N
        asm_i_P = cmps_asm.i_P
        adm_i_P = cmps_adm.i_P
        asm_cod = sum(asm_vals*asm_i_COD)
        asm_tkn = sum(asm_vals*asm_i_N) - sum(asm_vals[non_tkn_idx])
        asm_tp = sum(asm_vals*asm_i_P)
        cod_bl, cod_err, cod_tol, adm_cod = self.isbalanced(asm_cod, adm_vals, adm_i_COD)
        tkn_bl, tkn_err, tkn_tol, adm_tkn = self.isbalanced(asm_tkn, adm_vals, adm_i_N)
        tp_bl, tp_err, tp_tol, adm_tp = self.isbalanced(asm_tp, adm_vals, adm_i_P)
        
        if tkn_bl and tp_bl:
            if cod_bl:
                return adm_vals
            else:
                if cod_err > 0: dcod = -(cod_err - cod_tol)/adm_cod
                else: dcod = -(cod_err + cod_tol)/adm_cod
                _adm_vals = adm_vals * (1 + (adm_i_COD>0)*dcod)
                _tkn_bl, _tkn_err, _tkn_tol, _adm_tkn = self.isbalanced(asm_tkn, _adm_vals, adm_i_N)
                _tp_bl, _tp_err, _tp_tol, _adm_tp = self.isbalanced(asm_tp, _adm_vals, adm_i_P)
                if _tkn_bl and _tp_bl: return _adm_vals
                else: 
                    warn('cannot balance COD, TKN, and TP at the same \n'
                        f'time with rtol={self.rtol} and atol={self.atol}.\n '
                        f'influent (ASM) TKN is {asm_tkn}\n '
                        f'effluent (ADM) TKN is {adm_tkn} or {_adm_tkn}\n '
                        f'influent TP is {asm_tp}\n ' 
                        f'effluent TP is {adm_tp} or {_adm_tp}. '
                        f'influent COD is {asm_cod}\n ' 
                        f'effluent COD is {adm_cod} or {adm_cod*(1+dcod)}. ')
                    return adm_vals
        elif cod_bl and tp_bl:
            if tkn_bl:
                return adm_vals
            else:
                if tkn_err > 0: dtkn = -(tkn_err - tkn_tol)/adm_tkn
                else: dtkn = -(tkn_err + tkn_tol)/adm_tkn
                _adm_vals = adm_vals * (1 + (adm_i_N>0)*dtkn)
                _cod_bl, _cod_err, _cod_tol, _adm_cod = self.isbalanced(asm_cod, _adm_vals, adm_i_COD)
                _tp_bl, _tp_err, _tp_tol, _adm_tp = self.isbalanced(asm_tp, _adm_vals, adm_i_P)
                if _cod_bl and _tp_bl: return _adm_vals
                else: 
                    warn('cannot balance COD, TKN, and TP at the same time'
                        f'time with rtol={self.rtol} and atol={self.atol}.\n '
                        f'influent (ASM) COD is {asm_cod}\n '
                        f'effluent (ADM) COD is {adm_cod} or {_adm_cod}\n '
                        f'influent TP is {asm_tp}\n ' 
                        f'effluent TP is {adm_tp} or {_adm_tp}. '
                        f'influent TKN is {asm_tkn}\n ' 
                        f'effluent TKN is {adm_tkn} or {adm_tkn*(1+dtkn)}. '
                        'To balance TKN please ensure ASM2d(X_I.i_N) = ADM1(X_I.i_N)')
                    return adm_vals
        elif cod_bl and tkn_bl: 
            if tp_bl:
                return adm_vals
            else:
                if tp_err > 0: dtp = -(tp_err - tp_tol)/adm_tp
                else: dtp = -(tp_err + tp_tol)/adm_tp
                _adm_vals = adm_vals * (1 + (adm_i_P>0)*dtp)
                _cod_bl, _cod_err, _cod_tol, _adm_cod = self.isbalanced(asm_cod, _adm_vals, adm_i_COD)
                _tkn_bl, _tkn_err, _tkn_tol, _adm_tkn = self.isbalanced(asm_tkn, _adm_vals, adm_i_N)
                if _cod_bl and _tkn_bl: return _adm_vals
                else: 
                    warn('cannot balance COD, TKN, and TP at the same time'
                        f'time with rtol={self.rtol} and atol={self.atol}.\n '
                        f'influent (ASM) COD is {asm_cod}\n '
                        f'effluent (ADM) COD is {adm_cod} or {_adm_cod}\n '
                        f'influent TKN is {asm_tkn}\n ' 
                        f'effluent TKN is {adm_tkn} or {_adm_tkn}. '
                        f'influent TP is {asm_tp}\n ' 
                        f'effluent TP is {adm_tp} or {adm_tp*(1+dtp)}. ')
                    return adm_vals
        else:
            warn('cannot balance COD, TKN and TP at the same time. \n'
                 'Atleast two of the three COD, TKN, and TP are not balanced \n'
                f'time with rtol={self.rtol} and atol={self.atol}.\n '
                f'influent (ASM) COD is {asm_cod}\n '
                f'effluent (ADM) COD is {adm_cod}\n '
                f'influent TP is {asm_tp}\n ' 
                f'effluent TP is {adm_tp}'
                f'influent TKN is {asm_tkn}\n ' 
                f'effluent TKN is {adm_tkn}. ')
            return adm_vals
                
    def _compile_reactions(self):
        # Retrieve constants
        # rtol = self.rtol
        # atol = self.atol

        cmps_asm = self.ins[0].components
        cmps_adm = self.outs[0].components
        self.check_component_properties(cmps_asm, cmps_adm)
        
        # For COD balance 
        S_NO3_i_COD = cmps_asm.S_NO3.i_COD
        
        # For N balance 
        X_H_i_N = cmps_asm.X_H.i_N
        X_AUT_i_N = cmps_asm.X_AUT.i_N
        S_F_i_N = cmps_asm.S_F.i_N
        X_S_i_N = cmps_asm.X_S.i_N
        # S_A_i_N = cmps_asm.S_A.i_N
        
        # Due to issue with mapping of X_I across ASM2d and ADM1, making this user dependent is important
        
        # # PUT A CHECK ON USER DEFINED N AND P CONTENT IN SI and XI here.
        # if self.asm_X_I_i_N == None:
        #     asm_X_I_i_N = cmps_asm.X_I.i_N
        # else:
        #     asm_X_I_i_N = self.asm_X_I_i_N
            
        # if self.asm_X_I_i_P == None:
        #     asm_X_I_i_P = cmps_asm.X_I.i_P
        # else:
        #     asm_X_I_i_P = self.asm_X_I_i_P
        
        # if self.asm_S_I_i_N == None:
        #     asm_S_I_i_N = cmps_asm.S_I.i_N
        # else:
        #     asm_S_I_i_N = self.asm_X_I_i_N
            
        # if self.asm_S_I_i_P == None:
        #     asm_S_I_i_P = cmps_asm.S_I.i_P
        # else:
        #     asm_S_I_i_P = self.asm_S_I_i_P
        
        # For P balance
        X_H_i_P = cmps_asm.X_H.i_P
        X_AUT_i_P = cmps_asm.X_AUT.i_P
        S_F_i_P = cmps_asm.S_F.i_P
        X_S_i_P = cmps_asm.X_S.i_P
        # S_A_i_P = cmps_asm.S_A.i_P
        
        # if cmps_asm.S_A.i_N > 0: 
        #     warn(f'S_A in ASM has positive nitrogen content: {cmps_asm.S_S.i_N} gN/gCOD. '
        #          'These nitrogen will be ignored by the interface model '
        #          'and could lead to imbalance of TKN after conversion.')
            
        # if cmps_asm.S_A.i_P > 0: 
        #     warn(f'S_A in ASM has positive phosphorous content: {cmps_asm.S_S.i_P} gN/gCOD. '
        #          'These phosphorous will be ignored by the interface model '
        #          'and could lead to imbalance of TP after conversion.')
        
        # ------------------------------Rai Version---------------------------------------------
        # if cmps_asm.S_I.i_P > 0:
        #     warn(f'S_I in ASM has positive phosphorous content: {cmps_asm.S_I.i_P} gN/gCOD. '
        #          'These phosphorous will be ignored by the interface model '
        #          'and could lead to imbalance of TP after conversion.')
        # ------------------------------Rai Version---------------------------------------------
        
        # We do not need to check if X_S.i_N != 0 since we take care of it using X_ND_asm1
        # We do not need to check if S_F.i_N != 0 since we take care of it using S_ND_asm1
        
        # For nitrogen balance 
        # S_ac_i_N = cmps_adm.S_ac.i_N
        S_aa_i_N = cmps_adm.S_aa.i_N
        X_pr_i_N = cmps_adm.X_pr.i_N
        # adm_S_I_i_N = cmps_adm.S_I.i_N
        adm_X_I_i_N = cmps_adm.X_I.i_N
        
        # For phosphorous balance 
        # S_ac_i_P = cmps_adm.S_ac.i_P
        X_pr_i_P = cmps_adm.X_pr.i_P
        # adm_S_I_i_P = cmps_adm.S_I.i_P
        adm_X_I_i_P = cmps_adm.X_I.i_P
        
        # Checks for direct mapping of X_PAO, X_PP, X_PHA
        
        # Check for X_PAO (measured as COD so i_COD = 1 in both ASM2d and ADM1)
        # asm_X_PAO_i_N = cmps_asm.X_PAO.i_N
        
        # if self.adm_X_PAO_i_N == None:
        #     adm_X_PAO_i_N = cmps_adm.X_PAO.i_N
        # else:
        #     adm_X_PAO_i_N = self.adm_X_PAO_i_N
            
        # if asm_X_PAO_i_N != adm_X_PAO_i_N:
        #     raise RuntimeError('X_PAO cannot be directly mapped as N content'
        #                        f'in asm2d_X_PAO_i_N = {asm_X_PAO_i_N} is not equal to'
        #                        f'adm_X_PAO_i_N = {adm_X_PAO_i_N}')
            
        # asm_X_PAO_i_P = cmps_asm.X_PAO.i_P
        
        # if self.adm_X_PAO_i_P == None:
        #     adm_X_PAO_i_P = cmps_adm.X_PAO.i_P
        # else:
        #     adm_X_PAO_i_P = self.adm_X_PAO_i_P
        
        # if asm_X_PAO_i_P != adm_X_PAO_i_P:
        #     raise RuntimeError('X_PAO cannot be directly mapped as P content'
        #                        f'in asm2d_X_PAO_i_P = {asm_X_PAO_i_P} is not equal to'
        #                        f'adm_X_PAO_i_P = {adm_X_PAO_i_P}')
        
        # # Checks not required for X_PP as measured as P in both, with i_COD = i_N = 0
        # # Checks not required for X_PHA as measured as COD in both, with i_N = i_P = 0
        
        # if asm_X_I_i_N != adm_X_I_i_N:
        #     raise RuntimeError('X_I cannot be directly mapped as N content'
        #                        f'in asm2d_X_I_i_N = {asm_X_I_i_N} is not equal to'
        #                        f'adm_X_I_i_N = {adm_X_I_i_N}')
            
        # if asm_X_I_i_P != adm_X_I_i_P:
        #     raise RuntimeError('X_I cannot be directly mapped as P content'
        #                        f'in asm2d_X_I_i_P = {asm_X_I_i_P} is not equal to'
        #                        f'adm_X_I_i_P = {adm_X_I_i_P}')
        
        # if asm_S_I_i_N != adm_X_I_i_N:
        #     raise RuntimeError('S_I cannot be directly mapped as N content'
        #                        f'in asm2d_S_I_i_N = {asm_S_I_i_N} is not equal to'
        #                        f'adm_S_I_i_N = {adm_S_I_i_N}')
            
        # if asm_S_I_i_P != adm_S_I_i_P:
        #     raise RuntimeError('S_I cannot be directly mapped as P content'
        #                        f'in asm2d_S_I_i_P = {asm_S_I_i_P} is not equal to \n'
        #                        f'adm_S_I_i_P = {adm_S_I_i_P}')
        
        adm_ions_idx = cmps_adm.indices(['S_IN', 'S_IP', 'S_IC', 'S_cat', 'S_an'])
        
        frac_deg = self.frac_deg
        f_corr = self.balance_cod_tkn_tp

        # To convert components from ASM2d to mADM1 (asm2d-2-madm1)
        def asm2d2madm1(asm_vals):
            # S_I, S_S, X_I, X_S, X_BH, X_BA, X_P, S_O, S_NO, S_NH, S_ND, X_ND, S_ALK, S_N2, H2O = asm_vals
            
            S_O2, S_N2, S_NH4, S_NO3, S_PO4, S_F, S_A, S_I, S_ALK, X_I, X_S, X_H, \
                X_PAO, X_PP, X_PHA, X_AUT, X_MeOH, X_MeP, H2O = asm_vals

            # Step 0: charged component snapshot (# pg. 84 of IWA ASM textbook)
            _sno3 = S_NO3
            _snh4 = S_NH4
            _salk = S_ALK  
            _spo4 = S_PO4
            _sa = S_A 
            _xpp = X_PP
              
            # Step 1: remove any remaining COD demand
            O2_coddm = S_O2
            NO3_coddm = -S_NO3*S_NO3_i_COD
            
            # cod_spl = S_S + X_S + X_BH + X_BA
            # Replacing S_S with S_F + S_A (IWA ASM textbook)
            
            cod_spl = (S_A + S_F) + X_S + (X_H + X_AUT)
            
            # bioN = X_BH*X_BH_i_N + X_BA*X_BA_i_N
            
            bioN = X_H*X_H_i_N + X_AUT*X_AUT_i_N
            bioP = X_H*X_H_i_P + X_AUT*X_AUT_i_P
            
            # To be used in Step 2
            S_ND_asm1 = S_F*S_F_i_N   #S_ND (in asm1) equals the N content in S_F
            # To be used in Step 3
            X_ND_asm1 = X_S*X_S_i_N   #X_ND (in asm1) equals the N content in X_S
            # To be used in Step 5 (a)
            X_S_P = X_S*X_S_i_P
            # To be used in Step 5 (b)
            S_F_P = S_F*S_F_i_P 
            
            if cod_spl <= O2_coddm:
                S_O2 = O2_coddm - cod_spl
                S_F = S_A =  X_S = X_H = X_AUT = 0
            elif cod_spl <= O2_coddm + NO3_coddm:
                S_O2 = 0
                S_NO3 = -(O2_coddm + NO3_coddm - cod_spl)/S_NO3_i_COD
                S_A = S_F = X_S = X_H = X_AUT = 0
            else:
                S_A -= O2_coddm + NO3_coddm
                if S_A < 0:
                    S_F += S_A
                    S_A = 0
                    if S_F < 0:
                        X_S += S_F
                        S_F = 0
                        if X_S < 0:
                            X_H += X_S
                            X_S = 0
                            if X_H < 0:
                                X_AUT += X_H
                                X_H = 0
                S_O2 = S_NO3 = 0
            
            # Step 2: convert any readily biodegradable 
            # COD and TKN into amino acids and sugars
            
            # Directly map acetate
            S_ac = S_A
            # if S_ac_i_N == S_A_i_N and S_ac_i_P == S_A_i_P:
            #     S_ac = S_A
            # else:
            #     raise RuntimeError('N and P content should be the same in S_A and S_ac for direct translation')
            
            # First we calculate the amount of amino acid required in ADM1
            # if all available soluble organic N can be mapped to amino acid
            req_scod = S_ND_asm1 / S_aa_i_N
            
            # if available S_F is not enough to fulfill that amino acid requirement 
            if S_F < req_scod: 
                # then all available S_F is mapped to amino acids 
                S_aa = S_F
                # and no S_S would be available for conversion to sugars
                S_su = 0
                # This needs to be followed by a corresponding loss in soluble organic N 
                S_ND_asm1 -= S_aa * S_aa_i_N
            # if available S_F is more than enough to fulfill that amino acid requirement 
            else:
                # All soluble organic N will be mapped to amino acid
                S_aa = req_scod
                # The line above implies that a certain portion of S_S would also be consumed to form amino acid
                # The S_S which is left would form sugar 
                # In simpler terms; S_F = S_F - S_aa; S_su = S_F 
                S_su = S_F - S_aa
                # All soluble organic N would thus be consumed in amino acid formation
                S_ND_asm1 = 0

            # Step 3: convert slowly biodegradable COD and TKN
            # into proteins, lipids, and carbohydrates
            
            # First we calculate the amount of protein required in ADM1
            # if all available particulate organic N can be mapped to amino acid
            req_xcod = X_ND_asm1 / X_pr_i_N
            # Since X_pr_i_N >> X_pr_i_P there's no need to check req_xcod for N and P separately (CONFIRM LATER 05/16)
            
            # if available X_S is not enough to fulfill that protein requirement
            if X_S < req_xcod:
                # then all available X_S is mapped to amino acids
                X_pr = X_S
                # and no X_S would be available for conversion to lipid or carbohydrates 
                X_li = X_ch = 0
                # This needs to be followed by a corresponding loss in particulate organic N 
                X_ND_asm1 -= X_pr * X_pr_i_N
                
                # For P balance 
                # This needs to be followed by a corresponding loss in particulate organic P
                X_S_P -= X_pr * X_pr_i_P
                
            # if available X_S is more than enough to fulfill that protein requirement
            else:
                # All particulate organic N will be mapped to amino acid
                X_pr = req_xcod
                # The line above implies that a certain portion of X_S would also be consumed to form protein
                # The X_S which is left would form lipid and carbohydrates in a percentage define by the user  
                X_li = self.xs_to_li * (X_S - X_pr)
                X_ch = (X_S - X_pr) - X_li
                # All particulate organic N would thus be consumed in amino acid formation
                X_ND_asm1 = 0
                
                # For P balance 
                # This needs to be followed by a corresponding loss in particulate organic P
                X_S_P -= X_pr * X_pr_i_P
            
            # Step 4: convert active biomass into protein, lipids, 
            # carbohydrates and potentially particulate TKN
            
            # First the amount of biomass N/P available for protein, lipid etc is determined
            # For this calculation, from total biomass N available the amount 
            # of particulate inert N/P expected in ADM1 is subtracted 
            
            available_bioN = bioN - (X_H + X_AUT) * (1-frac_deg) * adm_X_I_i_N
            if available_bioN < 0:
                raise RuntimeError('Not enough N in X_H and X_AUT to fully convert '
                                   'the non-biodegradable portion into X_I in ADM1.')
                
            available_bioP = bioP - (X_H + X_AUT) * (1-frac_deg) * adm_X_I_i_P
            if available_bioP < 0:
                raise RuntimeError('Not enough P in X_H and X_AUT to fully convert '
                                   'the non-biodegradable portion into X_I in ADM1.')
                
            # Then the amount of biomass N/P required for biomass conversion to protein is determined
            req_bioN = (X_H + X_AUT) * frac_deg * X_pr_i_N
            req_bioP = (X_H + X_AUT) * frac_deg * X_pr_i_P
            
            # Case I: if both available biomass N/P and particulate organic N/P is greater than 
            # required biomass N/P for conversion to protein
            if available_bioN + X_ND_asm1 >= req_bioN and available_bioP + X_S_P >= req_bioP:
                # then all biodegradable biomass N/P (corrsponding to protein demand) is converted to protein
                X_pr += (X_H + X_AUT) * frac_deg
                # the remaining biomass N/P is transfered as organic N/P
                X_ND_asm1 += available_bioN - req_bioN 
                X_S_P += available_bioP - req_bioP 
                
            # Case II: if available biomass N and particulate organic N is less than 
            # required biomass N for conversion to protein, but available biomass P and  
            # particulate organic P is greater than required biomass P for conversion to protein
            
            # Case III: if available biomass P and particulate organic P is less than 
            # required biomass P for conversion to protein, but available biomass N and  
            # particulate organic N is greater than required biomass N for conversion to protein
            
            # Case IV: if both available biomass N/P and particulate organic N/P is less than 
            # required biomass N/P for conversion to protein
            else:
                
                if (available_bioP + X_S_P)/X_pr_i_P < (available_bioN + X_ND_asm1)/X_pr_i_N:
                    # all available P and particulate organic P is converted to protein
                    bio2pr = (available_bioP + X_S_P)/X_pr_i_P
                    X_pr += bio2pr
                    # Biodegradable biomass available after conversion to protein is calculated 
                    bio_to_split = (X_H + X_AUT) * frac_deg - bio2pr
                    # Part of the remaining biomass is mapped to lipid based on user defined value 
                    bio_split_to_li = bio_to_split * self.bio_to_li
                    X_li += bio_split_to_li
                    # The other portion of the remanining biomass is mapped to carbohydrates 
                    X_ch += (bio_to_split - bio_split_to_li)
                    # Since all organic P has been mapped to protein, none is left
                    X_S_P = 0
                    
                    # the remaining biomass N is transfered as organic N
                    X_ND_asm1 += available_bioN - (bio2pr*X_pr_i_N)
                
                else:
                    # all available N and particulate organic N is converted to protein
                    bio2pr = (available_bioN + X_ND_asm1)/X_pr_i_N
                    X_pr += bio2pr
                    # Biodegradable biomass available after conversion to protein is calculated 
                    bio_to_split = (X_H + X_AUT) * frac_deg - bio2pr
                    # Part of the remaining biomass is mapped to lipid based on user defined value 
                    bio_split_to_li = bio_to_split * self.bio_to_li
                    X_li += bio_split_to_li
                    # The other portion of the remanining biomass is mapped to carbohydrates 
                    X_ch += (bio_to_split - bio_split_to_li)
                    # Since all organic N has been mapped to protein, none is left
                    X_ND_asm1 = 0
                    
                    # the remaining biomass P is transfered as organic P
                    X_S_P += available_bioP - (bio2pr*X_pr_i_P)
            
            
            # Step 5: map particulate inerts
            
            # 5 (a)
            
            # --------------------------Rai version--------------------------------------
            # # First determine the amount of particulate inert N/P available from ASM2d
            # xi_nsp_asm2d = X_I * asm_X_I_i_N
            # xi_psp_asm2d = X_I * asm_X_I_i_P
            
            # # Then determine the amount of particulate inert N/P that could be produced 
            # # in ADM1 given the ASM1 X_I
            # xi_ndm = X_I * adm_X_I_i_N
            # xi_pdm = X_I * adm_X_I_i_P

            # # if particulate inert N available in ASM1 is greater than ADM1 demand
            # if xi_nsp_asm2d + X_ND_asm1 >= xi_ndm:
            #     deficit = xi_ndm - xi_nsp_asm2d
            #     # COD balance 
            #     X_I += (X_H+X_AUT) * (1-frac_deg)
            #     # N balance 
            #     X_ND_asm1 -= deficit
            #     # P balance 
            #     if xi_psp_asm2d + X_S_P >= xi_pdm:
            #         deficit = xi_pdm - xi_psp_asm2d
            #         X_S_P -= deficit
            #     elif isclose(xi_psp_asm2d+X_S_P, xi_pdm, rel_tol=rtol, abs_tol=atol):
            #         X_S_P  = 0
            #     else:
            #         raise RuntimeError('Not enough P in X_I, X_S to fully '
            #                            'convert X_I in ASM2d into X_I in ADM1.')
            # elif isclose(xi_nsp_asm2d+X_ND_asm1, xi_ndm, rel_tol=rtol, abs_tol=atol):
            #     # COD balance 
            #     X_I += (X_H+X_AUT) * (1-frac_deg)
            #     # N balance 
            #     X_ND_asm1 = 0
            #     # P balance 
            #     if xi_psp_asm2d + X_S_P >= xi_pdm:
            #         deficit = xi_pdm - xi_psp_asm2d
            #         X_S_P -= deficit
            #     elif isclose(xi_psp_asm2d+X_S_P, xi_pdm, rel_tol=rtol, abs_tol=atol):
            #         X_S_P  = 0
            #     else:
            #         raise RuntimeError('Not enough P in X_I, X_S to fully '
            #                            'convert X_I in ASM2d into X_I in ADM1.')
            # else:
            # # Since the N balance cannot hold, the P balance is not futher checked 
            #     raise RuntimeError('Not enough N in X_I, X_S to fully '
            #                        'convert X_I in ASM2d into X_I in ADM1.')
            # --------------------------Rai version--------------------------------------
            
            # --------------------------Flores Alsina et al. 2016 version--------------------------------------
            # X_I = X_I
            # COD balance 
            X_I += (X_H+X_AUT) * (1-frac_deg)
            # --------------------------Flores Alsina et al. 2016 version--------------------------------------
            
            # # 5(b)
            
            # --------------------------Rai version--------------------------------------
            # # Then determine the amount of soluble inert N/P that could be produced 
            # # in ADM1 given the ASM1 X_I
            # req_sn = S_I * adm_S_I_i_N
            # req_sp = S_I * adm_S_I_i_P
            
            # supply_inert_n_asm2d = S_I * asm_S_I_i_N
            
            # # N balance 
            # if req_sn <= S_ND_asm1 + supply_inert_n_asm2d:
            #     S_ND_asm1 -= (req_sn - supply_inert_n_asm2d)
            #     supply_inert_n_asm2d = 0 
            #     # P balance 
            #     if req_sp <= S_F_P:
            #         S_F_P -= req_sp
            #     elif req_sp <= S_F_P + X_S_P:
            #         X_S_P -= (req_sp - S_F_P)
            #         S_F_P = 0
            #     elif req_sp <= S_F_P + X_S_P + S_PO4:
            #         S_PO4 -= (req_sp - S_F_P - X_S_P)
            #         S_F_P = X_S_P = 0
            #     else:
            #         S_PO4 -= (req_sp - S_F_P - X_S_P)
            #         S_F_P =  X_S_P = 0
            # # N balance
            # elif req_sn <= S_ND_asm1 + X_ND_asm1 + supply_inert_n_asm2d:
            #     X_ND_asm1 -= (req_sn - S_ND_asm1 - supply_inert_n_asm2d)
            #     S_ND_asm1 = supply_inert_n_asm2d = 0
            #     # P balance
            #     if req_sp <= S_F_P:
            #         S_F_P -= req_sp
            #     elif req_sp <= S_F_P + X_S_P:
            #         X_S_P -= (req_sp - S_F_P)
            #         S_F_P = 0
            #     elif req_sp <= S_F_P + X_S_P + S_PO4:
            #         S_PO4 -= (req_sp - S_F_P - X_S_P)
            #         S_F_P = X_S_P = 0
            #     else:
            #         S_PO4 -= (req_sp - S_F_P - X_S_P)
            #         S_F_P =  X_S_P = 0
            # # N balance
            # elif req_sn <= S_ND_asm1 + X_ND_asm1 + S_NH4 + supply_inert_n_asm2d:
            #     S_NH4 -= (req_sn - S_ND_asm1 - X_ND_asm1 - supply_inert_n_asm2d)
            #     S_ND_asm1 = X_ND_asm1 = supply_inert_n_asm2d = 0
            #     # P balance 
            #     if req_sp <= S_F_P:
            #         S_F_P -= req_sp
            #     elif req_sp <= S_F_P + X_S_P:
            #         X_S_P -= (req_sp - S_F_P)
            #         S_F_P = 0
            #     elif req_sp <= S_F_P + X_S_P + S_PO4:
            #         S_PO4 -= (req_sp - S_F_P - X_S_P)
            #         S_F_P = X_S_P = 0
            #     else:
            #         S_PO4 -= (req_sp - S_F_P - X_S_P)
            #         S_F_P =  X_S_P = 0
            # elif req_sp <= S_F_P or req_sp <= S_F_P + X_S_P or req_sp <= S_F_P + X_S_P + S_PO4:
                
            #     S_NH4 -= (req_sn - S_ND_asm1 - X_ND_asm1 - supply_inert_n_asm2d)
            #     S_ND_asm1 = X_ND_asm1 = supply_inert_n_asm2d = 0
                
            #     if req_sp <= S_F_P:
            #         S_F_P -= req_sp
            #     elif req_sp <= S_F_P + X_S_P:
            #         X_S_P -= (req_sp - S_F_P)
            #         S_F_P = 0
            #     elif req_sp <= S_F_P + X_S_P + S_PO4:
            #         S_PO4 -= (req_sp - S_F_P - X_S_P)
            #         S_F_P = X_S_P = 0
            # else:
            #     if (S_ND_asm1 + X_ND_asm1 + S_NH4 + supply_inert_n_asm2d)/adm_S_I_i_N < (S_F_P + X_S_P + S_PO4)/adm_S_I_i_P:
            #         warn('Additional soluble inert COD is mapped to S_su.')
            #         SI_cod = (S_ND_asm1 + X_ND_asm1 + S_NH4 + supply_inert_n_asm2d)/adm_S_I_i_N
            #         S_su += S_I - SI_cod
            #         S_I = SI_cod
            #         S_ND_asm1 = X_ND_asm1 = S_NH4 = supply_inert_n_asm2d = 0
                    
            #         req_sp = S_I * adm_S_I_i_P
            #         S_PO4 -= (req_sp - S_F_P - X_S_P)
            #         S_F_P = X_S_P = 0
            #     else:
            #         warn('Additional soluble inert COD is mapped to S_su.')
            #         SI_cod = (S_F_P + X_S_P + S_PO4)/adm_S_I_i_P
            #         S_su += S_I - SI_cod
            #         S_I = SI_cod
            #         S_F_P = X_S_P = S_PO4 = 0
                    
            #         req_sn = S_I * adm_S_I_i_N
            #         S_NH4 -= (req_sn - S_ND_asm1 - X_ND_asm1 - supply_inert_n_asm2d)
            #         S_ND_asm1 = X_ND_asm1 = supply_inert_n_asm2d = 0
            
            # if S_PO4 < 0:
            #     raise RuntimeError(f'S_PO4 = {S_PO4}\n''Not enough P in S_F_P, X_S_P and S_PO4 to fully '
            #                        'convert S_I in ASM2d into S_I in ADM1. Consider '
            #                        'increasing the value of P content in S_I (ASM2d)')
                
            # if S_NH4 < 0:
            #     raise RuntimeError(f'S_NH4 = {S_NH4}\n''Not enough N in S_I, S_ND_asm1, X_ND_asm1, and S_NH4 to fully '
            #                         'convert S_I in ASM2d into S_I in ADM1. Consider '
            #                         'increasing the value of N content in S_I (ASM2d)')
            # --------------------------Rai version--------------------------------------
            
            # --------------------------Flores Alsina et al. 2016 version--------------------------------------
            # S_I = S_I
            # --------------------------Flores Alsina et al. 2016 version--------------------------------------
        
            # Step 6: Step map any remaining TKN/P
            
            # --------------------------Rai version--------------------------------------
            # S_IN = S_ND_asm1 + X_ND_asm1 + S_NH4 + supply_inert_n_asm2d
            # --------------------------Rai version--------------------------------------
            
            # --------------------------Flores Alsina et al. 2016 version--------------------------------------
            S_IN = S_ND_asm1 + X_ND_asm1 + S_NH4
            # --------------------------Flores Alsina et al. 2016 version--------------------------------------
            
            S_IP = S_F_P + X_S_P + S_PO4            
            
            # Step 8: check COD and TKN balance
            # has TKN: S_aa, S_IN, S_I, X_pr, X_I
            S_IC = S_cat = S_an = 0
            
            # When mapping components directly in Step 9 ensure the values of
            # cmps.i_N, cmps.i_P, and cmps.i_COD are same in both ASM2d and ADM1
            
            # Step 9: Mapping common state variables directly    
            # The next three commented lines are executed when outputting
            # array of ADM1 components 
            # X_PAO (ADM1) = X_PAO (ASM2d)
            # X_PP (ADM1) = X_PP (ASM2d)
            # X_PHA (ADM1) = X_PHA (ASM2d)
            # X_MeOH (ADM1) = X_MeOH (ASM2d)
            # X_MeP (ADM1) = X_MeP (ASM2d)
            
            adm_vals = np.array([
                S_su, S_aa, 
                0, 0, 0, 0, S_ac, # S_fa, S_va, S_bu, S_pro,
                0, 0, # S_h2, S_ch4,
                S_IC, S_IN, S_IP, S_I, 
                X_ch, X_pr, X_li, 
                0, 0, 0, 0, 0, 0, 0, # X_su, X_aa, X_fa, X_c4, X_pro, X_ac, X_h2,
                X_I, X_PHA, X_PP, X_PAO, 
                0, 0,  # S_K, S_Mg,
                X_MeOH, X_MeP,
                S_cat, S_an, H2O])
            
            adm_vals = f_corr(asm_vals, adm_vals)
            
            # Step 7: charge balance
            asm_charge = - _sa/64 + _snh4/14 - _sno3/14 - 1.5*_spo4/31 - _salk - _xpp/31 #Based on page 84 of IWA ASM handbook
            S_IN, S_IP = adm_vals[adm_ions_idx[:2]]
            
            #!!! charge balance should technically include VFAs, S_K, S_Mg,
            # but since their concentrations are assumed zero it is acceptable.
            S_IC = (asm_charge - S_IN*self.alpha_IN - S_IP*self.alpha_IP)/self.alpha_IC
            proton_charge = 10**(-self.pKa[0]+self.pH) - 10**(-self.pH) # self.pKa[0] is pKw
            # net_Scat = Scat - San
            net_Scat = asm_charge + proton_charge
            
            if net_Scat > 0:  
                S_cat = net_Scat
                S_an = 0
            else:
                S_cat = 0
                S_an = -net_Scat
            
            adm_vals[adm_ions_idx[2:]] = [S_IC, S_cat, S_an]
            
            return adm_vals
        
        self._reactions = asm2d2madm1
        
#%% A1junction

class A1junction(ADMjunction):
    '''
    An abstract superclass holding common properties of modified ADM interface classes.
    Users should use its subclasses (e.g., ``mASM2dtoADM1p``, ``ADM1ptomASM2d``) instead.
    
    See Also
    --------
    :class:`qsdsan.sanunits.ADMJunction`
    
    :class:`qsdsan.sanunits.ADM1ptomASM2d`
    
    :class:`qsdsan.sanunits.mASM2dtoADM1p`
    '''
    _parse_reactions = Junction._no_parse_reactions
    rtol = 1e-2
    atol = 1e-6
    cod_vfa = np.array([64, 112, 160, 208])

    def __init__(self, ID='', upstream=None, downstream=(), thermo=None,
                 init_with='WasteStream', F_BM_default=None, isdynamic=False,
                 adm1_model=None, asm2d_model=None):
        self.asm2d_model = asm2d_model
        super().__init__(ID=ID, upstream=upstream, downstream=downstream,
                         thermo=thermo, init_with=init_with, 
                         F_BM_default=F_BM_default, isdynamic=isdynamic,
                         adm1_model=adm1_model)
    
    @property
    def asm2d_model(self):
        '''[qsdsan.CompiledProcesses] ASM2d process model.'''
        return self._asm2d_model
    @asm2d_model.setter
    def asm2d_model(self, model):
        if not isinstance(model, (pc.ASM2d, pc.mASM2d)):
            raise ValueError('`asm2d_model` must be an `ASM2d` object, '
                              f'the given object is {type(model).__name__}.')
        self._asm2d_model = model

    @property
    def adm1_model(self):
        '''[qsdsan.CompiledProcesses] mADM1 process model.'''
        return self._adm1_model
    @adm1_model.setter
    def adm1_model(self, model):
        if not isinstance(model, (pc.ADM1_p_extension, pc.ADM1p)):
            raise ValueError('`adm1_model` must be an `ADM1_p_extension` object, '   #!!! update error message
                              f'the given object is {type(model).__name__}.')
        self._adm1_model = model

    def check_component_properties(self, cmps_asm, cmps_adm):
        get = getattr
        setv = setattr
        for name in ('X_PHA', 'X_PP', 'X_PAO'):
            casm = get(cmps_asm, name)
            cadm = get(cmps_adm, name)
            for attr in ('measured_as', 'i_COD', 'i_C', 'i_N', 'i_P'):
                vasm = get(casm, attr)
                if get(cadm, attr) != vasm:
                    setv(cadm, attr, vasm)
                    warn(f"ADM component {name}'s {attr} is changed to match "
                         "the corresponding ASM component")
        
        for name in ('S_I', 'X_I'):
            casm = get(cmps_asm, name)
            cadm = get(cmps_adm, name)
            for attr in ('measured_as', 'i_C', 'i_N', 'i_P'):
                vadm = get(cadm, attr)
                if get(casm, attr) != vadm:
                    setv(casm, attr, vadm)
                    warn(f"ASM component {name}'s {attr} is changed to match "
                         "the corresponding ADM component")        
        
        for attr in ('i_N', 'i_P'):
            vadm = get(cmps_adm.S_ac, attr)
            if get(cmps_asm.S_A, attr) != vadm:
                cmps_asm.S_A.i_N = vadm
                warn(f"ASM component S_A's {attr} is changed to match "
                     "the ADM component S_ac.") 
        
        if cmps_asm.S_IC.measured_as != cmps_adm.S_IC.measured_as:
            raise RuntimeError('S_ALK in ASM and S_IC in ADM must both be measured as "C".')
        if cmps_asm.S_NH4.measured_as != cmps_adm.S_IN.measured_as:
            raise RuntimeError('S_NH4 in ASM and S_IN in ADM must both be measured as "N".')
        if cmps_asm.S_PO4.measured_as != cmps_adm.S_IP.measured_as:
            raise RuntimeError('S_PO4 in ASM and S_IP in ADM must both be measured as "P".')
        cmps_asm.refresh_constants()
        cmps_adm.refresh_constants()


#%% ADM1ptomASM2d
class ADM1ptomASM2d(A1junction):
    '''
    Interface unit to convert ADM1 state variables
    to ASM2d components, following the A1 algorithm in [1]_.
    
    Parameters
    ----------
    upstream : stream or str
        Influent stream with ADM components.
    downstream : stream or str
        Effluent stream with ASM components.
    adm1_model : obj
        The anaerobic digestion process model (:class:`qsdsan.processes.ADM1_p_extension`).
    rtol : float
        Relative tolerance for COD and TKN balance.
    atol : float
        Absolute tolerance for COD and TKN balance.
    
    References
    ----------
    [1] Flores-Alsina, X., Solon, K., Kazadi Mbamba, C., Tait, S., Gernaey, K. V., 
    Jeppsson, U., & Batstone, D. J. (2016). Modelling phosphorus (P), sulfur (S) 
    and iron (FE) interactions for dynamic simulations of anaerobic digestion processes. 
    Water Research, 95, 370–382.
    
    See Also
    --------
    :class:`qsdsan.sanunits.A1junction`
    
    :class:`qsdsan.sanunits.mASM2dtoADM1p`
    '''
    
    def balance_cod_tkn(self, adm_vals, asm_vals):
        cmps_adm = self.ins[0].components
        cmps_asm = self.outs[0].components
        asm_i_COD = cmps_asm.i_COD
        adm_i_COD = cmps_adm.i_COD
        gas_idx = cmps_adm.indices(('S_h2', 'S_ch4'))
        
        asm_i_N = cmps_asm.i_N
        adm_i_N = cmps_adm.i_N    
        asm_i_P = cmps_asm.i_P
        adm_i_P = cmps_adm.i_P
        adm_cod = sum(adm_vals*adm_i_COD) - sum(adm_vals[gas_idx])
        adm_tkn = sum(adm_vals*adm_i_N)
        adm_tp = sum(adm_vals*adm_i_P)
        
        cod_bl, cod_err, cod_tol, asm_cod = self.isbalanced(adm_cod, asm_vals, asm_i_COD)
        tkn_bl, tkn_err, tkn_tol, asm_tkn = self.isbalanced(adm_tkn, asm_vals, asm_i_N)
        tp_bl, tp_err, tp_tol, asm_tp = self.isbalanced(adm_tp, asm_vals, asm_i_P)
        
        if tkn_bl and tp_bl:
            if cod_bl:
                return asm_vals
            else:
                if cod_err > 0: dcod = -(cod_err - cod_tol)/asm_cod
                else: dcod = -(cod_err + cod_tol)/asm_cod
                _asm_vals = asm_vals * (1 + (asm_i_COD>0)*dcod)
                _tkn_bl, _tkn_err, _tkn_tol, _asm_tkn = self.isbalanced(adm_tkn, _asm_vals, asm_i_N)
                _tp_bl, _tp_err, _tp_tol, _asm_tp = self.isbalanced(adm_tp, _asm_vals, asm_i_P)
                if _tkn_bl and _tp_bl: return _asm_vals
                else: 
                    warn('cannot balance COD, TKN, and TP at the same \n'
                        f'time with rtol={self.rtol} and atol={self.atol}.\n '
                        f'influent (ADM) TKN is {adm_tkn}\n '
                        f'effluent (ASM) TKN is {asm_tkn} or {_asm_tkn}\n '
                        f'influent (ADM) TP is {adm_tp}\n ' 
                        f'effluent (ASM) TP is {asm_tp} or {_asm_tp}. '
                        f'influent (ADM) COD is {adm_cod}\n ' 
                        f'effluent (ASM) COD is {asm_cod} or {asm_cod*(1+dcod)}. ')
                    return asm_vals
        elif cod_bl and tp_bl:
            if tkn_bl:
                return asm_vals
            else:
                if tkn_err > 0: dtkn = -(tkn_err - tkn_tol)/asm_tkn
                else: dtkn = -(tkn_err + tkn_tol)/asm_tkn
                _asm_vals = asm_vals * (1 + (asm_i_N>0)*dtkn)
                _cod_bl, _cod_err, _cod_tol, _asm_cod = self.isbalanced(adm_cod, _asm_vals, asm_i_COD)
                _tp_bl, _tp_err, _tp_tol, _asm_tp = self.isbalanced(adm_tp, _asm_vals, asm_i_P)
                if _cod_bl and _tp_bl: return _asm_vals
                else: 
                    warn('cannot balance COD, TKN, and TP at the same time'
                        f'time with rtol={self.rtol} and atol={self.atol}.\n '
                        f'influent (ADM) COD is {adm_cod}\n '
                        f'effluent (ASM) COD is {asm_cod} or {_asm_cod}\n '
                        f'influent (ADM) TP is {adm_tp}\n ' 
                        f'effluent (ASM) TP is {asm_tp} or {_asm_tp}. '
                        f'influent (ADM) TKN is {adm_tkn}\n ' 
                        f'effluent (ASM) TKN is {asm_tkn} or {asm_tkn*(1+dtkn)}. ')
                    return asm_vals
        elif cod_bl and tkn_bl:
            if tp_bl:
                return asm_vals
            else:
                if tp_err > 0: dtp = -(tp_err - tp_tol)/asm_tp
                else: dtp = -(tp_err + tp_tol)/asm_tp
                _asm_vals = asm_vals * (1 + (asm_i_P>0)*dtp)
                _cod_bl, _cod_err, _cod_tol, _asm_cod = self.isbalanced(adm_cod, _asm_vals, asm_i_COD)
                _tkn_bl, _tkn_err, _tkn_tol, _asm_tkn = self.isbalanced(adm_tkn, _asm_vals, asm_i_N)
                if _cod_bl and _tkn_bl: return _asm_vals
                else: 
                    warn('cannot balance COD, TKN, and TP at the same time'
                        f'time with rtol={self.rtol} and atol={self.atol}.\n '
                        f'influent (ADM) COD is {adm_cod}\n '
                        f'effluent (ASM) COD is {asm_cod} or {_asm_cod}\n '
                        f'influent (ADM) TKN is {adm_tkn}\n ' 
                        f'effluent (ASM) TKN is {asm_tkn} or {_asm_tkn}. '
                        f'influent (ADM) TP is {adm_tp}\n ' 
                        f'effluent (ASM) TP is {asm_tp} or {asm_tp*(1+dtp)}. ')
                    return asm_vals
        else:
            warn('cannot balance COD, TKN and TP at the same time. \n'
                 'Atleast two of the three COD, TKN, and TP are not balanced \n'
                f'time with rtol={self.rtol} and atol={self.atol}.\n '
                f'influent (ADM) COD is {adm_cod}\n '
                f'effluent (ASM) COD is {asm_cod}\n '
                f'influent (ADM) TP is {adm_tp}\n ' 
                f'effluent (ASM) TP is {asm_tp}'
                f'influent (ADM) TKN is {adm_tkn}\n ' 
                f'effluent (ASM) TKN is {asm_tkn}. ')
            return asm_vals
    
    def _compile_reactions(self):        
        cmps_adm = self.ins[0].components
        cmps_asm = self.outs[0].components
        self.check_component_properties(cmps_asm, cmps_adm)
        
        _asm_ids = cmps_asm.indices(['S_F', 'X_S', 'S_A'])
        _adm_ids = cmps_adm.indices(['S_su', 'S_aa', 'S_fa', 
                                     'S_va', 'S_bu', 'S_pro', 'S_ac',
                                     'X_pr', 'X_li', 'X_ch'])

        # For carbon balance
        C_SF, C_XS, C_SA = cmps_asm.i_C[_asm_ids]
        C_su, C_aa, C_fa, C_va, C_bu, C_pro, C_ac, C_pr, C_li, C_ch = cmps_adm.i_C[_adm_ids]
        
        # For nitrogen balance 
        N_SF, N_XS, N_SA = cmps_asm.i_N[_asm_ids]
        N_su, N_aa, N_fa, N_va, N_bu, N_pro, N_ac, N_pr, N_li, N_ch = cmps_adm.i_N[_adm_ids]
        
        # For phosphorous balance 
        P_SF, P_XS, P_SA = cmps_asm.i_P[_asm_ids]
        P_su, P_aa, P_fa, P_va, P_bu, P_pro, P_ac, P_pr, P_li, P_ch = cmps_adm.i_P[_adm_ids]

        adm = self.adm1_model
        adm_p1_idx = cmps_adm.indices(('X_su', 'X_aa', 'X_fa', 'X_c4', 
                                       'X_pro', 'X_ac', 'X_h2', 
                                       'X_PAO', 'X_PP', 'X_PHA'))
        decay_idx = [i for i in adm.IDs if i.startswith(('decay', 'lysis'))]
        decay_stoichio = np.asarray(adm.stoichiometry.loc[decay_idx])
        # f_corr = self.balance_cod_tkn

        # To convert components from ADM1p to ASM2d (A1)
        def adm1p2masm2d(adm_vals):    
            
            _adm_vals = adm_vals.copy()
                
            # PROCESS 1: decay of biomas, X_PP, X_PHA
            bio_xpp_pha = _adm_vals[adm_p1_idx]
            _adm_vals += np.dot(bio_xpp_pha, decay_stoichio)
            
            # PROCESS 2: strip biogas. Omitted because no S_ch4 or S_h2 in ASM2d components
            
            S_su, S_aa, S_fa, S_va, S_bu, S_pro, S_ac, S_h2, S_ch4, S_IC, S_IN, S_IP, S_I, \
                X_ch, X_pr, X_li, X_su, X_aa, X_fa, X_c4, X_pro, X_ac, X_h2, X_I, \
                    X_PHA, X_PP, X_PAO, S_K, S_Mg, S_Ca, X_CaCO3, X_struv, \
                        X_newb, X_ACP, X_MgCO3, X_AlOH, X_AlPO4, X_FeOH, X_FePO4, \
                            S_Na, S_Cl, H2O = _adm_vals     
            
            if S_h2 > 0 or S_ch4 > 0: warn('Ignored dissolved H2 or CH4.')

            S_NH4 = S_IN
            S_PO4 = S_IP

            # CONV 1: convert X_pr, X_li, X_ch to X_S
            X_S = X_pr + X_li + X_ch
            S_IC += X_pr*C_pr + X_li*C_li + X_ch*C_ch - X_S*C_XS
            S_NH4 += X_pr*N_pr + X_li*N_li + X_ch*N_ch - X_S*N_XS
            S_PO4 += X_pr*P_pr + X_li*P_li + X_ch*P_ch - X_S*P_XS

            # CONV 2: convert S_su, S_aa, S_fa to S_F
            S_F = S_su + S_aa + S_fa
            S_IC += S_su*C_su + S_aa*C_aa + S_fa*C_fa - S_F*C_SF
            S_NH4 += S_su*N_su + S_aa*N_aa + S_fa*N_fa - S_F*N_SF
            S_PO4 += S_su*P_su + S_aa*P_aa + S_fa*P_fa - S_F*P_SF
            
            # CONV 3: convert VFAs to S_A
            S_A = S_va + S_bu + S_pro + S_ac
            S_IC += S_va*C_va + S_bu*C_bu + S_pro*C_pro + S_ac*C_ac - S_A*C_SA
            # S_NH4 += S_va*N_va + S_bu*N_bu + S_pro*N_pro + S_ac*N_ac - S_A*N_SA
            # S_PO4 += S_va*P_va + S_bu*P_bu + S_pro*P_pro + S_ac*P_ac - S_A*P_SA

            asm_vals = np.array(([
                0, 0, # S_O2, S_N2,
                S_NH4, 
                0, # S_NO3
                S_PO4, S_F, S_A, S_I, 
                S_IC, S_K, S_Mg,
                X_I, X_S, 
                0,0,0,0,0, # X_H, X_PAO, X_PP, X_PHA, X_AUT,
                S_Ca, X_CaCO3, X_struv, X_newb, X_ACP, X_MgCO3, # directly mapped
                X_AlOH, X_AlPO4, X_FeOH, X_FePO4, S_Na, S_Cl, H2O])) 
                      
            # asm_vals = f_corr(adm_vals, asm_vals)
            return asm_vals
        
        self._reactions = adm1p2masm2d

#%% mASM2dtoADM1p

class mASM2dtoADM1p(A1junction):
    '''
    Interface unit to convert ASM2d state variables
    to ADM1 components, following the A1 scenario in [1]_.
    
    Parameters
    ----------
    upstream : stream or str
        Influent stream with ASM components.
    downstream : stream or str
        Effluent stream with ADM components.
    adm1_model : :class:`qsdsan.processes.ADM1_p_extension`
        The anaerobic digestion process model.
    xs_to_li : float
        Split of slowly biodegradable substrate COD to lipid, 
        after all N is mapped into protein.
    rtol : float
        Relative tolerance for COD and TKN balance.
    atol : float
        Absolute tolerance for COD and TKN balance.
    
    References
    ----------
    [1] Flores-Alsina, X., Solon, K., Kazadi Mbamba, C., Tait, S., Gernaey, K. V., 
    Jeppsson, U., & Batstone, D. J. (2016). Modelling phosphorus (P), sulfur (S) 
    and iron (FE) interactions for dynamic simulations of anaerobic digestion processes. 
    Water Research, 95, 370–382. 
    
    See Also
    --------
    :class:`qsdsan.sanunits.A1junction`
    
    :class:`qsdsan.sanunits.ADM1ptomASM2d` 
    
    '''    
    # User defined values
    xs_to_li = 0.6
    
    def balance_cod_tkn_tp(self, asm_vals, adm_vals):
        cmps_asm = self.ins[0].components
        cmps_adm = self.outs[0].components
        asm_i_COD = cmps_asm.i_COD
        adm_i_COD = cmps_adm.i_COD
        non_tkn_idx = cmps_asm.indices(('S_N2', 'S_NO3'))
        asm_i_N = cmps_asm.i_N
        
        adm_i_N = cmps_adm.i_N
        asm_i_P = cmps_asm.i_P
        adm_i_P = cmps_adm.i_P
        asm_cod = sum(asm_vals*asm_i_COD)
        asm_tkn = sum(asm_vals*asm_i_N) - sum(asm_vals[non_tkn_idx])
        asm_tp = sum(asm_vals*asm_i_P)
        cod_bl, cod_err, cod_tol, adm_cod = self.isbalanced(asm_cod, adm_vals, adm_i_COD)
        tkn_bl, tkn_err, tkn_tol, adm_tkn = self.isbalanced(asm_tkn, adm_vals, adm_i_N)
        tp_bl, tp_err, tp_tol, adm_tp = self.isbalanced(asm_tp, adm_vals, adm_i_P)
        
        if tkn_bl and tp_bl:
            if cod_bl:
                return adm_vals
            else:
                if cod_err > 0: dcod = -(cod_err - cod_tol)/adm_cod
                else: dcod = -(cod_err + cod_tol)/adm_cod
                _adm_vals = adm_vals * (1 + (adm_i_COD>0)*dcod)
                _tkn_bl, _tkn_err, _tkn_tol, _adm_tkn = self.isbalanced(asm_tkn, _adm_vals, adm_i_N)
                _tp_bl, _tp_err, _tp_tol, _adm_tp = self.isbalanced(asm_tp, _adm_vals, adm_i_P)
                if _tkn_bl and _tp_bl: return _adm_vals
                else: 
                    warn('cannot balance COD, TKN, and TP at the same \n'
                        f'time with rtol={self.rtol} and atol={self.atol}.\n '
                        f'influent (ASM) TKN is {asm_tkn}\n '
                        f'effluent (ADM) TKN is {adm_tkn} or {_adm_tkn}\n '
                        f'influent TP is {asm_tp}\n ' 
                        f'effluent TP is {adm_tp} or {_adm_tp}. '
                        f'influent COD is {asm_cod}\n ' 
                        f'effluent COD is {adm_cod} or {adm_cod*(1+dcod)}. ')
                    return adm_vals
        elif cod_bl and tp_bl:
            if tkn_bl:
                return adm_vals
            else:
                if tkn_err > 0: dtkn = -(tkn_err - tkn_tol)/adm_tkn
                else: dtkn = -(tkn_err + tkn_tol)/adm_tkn
                _adm_vals = adm_vals * (1 + (adm_i_N>0)*dtkn)
                _cod_bl, _cod_err, _cod_tol, _adm_cod = self.isbalanced(asm_cod, _adm_vals, adm_i_COD)
                _tp_bl, _tp_err, _tp_tol, _adm_tp = self.isbalanced(asm_tp, _adm_vals, adm_i_P)
                if _cod_bl and _tp_bl: return _adm_vals
                else: 
                    warn('cannot balance COD, TKN, and TP at the same time'
                        f'time with rtol={self.rtol} and atol={self.atol}.\n '
                        f'influent (ASM) COD is {asm_cod}\n '
                        f'effluent (ADM) COD is {adm_cod} or {_adm_cod}\n '
                        f'influent TP is {asm_tp}\n ' 
                        f'effluent TP is {adm_tp} or {_adm_tp}. '
                        f'influent TKN is {asm_tkn}\n ' 
                        f'effluent TKN is {adm_tkn} or {adm_tkn*(1+dtkn)}. '
                        'To balance TKN please ensure ASM2d(X_I.i_N) = ADM1(X_I.i_N)')
                    return adm_vals
        elif cod_bl and tkn_bl: 
            if tp_bl:
                return adm_vals
            else:
                if tp_err > 0: dtp = -(tp_err - tp_tol)/adm_tp
                else: dtp = -(tp_err + tp_tol)/adm_tp
                _adm_vals = adm_vals * (1 + (adm_i_P>0)*dtp)
                _cod_bl, _cod_err, _cod_tol, _adm_cod = self.isbalanced(asm_cod, _adm_vals, adm_i_COD)
                _tkn_bl, _tkn_err, _tkn_tol, _adm_tkn = self.isbalanced(asm_tkn, _adm_vals, adm_i_N)
                if _cod_bl and _tkn_bl: return _adm_vals
                else: 
                    warn('cannot balance COD, TKN, and TP at the same time'
                        f'time with rtol={self.rtol} and atol={self.atol}.\n '
                        f'influent (ASM) COD is {asm_cod}\n '
                        f'effluent (ADM) COD is {adm_cod} or {_adm_cod}\n '
                        f'influent TKN is {asm_tkn}\n ' 
                        f'effluent TKN is {adm_tkn} or {_adm_tkn}. '
                        f'influent TP is {asm_tp}\n ' 
                        f'effluent TP is {adm_tp} or {adm_tp*(1+dtp)}. ')
                    return adm_vals
        else:
            warn('cannot balance COD, TKN, and TP at the same time. \n'
                 'Atleast two of the three COD, TKN, and TP are not balanced \n'
                f'time with rtol={self.rtol} and atol={self.atol}.\n '
                f'influent (ASM) COD is {asm_cod}\n '
                f'effluent (ADM) COD is {adm_cod}\n '
                f'influent TP is {asm_tp}\n ' 
                f'effluent TP is {adm_tp}'
                f'influent TKN is {asm_tkn}\n ' 
                f'effluent TKN is {adm_tkn}. ')
            return adm_vals
    
    def _compile_reactions(self):
        cmps_asm = self.ins[0].components
        cmps_adm = self.outs[0].components
        self.check_component_properties(cmps_asm, cmps_adm)
        
        _asm_ids = cmps_asm.indices(['S_F', 'X_S', 'X_H', 'S_I', 'X_I'])
        _adm_ids = cmps_adm.indices(['S_aa', 'S_su', 'X_pr', 'X_li', 'X_ch'])
        # For carbon balance
        C_SF, C_XS, C_XB, C_SI, C_XI = cmps_asm.i_C[_asm_ids]
        C_aa, C_su, C_pr, C_li, C_ch = cmps_adm.i_C[_adm_ids]
        
        # For nitrogen balance 
        N_SF, N_XS, N_XB, N_SI, N_XI = cmps_asm.i_N[_asm_ids]
        N_aa, N_su, N_pr, N_li, N_ch = cmps_adm.i_N[_adm_ids]
        
        # For phosphorous balance 
        P_SF, P_XS, P_XB, P_SI, P_XI = cmps_asm.i_P[_asm_ids]
        P_aa, P_su, P_pr, P_li, P_ch = cmps_adm.i_P[_adm_ids]

        S_O2_idx, S_NO3_idx = cmps_asm.indices(['S_O2', 'S_NO3'])
        # f_corr = self.balance_cod_tkn_tp

        asm = self.asm2d_model
        adm = self.adm1_model
        p1_stoichio = np.asarray(asm.stoichiometry.loc['hetero_growth_S_A'])
        p1_stoichio /= abs(p1_stoichio[S_O2_idx])
        p2_stoichio = np.asarray(asm.stoichiometry.loc['denitri_S_A'])
        p2_stoichio /= abs(p2_stoichio[S_NO3_idx])
        p3_stoichio = np.array([adm.parameters[f'f_{k}_xb'] for k in ('sI', 'ch', 'pr', 'li', 'xI')])
        
        xs_to_li = self.xs_to_li
        
        # To convert components from ASM2d to mADM1 (asm2d-2-madm1)
        def masm2d2adm1p(asm_vals):           
            _asm_vals = asm_vals.copy()

            # PROCESS 1: remove S_O2 with S_A with associated X_H growth (aerobic growth of X_H on S_A)
            O2_coddm = _asm_vals[S_O2_idx]
            _asm_vals += O2_coddm * p1_stoichio # makes S_O2 = 0
            
            # PROCESS 2: remove S_NO3 with S_A with associated X_H growth (denitrification on S_A)
            NO3_coddm = _asm_vals[S_NO3_idx]
            _asm_vals += NO3_coddm * p2_stoichio # makes S_NO3 = 0
                      
            S_O2, S_N2, S_NH4, S_NO3, S_PO4, S_F, S_A, S_I, S_IC, S_K, S_Mg, \
                X_I, X_S, X_H, X_PAO, X_PP, X_PHA, X_AUT, S_Ca, X_CaCO3, \
                    X_struv, X_newb, X_ACP, X_MgCO3, X_AlOH, X_AlPO4, \
                        X_FeOH, X_FePO4, S_Na, S_Cl, H2O = _asm_vals
                
            S_IN = S_NH4
            S_IP = S_PO4
            
            # CONV 1: transform S_F into S_aa, S_su, S_fa
            S_ND = S_F*N_SF   # N in S_F
            req_scod = S_ND / N_aa
                       
            if S_F < req_scod: # if S_F cod is not enough to convert all organic soluble N into aa
                S_aa = S_F
                S_su = 0
            else: # if S_F cod is more than enough to convert all organic soluble N into aa
                S_aa = req_scod # All soluble organic N will be mapped to amino acid
                S_su = S_F - S_aa
            
            S_IN += S_ND - S_aa*N_aa
            S_IC += S_F*C_SF - (S_aa*C_aa + S_su*C_su) 
            S_IP += S_F*P_SF
            
            # PROCESS 3: biomass decay (X_H, X_AUT lysis) anaerobic
            bio = X_H + X_AUT
            _si, _ch, _pr, _li, _xi = bio * p3_stoichio
            S_IC += bio*C_XB - (_si*C_SI + _ch*C_ch + _pr*C_pr + _li*C_li + _xi*C_XI)
            S_IN += bio*N_XB - (_si*N_SI + _ch*N_ch + _pr*N_pr + _li*N_li + _xi*N_XI)
            S_IP += bio*P_XB - (_si*P_SI + _ch*P_ch + _pr*P_pr + _li*P_li + _xi*P_XI)
                        
            # CONV 2: transform asm X_S into X_pr, X_li, X_ch
            X_ND = X_S*N_XS
            req_xcod = X_ND / N_pr
            # Since X_pr_i_N >> X_pr_i_P there's no need to check req_xcod for N and P separately (CONFIRM LATER 05/16)
            
            # if available X_S is not enough to fulfill that protein requirement
            if X_S < req_xcod: # if X_S cod is not enough to convert all organic particulate N into pr
                X_pr = X_S
                X_li = X_ch = 0
            else:
                X_pr = req_xcod 
                X_li = xs_to_li * (X_S - X_pr)
                X_ch = (X_S - X_pr) - X_li

            S_IN += X_ND - X_pr*N_pr
            S_IC += X_S*C_XS - (X_pr*C_pr + X_li*C_li + X_ch*C_ch)
            S_IP += X_S*P_XS - (X_pr*P_pr + X_li*P_li + X_ch*P_ch)

            X_pr += _pr
            X_li += _li
            X_ch += _ch
            S_I += _si
            X_I += _xi

            # PROCESS 4-5: omitted, PAO related components mapped directly
            # CONV 3-5: convert S_A, S_I, X_I; conversion is immediate because identical component composition is enforced
            S_ac = S_A
            
            adm_vals = np.array([
                S_su, S_aa, 
                0, 0, 0, 0, S_ac, # S_fa, S_va, S_bu, S_pro,
                0, 0, # S_h2, S_ch4,
                S_IC, S_IN, S_IP, S_I, 
                X_ch, X_pr, X_li, 
                0, 0, 0, 0, 0, 0, 0, # X_su, X_aa, X_fa, X_c4, X_pro, X_ac, X_h2,
                X_I, X_PHA, X_PP, X_PAO, 
                S_K, S_Mg, S_Ca, X_CaCO3, X_struv, X_newb, X_ACP, X_MgCO3, 
                X_AlOH, X_AlPO4, X_FeOH, X_FePO4, S_Na, S_Cl, H2O])
            
            # adm_vals = f_corr(asm_vals, adm_vals)
            # adm_vals = f_corr(_asm_vals, adm_vals)
            return adm_vals
        
        self._reactions = masm2d2adm1p
