#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Joy Cheung <joycheung1994@gmail.com>
    Yalin Li <zoe.yalin.li@gmail.com>

Part of this module is based on the Thermosteam package:
https://github.com/BioSTEAMDevelopmentGroup/thermosteam

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

'''
TODO:
    Check compatibility for multiphase for phase-equilibrium
'''



# %%

import numpy as np
from thermosteam import utils, settings
from . import Components, Stream, MultiStream, SanStream, MissingSanStream, \
    set_thermo
from ._units_of_measure import auom

__all__ = ('WasteStream', 'MissingWasteStream')


_defined_composite_vars = ('COD', 'BOD5', 'BOD', 'uBOD', 'NOD', 'ThOD', 'cnBOD',
                           'C', 'N', 'P', 'K', 'Mg', 'Ca', 'solids', 'charge')

_common_composite_vars = ('_COD', '_BOD', '_uBOD', '_TC', '_TOC', '_TN', 
                          '_TKN', '_TP', '_TK', '_TMg', '_TCa', 
                          '_dry_mass', '_charge', '_ThOD', '_cnBOD')

_ws_specific_slots = (*_common_composite_vars,
                      '_pH', '_SAlk', '_ratios', '_impact_item')

_specific_groups = {'S_VFA': ('S_Ac', 'S_Prop'),
                    'X_Stor': ('X_OHO_PHA', 'X_GAO_PHA', 'X_PAO_PHA', 
                              'X_GAO_Gly', 'X_PAO_Gly'),
                    'X_ANO': ('X_AOO', 'X_NOO'),
                    'X_Bio': ('X_OHO', 'X_AOO', 'X_NOO', 'X_AMO', 'X_PAO', 
                             'X_MEOLO', 'X_ACO', 'X_HMO', 'X_PRO', 'X_FO'),
                    'S_NOx': ('S_NO2', 'S_NO3'),
                    'X_PAO_PP': ('X_PAO_PP_Lo', 'X_PAO_PP_Hi'),
                    'TKN': ()}

_default_ratios = {'iHi_XPAOPP': 0,
                   'iCB_XCB': 0.15,
                   'iBAP_CB': 0.,
                   'iUAP_CB': 0.,
                   'iCUInf_XCUInf': 0.,
                   'iSUInf_SU': 1.,
                   'iXUOHOE_XUE': None,}


vol_unit = auom('L/hr')
conc_unit = auom('mg/L')


#%% functions for calibrations of N, P contents and BOD:COD ratio of certain components
        
def _calib_SF_iN(components, concentrations, STKN):    
    cmp_c = np.asarray([v for v in concentrations.values()])
    SN = (cmp_c * components.i_N * (components.s + components.c)).sum()
    SF_N = concentrations['S_F'] * components.S_F.i_N
    SNOx_N = concentrations['S_NO2'] * components.S_NO2.i_N + concentrations['S_NO3'] * components.S_NO3.i_N
    other_stkn = SN - SF_N - SNOx_N
    SF_N = STKN - other_stkn           
    if SF_N < 0:
        raise ValueError("Negative N content for S_F was estimated.")                        
    return SF_N/concentrations['S_F']
    
def _calib_XBsub_iN(components, concentrations, XTKN):
    cmp_c = np.asarray([v for v in concentrations.values()])        
    other_xtkn = (cmp_c * components.i_N * components.x).sum() - concentrations['X_B_Subst'] * components.X_B_Subst.i_N                
    XB_Subst_N = XTKN - other_xtkn
    if XB_Subst_N < 0:
        raise ValueError("Negative N content for X_B_Subst was estimated.")            
    return XB_Subst_N/concentrations['X_B_Subst']
    

def _calib_XBsub_iP(components, concentrations, TP):
    cmp_c = np.asarray([v for v in concentrations.values()])    
    other_p = (cmp_c * components.i_P).sum() - concentrations['X_B_Subst'] * components.X_B_Subst.i_P
    XB_Subst_P = TP - other_p
    if XB_Subst_P < 0:
        raise ValueError("Negative P content for X_B_Subst was estimated.")    
    return XB_Subst_P/concentrations['X_B_Subst']
    
def _calib_XBsub_fBODCOD(components, concentrations, substrate_IDs, BOD):
    cmp_c = np.asarray([v for v in concentrations.values()])
    c_sub = np.asarray([v for k,v in concentrations.items() if k in substrate_IDs])
    XB_sub = components.subgroup(substrate_IDs)
    other_BOD = (cmp_c * (components.x + components.c + components.s) * components.f_BOD5_COD).sum() - (c_sub * XB_sub.f_BOD5_COD).sum()
    fbodtocod_sub = (BOD - other_BOD)/c_sub.sum()
    if fbodtocod_sub > 1 or fbodtocod_sub < 0:
        raise ValueError("BOD5-to-COD ratio for X_B_Subst and X_Stor was estimated out of range [0,1].")
    return fbodtocod_sub



# %%

# =============================================================================
# Define the WasteStream class
# =============================================================================

_load_components = settings.get_default_chemicals

class WasteStream(SanStream):
    '''
    A subclass of :class:`~.SanStream` with additional attributes and methods
    for wastewater.
    
    See Also
    --------
    `thermosteam.Stream <https://thermosteam.readthedocs.io/en/latest/Stream.html>`_
    '''
    
    __slots__ = SanStream.__slots__ + _ws_specific_slots
    _default_ratios = _default_ratios
    ticket_name = 'ws'
    
    def __init__(self, ID='', flow=(), phase='l', T=298.15, P=101325.,
                 units='kg/hr', price=0., thermo=None, 
                 pH=7., SAlk=2.5, COD=None, BOD=None, uBOD=None,
                 TC=None, TOC=None, TN=None, TKN=None, TP=None, TK=None,
                 TMg=None, TCa=None, dry_mass=None, charge=None, ratios=None,
                 ThOD=None, cnBOD=None, impact_item=None, **chemical_flows):
        
        SanStream.__init__(self=self, ID=ID, flow=flow, phase=phase, T=T, P=P,
                           units=units, price=price, thermo=thermo,
                           impact_item=impact_item, **chemical_flows)
        
        self._init_ws(pH, SAlk, COD, BOD, uBOD, TC, TOC, TN, TKN,
                      TP, TK, TMg, TCa, ThOD, cnBOD, dry_mass, charge, ratios)

    def _init_ws(self, pH=7., SAlk=None, COD=None, BOD=None,
                 uBOD=None, TC=None, TOC=None, TN=None, TKN=None,
                 TP=None, TK=None, TMg=None, TCa=None, ThOD=None, cnBOD=None,
                 dry_mass=None, charge=None, ratios=None):

        self._pH = pH
        self._SAlk = SAlk
        self._COD = COD
        self._BOD = BOD
        self._uBOD = uBOD
        self._TC = TC
        self._TOC = TOC
        self._TN = TN
        self._TKN = TKN
        self._TP = TP
        self._TK = TK
        self._TMg = TMg
        self._TCa = TCa
        self._ThOD = ThOD
        self._cnBOD = cnBOD
        self._dry_mass = dry_mass
        self._charge = charge
        self._ratios = ratios

    @staticmethod
    def from_stream(cls, stream, **kwargs):
        '''
        Cast a :class:`thermosteam.Stream` or :class:`biosteam.utils.MissingStream`
        to :class:`WasteStream` or :class:`MissingWasteStream`.
        '''
        
        new = SanStream.from_stream(cls, stream, **kwargs)
        if isinstance(new, MissingSanStream):
            missing_new = MissingWasteStream.__new__(MissingWasteStream)
            return missing_new
        
        new._init_ws()
        
        # if isinstance(new, MissingSanStream):
        #     setattr(new, '__class__', MissingWasteStream)
        # else:
        #     new._init_ws()
            
        return new

    
    def show(self, T='K', P='Pa', flow='g/hr', composition=False, N=15,
             stream_info=True, details=True):
        '''
        Print waste stream information.

        Parameters
        ----------
        T : str, optional
            The unit for temperature. The default is 'K'.
        P : float, optional
            The unit for pressure. The default is 'Pa'.
        flow : str, optional
            The unit for the flow. The default is 'kg/hr'.
        composition : bool, optional
            Whether to show flow information of different :class:`Component` objects in
            this waste stream as a percentage. The default is False.
        N : int, optional
            Number of components to print out, when left as None,
            the number depends on the default of :class:`thermosteam`. The default is 15.
        stream_info : bool, optional
            Whether to print stream-specific information. The default is True.
        details : bool, optional
            Whether to show the all composite variables of this waste stream. The default is True.

        '''

        info = ''
        # Stream-related specifications
        if stream_info:
            super().show(T, P, flow, composition, N)
        else:
            info += self._basic_info()
            display_units = self.display_units
            T_units = T or display_units.T
            P_units = P or display_units.P
            info += self._info_phaseTP(self.phase, T_units, P_units)
        info += self._wastestream_info(details=details)
        print(info)
        
    _ipython_display_ = show
    
    
    def _wastestream_info(self, details=True):
        _ws_info = ' WasteStream-specific properties:'
        # Wastewater-related properties are not relevant for gas or solids
        if self.phase != 'l':
            _ws_info += ' None for non-liquid waste streams\n'
        elif self.F_mass == 0:
            _ws_info += ' None for empty waste streams\n'
        else:
            _ws_info += '\n'
            # Only non-zero properties are shown
            _ws_info += int(bool(self.pH))*f'  pH         : {self.pH:.1f}\n'
            _ws_info += int(bool(self.SAlk))*f'  Alkalinity : {self.SAlk:.1f} mg/L\n'
            if details:
                _ws_info += int(bool(self.COD))   *f'  COD        : {self.COD:.1f} mg/L\n'
                _ws_info += int(bool(self.BOD))   *f'  BOD        : {self.BOD:.1f} mg/L\n'
                _ws_info += int(bool(self.TC))    *f'  TC         : {self.TC:.1f} mg/L\n'
                _ws_info += int(bool(self.TOC))   *f'  TOC        : {self.TOC:.1f} mg/L\n'
                _ws_info += int(bool(self.TN))    *f'  TN         : {self.TN:.1f} mg/L\n'
                _ws_info += int(bool(self.TKN))   *f'  TKN        : {self.TKN:.1f} mg/L\n'
                _ws_info += int(bool(self.TP))    *f'  TP         : {self.TP:.1f} mg/L\n'
                _ws_info += int(bool(self.TK))    *f'  TK         : {self.TK:.1f} mg/L\n'
                # _ws_info += int(bool(self.charge))*f'  charge     : {self.charge:.1f} mmol/L\n'
            else:
                _ws_info += '  ...\n'
            
        return _ws_info

    @property
    def ratios(self):
        '''
        The ratios used for estimating waste stream composition based on user input upon initialization.
        Only meaningful for creating a :class:`WasteStream` object from scratch.
        If not used or specified, default as None.
        '''
        return self._ratios
    @ratios.setter
    def ratios(self, ratios):
        r = self._ratios or WasteStream._default_ratios
        for name, ratio in ratios.items():
            if name not in r.keys():
                raise ValueError(f'Cannot identify ratio named "{name}".'
                                 f'Must be one of {r.keys()}.')
            elif isinstance(ratio, (int, float)) and (ratio > 1 or ratio < 0):
                raise ValueError(f"ratio {name}: {ratio} is out of range [0,1].")
            r[name] = ratio
        self._ratios = r

    def composite(self, variable, subgroup=None, particle_size=None, 
                  degradability=None, organic=None, volatile=None,
                  specification=None):
        
        """
        Calculate any composite variable by specifications.

        Parameters
        ----------
        variable : str
            The composite variable to calculate. One of the followings:
                ("COD", "BOD5", "BOD", "uBOD", "NOD", "ThOD", "cnBOD",
                "C", "N", "P", "K", "Mg", "Ca", 
                "solids", "charge").
        subgroup : CompiledComponents, optional
            A subgroup of :class:`CompiledComponents`. The default is None.
        particle_size : "g", "s", "c", or "x", optional 
            Dissolved gas ("g"), soluble ("s"), colloidal ("c"), particulate ("x"). 
            The default is None.
        degradability : "rb", "sb", "b" or "u", optional
            Readily biodegradable ("rb"), slowly biodegradable ("sb"), 
            biodegradable ("b"), or undegradable ("u"). The default is None.
        organic : bool, optional
            Organic (True) or inorganic (False). The default is None.
        volatile : bool, optional
            Volatile (True) or involatile (False). The default is None.
        specification : str, optional
            One of ("SVFA", "XStor", "XANO", "XBio", "SNOx", "XPAO_PP", "TKN"). 
            The default is None.

        Returns
        -------
        value : float
            The estimated value of the composite variable, in [mg/L] or [mmol/L] (for "Charge").

        """
        _get = getattr
        if self.F_vol == 0.:
            return 0.
        
        if variable not in _defined_composite_vars:
            raise KeyError(f"Undefined composite variable {variable},"
                           f"Must be one of {_defined_composite_vars}.")            
        
        #!!! assuming it's a liquid WasteStream
        #TODO: deal with units
        if subgroup:
            cmps = subgroup
            
        else:
            cmps = self.components

        IDs = list(cmps.IDs)            
        if 'H2O' in IDs: IDs.remove('H2O')

        if specification:
            if specification == 'TKN': 
                IDs = [ID for ID in IDs if ID not in ('S_N2','S_NO2','S_NO3')]
            elif specification not in _specific_groups.keys():
                raise KeyError(f"Undefined specification {specification}."
                               f"Must be one of {_specific_groups.keys()}."
                               "Or, try defining 'subgroup'.")
            else: 
                IDs = [ID for ID in IDs if ID in _specific_groups[specification]]
                
        IDs = tuple(IDs)
        cmps = cmps.subgroup(IDs)
        cmp_c = self.imass[IDs]/self.F_vol*1e3      #[mg/L]
        exclude_gas = _get(cmps, 's')+_get(cmps, 'c')+_get(cmps, 'x')        
        
        if variable == 'COD': 
            var = cmps.i_COD * cmp_c * exclude_gas * (cmps.i_COD >= 0)
        elif variable == 'uBOD': 
            var = cmps.i_COD * cmps.f_uBOD_COD * cmp_c * exclude_gas * (cmps.i_COD >= 0)
        elif variable in ('BOD5', 'BOD'): 
            var = cmps.i_COD * cmps.f_BOD5_COD * cmp_c * exclude_gas * (cmps.i_COD >= 0)
        elif variable == 'NOD':
            var = cmps.i_NOD * cmp_c * exclude_gas
        elif variable == 'ThOD':
            var = (cmps.i_NOD + cmps.i_COD * (cmps.i_COD >= 0)) * cmp_c
        elif variable == 'cnBOD':
            var = (cmps.i_NOD + cmps.i_COD * cmps.f_BOD5_COD * (cmps.i_COD >= 0)) * cmp_c * exclude_gas
        elif variable == 'C':
            var = cmps.i_C * cmp_c
        elif variable == 'N':
            var = cmps.i_N * cmp_c * exclude_gas
        elif variable == 'P':
            var = cmps.i_P * cmp_c
        elif variable == 'K':
            var = cmps.i_K * cmp_c
        elif variable == 'Mg':
            var = cmps.i_Mg * cmp_c
        elif variable == 'Ca':
            var = cmps.i_Ca * cmp_c
        elif variable == 'solids':
            var = cmps.i_mass * cmp_c * exclude_gas
            if volatile != None:
                if volatile: var *= cmps.f_Vmass_Totmass
                else: var *= 1-cmps.f_Vmass_Totmass
        else:
            var = cmps.i_charge * cmp_c
        
        dummy = np.ones(len(cmp_c))
        if particle_size:
            if particle_size == 'g': 
                dummy *= 1-exclude_gas
            else:
                dummy *= _get(cmps, particle_size)
        
        if degradability:
            if degradability == 'u': dummy *= 1-_get(cmps, 'b')
            elif degradability == 'b': dummy *= _get(cmps, 'b')
            elif degradability == 'rb': dummy *= _get(cmps, 'rb')
            else: dummy *= _get(cmps, 'b')-_get(cmps, 'rb')

        if organic != None:
            if organic: dummy *= _get(cmps, 'org')
            else: dummy *= 1-_get(cmps, 'org')
        
        return (dummy*var).sum()

    
    def _liq_sol_properties(self, prop, value):
        if self.phase != 'g':
            return getattr(self, '_'+prop) or value
        else:
            raise AttributeError(f'{self.phase} phase waste stream does not have {prop}.')
    
    @property
    def pH(self):
        '''[float] pH, unitless.'''
        return self._liq_sol_properties('pH', 7.)

    @property
    def SAlk(self):
        '''[float] Alkalinity in meq/L (or mmol HCO3-/L). Assumed to be mainly bicarbonate.'''
        return self._liq_sol_properties('SAlk', 0.)

    @property
    def COD(self):
        '''[float] Chemical oxygen demand in mg/L.'''
        return self._liq_sol_properties('COD', self.composite('COD'))

    @property    
    def BOD(self):
        '''[float] Biochemical oxygen demand in mg/L. Same as BOD5.'''
        return self._liq_sol_properties('BOD', self.composite('BOD'))

    @property    
    def BOD5(self):
        '''[float] 5-day biochemical oxygen demand, in mg/L. Same as BOD.'''
        return self.BOD

    @property    
    def uBOD(self):
        '''[float] Ultimate biochemical oxygen demand, in mg/L.'''
        return self._liq_sol_properties('uBOD', self.composite('uBOD'))

    @property    
    def cnBOD(self):
        '''[float] Carbonaceous nitrogenous BOD, in mg/L. Biochemical oxygen demand including nitrification.'''
        return self._liq_sol_properties('cnBOD', self.composite('cnBOD'))

    @property    
    def ThOD(self):
        '''[float] Theoretical oxygen demand, in mg/L.'''
        return self._liq_sol_properties('ThOD', self.composite('ThOD'))
    
    #!!! Maybe include C_frac, etc. to calculate C_mass/F_mass - valid for all phases
    # Or a function to calculate it?
    @property
    def TC(self):
        '''[float] Total carbon, in mg/L.'''
        return self._liq_sol_properties('TC', self.composite('C'))
    
    @property
    def TOC(self):
        '''[float] Total organic carbon, in mg/L.'''
        return self._liq_sol_properties('TOC', self.composite('C', organic=True))
        
    @property
    def TN(self):
        '''[float] Total nitrogen, in mg/L.'''
        return self._liq_sol_properties('TN', self.composite('N'))
    
    @property
    def TKN(self):
        '''[float] Total Kjeldahl nitrogen, in mg/L.'''
        return self._liq_sol_properties('TKN', self.composite('N', specification='TKN'))
    
    @property
    def TP(self):
        '''[float] Total phosphorus, in mg/L.'''
        return self._liq_sol_properties('TP', self.composite('P'))
    
    @property
    def TK(self):
        '''[float] Total potassium, in mg/L.'''
        return self._liq_sol_properties('TK', self.composite('K'))
    
    @property
    def TMg(self):
        '''[float] Total magnesium, in mg/L.'''
        return self._liq_sol_properties('TMg', self.composite('Mg'))
    
    @property
    def TCa(self):
        '''[float] Total calcium, in mg/L.'''
        return self._liq_sol_properties('TCa', self.composite('Ca'))
    
    @property
    def dry_mass(self):
        '''[float] Total solids, dry mass of dissolved and suspended solids, in mg/L.'''
        return self._liq_sol_properties('solids', self.composite('solids'))
    
    # TODO: calibrate Charge when weak acids are involved
    # @property
    # def charge(self):
    #     return self._liq_sol_properties('charge', self.composite('charge'))
    

    def copy(self, ID='', ws_properties=True):
        '''
        Copy the information of another stream.
        
        Parameters
        ----------
        ID : str
            ID of the new stream, a default ID will be assigned if not provided.
        ws_properties : bool
            Whether to copy wastewater-related properties to the new stream.
        
        
        .. note::
            
            [1] Price of the original stream is not copied.
            
            [2] If the original stream has an :class:`~.StreamImpactItem`,
            then a new :class:`~.StreamImpactItem` will be created for the new stream
            and the new impact item will be linked to the original impact item.
        '''
        new = Stream.copy(self)
        if ws_properties:
            new._init_ws()
            for slot in _ws_specific_slots:
                value = getattr(self, slot)
                if slot == '_impact_item' and value:
                    value.copy(stream=new)
                else:
                    setattr(new, slot, utils.copy_maybe(value))
        return new
    __copy__ = copy

    # TODO: add documents for these functions, differentiate copy, copy_like, and copy_flow
    def copy_like(self, other):
        Stream.copy_like(self, other)
        
        if not isinstance(other, WasteStream):
            return

        for slot in _ws_specific_slots:
            value = getattr(other, slot)
            setattr(self, slot, value)
    
    def copy_flow(self, other, IDs=..., *, remove=False, exclude=False, if_copy_ws=False):
        #!!! How to inherit the Stream copy_flow function?
        # Stream.copy_flow(self, other, IDs, remove, exclude)

        chemicals = self.chemicals
        mol = other.mol
        if exclude:
            IDs = chemicals.get_index(IDs)
            index = np.ones(chemicals.size, dtype=bool)
            index[IDs] = False
        else:
            index = chemicals.get_index(IDs)
        
        self.mol[index] = mol[index]
        if remove: 
            if isinstance(other, MultiStream):
                other.imol.data[:, index] = 0
            else:
                mol[index] = 0

        if if_copy_ws:
            for slot in _ws_specific_slots:
                value = getattr(other, slot)
                setattr(self, slot, value)

    def mix_from(self, others):
        others = [s for s in others if not 'Missing' in type(s).__name__]
        Stream.mix_from(self, others)


        for slot in _ws_specific_slots:
            #!!! This need reviewing, might not be good to calculate some
            # attributes like pH
            try:
                tot = sum(float(getattr(i, slot))*i.F_vol for i in others
                          if hasattr(i, slot))
            except: continue
            if tot == 0.:
                setattr(self, slot, None)
            else:
                setattr(self, slot, tot/self.F_vol)




    def get_TDS(self, include_colloidal=True):
        '''
        Total dissolved solids (TDS).

        Parameters
        ----------
        include_colloidal : bool, optional
            Whether to include colloidal components as TDS. The default is True.

        Returns
        -------
        TDS : float
            In mg/L.

        '''
        TDS = self.composite('solids', particle_size='s')
        if include_colloidal:
            TDS += self.composite('solids', particle_size='c')
        return TDS
    
    def get_TSS(self, include_colloidal=False):
        '''
        Total suspended solids (TSS).

        Parameters
        ----------
        include_colloidal : bool, optional
            Whether to include colloidal components as TSS. The default is False.

        Returns
        -------
        TSS : float
            In mg/L.

        '''
        TSS = self.composite('solids', particle_size='x')
        if include_colloidal:
            TSS += self.composite('solids', particle_size='c')        
        return TSS
    
    def get_VSS(self, include_colloidal=False):
        '''[float] Volatile suspended solids, in mg/L.'''
        VSS = self.composite('solids', particle_size='x', volatile=True)
        if include_colloidal:
            VSS += self.composite('solids', particle_size='c', volatile=True)        
        return VSS        
    
    def get_ISS(self):
        '''[float] Inorganic/involatile suspended solids, in mg/L.'''
        return self.composite('solids', particle_size='x', volatile=False)

    @classmethod
    def codstates_inf_model(cls, ID, flow_tot=0., units = ('L/hr', 'mg/L'), 
                            phase='l', T=298.15, P=101325., price=0., thermo=None, 
                            pH=7., SAlk=10., ratios=None, 
                            COD=430., TKN=40., TP=10., iVSS_TSS=0.75, iSNH_STKN=0.9,
                            S_NH4=25., S_NO2=0., S_NO3=0., S_PO4=8., 
                            S_Ca=140., S_Mg=50., S_K=28., S_CAT=3., S_AN=12., S_N2=18., 
                            frSUInf=0.05, frSF=0.2, frXCUInf=0.13, 
                            frSUE=0., frSCH3OH=0., frSAc=0., frSProp=0., 
                            frXOHO=0., frXAOO=0., frXNOO=0., frXAMO=0., frXPAO=0., 
                            frXPRO=0., frXACO=0., frXHMO=0., frXMEOLO=0., frXFO=0.,
                            frXOHO_PHA=0., frXGAO_PHA=0., frXPAO_PHA=0., 
                            frXGAO_Gly=0., frXPAO_Gly=0., frXU_OHO_E=0., frXU_PAO_E=0.,
                            X_FePO4=0., X_AlPO4=0., X_FeOH=0., X_AlOH=0., 
                            X_MAP=0., X_HAP=0., X_HDP=0., X_PAO_PP=0., 
                            X_MgCO3=0., X_CaCO3=0., DO=0., S_H2=0., S_CH4=0.):
        
        if thermo: cmps = thermo.chemicals
        else: cmps = _load_components(None)
                
        cmp_dct = dict.fromkeys(cmps.IDs, 0.)

        new = cls(ID=ID, phase=phase, T=T, P=P, units='kg/hr', price=price, 
                  thermo=thermo, pH=pH, SAlk=SAlk)

        if ratios: new.ratios = ratios
        else: new.ratios = WasteStream._default_ratios
        r = new._ratios

        #************ user-defined states **************        
        cmp_dct['S_H2'] = S_H2
        cmp_dct['S_CH4'] = S_CH4
        cmp_dct['S_N2'] = S_N2
        cmp_dct['S_O2'] = DO
        cmp_dct['S_NH4'] = S_NH4
        cmp_dct['S_NO2'] = S_NO2
        cmp_dct['S_NO3'] = S_NO3
        cmp_dct['S_PO4'] = S_PO4
        cmp_dct['S_CO3'] = SAlk * 12 * conc_unit.conversion_factor(units[1])             # 1 meq/L SAlk ~ 1 mmol/L HCO3- ~ 12 mg C/L (12 mg C/mmol HCO3-)
        cmp_dct['S_Ca'] = S_Ca
        cmp_dct['S_Mg'] = S_Mg
        cmp_dct['S_K'] = S_K
        cmp_dct['X_MAP'] = X_MAP
        cmp_dct['X_HAP'] = X_HAP
        cmp_dct['X_HDP'] = X_HDP
        cmp_dct['X_FePO4'] = X_FePO4
        cmp_dct['X_AlPO4'] = X_AlPO4
        cmp_dct['X_FeOH'] = X_FeOH
        cmp_dct['X_AlOH'] = X_AlOH
        cmp_dct['X_MgCO3'] = X_MgCO3
        cmp_dct['X_CaCO3'] = X_CaCO3
        cmp_dct['S_CAT'] = S_CAT
        cmp_dct['S_AN'] = S_AN
        
        #************ organic components **************
        cmp_dct['S_CH3OH'] = COD * frSCH3OH
        cmp_dct['S_Ac'] = COD * frSAc
        cmp_dct['S_Prop'] = COD * frSProp
        cmp_dct['S_F'] = COD * frSF
        cmp_dct['S_U_Inf'] = COD * frSUInf
        cmp_dct['S_U_E'] = COD * frSUE
        
        S_Org = sum([v for k,v in cmp_dct.items() if k in ('S_CH3OH','S_Ac','S_Prop','S_F','S_U_Inf','S_U_E')])
        
        XC_U_Inf = COD * frXCUInf
        cmp_dct['C_U_Inf'] = XC_U_Inf * r['iCUInf_XCUInf']
        cmp_dct['X_U_Inf'] = XC_U_Inf * (1 - r['iCUInf_XCUInf'])
        
        cmp_dct['X_OHO'] = COD * frXOHO
        cmp_dct['X_AOO'] = COD * frXAOO
        cmp_dct['X_NOO'] = COD * frXNOO
        cmp_dct['X_AMO'] = COD * frXAMO
        cmp_dct['X_PAO'] = COD * frXPAO
        cmp_dct['X_ACO'] = COD * frXACO
        cmp_dct['X_HMO'] = COD * frXHMO
        cmp_dct['X_PRO'] = COD * frXPRO
        cmp_dct['X_MEOLO'] = COD * frXMEOLO
        cmp_dct['X_FO'] = COD * frXFO
        
        X_Bio = sum([v for k,v in cmp_dct.items() if k.startswith('X_') and k.endswith('O')])
        
        cmp_dct['X_OHO_PHA'] = COD * frXOHO_PHA
        cmp_dct['X_GAO_PHA'] = COD * frXGAO_PHA
        cmp_dct['X_PAO_PHA'] = COD * frXPAO_PHA
        cmp_dct['X_GAO_Gly'] = COD * frXGAO_Gly
        cmp_dct['X_PAO_Gly'] = COD * frXPAO_Gly
        
        X_Stor = sum([v for k,v in cmp_dct.items() if k.endswith(('PHA','Gly'))])
        
        cmp_dct['X_U_OHO_E'] = COD * frXU_OHO_E
        cmp_dct['X_U_PAO_E'] = COD * frXU_PAO_E
        
        X_U_E = cmp_dct['X_U_OHO_E'] + cmp_dct['X_U_PAO_E']
        
        XC_B = COD - S_Org - XC_U_Inf - X_U_E
        C_B = XC_B * r['iCB_XCB']
        cmp_dct['C_B_BAP'] = C_B * r['iBAP_CB']
        cmp_dct['C_B_UAP'] = C_B * r['iUAP_CB']
        cmp_dct['C_B_Subst'] = C_B - cmp_dct['C_B_BAP'] - cmp_dct['C_B_UAP']
        
        cmp_dct['X_B_Subst'] = XC_B - C_B - X_Bio - X_Stor
        
        cmp_c = np.asarray([v for v in cmp_dct.values()])
        VSS = (cmp_c * cmps.i_mass * cmps.f_Vmass_Totmass * cmps.x * cmps.org).sum()
        TSS = VSS/iVSS_TSS
        X_Org_ISS = (cmp_c * cmps.i_mass * (1-cmps.f_Vmass_Totmass) * cmps.x * cmps.org).sum()

        del S_Org, XC_U_Inf, X_Bio, X_Stor, X_U_E, XC_B, C_B
        
        #************ inorganic components **************
        cmp_dct['X_PAO_PP_Hi'] = X_PAO_PP * r['iHi_XPAOPP']
        cmp_dct['X_PAO_PP_Lo'] = X_PAO_PP * (1 - r['iHi_XPAOPP'])
        
        ISS = TSS - VSS
        cmp_c = np.asarray([v for v in cmp_dct.values()])
        other_ig_iss = (cmp_c * cmps.i_mass * cmps.x * (1-cmps.org)).sum()
        cmp_dct['X_Ig_ISS'] = ISS - X_Org_ISS - other_ig_iss
        
        del ISS, VSS, TSS, X_Org_ISS, other_ig_iss, cmp_c
        
        # TODO: calibrate pH, SAlk, SCAT, SAN
        bad_vars = {k:v for k,v in cmp_dct.items() if v<0}
        if len(bad_vars) > 0:
            raise ValueError(f"The following state variable(s) was found negative: {bad_vars}.")
        
        del bad_vars
        
        #************ calibrate XB_subst, SF's N, P content *************
        if S_NH4 > 0 and cmp_dct['S_F'] > 0:
            SFi_N = _calib_SF_iN(cmps, cmp_dct, S_NH4/iSNH_STKN)
        XB_Substi_N = _calib_XBsub_iN(cmps, cmp_dct, TKN - S_NH4/iSNH_STKN)
        XB_Substi_P = _calib_XBsub_iP(cmps, cmp_dct, TP)
        
        cmps = Components(cmps)
        if S_NH4 > 0 and cmp_dct['S_F'] > 0: cmps.S_F.i_N = SFi_N
        cmps.X_B_Subst.i_N = XB_Substi_N
        cmps.X_B_Subst.i_P = XB_Substi_P
        cmps.compile()
        set_thermo(cmps)
        
        #************ convert concentrations to flow rates *************
        flow_tot /= vol_unit.conversion_factor(units[0])
        factor = conc_unit.conversion_factor(units[1])
            
        cmp_dct = {k:v/factor*flow_tot*1e-6 for k,v in cmp_dct.items()}       # [mg/L]*[L/hr]*1e-6[kg/mg] = [kg/hr]
        dwt = sum(cmp_dct.values())         # dry weight
        
        den = 1
        i = 0
        while True:
            den0 = den
            cmp_dct['H2O'] = flow_tot*den0 - dwt

            #---------- Yalin's comments ----------
            # This is to avoid BioSTEAM showing the warning that this ID exists
            # in the registry
            #!!! But maybe there's a better way to deal with this?
            # Also this appears in all of the four models, consider making this into a function
            temp = cls(ID=ID+str(i), phase=phase, T=T, P=P, units='kg/hr', price=price, 
                      thermo=thermo, pH=pH, SAlk=SAlk, **cmp_dct)
            den = flow_tot*den0/(temp.F_vol*1e3)            
            i += 1
            if abs(den-den0) <= 1e-3: break
            if i > 50: raise ValueError('Density calculation failed to converge within 50 iterations.')
        
        new.mass = temp.mass
        new.ratios = r
        
        return new


    @classmethod
    def codbased_inf_model(cls, ID, flow_tot=0., units = ('L/hr', 'mg/L'), 
                           phase='l', T=298.15, P=101325., price=0., thermo=None, 
                           pH=7., SAlk=10., ratios=None, 
                           COD=430., TKN=40., TP=10., iVSS_TSS=0.75, iSNH_STKN=0.9,
                           iSCOD_COD=0.25, iSBOD_SCOD=0.50, iBOD_COD=0.58,
                           S_NH4=25., S_NO2=0., S_NO3=0., S_PO4=8., 
                           S_Ca=140., S_Mg=50., S_K=28., S_CAT=3., S_AN=12., S_N2=18., 
                           C_B=40., S_CH3OH=0., S_Ac=0., S_Prop=0., 
                           X_OHO_PHA=0., X_GAO_PHA=0., X_PAO_PHA=0., 
                           X_GAO_Gly=0., X_PAO_Gly=0., X_U_OHO_E=0., X_U_PAO_E=0.,
                           X_OHO=0., X_AOO=0., X_NOO=0., X_AMO=0., X_PAO=0., 
                           X_PRO=0., X_ACO=0., X_HMO=0., X_MEOLO=0., X_FO=0.,
                           X_FePO4=0., X_AlPO4=0., X_FeOH=0., X_AlOH=0., 
                           X_MAP=0., X_HAP=0., X_HDP=0., X_PAO_PP=0., 
                           X_MgCO3=0., X_CaCO3=0., DO=0., S_H2=0., S_CH4=0.):
        
           
        if thermo: cmps = thermo.chemicals
        else: cmps = _load_components(None)
        
        cmp_dct = dict.fromkeys(cmps.IDs, 0.)

        new = cls(ID=ID, phase=phase, T=T, P=P, units='kg/hr', price=price, 
                  thermo=thermo, pH=pH, SAlk=SAlk)

        if ratios: new.ratios = ratios
        else: new.ratios = WasteStream._default_ratios
        r = new._ratios

        #************ user-defined inorganic states **************        
        cmp_dct['S_H2'] = S_H2
        cmp_dct['S_CH4'] = S_CH4
        cmp_dct['S_N2'] = S_N2
        cmp_dct['S_O2'] = DO
        cmp_dct['S_NH4'] = S_NH4
        cmp_dct['S_NO2'] = S_NO2
        cmp_dct['S_NO3'] = S_NO3
        cmp_dct['S_PO4'] = S_PO4
        cmp_dct['S_CO3'] = SAlk * 12 * conc_unit.conversion_factor(units[1])       # 1 meq/L SAlk ~ 1 mmol/L HCO3- ~ 12 mg C/L (12 mg C/mmol HCO3-)
        cmp_dct['S_Ca'] = S_Ca
        cmp_dct['S_Mg'] = S_Mg
        cmp_dct['S_K'] = S_K
        cmp_dct['X_MAP'] = X_MAP
        cmp_dct['X_HAP'] = X_HAP
        cmp_dct['X_HDP'] = X_HDP
        cmp_dct['X_FePO4'] = X_FePO4
        cmp_dct['X_AlPO4'] = X_AlPO4
        cmp_dct['X_FeOH'] = X_FeOH
        cmp_dct['X_AlOH'] = X_AlOH
        cmp_dct['X_MgCO3'] = X_MgCO3
        cmp_dct['X_CaCO3'] = X_CaCO3
        cmp_dct['S_CAT'] = S_CAT
        cmp_dct['S_AN'] = S_AN
        
        #************ organic components **************
        sCOD = COD * iSCOD_COD
        sBOD = sCOD * iSBOD_SCOD
        cmp_dct['S_CH3OH'] = S_CH3OH
        cmp_dct['S_Ac'] = S_Ac
        cmp_dct['S_Prop'] = S_Prop
        
        cmp_c = np.asarray([v for v in cmp_dct.values()])
        other_sBOD = (cmp_c * cmps.s * cmps.org * cmps.f_BOD5_COD).sum()
        cmp_dct['S_F'] = (sBOD - other_sBOD)/cmps.S_F.f_BOD5_COD
        
        cmp_c = np.asarray([v for v in cmp_dct.values()])        
        S_U = sCOD - (cmp_c * cmps.s * cmps.org).sum()
        cmp_dct['S_U_Inf'] = S_U * r['iSUInf_SU']
        cmp_dct['S_U_E'] = S_U - cmp_dct['SU_Inf']
        
        cmp_dct['C_B_BAP'] = C_B * r['iBAP_CB']
        cmp_dct['C_B_UAP'] = C_B * r['iUAP_CB']
        cmp_dct['C_B_Subst'] = C_B - cmp_dct['C_B_BAP'] - cmp_dct['C_B_UAP']

        XC_B = C_B/r['iCB_XCB']
        XC_U = COD - sCOD - XC_B
        cmp_dct['X_U_OHO_E'] = X_U_OHO_E
        cmp_dct['X_U_PAO_E'] = X_U_PAO_E
        
        XC_U_Inf = XC_U - X_U_OHO_E - X_U_PAO_E
        cmp_dct['C_U_Inf'] = XC_U_Inf * r['iCUInf_XCUInf']
        cmp_dct['X_U_Inf'] = XC_U_Inf * (1 - r['iCUInf_XCUInf'])
        
        cmp_dct['X_OHO'] = X_OHO
        cmp_dct['X_AOO'] = X_AOO
        cmp_dct['X_NOO'] = X_NOO
        cmp_dct['X_AMO'] = X_AMO
        cmp_dct['X_PAO'] = X_PAO
        cmp_dct['X_ACO'] = X_ACO
        cmp_dct['X_HMO'] = X_HMO
        cmp_dct['X_PRO'] = X_PRO
        cmp_dct['X_MEOLO'] = X_MEOLO
        cmp_dct['X_FO'] = X_FO
        cmp_dct['X_OHO_PHA'] = X_OHO_PHA
        cmp_dct['X_GAO_PHA'] = X_GAO_PHA
        cmp_dct['X_PAO_PHA'] = X_PAO_PHA
        cmp_dct['X_GAO_Gly'] = X_GAO_Gly
        cmp_dct['X_PAO_Gly'] = X_PAO_Gly
                
        cmp_c = np.asarray([v for v in cmp_dct.values()])
        cmp_dct['X_B_Subst'] = COD - (cmp_c * cmps.org * (cmps.s + cmps.c + cmps.x)).sum()
        
        cmp_c = np.asarray([v for v in cmp_dct.values()])
        VSS = (cmp_c * cmps.i_mass * cmps.f_Vmass_Totmass * cmps.x * cmps.org).sum()
        TSS = VSS/iVSS_TSS
        X_Org_ISS = (cmp_c * cmps.i_mass * (1-cmps.f_Vmass_Totmass) * cmps.x * cmps.org).sum()

        del sCOD, sBOD, other_sBOD, S_U, XC_U, XC_B, XC_U_Inf
        
        #************ inorganic components **************
        cmp_dct['X_PAO_PP_Hi'] = X_PAO_PP * r['iHi_XPAOPP']
        cmp_dct['X_PAO_PP_Lo'] = X_PAO_PP * (1 - r['iHi_XPAOPP'])
        
        ISS = TSS - VSS
        cmp_c = np.asarray([v for v in cmp_dct.values()])
        other_ig_iss = (cmp_c * cmps.i_mass * cmps.x * (1-cmps.org)).sum()
        cmp_dct['X_Ig_ISS'] = ISS - X_Org_ISS - other_ig_iss
        
        del ISS, VSS, TSS, X_Org_ISS, other_ig_iss, cmp_c
        
        # TODO: calibrate pH, SAlk, SCAT, SAN
        bad_vars = {k:v for k,v in cmp_dct.items() if v<0}
        if len(bad_vars) > 0:
            raise ValueError(f"The following state variable(s) was found negative: {bad_vars}.")
        
        del bad_vars
        
        #************ calibrate XB_subst, SF's N, P content *************
        if S_NH4 > 0 and cmp_dct['S_F'] > 0:
            SFi_N = _calib_SF_iN(cmps, cmp_dct, S_NH4/iSNH_STKN)
        XB_Substi_N = _calib_XBsub_iN(cmps, cmp_dct, TKN - S_NH4/iSNH_STKN)
        XB_Substi_P = _calib_XBsub_iP(cmps, cmp_dct, TP)
        
        BOD = COD * iBOD_COD
        sub_IDs = ('X_B_Subst', 'X_OHO_PHA', 'X_GAO_PHA', 'X_PAO_PHA', 'X_GAO_Gly', 'X_PAO_Gly')
        fbodtocod_sub = _calib_XBsub_fBODCOD(cmps, cmp_dct, sub_IDs, BOD)
        
        cmps = Components(cmps)
        if S_NH4 > 0 and cmp_dct['S_F'] > 0: cmps.S_F.i_N = SFi_N
        cmps.X_B_Subst.i_N = XB_Substi_N
        cmps.X_B_Subst.i_P = XB_Substi_P
        for i in sub_IDs: cmps[i].f_BOD5_COD = fbodtocod_sub
        cmps.compile()
        set_thermo(cmps)
        
        #************ convert concentrations to flow rates *************
        flow_tot /= vol_unit.conversion_factor(units[0])
        factor = conc_unit.conversion_factor(units[1])
            
        cmp_dct = {k:v/factor*flow_tot*1e-6 for k,v in cmp_dct.items()}       # [mg/L]*[L/hr]*1e-6[kg/mg] = [kg/hr]
        dwt = sum(cmp_dct.values())
        
        den = 1
        i = 0
        while True:
            den0 = den
            cmp_dct['H2O'] = flow_tot*den0 - dwt
            temp = cls(ID=ID+str(i), phase=phase, T=T, P=P, units='kg/hr', price=price, 
                      thermo=thermo, pH=pH, SAlk=SAlk, **cmp_dct)
            den = flow_tot*den0/(temp.F_vol*1e3)            
            i += 1
            if abs(den-den0) <= 1e-3: break
            if i > 50: raise ValueError('Density calculation failed to converge within 50 iterations.')

        new.mass = temp.mass
        new.ratios = r
        
        return new


    @classmethod
    def bodbased_inf_model(cls, ID, flow_tot=0., units = ('L/hr', 'mg/L'), 
                           phase='l', T=298.15, P=101325., price=0., thermo=None, 
                           pH=7., SAlk=10., ratios=None, 
                           BOD=250., TKN=40., TP=10., iVSS_TSS=0.75, iSNH_STKN=0.9,
                           iSBOD_BOD=0.25, iSBOD_SCOD=0.50, iBOD_COD=0.58,
                           S_N2=18., S_NH4=25., S_NO2=0., S_NO3=0., S_PO4=8., 
                           S_Ca=140., S_Mg=50., S_K=28., S_CAT=3., S_AN=12.,  
                           C_B=40., S_CH3OH=0., S_Ac=0., S_Prop=0., 
                           X_OHO_PHA=0., X_GAO_PHA=0., X_PAO_PHA=0., 
                           X_GAO_Gly=0., X_PAO_Gly=0., X_U_OHO_E=0., X_U_PAO_E=0.,
                           X_OHO=0., X_AOO=0., X_NOO=0., X_AMO=0., X_PAO=0., 
                           X_PRO=0., X_ACO=0., X_HMO=0., X_MEOLO=0., X_FO=0.,
                           X_FePO4=0., X_AlPO4=0., X_FeOH=0., X_AlOH=0., 
                           X_MAP=0., X_HAP=0., X_HDP=0., X_PAO_PP=0., 
                           X_MgCO3=0., X_CaCO3=0., DO=0., S_H2=0., S_CH4=0.):
        
           
        if thermo: cmps = thermo.chemicals
        else: cmps = _load_components(None)
        
        cmp_dct = dict.fromkeys(cmps.IDs, 0.)

        new = cls(ID=ID, phase=phase, T=T, P=P, units='kg/hr', price=price, 
                  thermo=thermo, pH=pH, SAlk=SAlk)

        if ratios: new.ratios = ratios
        else: new.ratios = WasteStream._default_ratios
        r = new._ratios

        #************ user-defined inorganic states **************        
        cmp_dct['S_H2'] = S_H2
        cmp_dct['S_CH4'] = S_CH4
        cmp_dct['S_N2'] = S_N2
        cmp_dct['S_O2'] = DO
        cmp_dct['S_NH4'] = S_NH4
        cmp_dct['S_NO2'] = S_NO2
        cmp_dct['S_NO3'] = S_NO3
        cmp_dct['S_PO4'] = S_PO4
        cmp_dct['S_CO3'] = SAlk * 12 * conc_unit.conversion_factor(units[1])       # 1 meq/L SAlk ~ 1 mmol/L HCO3- ~ 12 mg C/L (12 mg C/mmol HCO3-)
        cmp_dct['S_Ca'] = S_Ca
        cmp_dct['S_Mg'] = S_Mg
        cmp_dct['S_K'] = S_K
        cmp_dct['X_MAP'] = X_MAP
        cmp_dct['X_HAP'] = X_HAP
        cmp_dct['X_HDP'] = X_HDP
        cmp_dct['X_FePO4'] = X_FePO4
        cmp_dct['X_AlPO4'] = X_AlPO4
        cmp_dct['X_FeOH'] = X_FeOH
        cmp_dct['X_AlOH'] = X_AlOH
        cmp_dct['X_MgCO3'] = X_MgCO3
        cmp_dct['X_CaCO3'] = X_CaCO3
        cmp_dct['S_CAT'] = S_CAT
        cmp_dct['S_AN'] = S_AN
        
        #************ organic components **************
        COD = BOD / iBOD_COD
        sBOD = BOD * iSBOD_BOD
        cmp_dct['S_CH3OH'] = S_CH3OH
        cmp_dct['S_Ac'] = S_Ac
        cmp_dct['S_Prop'] = S_Prop
        
        cmp_c = np.asarray([v for v in cmp_dct.values()])
        other_sBOD = (cmp_c * cmps.s * cmps.org * cmps.f_BOD5_COD).sum()
        cmp_dct['S_F'] = (sBOD - other_sBOD)/cmps.S_F.f_BOD5_COD
        
        sCOD = sBOD / iSBOD_SCOD
        cmp_c = np.asarray([v for v in cmp_dct.values()])        
        S_U = sCOD - (cmp_c * cmps.s * cmps.org).sum()
        cmp_dct['S_U_Inf'] = S_U * r['iSUInf_SU']
        cmp_dct['S_U_E'] = S_U - cmp_dct['SU_Inf']
        
        cmp_dct['C_B_BAP'] = C_B * r['iBAP_CB']
        cmp_dct['C_B_UAP'] = C_B * r['iUAP_CB']
        cmp_dct['C_B_Subst'] = C_B - cmp_dct['C_B_BAP'] - cmp_dct['C_B_UAP']

        XC_B = C_B/r['iCB_XCB']
        XC_U = COD - sCOD - XC_B
        cmp_dct['X_U_OHO_E'] = X_U_OHO_E
        cmp_dct['X_U_PAO_E'] = X_U_PAO_E
        
        XC_U_Inf = XC_U - X_U_OHO_E - X_U_PAO_E
        cmp_dct['C_U_Inf'] = XC_U_Inf * r['iCUInf_XCUInf']
        cmp_dct['X_U_Inf'] = XC_U_Inf * (1 - r['iCUInf_XCUInf'])
        
        cmp_dct['X_OHO'] = X_OHO
        cmp_dct['X_AOO'] = X_AOO
        cmp_dct['X_NOO'] = X_NOO
        cmp_dct['X_AMO'] = X_AMO
        cmp_dct['X_PAO'] = X_PAO
        cmp_dct['X_ACO'] = X_ACO
        cmp_dct['X_HMO'] = X_HMO
        cmp_dct['X_PRO'] = X_PRO
        cmp_dct['X_MEOLO'] = X_MEOLO
        cmp_dct['X_FO'] = X_FO
        cmp_dct['X_OHO_PHA'] = X_OHO_PHA
        cmp_dct['X_GAO_PHA'] = X_GAO_PHA
        cmp_dct['X_PAO_PHA'] = X_PAO_PHA
        cmp_dct['X_GAO_Gly'] = X_GAO_Gly
        cmp_dct['X_PAO_Gly'] = X_PAO_Gly
                
        cmp_c = np.asarray([v for v in cmp_dct.values()])
        cmp_dct['X_B_Subst'] = COD - (cmp_c * cmps.org * (cmps.s + cmps.c + cmps.x)).sum()
        
        cmp_c = np.asarray([v for v in cmp_dct.values()])
        VSS = (cmp_c * cmps.i_mass * cmps.f_Vmass_Totmass * cmps.x * cmps.org).sum()
        TSS = VSS/iVSS_TSS
        X_Org_ISS = (cmp_c * cmps.i_mass * (1-cmps.f_Vmass_Totmass) * cmps.x * cmps.org).sum()

        del sCOD, sBOD, other_sBOD, S_U, XC_U, XC_B, XC_U_Inf
        
        #************ inorganic components **************
        cmp_dct['X_PAO_PP_Hi'] = X_PAO_PP * r['iHi_XPAOPP']
        cmp_dct['X_PAO_PP_Lo'] = X_PAO_PP * (1 - r['iHi_XPAOPP'])
        
        ISS = TSS - VSS
        cmp_c = np.asarray([v for v in cmp_dct.values()])
        other_ig_iss = (cmp_c * cmps.i_mass * cmps.x * (1-cmps.org)).sum()
        cmp_dct['X_Ig_ISS'] = ISS - X_Org_ISS - other_ig_iss
        
        del ISS, VSS, TSS, X_Org_ISS, other_ig_iss, cmp_c
        
        # TODO: calibrate pH, SAlk, SCAT, SAN
        bad_vars = {k:v for k,v in cmp_dct.items() if v<0}
        if len(bad_vars) > 0:
            raise ValueError(f"The following state variable(s) was found negative: {bad_vars}.")
        
        del bad_vars
        
        #************ calibrate XB_subst, SF's N, P content *************
        if S_NH4 > 0 and cmp_dct['S_F'] > 0:
            SFi_N = _calib_SF_iN(cmps, cmp_dct, S_NH4/iSNH_STKN)
        XB_Substi_N = _calib_XBsub_iN(cmps, cmp_dct, TKN - S_NH4/iSNH_STKN)
        XB_Substi_P = _calib_XBsub_iP(cmps, cmp_dct, TP)
        
        BOD = COD * iBOD_COD
        sub_IDs = ('X_B_Subst', 'X_OHO_PHA', 'X_GAO_PHA', 'X_PAO_PHA', 'X_GAO_Gly', 'X_PAO_Gly')
        fbodtocod_sub = _calib_XBsub_fBODCOD(cmps, cmp_dct, sub_IDs, BOD)
        
        cmps = Components(cmps)
        if S_NH4 > 0 and cmp_dct['S_F'] > 0: cmps.S_F.i_N = SFi_N
        cmps.X_B_Subst.i_N = XB_Substi_N
        cmps.X_B_Subst.i_P = XB_Substi_P
        for i in sub_IDs: cmps[i].f_BOD5_COD = fbodtocod_sub
        cmps.compile()
        set_thermo(cmps)

        #************ convert concentrations to flow rates *************
        flow_tot /= vol_unit.conversion_factor(units[0])
        factor = conc_unit.conversion_factor(units[1])
            
        cmp_dct = {k:v/factor*flow_tot*1e-6 for k,v in cmp_dct.items()}       # [mg/L]*[L/hr]*1e-6[kg/mg] = [kg/hr]
        dwt = sum(cmp_dct.values())
        
        den = 1
        i = 0
        while True:
            den0 = den
            cmp_dct['H2O'] = flow_tot*den0 - dwt
            temp = cls(ID=ID+str(i), phase=phase, T=T, P=P, units='kg/hr', price=price, 
                      thermo=thermo, pH=pH, SAlk=SAlk, **cmp_dct)
            den = flow_tot*den0/(temp.F_vol*1e3)            
            i += 1
            if abs(den-den0) <= 1e-3: break
            if i > 50: raise ValueError('Density calculation failed to converge within 50 iterations.')
        
        new.mass = temp.mass
        new.ratios = r
        
        return new


    @classmethod
    def sludge_inf_model(cls, ID, flow_tot=0., units = ('L/hr', 'mg/L'), 
                         phase='l', T=298.15, P=101325., price=0., thermo=None, 
                         pH=7., SAlk=10., ratios=None, 
                         TSS=1e4, TKN=750., TP=250., S_NH4=100., S_PO4=50., 
                         iVSS_TSS=0.65, iscCOD_COD=0.01, iSNH_STKN=0.9,
                         frXUInf_VSS=0.4, frXUE_VSS=0.3, frXOHO_VSS=0.2,  
                         frXAOO_VSS=0.01, frXNOO_VSS=0.01,frXPAO_VSS=0.01, 
                         frCB_scCOD=0.1, frSU_scCOD=0.8, 
                         S_Ca=140., S_Mg=50., S_K=28., S_CAT=3., S_AN=12., S_N2=18.,
                         frXACO_VSS=0., frXHMO_VSS=0., frXPRO_VSS=0., frXFO_VSS=0.,
                         frXMEOLO_VSS=0., frXAMO_VSS=0., 
                         frXOHO_PHA_VSS=0., frXGAO_PHA_VSS=0., frXPAO_PHA_VSS=0., 
                         frXGAO_Gly_VSS=0., frXPAO_Gly_VSS=0., 
                         frSCH3OH_scCOD=0., frSAc_scCOD=0., frSProp_scCOD=0.,
                         S_NO2=0., S_NO3=0., X_PAO_PP=0., X_FeOH=0., X_AlOH=0., 
                         X_FePO4=0., X_AlPO4=0., X_MAP=0., X_HAP=0., X_HDP=0., 
                         X_MgCO3=0., X_CaCO3=0., DO=0., S_H2=0., S_CH4=0.):
        
        
        if thermo: cmps = thermo.chemicals
        else: cmps = _load_components(None)
        
        cmp_dct = dict.fromkeys(cmps.IDs, 0.)

        new = cls(ID=ID, phase=phase, T=T, P=P, units='kg/hr', price=price, 
                  thermo=thermo, pH=pH, SAlk=SAlk)

        if ratios: new.ratios = ratios
        else: new.ratios = WasteStream._default_ratios
        r = new._ratios

        #************ user-defined inorganic states **************        
        cmp_dct['S_H2'] = S_H2
        cmp_dct['S_CH4'] = S_CH4
        cmp_dct['S_N2'] = S_N2
        cmp_dct['S_O2'] = DO
        cmp_dct['S_NH4'] = S_NH4
        cmp_dct['S_NO2'] = S_NO2
        cmp_dct['S_NO3'] = S_NO3
        cmp_dct['S_PO4'] = S_PO4
        cmp_dct['S_CO3'] = SAlk * 12 * conc_unit.conversion_factor(units[1])       # 1 meq/L SAlk ~ 1 mmol/L HCO3- ~ 12 mg C/L (12 mg C/mmol HCO3-)
        cmp_dct['S_Ca'] = S_Ca
        cmp_dct['S_Mg'] = S_Mg
        cmp_dct['S_K'] = S_K
        cmp_dct['X_MAP'] = X_MAP
        cmp_dct['X_HAP'] = X_HAP
        cmp_dct['X_HDP'] = X_HDP
        cmp_dct['X_FePO4'] = X_FePO4
        cmp_dct['X_AlPO4'] = X_AlPO4
        cmp_dct['X_FeOH'] = X_FeOH
        cmp_dct['X_AlOH'] = X_AlOH
        cmp_dct['X_MgCO3'] = X_MgCO3
        cmp_dct['X_CaCO3'] = X_CaCO3
        cmp_dct['S_CAT'] = S_CAT
        cmp_dct['S_AN'] = S_AN
        
        #************ particulate components **************
        VSS = TSS * iVSS_TSS
        cmp_dct['X_U_Inf'] = frXUInf_VSS * VSS
        
        if r['iXUOHOE_XUE']: frOHO = r['iXUOHOE_XUE']
        else: 
            try: frOHO = frXOHO_VSS/(frXOHO_VSS + frXPAO_VSS)
            except ZeroDivisionError: frOHO = 0.5
        cmp_dct['X_U_OHO_E'] = frXUE_VSS * VSS * frOHO
        cmp_dct['X_U_PAO_E'] = frXUE_VSS * VSS * (1-frOHO)
        
        cmp_dct['X_OHO'] = frXOHO_VSS * VSS
        cmp_dct['X_AOO'] = frXAOO_VSS * VSS
        cmp_dct['X_NOO'] = frXNOO_VSS * VSS
        cmp_dct['X_AMO'] = frXAMO_VSS * VSS
        cmp_dct['X_PAO'] = frXPAO_VSS * VSS
        cmp_dct['X_ACO'] = frXACO_VSS * VSS
        cmp_dct['X_HMO'] = frXHMO_VSS * VSS
        cmp_dct['X_PRO'] = frXPRO_VSS * VSS
        cmp_dct['X_MEOLO'] = frXMEOLO_VSS * VSS
        cmp_dct['X_FO'] = frXFO_VSS * VSS
        cmp_dct['X_OHO_PHA'] = frXOHO_PHA_VSS * VSS
        cmp_dct['X_GAO_PHA'] = frXGAO_PHA_VSS * VSS
        cmp_dct['X_PAO_PHA'] = frXPAO_PHA_VSS * VSS
        cmp_dct['X_GAO_Gly'] = frXGAO_Gly_VSS * VSS
        cmp_dct['X_PAO_Gly'] = frXPAO_Gly_VSS * VSS
        
        cmp_c = np.asarray([v for v in cmp_dct.values()])
        cmp_dct['X_B_Subst'] = VSS - (cmp_c * cmps.org * cmps.x).sum()
        
        # convert gVSS to gCOD
        for cmp in cmps:
            if cmp.organic and cmp.particle_size == 'Particulate':
                cmp_dct[cmp.ID] /= cmp.i_mass * cmp.f_Vmass_Totmass
        
        cmp_dct['X_PAO_PP_Hi'] = X_PAO_PP * r['iHi_XPAOPP']
        cmp_dct['X_PAO_PP_Lo'] = X_PAO_PP * (1 - r['iHi_XPAOPP'])
        
        cmp_c = np.asarray([v for v in cmp_dct.values()])
        ig_ISS = TSS - (cmp_c * cmps.i_mass * cmps.x * cmps.org).sum()
        other_ig_iss = (cmp_c * cmps.i_mass * cmps.x * (1-cmps.org)).sum()
        cmp_dct['X_Ig_ISS'] = ig_ISS - other_ig_iss
        
        del other_ig_iss, cmp_c
        
        #*********** soluble and colloidal components *************
        cmp_c = np.asarray([v for v in cmp_dct.values()])
        xCOD = (cmp_c * cmps.x * cmps.org).sum()
        scCOD = xCOD * iscCOD_COD / (1-iscCOD_COD)
        
        cmp_dct['S_CH3OH'] = frSCH3OH_scCOD * scCOD
        cmp_dct['S_Ac'] = frSAc_scCOD * scCOD
        cmp_dct['S_Prop'] = frSProp_scCOD * scCOD
        cmp_dct['S_U_Inf'] = frSU_scCOD * scCOD * r['iSUInf_SU']
        cmp_dct['S_U_E'] = frSU_scCOD * scCOD - cmp_dct['SU_Inf']
        
        C_B = frCB_scCOD * scCOD
        cmp_dct['C_B_BAP'] = C_B * r['iBAP_CB']
        cmp_dct['C_B_UAP'] = C_B * r['iUAP_CB']
        cmp_dct['C_B_Subst'] = C_B - cmp_dct['C_B_BAP'] - cmp_dct['C_B_UAP']        
        cmp_dct['C_U_Inf'] = cmp_dct['X_U_Inf'] * r['iCUInf_XCUInf'] / (1-r['iCUInf_XCUInf'])
        
        cmp_c = np.asarray([v for v in cmp_dct.values()])
        cmp_dct['S_F'] = scCOD - (cmp_c * (cmps.s + cmps.c) * cmps.org).sum()
                
        # TODO: calibrate pH, SAlk, SCAT, SAN
        bad_vars = {k:v for k,v in cmp_dct.items() if v<0}
        if len(bad_vars) > 0:
            raise ValueError(f"The following state variable(s) was found negative: {bad_vars}.")
        
        del bad_vars
        
        #************ calibrate XB_subst, SF's N, P content *************
        if S_NH4 > 0 and cmp_dct['S_F'] > 0:
            SFi_N = _calib_SF_iN(cmps, cmp_dct, S_NH4/iSNH_STKN)
        XB_Substi_N = _calib_XBsub_iN(cmps, cmp_dct, TKN - S_NH4/iSNH_STKN)
        XB_Substi_P = _calib_XBsub_iP(cmps, cmp_dct, TP)
        
        cmps = Components(cmps)
        if S_NH4 > 0 and cmp_dct['S_F'] > 0: cmps.S_F.i_N = SFi_N
        cmps.X_B_Subst.i_N = XB_Substi_N
        cmps.X_B_Subst.i_P = XB_Substi_P
        cmps.compile()
        set_thermo(cmps)

        #************ convert concentrations to flow rates *************
        flow_tot /= vol_unit.conversion_factor(units[0])
        factor = conc_unit.conversion_factor(units[1])
            
        cmp_dct = {k:v/factor*flow_tot*1e-6 for k,v in cmp_dct.items()}       # [mg/L]*[L/hr]*1e-6[kg/mg] = [kg/hr]
        dwt = sum(cmp_dct.values())
        
        den = 1
        i = 0
        while True:
            den0 = den
            cmp_dct['H2O'] = flow_tot*den0 - dwt
            temp = cls(ID=ID+str(i), phase=phase, T=T, P=P, units='kg/hr', price=price, 
                      thermo=thermo, pH=pH, SAlk=SAlk, **cmp_dct)
            den = flow_tot*den0/(temp.F_vol*1e3)            
            i += 1
            if abs(den-den0) <= 1e-3: break
            if i > 50: raise ValueError('Density calculation failed to converge within 50 iterations.')
        
        new.mass = temp.mass
        new.ratios = r

        return new


# %%

class MissingWasteStream(MissingSanStream):
    '''
    A subclass of :class:`MissingSanStream`, create a special object
    that acts as a dummy until replaced by an actual :class:`WasteStream`.
    '''

    # TODO: add others
    @property
    def pH(self):
        '''[float] pH, unitless.'''
        return 7.
    
    @property
    def SAlk(self):
        '''[float] Alkalinity.'''
        return 0.
    
    @property
    def COD(self):
        '''[float] Chemical oxygen demand.'''
        return 0.
    
    @property
    def BOD(self):
        '''[float] Biochemical oxygen demand, same as BOD5.'''
        return 0.
    
    @property    
    def BOD5(self):
        '''[float] 5-day biochemical oxygen demand, same as BOD.'''
        return self.BOD
    
    @property    
    def uBOD(self):
        '''[float] Ultimate biochemical oxygen demand.'''
        return 0.
    
    @property    
    def cnBOD(self):
        '''[float] Carbonaceous nitrogenous BOD.'''
        return 0.
    
    @property    
    def ThOD(self):
        '''[float] Theoretical oxygen demand.'''
        return 0.
    
    @property
    def TC(self):
        '''[float] Total carbon.'''
        return 0.
    
    @property
    def TOC(self):
        '''[float] Total organic carbon.'''
        return 0.
    
    @property
    def TN(self):
        '''[float] Total nitrogen.'''
        return 0.
    
    @property
    def TKN(self):
        '''[float] Total Kjeldahl nitrogen.'''
        return 0.
    
    @property
    def TP(self):
        '''[float] Total phosphorus.'''
        return 0.
    
    @property
    def TK(self):
        '''[float] Total potassium.'''
        return 0.
    
    @property
    def TMg(self):
        '''[float] Total magnesium.'''
        return 0.
    
    @property
    def TCa(self):
        '''[float] Total calcium.'''
        return 0.
    
    @property
    def dry_mass(self):
        '''[float] Total solids.'''
        return 0.
    
    #!!! Keep this up-to-date with WasteStream
    # @property
    # def charge(self):
    #     return 0.

    def __repr__(self):
        return '<MissingWasteStream>'
    
    def __str__(self):
        return 'missing waste stream'



