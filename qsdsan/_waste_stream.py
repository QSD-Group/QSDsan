#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Joy Cheung <joycheung1994@gmail.com>
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/master/LICENSE.txt
for license details.
'''

# %%
import numpy as np
import biosteam as bst
from thermosteam import Stream, MultiStream, utils
from . import Components
from ._units_of_measure import auom


__all__ = ('WasteStream',)


_defined_composite_vars = ('COD', 'BOD5', 'BOD', 'uBOD', 'NOD', 'ThOD', 'cnBOD',
                           'C', 'N', 'P', 'K', 'Mg', 'Ca', 'solids', 'charge')

_common_composite_vars = ('_COD', '_BOD', '_uBOD', '_TC', '_TOC', '_TN', 
                          '_TKN', '_TP', '_TK', '_TMg', '_TCa', 
                          '_dry_mass', '_charge', '_ThOD', '_cnBOD')

_ws_specific_slots = (*_common_composite_vars,
                      '_pH', '_SAlk', '_ratios', '_impact_item')

_specific_groups = {'SVFA': ('SAc', 'SProp'),
                    'XStor': ('XOHO_PHA', 'XGAO_PHA', 'XPAO_PHA', 
                              'XGAO_Gly', 'XPAO_Gly'),
                    'XANO': ('XAOO', 'XNOO'),
                    'XBio': ('XOHO', 'XAOO', 'XNOO', 'XAMO', 'XPAO', 
                             'XMEOLO', 'XACO', 'XHMO', 'XPRO', 'XFO'),
                    'SNOx': ('SNO2', 'SNO3'),
                    'XPAO_PP': ('XPAO_PP_Lo', 'XPAO_PP_Hi'),
                    'TKN': ()}

_default_ratios = {'iHi_XPAOPP': 0.5,
                   'iCB_XCB': 0.15,
                   'iBAP_CB': 0.,
                   'iUAP_CB': 0.,
                   'iCUInf_XCUInf': 0.,
                   'iSUInf_SU': 1.,
                   'iXUOHOE_XUE': None,}


vol_unit = auom('L/hr')
conc_unit = auom('mg/L')


# %%

# =============================================================================
# Define the WasteStream class
# =============================================================================

class WasteStream(Stream):
    '''
    A subclass of :class:`thermosteam.Stream` with additional attributes
    and methods for waste treatment.
    
    See Also
    --------
    `thermosteam.Stream <https://thermosteam.readthedocs.io/en/latest/Stream.html>`_
    
    '''
    
    # Child class will inherit parent class's slots
    __slots__ = _ws_specific_slots
    _default_ratios = _default_ratios
    ticket_name = 'ws'
    
    def __init__(self, ID='', flow=(), phase='l', T=298.15, P=101325.,
                 units='kg/hr', price=0., thermo=None, 
                 pH=7., SAlk=2.5, COD=None, BOD=None, uBOD=None,
                 TC=None, TOC=None, TN=None, TKN=None, TP=None, TK=None,
                 TMg=None, TCa=None, dry_mass=None, charge=None, ratios=None,
                 ThOD=None, cnBOD=None, impact_item=None, **chemical_flows):
        
        super().__init__(ID=ID, flow=flow, phase=phase, T=T, P=P,
                         units=units, price=price, thermo=thermo, **chemical_flows)
        self._init_ws(pH, SAlk, COD, BOD, uBOD, TC, TOC, TN, TKN,
                      TP, TK, TMg, TCa, ThOD, cnBOD, dry_mass, charge, ratios,
                      impact_item)

    def _init_ws(self, pH=7., SAlk=None, COD=None, BOD=None,
                  uBOD=None, TC=None, TOC=None, TN=None, TKN=None,
                  TP=None, TK=None, TMg=None, TCa=None, ThOD=None, cnBOD=None,
                  dry_mass=None, charge=None, ratios=None, impact_item=None):

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
        if impact_item:
            impact_item._linked_stream = self
        self._impact_item = impact_item

    
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
        _ws_info = '\n WasteStream-specific properties:'
        # Wastewater-related properties are not relevant for gas or solids
        if self.phase != 'l':
            _ws_info += ' None for non-liquid WasteStreams'
        elif self.F_mass == 0:
            _ws_info += ' None for empty WasteStreams'
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
    def components(self):
        return self.chemicals

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
                IDs = [ID for ID in IDs if ID not in ('SN2','SNO2','SNO3')]
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

    
    @property
    def impact_item(self):
        '''[StreamImpactItem] The :class:`StreamImpactItem` this waste stream is linked to.'''
        return self._impact_item
    @impact_item.setter
    def impact_item(self, i):
        self._impact_item = i
        i.linked_stream = self

    
    def _liq_sol_properties(self, prop, value):
        if self.phase != 'g':
            return getattr(self, '_'+prop) or value
        else:
            raise AttributeError(f'{self.phase} phase waste stream does not have {prop}.')
    
    #!!! Add some document
    @property
    def pH(self):
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
    

    def copy(self, ID=None):
        new = super().copy()
        new._init_ws()
        for slot in _ws_specific_slots:
            value = getattr(self, slot)
            if slot == '_impact_item' and value:
                value.copy(stream=new)
            else:
                setattr(new, slot, utils.copy_maybe(value))
        return new
    __copy__ = copy

    def copy_like(self, other):
        Stream.copy_like(self, other)
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
        Stream.mix_from(self, others)
        for slot in _ws_specific_slots:
            #!!! This need reviewing, might not be good to calculate some
            # attributes like pH
            try: tot = sum(float(getattr(i, slot))*i.F_vol for i in others)
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
                            SNH4=25., SNO2=0., SNO3=0., SPO4=8., 
                            SCa=140., SMg=50., SK=28., SCAT=3., SAN=12., SN2=18., 
                            frSUInf=0.05, frSF=0.2, frXCUInf=0.13, 
                            frSUE=0., frSCH3OH=0., frSAc=0., frSProp=0., 
                            frXOHO=0., frXAOO=0., frXNOO=0., frXAMO=0., frXPAO=0., 
                            frXPRO=0., frXACO=0., frXHMO=0., frXMEOLO=0., frXFO=0.,
                            frXOHO_PHA=0., frXGAO_PHA=0., frXPAO_PHA=0., 
                            frXGAO_Gly=0., frXPAO_Gly=0., frXU_OHO_E=0., frXU_PAO_E=0.,
                            XFePO4=0., XAlPO4=0., XFeOH=0., XAlOH=0., 
                            XMAP=0., XHAP=0., XHDP=0., XPAO_PP=0., 
                            XMgCO3=0., XCaCO3=0., DO=0., SH2=0., SCH4=0.):
        
           
        cmps = Components.load_default(default_compile=True)
        bst.settings.set_thermo(cmps)
        
        cmp_dct = dict.fromkeys(cmps.IDs, 0.)

        new = cls(ID=ID, phase=phase, T=T, P=P, units='kg/hr', price=price, 
                  thermo=thermo, pH=pH, SAlk=SAlk)

        if ratios: new.ratios = ratios
        else: new.ratios = WasteStream._default_ratios
        r = new._ratios

        #************ user-defined states **************        
        cmp_dct['SH2'] = SH2
        cmp_dct['SCH4'] = SCH4
        cmp_dct['SN2'] = SN2
        cmp_dct['SO2'] = DO
        cmp_dct['SNH4'] = SNH4
        cmp_dct['SNO2'] = SNO2
        cmp_dct['SNO3'] = SNO3
        cmp_dct['SPO4'] = SPO4
        cmp_dct['SCO3'] = SAlk * 12 * conc_unit.conversion_factor(units[1])             # 1 meq/L SAlk ~ 1 mmol/L HCO3- ~ 12 mg C/L (12 mg C/mmol HCO3-)
        cmp_dct['SCa'] = SCa
        cmp_dct['SMg'] = SMg
        cmp_dct['SK'] = SK
        cmp_dct['XMAP'] = XMAP
        cmp_dct['XHAP'] = XHAP
        cmp_dct['XHDP'] = XHDP
        cmp_dct['XFePO4'] = XFePO4
        cmp_dct['XAlPO4'] = XAlPO4
        cmp_dct['XFeOH'] = XFeOH
        cmp_dct['XAlOH'] = XAlOH
        cmp_dct['XMgCO3'] = XMgCO3
        cmp_dct['XCaCO3'] = XCaCO3
        cmp_dct['SCAT'] = SCAT
        cmp_dct['SAN'] = SAN
        
        #************ organic components **************
        cmp_dct['SCH3OH'] = COD * frSCH3OH
        cmp_dct['SAc'] = COD * frSAc
        cmp_dct['SProp'] = COD * frSProp
        cmp_dct['SF'] = COD * frSF
        cmp_dct['SU_Inf'] = COD * frSUInf
        cmp_dct['SU_E'] = COD * frSUE
        
        SOrg = sum([v for k,v in cmp_dct.items() if k in ('SCH3OH','SAc','SProp','SF','SU_Inf','SU_E')])
        
        XCU_Inf = COD * frXCUInf
        cmp_dct['CU_Inf'] = XCU_Inf * r['iCUInf_XCUInf']
        cmp_dct['XU_Inf'] = XCU_Inf * (1 - r['iCUInf_XCUInf'])
        
        cmp_dct['XOHO'] = COD * frXOHO
        cmp_dct['XAOO'] = COD * frXAOO
        cmp_dct['XNOO'] = COD * frXNOO
        cmp_dct['XAMO'] = COD * frXAMO
        cmp_dct['XPAO'] = COD * frXPAO
        cmp_dct['XACO'] = COD * frXACO
        cmp_dct['XHMO'] = COD * frXHMO
        cmp_dct['XPRO'] = COD * frXPRO
        cmp_dct['XMEOLO'] = COD * frXMEOLO
        cmp_dct['XFO'] = COD * frXFO
        
        XBio = sum([v for k,v in cmp_dct.items() if k.startswith('x') and k.endswith('O')])
        
        cmp_dct['XOHO_PHA'] = COD * frXOHO_PHA
        cmp_dct['XGAO_PHA'] = COD * frXGAO_PHA
        cmp_dct['XPAO_PHA'] = COD * frXPAO_PHA
        cmp_dct['XGAO_Gly'] = COD * frXGAO_Gly
        cmp_dct['XPAO_Gly'] = COD * frXPAO_Gly
        
        XStor = sum([v for k,v in cmp_dct.items() if k.endswith(('PHA','Gly'))])
        
        cmp_dct['XU_OHO_E'] = COD * frXU_OHO_E
        cmp_dct['XU_PAO_E'] = COD * frXU_PAO_E
        
        XU_E = cmp_dct['XU_OHO_E'] + cmp_dct['XU_PAO_E']
        
        XCB = COD - SOrg - XCU_Inf - XU_E
        CB = XCB * r['iCB_XCB']
        cmp_dct['CB_BAP'] = CB * r['iBAP_CB']
        cmp_dct['CB_UAP'] = CB * r['iUAP_CB']
        cmp_dct['CB_Subst'] = CB - cmp_dct['CB_BAP'] - cmp_dct['CB_UAP']
        
        cmp_dct['XB_Subst'] = XCB - CB - XBio - XStor
        
        cmp_c = np.asarray([v for v in cmp_dct.values()])
        VSS = (cmp_c * cmps.i_mass * cmps.f_Vmass_Totmass * cmps.x * cmps.org).sum()
        TSS = VSS/iVSS_TSS
        XOrg_ISS = (cmp_c * cmps.i_mass * (1-cmps.f_Vmass_Totmass) * cmps.x * cmps.org).sum()

        del SOrg, XCU_Inf, XBio, XStor, XU_E, XCB, CB
        
        #************ inorganic components **************
        cmp_dct['XPAO_PP_Hi'] = XPAO_PP * r['iHi_XPAOPP']
        cmp_dct['XPAO_PP_Lo'] = XPAO_PP * (1 - r['iHi_XPAOPP'])
        
        ISS = TSS - VSS
        cmp_c = np.asarray([v for v in cmp_dct.values()])
        other_ig_iss = (cmp_c * cmps.i_mass * cmps.x * (1-cmps.org)).sum()
        cmp_dct['XIg_ISS'] = ISS - XOrg_ISS - other_ig_iss
        
        del ISS, VSS, TSS, XOrg_ISS, other_ig_iss, cmp_c
        
        # TODO: calibrate pH, SAlk, SCAT, SAN
        bad_vars = {k:v for k,v in cmp_dct.items() if v<0}
        if len(bad_vars) > 0:
            raise ValueError(f"The following state variable(s) was found negative: {bad_vars}.")
        
        del bad_vars
        
        #************ calibrate XB_subst, SF's N, P content *************
        if SNH4 > 0 and cmp_dct['SF'] > 0:
            SFi_N = _calib_SF_iN(cmps, cmp_dct, SNH4/iSNH_STKN)
        XB_Substi_N = _calib_XBsub_iN(cmps, cmp_dct, TKN - SNH4/iSNH_STKN)
        XB_Substi_P = _calib_XBsub_iP(cmps, cmp_dct, TP)
        
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
            new = cls(ID=ID, phase=phase, T=T, P=P, units='kg/hr', price=price, 
                      thermo=thermo, pH=pH, SAlk=SAlk, **cmp_dct)
            den = flow_tot*den0/(new.F_vol*1e3)            
            i += 1
            if abs(den-den0) <= 1e-3: break
            if i > 50: raise ValueError('Density calculation failed to converge within 50 iterations.')
        
        new.components.SF.i_N = SFi_N
        new.components.XB_Subst.i_N = XB_Substi_N
        new.components.XB_Subst.i_P = XB_Substi_P
        
        return new


    @classmethod
    def codbased_inf_model(cls, ID, flow_tot=0., units = ('L/hr', 'mg/L'), 
                           phase='l', T=298.15, P=101325., price=0., thermo=None, 
                           pH=7., SAlk=10., ratios=None, 
                           COD=430., TKN=40., TP=10., iVSS_TSS=0.75, iSNH_STKN=0.9,
                           iSCOD_COD=0.25, iSBOD_SCOD=0.50, iBOD_COD=0.58,
                           SNH4=25., SNO2=0., SNO3=0., SPO4=8., 
                           SCa=140., SMg=50., SK=28., SCAT=3., SAN=12., SN2=18., 
                           CB=40., SCH3OH=0., SAc=0., SProp=0., 
                           XOHO_PHA=0., XGAO_PHA=0., XPAO_PHA=0., 
                           XGAO_Gly=0., XPAO_Gly=0., XU_OHO_E=0., XU_PAO_E=0.,
                           XOHO=0., XAOO=0., XNOO=0., XAMO=0., XPAO=0., 
                           XPRO=0., XACO=0., XHMO=0., XMEOLO=0., XFO=0.,
                           XFePO4=0., XAlPO4=0., XFeOH=0., XAlOH=0., 
                           XMAP=0., XHAP=0., XHDP=0., XPAO_PP=0., 
                           XMgCO3=0., XCaCO3=0., DO=0., SH2=0., SCH4=0.):
        
           
        cmps = Components.load_default(default_compile=True)
        bst.settings.set_thermo(cmps)
        
        cmp_dct = dict.fromkeys(cmps.IDs, 0.)

        new = cls(ID=ID, phase=phase, T=T, P=P, units='kg/hr', price=price, 
                  thermo=thermo, pH=pH, SAlk=SAlk)

        if ratios: new.ratios = ratios
        else: new.ratios = WasteStream._default_ratios
        r = new._ratios

        #************ user-defined inorganic states **************        
        cmp_dct['SH2'] = SH2
        cmp_dct['SCH4'] = SCH4
        cmp_dct['SN2'] = SN2
        cmp_dct['SO2'] = DO
        cmp_dct['SNH4'] = SNH4
        cmp_dct['SNO2'] = SNO2
        cmp_dct['SNO3'] = SNO3
        cmp_dct['SPO4'] = SPO4
        cmp_dct['SCO3'] = SAlk * 12 * conc_unit.conversion_factor(units[1])       # 1 meq/L SAlk ~ 1 mmol/L HCO3- ~ 12 mg C/L (12 mg C/mmol HCO3-)
        cmp_dct['SCa'] = SCa
        cmp_dct['SMg'] = SMg
        cmp_dct['SK'] = SK
        cmp_dct['XMAP'] = XMAP
        cmp_dct['XHAP'] = XHAP
        cmp_dct['XHDP'] = XHDP
        cmp_dct['XFePO4'] = XFePO4
        cmp_dct['XAlPO4'] = XAlPO4
        cmp_dct['XFeOH'] = XFeOH
        cmp_dct['XAlOH'] = XAlOH
        cmp_dct['XMgCO3'] = XMgCO3
        cmp_dct['XCaCO3'] = XCaCO3
        cmp_dct['SCAT'] = SCAT
        cmp_dct['SAN'] = SAN
        
        #************ organic components **************
        sCOD = COD * iSCOD_COD
        sBOD = sCOD * iSBOD_SCOD
        cmp_dct['SCH3OH'] = SCH3OH
        cmp_dct['SAc'] = SAc
        cmp_dct['SProp'] = SProp
        
        cmp_c = np.asarray([v for v in cmp_dct.values()])
        other_sBOD = (cmp_c * cmps.s * cmps.org * cmps.f_BOD5_COD).sum()
        cmp_dct['SF'] = (sBOD - other_sBOD)/cmps.SF.f_BOD5_COD
        
        cmp_c = np.asarray([v for v in cmp_dct.values()])        
        SU = sCOD - (cmp_c * cmps.s * cmps.org).sum()
        cmp_dct['SU_Inf'] = SU * r['iSUInf_SU']
        cmp_dct['SU_E'] = SU - cmp_dct['SU_Inf']
        
        cmp_dct['CB_BAP'] = CB * r['iBAP_CB']
        cmp_dct['CB_UAP'] = CB * r['iUAP_CB']
        cmp_dct['CB_Subst'] = CB - cmp_dct['CB_BAP'] - cmp_dct['CB_UAP']

        XCB = CB/r['iCB_XCB']
        XCU = COD - sCOD - XCB
        cmp_dct['XU_OHO_E'] = XU_OHO_E
        cmp_dct['XU_PAO_E'] = XU_PAO_E
        
        XCU_Inf = XCU - XU_OHO_E - XU_PAO_E
        cmp_dct['CU_Inf'] = XCU_Inf * r['iCUInf_XCUInf']
        cmp_dct['XU_Inf'] = XCU_Inf * (1 - r['iCUInf_XCUInf'])
        
        cmp_dct['XOHO'] = XOHO
        cmp_dct['XAOO'] = XAOO
        cmp_dct['XNOO'] = XNOO
        cmp_dct['XAMO'] = XAMO
        cmp_dct['XPAO'] = XPAO
        cmp_dct['XACO'] = XACO
        cmp_dct['XHMO'] = XHMO
        cmp_dct['XPRO'] = XPRO
        cmp_dct['XMEOLO'] = XMEOLO
        cmp_dct['XFO'] = XFO
        cmp_dct['XOHO_PHA'] = XOHO_PHA
        cmp_dct['XGAO_PHA'] = XGAO_PHA
        cmp_dct['XPAO_PHA'] = XPAO_PHA
        cmp_dct['XGAO_Gly'] = XGAO_Gly
        cmp_dct['XPAO_Gly'] = XPAO_Gly
                
        cmp_c = np.asarray([v for v in cmp_dct.values()])
        cmp_dct['XB_Subst'] = COD - (cmp_c * cmps.org * (cmps.s + cmps.c + cmps.x)).sum()
        
        cmp_c = np.asarray([v for v in cmp_dct.values()])
        VSS = (cmp_c * cmps.i_mass * cmps.f_Vmass_Totmass * cmps.x * cmps.org).sum()
        TSS = VSS/iVSS_TSS
        XOrg_ISS = (cmp_c * cmps.i_mass * (1-cmps.f_Vmass_Totmass) * cmps.x * cmps.org).sum()

        del sCOD, sBOD, other_sBOD, SU, XCU, XCB, XCU_Inf
        
        #************ inorganic components **************
        cmp_dct['XPAO_PP_Hi'] = XPAO_PP * r['iHi_XPAOPP']
        cmp_dct['XPAO_PP_Lo'] = XPAO_PP * (1 - r['iHi_XPAOPP'])
        
        ISS = TSS - VSS
        cmp_c = np.asarray([v for v in cmp_dct.values()])
        other_ig_iss = (cmp_c * cmps.i_mass * cmps.x * (1-cmps.org)).sum()
        cmp_dct['XIg_ISS'] = ISS - XOrg_ISS - other_ig_iss
        
        del ISS, VSS, TSS, XOrg_ISS, other_ig_iss, cmp_c
        
        # TODO: calibrate pH, SAlk, SCAT, SAN
        bad_vars = {k:v for k,v in cmp_dct.items() if v<0}
        if len(bad_vars) > 0:
            raise ValueError(f"The following state variable(s) was found negative: {bad_vars}.")
        
        del bad_vars
        
        #************ calibrate XB_subst, SF's N, P content *************
        if SNH4 > 0 and cmp_dct['SF'] > 0:
            SFi_N = _calib_SF_iN(cmps, cmp_dct, SNH4/iSNH_STKN)
        XB_Substi_N = _calib_XBsub_iN(cmps, cmp_dct, TKN - SNH4/iSNH_STKN)
        XB_Substi_P = _calib_XBsub_iP(cmps, cmp_dct, TP)
        
        BOD = COD * iBOD_COD
        sub_IDs = ('XB_Subst', 'XOHO_PHA', 'XGAO_PHA', 'XPAO_PHA', 'XGAO_Gly', 'XPAO_Gly')
        fbodtocod_sub = _calib_XBsub_fBODCOD(cmps, cmp_dct, sub_IDs, BOD)
        
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
            new = cls(ID=ID, phase=phase, T=T, P=P, units='kg/hr', price=price, 
                      thermo=thermo, pH=pH, SAlk=SAlk, **cmp_dct)
            den = flow_tot*den0/(new.F_vol*1e3)            
            i += 1
            if abs(den-den0) <= 1e-3: break
            if i > 50: raise ValueError('Density calculation failed to converge within 50 iterations.')

        new.components.SF.i_N = SFi_N
        new.components.XB_Subst.i_N = XB_Substi_N
        new.components.XB_Subst.i_P = XB_Substi_P
        for i in sub_IDs: new.components[i].f_BOD5_COD = fbodtocod_sub
        
        return new


    @classmethod
    def bodbased_inf_model(cls, ID, flow_tot=0., units = ('L/hr', 'mg/L'), 
                           phase='l', T=298.15, P=101325., price=0., thermo=None, 
                           pH=7., SAlk=10., ratios=None, 
                           BOD=250., TKN=40., TP=10., iVSS_TSS=0.75, iSNH_STKN=0.9,
                           iSBOD_BOD=0.25, iSBOD_SCOD=0.50, iBOD_COD=0.58,
                           SN2=18., SNH4=25., SNO2=0., SNO3=0., SPO4=8., 
                           SCa=140., SMg=50., SK=28., SCAT=3., SAN=12.,  
                           CB=40., SCH3OH=0., SAc=0., SProp=0., 
                           XOHO_PHA=0., XGAO_PHA=0., XPAO_PHA=0., 
                           XGAO_Gly=0., XPAO_Gly=0., XU_OHO_E=0., XU_PAO_E=0.,
                           XOHO=0., XAOO=0., XNOO=0., XAMO=0., XPAO=0., 
                           XPRO=0., XACO=0., XHMO=0., XMEOLO=0., XFO=0.,
                           XFePO4=0., XAlPO4=0., XFeOH=0., XAlOH=0., 
                           XMAP=0., XHAP=0., XHDP=0., XPAO_PP=0., 
                           XMgCO3=0., XCaCO3=0., DO=0., SH2=0., SCH4=0.):
        
           
        cmps = Components.load_default(default_compile=True)
        bst.settings.set_thermo(cmps)
        
        cmp_dct = dict.fromkeys(cmps.IDs, 0.)

        new = cls(ID=ID, phase=phase, T=T, P=P, units='kg/hr', price=price, 
                  thermo=thermo, pH=pH, SAlk=SAlk)

        if ratios: new.ratios = ratios
        else: new.ratios = WasteStream._default_ratios
        r = new._ratios

        #************ user-defined inorganic states **************        
        cmp_dct['SH2'] = SH2
        cmp_dct['SCH4'] = SCH4
        cmp_dct['SN2'] = SN2
        cmp_dct['SO2'] = DO
        cmp_dct['SNH4'] = SNH4
        cmp_dct['SNO2'] = SNO2
        cmp_dct['SNO3'] = SNO3
        cmp_dct['SPO4'] = SPO4
        cmp_dct['SCO3'] = SAlk * 12 * conc_unit.conversion_factor(units[1])       # 1 meq/L SAlk ~ 1 mmol/L HCO3- ~ 12 mg C/L (12 mg C/mmol HCO3-)
        cmp_dct['SCa'] = SCa
        cmp_dct['SMg'] = SMg
        cmp_dct['SK'] = SK
        cmp_dct['XMAP'] = XMAP
        cmp_dct['XHAP'] = XHAP
        cmp_dct['XHDP'] = XHDP
        cmp_dct['XFePO4'] = XFePO4
        cmp_dct['XAlPO4'] = XAlPO4
        cmp_dct['XFeOH'] = XFeOH
        cmp_dct['XAlOH'] = XAlOH
        cmp_dct['XMgCO3'] = XMgCO3
        cmp_dct['XCaCO3'] = XCaCO3
        cmp_dct['SCAT'] = SCAT
        cmp_dct['SAN'] = SAN
        
        #************ organic components **************
        COD = BOD / iBOD_COD
        sBOD = BOD * iSBOD_BOD
        cmp_dct['SCH3OH'] = SCH3OH
        cmp_dct['SAc'] = SAc
        cmp_dct['SProp'] = SProp
        
        cmp_c = np.asarray([v for v in cmp_dct.values()])
        other_sBOD = (cmp_c * cmps.s * cmps.org * cmps.f_BOD5_COD).sum()
        cmp_dct['SF'] = (sBOD - other_sBOD)/cmps.SF.f_BOD5_COD
        
        sCOD = sBOD / iSBOD_SCOD
        cmp_c = np.asarray([v for v in cmp_dct.values()])        
        SU = sCOD - (cmp_c * cmps.s * cmps.org).sum()
        cmp_dct['SU_Inf'] = SU * r['iSUInf_SU']
        cmp_dct['SU_E'] = SU - cmp_dct['SU_Inf']
        
        cmp_dct['CB_BAP'] = CB * r['iBAP_CB']
        cmp_dct['CB_UAP'] = CB * r['iUAP_CB']
        cmp_dct['CB_Subst'] = CB - cmp_dct['CB_BAP'] - cmp_dct['CB_UAP']

        XCB = CB/r['iCB_XCB']
        XCU = COD - sCOD - XCB
        cmp_dct['XU_OHO_E'] = XU_OHO_E
        cmp_dct['XU_PAO_E'] = XU_PAO_E
        
        XCU_Inf = XCU - XU_OHO_E - XU_PAO_E
        cmp_dct['CU_Inf'] = XCU_Inf * r['iCUInf_XCUInf']
        cmp_dct['XU_Inf'] = XCU_Inf * (1 - r['iCUInf_XCUInf'])
        
        cmp_dct['XOHO'] = XOHO
        cmp_dct['XAOO'] = XAOO
        cmp_dct['XNOO'] = XNOO
        cmp_dct['XAMO'] = XAMO
        cmp_dct['XPAO'] = XPAO
        cmp_dct['XACO'] = XACO
        cmp_dct['XHMO'] = XHMO
        cmp_dct['XPRO'] = XPRO
        cmp_dct['XMEOLO'] = XMEOLO
        cmp_dct['XFO'] = XFO
        cmp_dct['XOHO_PHA'] = XOHO_PHA
        cmp_dct['XGAO_PHA'] = XGAO_PHA
        cmp_dct['XPAO_PHA'] = XPAO_PHA
        cmp_dct['XGAO_Gly'] = XGAO_Gly
        cmp_dct['XPAO_Gly'] = XPAO_Gly
                
        cmp_c = np.asarray([v for v in cmp_dct.values()])
        cmp_dct['XB_Subst'] = COD - (cmp_c * cmps.org * (cmps.s + cmps.c + cmps.x)).sum()
        
        cmp_c = np.asarray([v for v in cmp_dct.values()])
        VSS = (cmp_c * cmps.i_mass * cmps.f_Vmass_Totmass * cmps.x * cmps.org).sum()
        TSS = VSS/iVSS_TSS
        XOrg_ISS = (cmp_c * cmps.i_mass * (1-cmps.f_Vmass_Totmass) * cmps.x * cmps.org).sum()

        del sCOD, sBOD, other_sBOD, SU, XCU, XCB, XCU_Inf
        
        #************ inorganic components **************
        cmp_dct['XPAO_PP_Hi'] = XPAO_PP * r['iHi_XPAOPP']
        cmp_dct['XPAO_PP_Lo'] = XPAO_PP * (1 - r['iHi_XPAOPP'])
        
        ISS = TSS - VSS
        cmp_c = np.asarray([v for v in cmp_dct.values()])
        other_ig_iss = (cmp_c * cmps.i_mass * cmps.x * (1-cmps.org)).sum()
        cmp_dct['XIg_ISS'] = ISS - XOrg_ISS - other_ig_iss
        
        del ISS, VSS, TSS, XOrg_ISS, other_ig_iss, cmp_c
        
        # TODO: calibrate pH, SAlk, SCAT, SAN
        bad_vars = {k:v for k,v in cmp_dct.items() if v<0}
        if len(bad_vars) > 0:
            raise ValueError(f"The following state variable(s) was found negative: {bad_vars}.")
        
        del bad_vars
        
        #************ calibrate XB_subst, SF's N, P content *************
        if SNH4 > 0 and cmp_dct['SF'] > 0:
            SFi_N = _calib_SF_iN(cmps, cmp_dct, SNH4/iSNH_STKN)
        XB_Substi_N = _calib_XBsub_iN(cmps, cmp_dct, TKN - SNH4/iSNH_STKN)
        XB_Substi_P = _calib_XBsub_iP(cmps, cmp_dct, TP)
        
        BOD = COD * iBOD_COD
        sub_IDs = ('XB_Subst', 'XOHO_PHA', 'XGAO_PHA', 'XPAO_PHA', 'XGAO_Gly', 'XPAO_Gly')
        fbodtocod_sub = _calib_XBsub_fBODCOD(cmps, cmp_dct, sub_IDs, BOD)
        
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
            new = cls(ID=ID, phase=phase, T=T, P=P, units='kg/hr', price=price, 
                      thermo=thermo, pH=pH, SAlk=SAlk, **cmp_dct)
            den = flow_tot*den0/(new.F_vol*1e3)            
            i += 1
            if abs(den-den0) <= 1e-3: break
            if i > 50: raise ValueError('Density calculation failed to converge within 50 iterations.')
        
        new.components.SF.i_N = SFi_N
        new.components.XB_Subst.i_N = XB_Substi_N
        new.components.XB_Subst.i_P = XB_Substi_P        
        for i in sub_IDs: new.components[i].f_BOD5_COD = fbodtocod_sub

        return new


    @classmethod
    def sludge_inf_model(cls, ID, flow_tot=0., units = ('L/hr', 'mg/L'), 
                         phase='l', T=298.15, P=101325., price=0., thermo=None, 
                         pH=7., SAlk=10., ratios=None, 
                         TSS=1e4, TKN=750., TP=250., SNH4=100., SPO4=50., 
                         iVSS_TSS=0.65, iscCOD_COD=0.01, iSNH_STKN=0.9,
                         frXUInf_VSS=0.4, frXUE_VSS=0.3, frXOHO_VSS=0.2,  
                         frXAOO_VSS=0.01, frXNOO_VSS=0.01,frXPAO_VSS=0.01, 
                         frCB_scCOD=0.1, frSU_scCOD=0.8, 
                         SCa=140., SMg=50., SK=28., SCAT=3., SAN=12., SN2=18.,
                         frXACO_VSS=0., frXHMO_VSS=0., frXPRO_VSS=0., frXFO_VSS=0.,
                         frXMEOLO_VSS=0., frXAMO_VSS=0., 
                         frXOHO_PHA_VSS=0., frXGAO_PHA_VSS=0., frXPAO_PHA_VSS=0., 
                         frXGAO_Gly_VSS=0., frXPAO_Gly_VSS=0., 
                         frSCH3OH_scCOD=0., frSAc_scCOD=0., frSProp_scCOD=0.,
                         SNO2=0., SNO3=0., XPAO_PP=0., XFeOH=0., XAlOH=0., 
                         XFePO4=0., XAlPO4=0., XMAP=0., XHAP=0., XHDP=0., 
                         XMgCO3=0., XCaCO3=0., DO=0., SH2=0., SCH4=0.):
        
        
        cmps = Components.load_default(default_compile=True)
        bst.settings.set_thermo(cmps)
        
        cmp_dct = dict.fromkeys(cmps.IDs, 0.)

        new = cls(ID=ID, phase=phase, T=T, P=P, units='kg/hr', price=price, 
                  thermo=thermo, pH=pH, SAlk=SAlk)

        if ratios: new.ratios = ratios
        else: new.ratios = WasteStream._default_ratios
        r = new._ratios

        #************ user-defined inorganic states **************        
        cmp_dct['SH2'] = SH2
        cmp_dct['SCH4'] = SCH4
        cmp_dct['SN2'] = SN2
        cmp_dct['SO2'] = DO
        cmp_dct['SNH4'] = SNH4
        cmp_dct['SNO2'] = SNO2
        cmp_dct['SNO3'] = SNO3
        cmp_dct['SPO4'] = SPO4
        cmp_dct['SCO3'] = SAlk * 12 * conc_unit.conversion_factor(units[1])       # 1 meq/L SAlk ~ 1 mmol/L HCO3- ~ 12 mg C/L (12 mg C/mmol HCO3-)
        cmp_dct['SCa'] = SCa
        cmp_dct['SMg'] = SMg
        cmp_dct['SK'] = SK
        cmp_dct['XMAP'] = XMAP
        cmp_dct['XHAP'] = XHAP
        cmp_dct['XHDP'] = XHDP
        cmp_dct['XFePO4'] = XFePO4
        cmp_dct['XAlPO4'] = XAlPO4
        cmp_dct['XFeOH'] = XFeOH
        cmp_dct['XAlOH'] = XAlOH
        cmp_dct['XMgCO3'] = XMgCO3
        cmp_dct['XCaCO3'] = XCaCO3
        cmp_dct['SCAT'] = SCAT
        cmp_dct['SAN'] = SAN
        
        #************ particulate components **************
        VSS = TSS * iVSS_TSS
        cmp_dct['XU_Inf'] = frXUInf_VSS * VSS
        
        if r['iXUOHOE_XUE']: frOHO = r['iXUOHOE_XUE']
        else: 
            try: frOHO = frXOHO_VSS/(frXOHO_VSS + frXPAO_VSS)
            except ZeroDivisionError: frOHO = 0.5
        cmp_dct['XU_OHO_E'] = frXUE_VSS * VSS * frOHO
        cmp_dct['XU_PAO_E'] = frXUE_VSS * VSS * (1-frOHO)
        
        cmp_dct['XOHO'] = frXOHO_VSS * VSS
        cmp_dct['XAOO'] = frXAOO_VSS * VSS
        cmp_dct['XNOO'] = frXNOO_VSS * VSS
        cmp_dct['XAMO'] = frXAMO_VSS * VSS
        cmp_dct['XPAO'] = frXPAO_VSS * VSS
        cmp_dct['XACO'] = frXACO_VSS * VSS
        cmp_dct['XHMO'] = frXHMO_VSS * VSS
        cmp_dct['XPRO'] = frXPRO_VSS * VSS
        cmp_dct['XMEOLO'] = frXMEOLO_VSS * VSS
        cmp_dct['XFO'] = frXFO_VSS * VSS
        cmp_dct['XOHO_PHA'] = frXOHO_PHA_VSS * VSS
        cmp_dct['XGAO_PHA'] = frXGAO_PHA_VSS * VSS
        cmp_dct['XPAO_PHA'] = frXPAO_PHA_VSS * VSS
        cmp_dct['XGAO_Gly'] = frXGAO_Gly_VSS * VSS
        cmp_dct['XPAO_Gly'] = frXPAO_Gly_VSS * VSS
        
        cmp_c = np.asarray([v for v in cmp_dct.values()])
        cmp_dct['XB_Subst'] = VSS - (cmp_c * cmps.org * cmps.x).sum()
        
        # convert gVSS to gCOD
        for cmp in cmps:
            if cmp.organic and cmp.particle_size == 'Particulate':
                cmp_dct[cmp.ID] /= cmp.i_mass * cmp.f_Vmass_Totmass
        
        cmp_dct['XPAO_PP_Hi'] = XPAO_PP * r['iHi_XPAOPP']
        cmp_dct['XPAO_PP_Lo'] = XPAO_PP * (1 - r['iHi_XPAOPP'])
        
        cmp_c = np.asarray([v for v in cmp_dct.values()])
        ig_ISS = TSS - (cmp_c * cmps.i_mass * cmps.x * cmps.org).sum()
        other_ig_iss = (cmp_c * cmps.i_mass * cmps.x * (1-cmps.org)).sum()
        cmp_dct['XIg_ISS'] = ig_ISS - other_ig_iss
        
        del other_ig_iss, cmp_c
        
        #*********** soluble and colloidal components *************
        cmp_c = np.asarray([v for v in cmp_dct.values()])
        xCOD = (cmp_c * cmps.x * cmps.org).sum()
        scCOD = xCOD * iscCOD_COD / (1-iscCOD_COD)
        
        cmp_dct['SCH3OH'] = frSCH3OH_scCOD * scCOD
        cmp_dct['SAc'] = frSAc_scCOD * scCOD
        cmp_dct['SProp'] = frSProp_scCOD * scCOD
        cmp_dct['SU_Inf'] = frSU_scCOD * scCOD * r['iSUInf_SU']
        cmp_dct['SU_E'] = frSU_scCOD * scCOD - cmp_dct['SU_Inf']
        
        CB = frCB_scCOD * scCOD
        cmp_dct['CB_BAP'] = CB * r['iBAP_CB']
        cmp_dct['CB_UAP'] = CB * r['iUAP_CB']
        cmp_dct['CB_Subst'] = CB - cmp_dct['CB_BAP'] - cmp_dct['CB_UAP']        
        cmp_dct['CU_Inf'] = cmp_dct['XU_Inf'] * r['iCUInf_XCUInf'] / (1-r['iCUInf_XCUInf'])
        
        cmp_c = np.asarray([v for v in cmp_dct.values()])
        cmp_dct['SF'] = scCOD - (cmp_c * (cmps.s + cmps.c) * cmps.org).sum()
                
        # TODO: calibrate pH, SAlk, SCAT, SAN
        bad_vars = {k:v for k,v in cmp_dct.items() if v<0}
        if len(bad_vars) > 0:
            raise ValueError(f"The following state variable(s) was found negative: {bad_vars}.")
        
        del bad_vars
        
        #************ calibrate XB_subst, SF's N, P content *************
        if SNH4 > 0 and cmp_dct['SF'] > 0:
            SFi_N = _calib_SF_iN(cmps, cmp_dct, SNH4/iSNH_STKN)
        XB_Substi_N = _calib_XBsub_iN(cmps, cmp_dct, TKN - SNH4/iSNH_STKN)
        XB_Substi_P = _calib_XBsub_iP(cmps, cmp_dct, TP)

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
            new = cls(ID=ID, phase=phase, T=T, P=P, units='kg/hr', price=price, 
                      thermo=thermo, pH=pH, SAlk=SAlk, **cmp_dct)
            den = flow_tot*den0/(new.F_vol*1e3)            
            i += 1
            if abs(den-den0) <= 1e-3: break
            if i > 50: raise ValueError('Density calculation failed to converge within 50 iterations.')
        

        new.components.SF.i_N = SFi_N
        new.components.XB_Subst.i_N = XB_Substi_N
        new.components.XB_Subst.i_P = XB_Substi_P        

        return new

#%% functions for calibrations of N, P contents and BOD:COD ratio of certain components
        
def _calib_SF_iN(components, concentrations, STKN):    
    cmp_c = np.asarray([v for v in concentrations.values()])
    SN = (cmp_c * components.i_N * (components.s + components.c)).sum()
    SF_N = concentrations['SF'] * components.SF.i_N
    SNOx_N = concentrations['SNO2'] * components.SNO2.i_N + concentrations['SNO3'] * components.SNO3.i_N
    other_stkn = SN - SF_N - SNOx_N
    SF_N = STKN - other_stkn           
    if SF_N < 0:
        raise ValueError("Negative N content for `SF` was estimated.")                        
    return SF_N/concentrations['SF']
    
def _calib_XBsub_iN(components, concentrations, XTKN):
    cmp_c = np.asarray([v for v in concentrations.values()])        
    other_xtkn = (cmp_c * components.i_N * components.x).sum() - concentrations['XB_Subst'] * components.XB_Subst.i_N                
    XB_Subst_N = XTKN - other_xtkn
    if XB_Subst_N < 0:
        raise ValueError("Negative N content for `XB_Subst` was estimated.")            
    return XB_Subst_N/concentrations['XB_Subst']
    

def _calib_XBsub_iP(components, concentrations, TP):
    cmp_c = np.asarray([v for v in concentrations.values()])    
    other_p = (cmp_c * components.i_P).sum() - concentrations['XB_Subst'] * components.XB_Subst.i_P
    XB_Subst_P = TP - other_p
    if XB_Subst_P < 0:
        raise ValueError("Negative P content for `XB_Subst` was estimated.")    
    return XB_Subst_P/concentrations['XB_Subst']
    
def _calib_XBsub_fBODCOD(components, concentrations, substrate_IDs, BOD):
    cmp_c = np.asarray([v for v in concentrations.values()])
    c_sub = np.asarray([v for k,v in concentrations.items() if k in substrate_IDs])
    XB_sub = components.subgroup(substrate_IDs)
    other_BOD = (cmp_c * (components.x + components.c + components.s) * components.f_BOD5_COD).sum() - (c_sub * XB_sub.f_BOD5_COD).sum()
    fbodtocod_sub = (BOD - other_BOD)/c_sub.sum()
    if fbodtocod_sub > 1 or fbodtocod_sub < 0:
        raise ValueError("BOD5-to-COD ratio for `XB_Subst` and `XStor` was estimated out of range [0,1].")
    return fbodtocod_sub

