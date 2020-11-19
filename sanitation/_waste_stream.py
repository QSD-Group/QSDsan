#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Sanitation Explorer: Sustainable design of non-sewered sanitation technologies
Copyright (C) 2020, Sanitation Explorer Development Group

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>
    Joy Cheung

This module is under the UIUC open-source license. Please refer to 
https://github.com/QSD-for-WaSH/sanitation/blob/master/LICENSE.txt
for license details.
'''

# %%
import numpy as np
import biosteam as bst
from thermosteam import Stream, utils, units_of_measure
from . import Components



__all__ = ('WasteStream',)

_defined_composite_vars = ('COD', 'BOD5', 'BOD', 'uBOD', 'C',
                           'N', 'P', 'K', 'Solids', 'Charge')

_specific_groups = {'SVFA': ('SAc', 'SProp'),
                    'XStor': ('XOHO_PHA', 'XGAO_PHA', 'XPAO_PHA', 
                              'XGAO_Gly', 'XPAO_Gly'),
                    'XANO': ('XAOO', 'XNOO'),
                    'XBio': ('XOHO', 'XAOO', 'XNOO', 'XAMO', 'XPAO', 
                             'XMEOLO', 'XACO', 'XHMO', 'XPRO'),
                    'SNOx': ('SNO2', 'SNO3'),
                    'XPAO_PP': ('XPAO_PP_Lo', 'XPAO_PP_Hi'),
                    'TKN': ()}

# path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'default_data/_ratios.csv')
# _default_ratios = pd.read_csv(path)
_default_ratios = {'iHi_XPAOPP': 0.5,
                   'iCB_XCB': 0.15,
                   'iBAP_CB': 0.,
                   'iUAP_CB': 0.,
                   'iCUInf_XCUInf': 0.,}

# del os, path
vol_unit = units_of_measure.AbsoluteUnitsOfMeasure('L/hr')
conc_unit = units_of_measure.AbsoluteUnitsOfMeasure('mg/L')

# %%

# =============================================================================
# Define the WasteStream class
# =============================================================================

@utils.registered(ticket_name='ws')
class WasteStream(Stream):
    '''
    A subclass of the Stream object in the thermosteam package with additional
    attributes and methods for waste treatment.
    '''
    
    _default_ratios = _default_ratios
    
    def __init__(self, ID='', flow=(), phase='l', T=298.15, P=101325.,
                 units='kg/hr', price=0., thermo=None, pH=7., SAlk=2.5,
                 **chemical_flows):
        
        super().__init__(ID=ID, flow=flow, phase=phase, T=T, P=P,
                         units=units, price=price, thermo=thermo, **chemical_flows)
        self._pH = pH
        self._SAlk = SAlk
        self._ratios = None
    
    def show(self, T='K', P='Pa', flow='kg/hr', composition=False, N=7,
             stream_info=True, details=True):
        '''
        Print WasteStream information.

        Parameters
        ----------
        T : [str], optional
            The unit for temperature. The default is 'K'.
        P : [float], optional
            The unit for pressure. The default is 'Pa'.
        flow : [str], optional
            The unit for the flow. The default is 'kg/hr'.
        composition : [bool], optional
            Whether to show flow information of different Component objects in
            the WasteStream as a percentage. The default is False.
        N : [int], optional
            Number of Component objects to print out, when left as None,
            the number depends on the default of thermosteam. The default is None.
        stream_info : [bool], optional
            Whether to print Stream-specific information. The default is True.
        details : [bool], optional
            Whether to show the all composite variables of the WasteStream. The default is True.

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
    
    #!!! Need reviewing
    def _wastestream_info(self, details=True):
        _ws_info = '\n WasteStream-specific properties:'
        if self.F_mass == 0:
            _ws_info += ' None'
        else:
            _ws_info += '\n'
            _ws_info += f'  pH         : {self.pH:.1f}\n'
            #TODO: unit and definition for following properties
            _ws_info += f'  Alkalinity : {self.SAlk:.1f} mmol/L\n'
            # Or we can just print all...
            if details:
                _ws_info += f'  COD        : {self.COD:.1f} mg/L\n'
                _ws_info += f'  BOD        : {self.BOD:.1f} mg/L\n'
                _ws_info += f'  TC         : {self.TC:.1f} mg/L\n'
                _ws_info += f'  TOC        : {self.TOC:.1f} mg/L\n'
                _ws_info += f'  TN         : {self.TN:.1f} mg/L\n'
                _ws_info += f'  TKN        : {self.TKN:.1f} mg/L\n'
                _ws_info += f'  TP         : {self.TP:.1f} mg/L\n'
                _ws_info += f'  TK         : {self.TK:.1f} mg/L\n'
                # _ws_info += f'  charge     : {self.charge:.1f} mmol/L\n'
            else:
                _ws_info += '  ...\n'
            
        return _ws_info

    @property
    def components(self):
        return self.chemicals

    @property
    def ratios(self):
        '''
        The ratios used for estimating WasteStream composition based on user input upon initialization.
        Only meaningful for creating a WasteStream object from scratch.
        If not used or specified, default as None.
        '''
        return self._ratios
    @ratios.setter
    def ratios(self, ratios):
        r = self._ratios or WasteStream._default_ratios
        for name, ratio in ratios.items():
            if name not in r.keys():
                raise ValueError(f"Cannot identify ratio named '{name}'."
                                 f"Must be one of {r.keys()}")
            elif ratio > 1 or ratio < 0:
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
        variable : [str]
            The composite variable to calculate. 
            One of ('COD', 'BOD5', 'BOD', 'uBOD', 'C', 'N', 'P', 'K', 'Solids', 'Charge').
        subgroup : CompiledComponents, optional
            A subgroup of CompiledComponents. The default is None.
        particle_size : 'g', 's', 'c', or 'x', optional 
            Dissolved gas ('g'), soluble ('s'), colloidal ('c'), particulate ('x'). 
            The default is None.
        degradability : 'rb', 'sb', or 'u', optional
            Readily biodegradable ('rb'), slowly biodegradable ('sb'), 
            biodegradable ('b'), or undegradable ('u'). The default is None.
        organic : [bool], optional
            Organic (True) or inorganic (False). The default is None.
        volatile : [bool], optional
            Volatile (True) or involatile (False). The default is None.
        specification : [str], optional
            One of ('SVFA', 'XStor', 'XANO', 'XBio', 'SNOx', 'XPAO_PP','TKN'). 
            The default is None.

        Returns
        -------
        [float]
            The estimated value of the composite variable, in [mg/L] or [mmol/L] (for "Charge").

        """
        
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
        
        if variable in ('COD', 'BOD5', 'BOD', 'uBOD'):            
            if organic == False: var = 0.
            else: 
                organic = True
                exclude_gas = getattr(cmps, 's')+getattr(cmps, 'c')+getattr(cmps, 'x')
                if variable == 'COD': var = cmp_c * exclude_gas
                elif variable == 'uBOD': var = cmps.f_uBOD_COD * cmp_c * exclude_gas
                else: var = cmps.f_BOD5_COD * cmp_c * exclude_gas
        elif variable == 'C':
            var = cmps.i_C * cmp_c
        elif variable == 'N':
            var = cmps.i_N * cmp_c
        elif variable == 'P':
            var = cmps.i_P * cmp_c
        elif variable == 'K':
            var = cmps.i_K * cmp_c
        elif variable == 'Solids':
            var = cmps.i_mass * cmp_c
            if volatile != None:
                if volatile: var *= cmps.f_Vmass_Totmass
                else: var *= 1-cmps.f_Vmass_Totmass
        else:
            var = cmps.i_charge * cmp_c
        
        dummy = np.ones(len(cmp_c))
        if particle_size:
            if particle_size == 'g': 
                dummy *= 1-getattr(cmps, 's')-getattr(cmps, 'c')-getattr(cmps, 'x')
            else:
                dummy *= getattr(cmps, particle_size)
        
        if degradability:
            if degradability == 'u': dummy *= 1-getattr(cmps, 'b')
            elif degradability == 'b': dummy *= getattr(cmps, 'b')
            elif degradability == 'rb': dummy *= getattr(cmps, 'rb')
            else: dummy *= getattr(cmps, 'b')-getattr(cmps, 'rb')
        
        if organic != None:
            if organic: dummy *= getattr(cmps, 'org')
            else: dummy *= 1-getattr(cmps, 'org')
        
        return sum(dummy*var)
    
    @property
    def pH(self):
        return self._pH

    @property
    def SAlk(self):
        return self._SAlk

    @property
    def COD(self):
        return self.composite('COD')

    @property    
    def BOD(self):
        return self.composite('BOD')
    
    @property
    def TC(self):
        return self.composite('C')
    
    @property
    def TOC(self):
        return self.composite('C', organic=True)
        
    @property
    def TN(self):
        return self.composite('N')
    
    @property
    def TKN(self):
        return self.composite('N', specification='TKN')
    
    @property
    def TP(self):
        return self.composite('P')
    
    @property
    def TK(self):
        return self.composite('K')
    
    # TODO: calibrate Charge when weak acids are involved
    # @property
    # def charge(self):
    #     return self.composite('Charge')
    
    # Below are funtions, not properties (i.e., need to be called), so changed names accordingly
    def get_TDS(self, include_colloidal=True):
        TDS = self.composite('Solids', particle_size='s')
        if include_colloidal:
            TDS += self.composite('Solids', particle_size='c')
        return TDS
    
    def get_TSS(self, include_colloidal=False):
        TSS = self.composite('Solids', particle_size='x')
        if include_colloidal:
            TSS += self.composite('Solids', particle_size='c')        
        return TSS
    
    def get_VSS(self, include_colloidal=False):
        VSS = self.composite('Solids', particle_size='x', volatile=True)
        if include_colloidal:
            VSS += self.composite('Solids', particle_size='c', volatile=True)        
        return VSS        
    
    def get_ISS(self):
        return self.composite('Solids', particle_size='x', volatile=False)
    
    @classmethod
    def from_composite_measures(cls, ID, flow_tot=0., phase='l', T=298.15, P=101325., 
                                units=('L/hr', 'mg/L'), price=0., thermo=None, 
                                pH=7., SAlk=150., ratios=None, COD=430., TKN=40., TP=10., 
                                SNH4=25., SNO2=0., SNO3=0., SPO4=8., SCa=140., SMg=50., 
                                SK=28., XMeP=0., XMeOH=0., XMAP=0., XHAP=0., XHDP=0., 
                                XPAO_PP=0., DO=0., SH2=0., SN2=18., SCH4=0., SCAT=3., SAN=12.):
                
        cmps = Components.load_default(default_compile=True)
        bst.settings.set_thermo(cmps)
        
        cmp_dct = dict.fromkeys(cmps.IDs, 0.)

        new = cls(ID=ID, phase=phase, T=T, P=P, units='kg/hr', price=price, 
                  thermo=thermo, pH=pH, SAlk=SAlk)

        if ratios:
            new.ratios = ratios
            r = new._ratios
        else: r = WasteStream._default_ratios

        #************ user-defined values **************        
        cmp_dct['SH2'] = SH2
        cmp_dct['SCH4'] = SCH4
        cmp_dct['SN2'] = SN2
        cmp_dct['SO2'] = DO
        cmp_dct['SNH4'] = SNH4
        cmp_dct['SNO2'] = SNO2
        cmp_dct['SNO3'] = SNO3
        cmp_dct['SPO4'] = SPO4
        cmp_dct['SCa'] = SCa
        cmp_dct['SMg'] = SMg
        cmp_dct['SK'] = SK
        cmp_dct['XMAP'] = XMAP
        cmp_dct['XHAP'] = XHAP
        cmp_dct['XHDP'] = XHDP
        # cmp_dct['XMeP'] = XMeP
        # cmp_dct['XMeOH'] = XMeOH
        cmp_dct['SCAT'] = SCAT
        cmp_dct['SAN'] = SAN
        
        #************ organic components **************
        cmp_dct['SCH3OH'] = COD * r['fSCH3OH_TotCOD']
        cmp_dct['SAc'] = COD * r['fSAc_TotCOD']
        cmp_dct['SProp'] = COD * r['fSProp_TotCOD']
        cmp_dct['SF'] = COD * r['fSF_TotCOD']
        cmp_dct['SU_Inf'] = COD * r['fSUInf_TotCOD']
        cmp_dct['SU_E'] = COD * r['fSUE_TotCOD']
        
        SOrg = sum([v for k,v in cmp_dct.items() if k in ('SCH3OH','SAc','SProp','SF','SU_Inf','SU_E')])
        
        XCU_Inf = COD * r['fXCUInf_TotCOD']
        cmp_dct['CU_Inf'] = XCU_Inf * r['fCUInf_XCUInf']
        cmp_dct['XU_Inf'] = XCU_Inf * (1 - r['fCUInf_XCUInf'])
        
        cmp_dct['XOHO'] = COD * r['fXOHO_TotCOD']
        cmp_dct['XAOO'] = COD * r['fXAOO_TotCOD']
        cmp_dct['XNOO'] = COD * r['fXNOO_TotCOD']
        cmp_dct['XAMO'] = COD * r['fXAMO_TotCOD']
        cmp_dct['XPAO'] = COD * r['fXPAO_TotCOD']
        cmp_dct['XACO'] = COD * r['fXACO_TotCOD']
        cmp_dct['XHMO'] = COD * r['fXHMO_TotCOD']
        cmp_dct['XPRO'] = COD * r['fXPRO_TotCOD']
        cmp_dct['XMEOLO'] = COD * r['fXMEOLO_TotCOD']
        
        XBio = sum([v for k,v in cmp_dct.items() if k.startswith('x') and k.endswith('O')])
        
        cmp_dct['XOHO_PHA'] = COD * r['fXOHOPHA_TotCOD']
        cmp_dct['XGAO_PHA'] = COD * r['fXGAOPHA_TotCOD']
        cmp_dct['XPAO_PHA'] = COD * r['fXPAOPHA_TotCOD']
        cmp_dct['XGAO_Gly'] = COD * r['fXGAOGly_TotCOD']
        cmp_dct['XPAO_Gly'] = COD * r['fXPAOGly_TotCOD']
        
        XStor = sum([v for k,v in cmp_dct.items() if k.endswith(('PHA','Gly'))])
        
        cmp_dct['XU_OHO_E'] = COD * r['fXUOHOE_TotCOD']
        cmp_dct['XU_PAO_E'] = COD * r['fXUPAOE_TotCOD']
        
        XU_E = cmp_dct['XU_OHO_E'] + cmp_dct['XU_PAO_E']
        
        XCB = COD - SOrg - XCU_Inf - XU_E
        CB = XCB * r['fCB_XCB']
        cmp_dct['CB_BAP'] = CB * r['fBAP_CB']
        cmp_dct['CB_UAP'] = CB * r['fUAP_CB']
        cmp_dct['CB_Subst'] = CB - cmp_dct['CB_BAP'] - cmp_dct['CB_UAP']
        
        cmp_dct['XB_Subst'] = XCB - CB - XBio - XStor
        
        cmp_c = np.asarray([v for v in cmp_dct.values()])
        VSS = (cmp_c * cmps.i_mass * cmps.f_Vmass_Totmass * cmps.x * cmps.org).sum()
        TSS = VSS/r['iVSS_TSS']
        XOrg_ISS = (cmp_c * cmps.i_mass * (1-cmps.f_Vmass_Totmass) * cmps.x * cmps.org).sum()

        del SOrg, XCU_Inf, XBio, XStor, XU_E, XCB, CB
        
        #************ inorganic components **************
        cmp_dct['XPAO_PP_Hi'] = XPAO_PP * r['fHi_XPAOPP']
        cmp_dct['XPAO_PP_Lo'] = XPAO_PP * (1 - r['fHi_XPAOPP'])
        
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
            STKN = SNH4/r['fSNH4_STKN']
            cmp_c = np.asarray([v for v in cmp_dct.values()])
            SN = (cmp_c * cmps.i_N * cmps.s).sum()
            SF_N = cmp_dct['SF'] * cmps.SF.i_N
            SNOx_N = SNO2 * cmps.SNO2.i_N + SNO3 * cmps.SNO3.i_N
            other_stkn = SN - SF_N - SNOx_N - SN2*cmps.SN2.i_N
            SF_N = STKN - other_stkn
            
            if SF_N < 0:
                raise ValueError("Negative N content for SF was estimated.")            
            
            cmps.SF.i_N = SF_N/cmp_dct['SF']
            
            del STKN, SN, SF_N, other_stkn
            
        other_tkn = (cmp_c*cmps.i_N).sum() - SNOx_N - SN2*cmps.SN2.i_N - cmp_dct['XB_Subst']*cmps.XB_Subst.i_N                
        XB_Subst_N = TKN - other_tkn
        if XB_Subst_N < 0:
            raise ValueError("Negative N content for XB_Subst was estimated.")            
        cmps.XB_Subst.i_N = XB_Subst_N/cmp_dct['XB_Subst']
        
        other_p = (cmp_c*cmps.i_P).sum() - cmp_dct['XB_Subst']*cmps.XB_Subst.i_P
        XB_Subst_P = TP - other_p
        if XB_Subst_P < 0:
            raise ValueError("Negative P content for XB_Subst was estimated.")    
        cmps.XB_Subst.i_P = XB_Subst_P/cmp_dct['XB_Subst']
        
        del other_tkn, XB_Subst_N, other_p, XB_Subst_P, cmp_c
        
        #************ convert concentrations to flow rates *************
        # TODO: other unit options
        cmp_dct = {k:v*flow_tot*1e-6 for k,v in cmp_dct.items()}       # [mg/L]*[L/hr]*1e-6[kg/mg] = [kg/hr]
        cmp_dct['H2O'] = flow_tot - sum(cmp_dct.values())              # [L/hr]*1[kg/L] = [kg/hr], assuming raw WW density = 1kg/L
        
        new = cls(ID=ID, phase=phase, T=T, P=P, units='kg/hr', price=price, 
                  thermo=thermo, pH=pH, SAlk=SAlk, **cmp_dct)
        # for i, j in cmp_dct.items():
        #     setattr(new, i, j)
        return new


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
        cmp_dct['SCO3'] = SAlk * 12              # 1 meq/L SAlk ~ 1 mmol/L HCO3- ~ 12 mg C/L (12 mg C/mmol HCO3-)
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
            STKN = SNH4/iSNH_STKN
            cmp_c = np.asarray([v for v in cmp_dct.values()])
            SN = (cmp_c * cmps.i_N * cmps.s).sum()
            SF_N = cmp_dct['SF'] * cmps.SF.i_N
            SNOx_N = SNO2 * cmps.SNO2.i_N + SNO3 * cmps.SNO3.i_N
            other_stkn = SN - SF_N - SNOx_N - SN2*cmps.SN2.i_N
            SF_N = STKN - other_stkn
            
            if SF_N < 0:
                raise ValueError("Negative N content for SF was estimated.")            
            
            cmps.SF.i_N = SF_N/cmp_dct['SF']
            
            del STKN, SN, SF_N, other_stkn
            
        other_tkn = (cmp_c*cmps.i_N).sum() - SNOx_N - SN2*cmps.SN2.i_N - cmp_dct['XB_Subst']*cmps.XB_Subst.i_N                
        XB_Subst_N = TKN - other_tkn
        if XB_Subst_N < 0:
            raise ValueError("Negative N content for XB_Subst was estimated.")            
        cmps.XB_Subst.i_N = XB_Subst_N/cmp_dct['XB_Subst']
        
        other_p = (cmp_c*cmps.i_P).sum() - cmp_dct['XB_Subst']*cmps.XB_Subst.i_P
        XB_Subst_P = TP - other_p
        if XB_Subst_P < 0:
            raise ValueError("Negative P content for XB_Subst was estimated.")    
        cmps.XB_Subst.i_P = XB_Subst_P/cmp_dct['XB_Subst']
        
        del other_tkn, XB_Subst_N, other_p, XB_Subst_P, cmp_c
        
        #************ convert concentrations to flow rates *************
        f1 = vol_unit.conversion_factor(units[0])
        f2 = conc_unit.conversion_factor(units[1])
            
        cmp_dct = {k:v/f2*flow_tot/f1*1e-6 for k,v in cmp_dct.items()}       # [mg/L]*[L/hr]*1e-6[kg/mg] = [kg/hr]
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
        
        return new

