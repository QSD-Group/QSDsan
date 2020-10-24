#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 17 07:48:41 2020

@author: yalinli_cabbi, joy_c
"""

import pandas as pd
import numpy as np
import thermosteam as tmo
from biosteam.settings import set_thermo
from thermosteam import Stream, utils
from sanitation import Component
from sanitation.utils import load_components



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
                    'TKN': (-'SNO2', -'SNO3')}
# %%

# =============================================================================
# Define the WasteStream class
# =============================================================================

@utils.registered(ticket_name='ws')
class WasteStream(Stream):
    '''A subclass of the Stream object in the thermosteam package with additional attributes and methods for waste treatment    '''
    
    # TODO: add class attribute - default ratios
    _ratios = pd.read_csv("_ratios.csv")
    _ratios = dict(zip(_ratios.Variable, _ratios.Default))

    # TODO: add other state variables to initiation (pH, SAlk, SCAT, SAN)
    
    def show(self, T=None, P=None, flow='kg/hr', composition=None, N=None,
             stream_info=True):
        '''Show WasteStream information'''        
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
        
        # # Component-related properties
        # info += '\n Component-specific properties:\n'
        # info += f'  charge: {self.charge} mol/hr\n'
        
        print(info)
        
    _ipython_display_ = show
    
    @property
    def components(self):
        return self._thermo.chemicals

    #!!! Suspect many of the units below aren't correct
    # double-check using self.mass or self.mol
    # @property
    # def TC(self):
    #     '''[float] Total carbon content of the stream in + g C/hr'''
    #     return (self._thermo.chemicals.i_C * self.mass).sum()

    # @property
    # def charge(self):
    #     '''[float] Total charge of the stream in + mol/hr'''
    #     return (self._thermo.chemicals.i_charge * self.mol).sum()

    def composite(self, variable, subgroup=None, particle_size=None, 
                  degradability=None, organic=None, volatile=None,
                  specification=None):
        """
        calculate any composite variable by specifications

        Parameters
        ----------
        variable : str
            The composite variable to calculate. 
            One of ('COD', 'BOD5', 'BOD', 'uBOD', 'C', 'N', 'P', 'K', 'Solids', 'Charge').
        subgroup=None: CompiledComponents, optional
            A subgroup of CompiledComponents. The default is None.
        particle_size=None: 'g', 's', 'c', or 'x', optional 
            Dissolved gas ('g'), soluble ('s'), colloidal ('c'), particulate ('x'). 
            The default is None.
        degradability=None: 'b' or 'u', optional
            Either degradable ('b') or undegradable ('u'). The default is None.
        organic=None: bool, optional
            Organic (True) or inorganic (False). The default is None.
        volatile=None: bool, optional
            Volatile (True) or involatile (False). The default is None.
        specification=None: str, optional
            One of ('SVFA', 'XStor', 'XANO', 'XBio', 'SNOx', 'XPAO_PP','TKN'). 
            The default is None.


        Returns
        -------
        float
            The estimated value of the composite variable.

        """
        
        if variable not in _defined_composite_vars:
            raise KeyError(f"Undefined composite variable {variable}."
                           f"Must be one of {_defined_composite_vars}")            
        
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
                IDs = [ID for ID in IDs if ID not in _specific_groups['SNOx']]
            elif specification not in _specific_groups.keys():
                raise KeyError(f"Undefined specification {specification}."
                               f"Must be one of _specific_groups.keys()."
                               f"Or, try defining 'subgroup'.")
            else: 
                IDs = [ID for ID in IDs if ID in _specific_groups[specification]]
                
        cmps = cmps.subgroup(IDs)
        cmp_c = self.imass[IDs]/self.ivol['H2O']*1e3      #[mg/L]
        
        if variable == 'COD':
            var = cmp_c
            if organic == None: organic = True
        elif variable in ('BOD5', 'BOD'):
            var = cmps.f_BOD5_COD * cmp_c
            if organic == None: organic = True
        elif variable == 'uBOD':
            var = cmps.f_uBOD_COD * cmp_c
            if organic == None: organic = True
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
            else: dummy *= getattr(cmps, 'b')
        
        if organic != None:
            if organic: dummy *= getattr(cmps, 'org')
            else: dummy *= 1-getattr(cmps, 'org')
        
        return sum(dummy*var)
    
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
    
    def TDS(self, include_colloidal=False):
        TDS = self.composite('Solids', particle_size='s')
        if include_colloidal:
            TDS += self.composite('Solids', particle_size='c')
        return TDS
    
    def TSS(self, include_colloidal=True):
        TSS = self.composite('Solids', particle_size='x')
        if include_colloidal:
            TSS += self.composite('Solids', particle_size='c')        
        return TSS
    
    def VSS(self, include_colloidal=True):
        VSS = self.composite('Solids', particle_size='x', volatile=True)
        if include_colloidal:
            VSS += self.composite('Solids', particle_size='c', volatile=True)        
        return VSS        
    
    def ISS(self):
        return self.composite('Solids', particle_size='x', volatile=False)
    
    @classmethod
    def from_composite_measures(cls, ID= '', flow_tot=0., phase='l', T=298.15, P=101325.,
                                units=('L/hr', 'mg/L'), price=0., thermo=None, pH=7., C_Alk=150., COD=430.,
                                TKN=40., TP=10., SNH4=25., SNO2=0., SNO3=0., SPO4=8., XPAO_PP=0.,
                                SCa=140., SMg=50., SK=28., XMeP=0., XMeOH=0., XMAP=0., XHAP=0., 
                                XHDP=0., DO=0., SH2=0., SN2=18., SCH4=0., SCAT=3., SAN=12.):
        
        cmps = load_components("sanitation/utils/default_components.csv")
        H2O = Component.from_chemical('H2O', tmo.Chemical('H2O'),
                                      i_C=0, i_N=0, i_P=0, i_K=0, i_mass=1,
                                      i_charge=0, f_BOD5_COD=0, f_uBOD_COD=0,
                                      f_Vmass_Totmass=0,
                                      particle_size='Soluble',
                                      degradability='Undegradable', organic=False)
        cmps.append(H2O)
        TMH = tmo.base.thermo_model_handle.ThermoModelHandle
        for i in cmps:
            i.default()
            for j in ('sigma', 'epsilon', 'kappa', 'V', 'Cn', 'mu'):
                if isinstance(getattr(i, j), TMH) and len(getattr(i, j).models) > 0: continue
                i.copy_models_from(cmps.H2O, names=(j,))
        
        cmps.compile()
        set_thermo(cmps)
        
        cmp_dct = dict.fromkeys(cmps.IDs, 0.)
        r = WasteStream._ratios

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
        cmp_dct['XMeP'] = XMeP
        cmp_dct['XMeOH'] = XMeOH
        
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
                raise ValueError(f"negative N content for SF was estimated.")            
            
            cmps.SF.i_N = SF_N/cmp_dct['SF']
            
            del STKN, SN, SF_N, other_stkn
            
        other_tkn = (cmp_c*cmps.i_N).sum() - SNOx_N - SN2*cmps.SN2.i_N - cmp_dct['XB_Subst']*cmps.XB_Subst.i_N                
        XB_Subst_N = TKN - other_tkn
        if XB_Subst_N < 0:
            raise ValueError(f"negative N content for XB_Subst was estimated.")            
        cmps.XB_Subst.i_N = XB_Subst_N/cmp_dct['XB_Subst']
        
        other_p = (cmp_c*cmps.i_P).sum() - cmp_dct['XB_Subst']*cmps.XB_Subst.i_P
        XB_Subst_P = TP - other_p
        if XB_Subst_P < 0:
            raise ValueError(f"negative P content for XB_Subst was estimated.")    
        cmps.XB_Subst.i_P = XB_Subst_P/cmp_dct['XB_Subst']
        
        del other_tkn, XB_Subst_N, other_p, XB_Subst_P, cmp_c
        
        #************ convert concentrations to flow rates *************
        # TODO: other unit options
        cmp_dct = {k:v*flow_tot*1e-6 for k,v in cmp_dct.items()}       # [mg/L]*[L/hr]*1e-6[kg/mg] = [kg/hr]
        cmp_dct['H2O'] = flow_tot                                      # [L/hr]*1[kg/L] = [kg/hr]
        
        new = cls.__init__(ID=ID, phase=phase, T=T, P=P, units='kg/hr', 
                           price=price, thermo=thermo, pH=pH, SAlk=C_Alk, 
                           SCAT=SCAT, SAN=SAN, **cmp_dct)
        
        return new