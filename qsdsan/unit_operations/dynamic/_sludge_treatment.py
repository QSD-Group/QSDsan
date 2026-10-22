# -*- coding: utf-8 -*-
'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:

    Saumitra Rai <raisaumitra9@gmail.com>

    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

from ... import SanUnit, WasteStream
import numpy as np
import thermosteam as tmo
from chemicals.elements import periodic_table, molecular_weight
from ..static._pumping import WWTpump
from warnings import warn
from ..static._pumping import (
    default_F_BM as default_WWTpump_F_BM,
    default_equipment_lifetime as default_WWTpump_equipment_lifetime,
    )

__all__ = ('Thickener', 'Centrifuge', 'Incinerator')

# Assign a bare module of 1 to all
default_F_BM = {
        'Wall concrete': 1.,
        'Slab concrete': 1.,
        'Wall stainless steel': 1.,
        'Scraper': 1,
        'v notch weir': 1,
        'Pumps': 1,
        
        # Centrifuge
        'Bowl stainless steel': 1, 
        'Conveyor': 1
        }
default_F_BM.update(default_WWTpump_F_BM)

#%% Thickener

def calc_f_thick(thickener_perc, TSS_in):
    """Returns thickening factor, i.e., thickened sludge solid concentration to influent solids concentration"""
    if TSS_in > 0:
        thickener_factor = thickener_perc*10000/TSS_in   # underlying assumption is density of mixed liquor = 1 kg/L = 1e6 mg/L
        if thickener_factor < 1: thickener_factor = 1
        return thickener_factor
    else: 
        raise ValueError(f'Influent TSS is not valid: ({TSS_in:.2f} mg/L).')
        
def calc_f_Qu_thin(TSS_removal_perc, thickener_factor):
    """Returns Qu factor (i.e., underflow flowrate to influent flowrate) and 
    thinning factor (i.e., overflow solids concentration to influent solids concentration)"""
    if thickener_factor <= 1:
        Qu_factor = 1
        thinning_factor=0
    else:
        Qu_factor = TSS_removal_perc/(100*thickener_factor)
        thinning_factor = (1 - TSS_removal_perc/100)/(1 - Qu_factor)
    return Qu_factor, thinning_factor
    
class Thickener(SanUnit):
    
    """
    Thickener based on BSM2 Layout.
    
    Parameters
    ----------
    ID : str
        ID for the Thickener. The default is ''.
    ins : class:`WasteStream`
        Influent to the clarifier. Expected number of influent is 1. 
    outs : class:`WasteStream`
        Thickened sludge and effluent.
    thickener_perc : float
        The percentage of solids in the underflow of the thickener.[1]
    TSS_removal_perc : float
        The percentage of suspended solids removed in the thickener.[1]
    solids_loading_rate : float
        Solid loading rate in the thickener in [(kg/hr)/m2]. Default is 4 kg/(m2*hr) [2]
        If the thickener is treating:
        Only Primary clarifier sludge, then expected range: 4-6 kg/(m2*hr)
        Only WAS (treated with air or oxygen): 0.5-1.5 kg/(m2*hr)
        Primary clarifier sludge + WAS: 1.5-3.5 kg/(m2/hr)
    h_thickener = float
        Side water depth of the thickener. Typically lies between 3-4 m. [2]
        Height of tank forming the thickener.
    downward_flow_velocity : float, optional
        Speed on the basis of which center feed diameter is designed [m/hr]. [3]
        The default is 36 m/hr. (10 mm/sec)
    F_BM : dict
        Equipment bare modules.
        
    Examples
    --------
    >>> from qsdsan import set_thermo, Components, WasteStream
    >>> cmps = Components.load_default()
    >>> cmps_test = cmps.subgroup(['S_F', 'S_NH4', 'X_OHO', 'H2O'])
    >>> set_thermo(cmps_test)
    >>> ws = WasteStream('ws', S_F = 10, S_NH4 = 20, X_OHO = 15, H2O=1000)
    >>> from qsdsan.unit_operations import Thickener
    >>> TC = Thickener(ID='TC', ins= (ws), outs=('sludge', 'effluent'))
    >>> TC.simulate()
    >>> sludge, effluent = TC.outs
    >>> sludge.imass['X_OHO']/ws.imass['X_OHO']
    0.98
    >>> TC.show() # doctest: +ELLIPSIS
    Thickener: TC
    ins...
    [0] ws
    phase: 'l', T: 298.15 K, P: 101325 Pa
    flow (g/hr): S_F    1e+04
                    S_NH4  2e+04
                    X_OHO  1.5e+04
                    H2O    1e+06
        WasteStream-specific properties:
         pH         : 7.0
         Alkalinity : 2.5 mmol/L
         COD        : 23916.5 mg/L
         BOD        : 14990.5 mg/L
         TC         : 8313.4 mg/L
         TOC        : 8313.4 mg/L
         TN         : 20400.3 mg/L
         TP         : 368.3 mg/L
         TK         : 68.4 mg/L
         TSS        : 11144.7 mg/L
    outs...
    [0] sludge
    phase: 'l', T: 298.15 K, P: 101325 Pa
    flow (g/hr): S_F    1.56e+03
                    S_NH4  3.12e+03
                    X_OHO  1.47e+04
                    H2O    1.56e+05
        WasteStream-specific properties:
         pH         : 7.0
         COD        : 95948.4 mg/L
         BOD        : 55752.4 mg/L
         TC         : 34693.6 mg/L
         TOC        : 34693.6 mg/L
         TN         : 24614.1 mg/L
         TP         : 1740.2 mg/L
         TK         : 413.6 mg/L
         TSS        : 67366.9 mg/L
    [1] effluent
    phase: 'l', T: 298.15 K, P: 101325 Pa
    flow (g/hr): S_F    8.44e+03
                    S_NH4  1.69e+04
                    X_OHO  300
                    H2O    8.44e+05
        WasteStream-specific properties:
         pH         : 7.0
         COD        : 9978.8 mg/L
         BOD        : 7103.3 mg/L
         TC         : 3209.0 mg/L
         TOC        : 3209.0 mg/L
         TN         : 19584.9 mg/L
         TP         : 102.9 mg/L
         TK         : 1.6 mg/L
         TSS        : 266.0 mg/L

    References
    ----------
    [1] Gernaey, Krist V., Ulf Jeppsson, Peter A. Vanrolleghem, and John B. Copp.
    Benchmarking of control strategies for wastewater treatment plants. IWA publishing, 2014.

    [2] Chapter-21: Solids Thicknening (Table 21.3). WEF Manual of Practice No. 8. 
    6th Edition. Virginia: McGraw-Hill, 2018. 
    
    [3] Introduction to Wastewater Clarifier Design by Nikolay Voutchkov, PE, BCEE.
    
    [4] Metcalf, Leonard, Harrison P. Eddy, and Georg Tchobanoglous. Wastewater 
    engineering: treatment, disposal, and reuse. Vol. 4. New York: McGraw-Hill, 1991.
    """
    
    _N_ins = 1
    _N_outs = 2   # [0] thickened sludge, [1] reject water
    _ins_size_is_fixed = False
    _outs_size_is_fixed = False
        
    def __init__(self, ID='', ins=None, outs=(), thermo=None, isdynamic=False, 
                  init_with='WasteStream', F_BM_default=default_F_BM, thickener_perc=7, 
                  TSS_removal_perc=98, solids_loading_rate=4, h_thickener=4, 
                  downward_flow_velocity= 36, F_BM=default_F_BM, **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo, isdynamic=isdynamic, 
                         init_with=init_with)
        self.thickener_perc = thickener_perc 
        self.TSS_removal_perc = TSS_removal_perc
        self.solids_loading_rate = solids_loading_rate 
        self.h_thickener = h_thickener
        self.downward_flow_velocity = downward_flow_velocity
        self.F_BM.update(F_BM)
        self._mixed = WasteStream(f'{ID}_mixed', thermo = thermo)        
        self._sludge = self.outs[0].copy(f'{ID}_sludge')
        self._thickener_factor = None
        self._thinning_factor = None
        self._Qu_factor = None
        
    @property
    def thickener_perc(self):
        '''The percentage of suspended solids in the thickened sludge, in %.'''
        return self._tp

    @thickener_perc.setter
    def thickener_perc(self, tp):
        if tp is not None:
            if tp>=100 or tp<=0:
                raise ValueError(f'should be between 0 and 100 not {tp}')
            self._tp = tp
        else: 
            raise ValueError('percentage of SS in the underflow of the thickener expected from user')
            
    @property
    def solids_loading_rate(self):
        '''solids_loading_rate is the loading in the thickener'''
        return self._slr

    @solids_loading_rate.setter
    def solids_loading_rate(self, slr):
        if slr is not None:
            self._slr = slr
        else: 
            raise ValueError('solids_loading_rate of the thickener expected from user')
            
    @property
    def TSS_removal_perc(self):
        '''The percentage of suspended solids removed in the thickener'''
        return self._TSS_rmv

    @TSS_removal_perc.setter
    def TSS_removal_perc(self, TSS_rmv):
        if TSS_rmv is not None:
            if TSS_rmv>=100 or TSS_rmv<=0:
                raise ValueError(f'should be between 0 and 100 not {TSS_rmv}')
            self._TSS_rmv = TSS_rmv
        else: 
            raise ValueError('percentage of suspended solids removed in the thickener expected from user')
            
    @property
    def thickener_factor(self):
        inf = self._mixed
        inf.mix_from(self.ins)
        if not self.ins: return
        elif inf.isempty(): return
        else: 
            TSS_in = inf.get_TSS()
            self._Qu_factor = None
            return calc_f_thick(self._tp, TSS_in)
    
    @property
    def thinning_factor(self):
        f_Qu, f_thin = calc_f_Qu_thin(self.TSS_removal_perc, self.thickener_factor)
        return f_thin
    
    @property
    def Qu_factor(self):
        f_Qu, f_thin = calc_f_Qu_thin(self.TSS_removal_perc, self.thickener_factor)
        return f_Qu

    def _update_parameters(self):
        cmps = self.components 
        TSS_in = np.sum(self._state[:-1]*cmps.i_mass*cmps.x)
        self._f_thick = f_thick = calc_f_thick(self._tp, TSS_in)
        self._f_Qu, self._f_thin = calc_f_Qu_thin(self._TSS_rmv, f_thick)
        
    def _run(self):
        mixed = self._mixed
        mixed.mix_from(self.ins)
        x = self.components.x
        uf, of = self.outs
        
        TSS_rmv = self._TSS_rmv
        TSS_in = mixed.get_TSS()
        f_thick = calc_f_thick(self._tp, TSS_in)
        f_Qu, f_thin = calc_f_Qu_thin(TSS_rmv, f_thick)
        
        if f_thick > 1: split_to_uf = (1-x)*f_Qu + x*TSS_rmv/100
        else: split_to_uf = 1
        mixed.split_to(uf, of, split_to_uf)
       
    def _init_state(self):
        Qs = self._ins_QC[:,-1]
        Cs = self._ins_QC[:,:-1]
        self._state = np.append(Qs @ Cs / Qs.sum(), Qs.sum())
        self._dstate = self._state * 0.
        self._update_parameters()
        
    def _update_state(self):
        '''updates conditions of output stream based on conditions of the Thickener''' 

        thickener_factor = self._f_thick
        thinning_factor = self._f_thin
        Qu_factor = self._f_Qu
        x = self.components.x
        
        uf, of = self.outs
        if uf.state is None: uf.state = np.zeros(len(self.components)+1)
        if of.state is None: of.state = np.zeros(len(self.components)+1)

        arr = self._state
        if thickener_factor <= 1: 
            uf.state[:] = arr
            of.state[:] = 0.
        else:
            # For sludge, the particulate concentrations (x) are multipled by thickener factor, and
            # flowrate is multiplied by Qu_factor. The soluble concentrations (1-x) remains same. 
            uf.state[:-1] = arr[:-1] * ((1-x) + x*thickener_factor)
            uf.state[-1] = arr[-1] * Qu_factor            
            # For effluent, the particulate concentrations (x) are multipled by thinning factor, and
            # flowrate is multiplied by Qu_factor. The soluble concentrations (1-x) remains same. 
            of.state[:-1] = arr[:-1] * ((1-x) + x*thinning_factor)
            of.state[-1] = arr[-1] * (1 - Qu_factor)

    def _update_dstate(self):
        '''updates rates of change of output stream from rates of change of the Thickener'''
        
        thickener_factor = self._f_thick
        thinning_factor = self._f_thin
        Qu_factor = self._f_Qu
        x = self.components.x

        uf, of = self.outs
        if uf.dstate is None: uf.dstate = np.zeros(len(self.components)+1)
        if of.dstate is None: of.dstate = np.zeros(len(self.components)+1)
        arr = self._dstate
        if thickener_factor <= 1:
            uf.dstate[:] = arr
            of.dstate[:] = 0.
        else:
            # For sludge, the particulate concentrations are multipled by thickener factor, and
            # flowrate is multiplied by Qu_factor. The soluble concentrations remains same. 
            uf.dstate[:-1] = arr[:-1] * ((1-x) + x*thickener_factor)
            uf.dstate[-1] = arr[-1] * Qu_factor
            
            # For effluent, the particulate concentrations are multipled by thinning factor, and
            # flowrate is multiplied by Qu_factor. The soluble concentrations remains same.
            of.dstate[:-1] = arr[:-1] * ((1-x) + x*thinning_factor)
            of.dstate[-1] = arr[-1] * (1 - Qu_factor)
    
    @property
    def AE(self):
        if self._AE is None:
            self._compile_AE()
        return self._AE

    def _compile_AE(self):
        
        _state = self._state
        _dstate = self._dstate
        _update_state = self._update_state
        _update_dstate = self._update_dstate
        _update_parameters = self._update_parameters
        def yt(t, QC_ins, dQC_ins):
            Q_ins = QC_ins[:, -1]
            C_ins = QC_ins[:, :-1]
            dQ_ins = dQC_ins[:, -1]
            dC_ins = dQC_ins[:, :-1]
            Q = Q_ins.sum()
            C = Q_ins @ C_ins / Q
            _state[-1] = Q
            _state[:-1] = C
            Q_dot = dQ_ins.sum()
            C_dot = (dQ_ins @ C_ins + Q_ins @ dC_ins - Q_dot * C)/Q
            _dstate[-1] = Q_dot
            _dstate[:-1] = C_dot
            _update_parameters()
            _update_state()
            _update_dstate()
        self._AE = yt
         

#%% Centrifuge

class Centrifuge(Thickener):

    """
    Centrifuge based on BSM2 Layout.
    
    Parameters
    ----------
    ID : str
        ID for the Dewatering Unit. The default is ''.
    ins : class:`WasteStream`
        Influent to the Dewatering Unit. Expected number of influent is 1. 
    outs : class:`WasteStream`
        Treated effluent and sludge.
    thickener_perc : float
        The percentage of Suspended Sludge in the underflow of the dewatering unit.[1]
    TSS_removal_perc : float
        The percentage of suspended solids removed in the dewatering unit.[1]
    solids_feed_rate : float
        Rate of solids processed by one centrifuge in dry tonne per day (dtpd).
        Default value is 150 dtpd. [6]  
    g_factor : float
        Factor by which g (9.81 m/s2) is multiplied to obtain centrifugal acceleration. 
        g_factor typically lies between 1500 and 3000. 
        centrifugal acceleration = g * g_factor = k * (RPM)^2 * diameter [3]
    rotational_speed : float
        rotational speed of the centrifuge in rpm. Typical rpm is between 2000-3000 rpm [MOP-8, PAGE 1733]
    LtoD: The ratio of length to diameter of the centrifuge.
        The value typically lies between 3-4. [4]
    polymer_dosage : float
        mass of polymer utilised (lb) per tonne of dry solid waste (lbs/tonne).[5]
        Depends on the type of influents, please refer to [5] for appropriate values. 
        Default value of 20 lbs/tonne is taken from [5], based on Primary + WAS aerated undigested value.  
    h_cylindrical: float
        length of cylindrical portion of dewatering unit.
    h_conical: float
        length of conical portion of dewatering unit.
    
    
    References
    ----------
    [1] Gernaey, Krist V., Ulf Jeppsson, Peter A. Vanrolleghem, and John B. Copp.
    Benchmarking of control strategies for wastewater treatment plants. IWA publishing, 2014.

    [2] Metcalf, Leonard, Harrison P. Eddy, and Georg Tchobanoglous. Wastewater 
    engineering: treatment, disposal, and reuse. Vol. 4. New York: McGraw-Hill, 1991.
    
    [3] Design of Municipal Wastewater Treatment Plants: WEF Manual of Practice 
    No. 8 ASCE Manuals and Reports on Engineering Practice No. 76, Fifth Edition. 
    
    [4] https://www.alibaba.com/product-detail/Multifunctional-Sludge-Dewatering-Decanter-Centrifuge_1600285055254.html?spm=a2700.galleryofferlist.normal_offer.d_title.1cd75229sPf1UW&s=p
    
    [5] United States Environmental Protection Agency (EPA) 'Biosolids Technology Fact Sheet Centrifuge Thickening and Dewatering'  
    
    [6] San Diego (.gov) Chapter - 3 'Solids Treatment Facility' 
    (https://www.sandiego.gov/sites/default/files/legacy/mwwd/pdf/mbc/chapterthree.pdf)
    """
    
    _N_ins = 1
    _N_outs = 2
    _ins_size_is_fixed = False
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, isdynamic=False, 
                  init_with='WasteStream', F_BM_default=default_F_BM, 
                  thickener_perc=28, TSS_removal_perc=98, 
                  solids_feed_rate=70, g_factor=2500, rotational_speed=2500, 
                  LtoD=4, F_BM=default_F_BM,
                  polymer_dosage=20, h_cylindrical=2, h_conical=1, **kwargs):
        
        Thickener.__init__(self, ID=ID, ins=ins, outs=outs, thermo=thermo, isdynamic=isdynamic,
                      init_with=init_with, F_BM_default=F_BM_default, 
                      thickener_perc=thickener_perc, 
                      TSS_removal_perc=TSS_removal_perc, **kwargs)
        
        self.solids_feed_rate = solids_feed_rate
        self.g_factor = g_factor #unitless, centrifugal acceleration = g_factor*9.81
        self.rotational_speed = rotational_speed #in revolution/min
        self.LtoD = LtoD 
        self.polymer_dosage = polymer_dosage #in (lbs,polymer/tonne,solids)
        self.h_cylindrical = h_cylindrical
        self.h_conical = h_conical
        
        
#%% Incinerator

class Incinerator(SanUnit):
    
    """
    Fluidized bed incinerator.
    
    Parameters
    ----------
    ID : str
        ID for the Incinerator Unit. The default is ''.
    ins : class:`WasteStream`
        Influent to the Incinerator Unit. Expected number of influent streams are 3. 
        Please remember the order of influents as {wastestream, air, fuel} 
    outs : class:`WasteStream`
        Flue gas and ash. 
    process_efficiency : float
        The process efficiency of the incinerator unit. Expected value between 0 and 1.
    calorific_value_sludge : float, optional
        The calorific value of the dry sludge in KJ/kg. If None (default),
        calculated from the properties of the components.
    calorific_value_fuel : float, optional
        The calorific value of the fuel in KJ/kg. If None (default),
        calculated from the properties of the components.
    heating_value_basis : str
        Basis of the heats of combustion, 'LHV' (lower heating value, default)
        or 'HHV' (higher heating value).
    ash_component_ID, nitrogen_ID, water_ID, carbon_dioxide_ID, oxygen_ID : str
        IDs of the ash, N2, H2O, CO2, and O2 components.

    Notes
    -----
    Complete combustion is assumed and calculated with BioSTEAM reactions
    (:class:`thermosteam.reaction.ParallelReaction`, see `combustion_reactions`)
    on the mole flows of all influents, i.e., C becomes CO2, H becomes H2O,
    and N becomes N2. As the flows of components with `measured_as` are in the
    measured-as unit, moles are from their real mass (`i_mass`), which is correct
    with or without adjusting the components' MW to the measured-as unit. The O2 needed is taken from the air and the excess is
    in the flue gas (a warning is raised if the air is not enough). Mineral
    solids and other elements go to the ash, and water, N2, CO2, and O2 pass through.

    For each component, BioSTEAM's combustion reaction (from its formula) is used
    when it is defined with the components of the system (and, for organic components,
    its O2 demand agrees with the component's COD). Otherwise (e.g., ASM
    state variables with no formula, or components with P and metals), the
    reaction is from the component's composition: the formula, or for components
    with no formula, the C and N contents together with the mass and COD (to solve
    for H and O); non-organic components with no formula go to the ash.

    The energy balance (in kJ/hr, after simulation), similar to BioSTEAM's
    :class:`biosteam.facilities.BoilerTurbogenerator`, uses the heating values of
    the components (their LHV/HHV if neutral and available, otherwise estimated
    from the O2 consumed by their combustion with Thornton's rule, i.e., 13.1 MJ of
    net heat per kg of O2 [3, 4]):
    `heat_sludge`, `heat_fuel`, and their sum `heat_combustion`; `heat_air` is the
    sensible heat of the air above 298.15 K; `heat_flue_gas` is `process_efficiency`
    times the heat input, and `heat_loss` is the rest. `heat_evaporation` is the
    energy to heat the water in the sludge to boiling and evaporate it, for reference.

    Examples
    --------
    >>> import qsdsan as qs
    >>> cmps = qs.Components.load_default()
    >>> CO2 = qs.Component.from_chemical('S_CO2', chemical='CO2',
    ...                                  particle_size='Soluble',
    ...                                  degradability='Undegradable', 
    ...                                  organic=False)
    >>> cmps_test = qs.Components([cmps.S_F, cmps.S_NH4, cmps.X_OHO, cmps.H2O, 
    ...                            cmps.S_CH4, cmps.S_O2, cmps.S_N2, cmps.S_H2, 
    ...                            cmps.X_Ig_ISS, CO2])
    >>> cmps_test.default_compile(ignore_inaccurate_molar_weight=True)
    >>> qs.set_thermo(cmps_test)
    >>> ws = qs.WasteStream('ws', S_F=10, S_NH4=20, X_OHO=15, H2O=1000)
    >>> natural_gas = qs.WasteStream('nat_gas', phase='g', S_CH4=1000)
    >>> air = qs.WasteStream('air', phase='g', S_O2=1500, S_N2=4900)
    >>> from qsdsan.unit_operations import Incinerator
    >>> Inc = Incinerator(ID='Inc', ins= (ws, air, natural_gas),
    ...                   outs=('flu_gas', 'ash'),
    ...                   isdynamic=False)
    >>> Inc.simulate()
    >>> Inc.show() # doctest: +ELLIPSIS
    Incinerator: Inc
    ins...
    [0] ws
    phase: 'l', T: 298.15 K, P: 101325 Pa
    flow (g/hr): S_F    1e+04
                    S_NH4  2e+04
                    X_OHO  1.5e+04
                    H2O    1e+06
        WasteStream-specific properties:
         pH         : 7.0
         Alkalinity : 2.5 mmol/L
         COD        : 23916.5 mg/L
         BOD        : 14990.5 mg/L
         TC         : 8313.4 mg/L
         TOC        : 8313.4 mg/L
         TN         : 20400.3 mg/L
         TP         : 368.3 mg/L
         TK         : 68.4 mg/L
         TSS        : 11144.7 mg/L
    [1] air
    phase: 'g', T: 298.15 K, P: 101325 Pa
    flow (g/hr): S_O2  1.5e+06
                    S_N2  4.9e+06
        WasteStream-specific properties: None for non-liquid waste streams
    [2] nat_gas
    phase: 'g', T: 298.15 K, P: 101325 Pa
    flow (g/hr): S_CH4  1e+06
        WasteStream-specific properties: None for non-liquid waste streams
    outs...
    [0] flu_gas
    phase: 'g', T: 298.15 K, P: 101325 Pa
    flow (g/hr): H2O    1.62e+06
                    S_O2   4.29e+05
                    S_N2   4.92e+06
                    S_CO2  7.2e+05
        WasteStream-specific properties: None for non-liquid waste streams
    [1] ash
    phase: 's', T: 298.15 K, P: 101325 Pa
    flow (g/hr): X_Ig_ISS  2.05e+03
        WasteStream-specific properties: None for non-liquid waste streams

    The energy balance [kJ/hr] (rounded here to thousands):

    >>> [round(i, -3) for i in (Inc.heat_sludge, Inc.heat_fuel, Inc.heat_flue_gas, Inc.heat_loss)]
    [926000.0, 12541000.0, 12120000.0, 1347000.0]
    >>> round(Inc.heat_evaporation, -3)
    2570000.0

    In a dynamic system, the unit's algebraic equations give the same flows as static simulation:

    >>> Inc_dyn = Incinerator(ID='Inc_dyn', outs=('flu_gas_dyn', 'ash_dyn'), isdynamic=True,
    ...                       ins=(ws.copy('ws_dyn'), air.copy('air_dyn'), natural_gas.copy('gas_dyn')))
    >>> sys = qs.System('sys', path=(Inc_dyn,))
    >>> sys.simulate(t_span=(0, 1), method='BDF')
    >>> flue_gas, ash = Inc_dyn.outs
    >>> [round(flue_gas.imass[i], 1) for i in ('S_CO2', 'S_N2', 'H2O', 'S_O2')]
    [719.5, 4921.3, 1620.9, 429.3]
    >>> round(ash.imass['X_Ig_ISS'], 1)
    2.1

    References
    ----------
    [1] Khuriati, A., P. Purwanto, H. S. Huboyo, Suryono Sumariyah, S. Suryono, and A. B. Putranto. 
    "Numerical calculation based on mass and energy balance of waste incineration in the fixed bed reactor." 
    In Journal of Physics: Conference Series, vol. 1524, no. 1, p. 012002. IOP Publishing, 2020.
    
    [2] Omari, Arthur, Karoli N. Njau, Geoffrey R. John, Joseph H. Kihedu, and Peter L. Mtui. 
    "Mass And Energy Balance For Fixed Bed Incinerators A case of a locally designed incinerator in Tanzania."

    [3] Thornton, W. M. "The relation of oxygen to the heat of combustion of organic compounds."
    Philosophical Magazine and Journal of Science, 6th series, vol. 33, pp. 196-203, 1917.

    [4] Huggett, C. "Estimation of rate of heat release by means of oxygen consumption measurements."
    Fire and Materials, vol. 4, no. 2, pp. 61-65, 1980 (13.1 MJ/kg O2, +/-5%).
    """

    #These are class attributes
    _N_ins = 3
    _N_outs = 2   
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, isdynamic=False, 
                  init_with='WasteStream', F_BM_default=None, process_efficiency=0.90, 
                  calorific_value_sludge=None, calorific_value_fuel=None, 
                  heating_value_basis='LHV', 
                  ash_component_ID='X_Ig_ISS', nitrogen_ID='S_N2', water_ID='H2O',
                  carbon_dioxide_ID='S_CO2', oxygen_ID='S_O2', **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo, isdynamic=isdynamic, 
                         init_with=init_with, F_BM_default=F_BM_default)
        
        self.calorific_value_sludge = calorific_value_sludge #in KJ/kg
        self.calorific_value_fuel  = calorific_value_fuel #in KJ/kg (here the considered fuel is natural gas) 
        self.process_efficiency = process_efficiency
        self.ash_component_ID = ash_component_ID
        self.nitrogen_ID = nitrogen_ID
        self.water_ID = water_ID
        self.carbon_dioxide_ID = carbon_dioxide_ID
        self.oxygen_ID = oxygen_ID
        self.heating_value_basis = heating_value_basis
        self._combustion_data = None
        self.combustion_reactions = None
        self.heat_sludge = None
        self.heat_fuel = None
        self.heat_combustion = None
        self.heat_air = None
        self.heat_evaporation = None
        self.heat_flue_gas = None
        self.heat_loss = None

    @property
    def process_efficiency(self):
        '''Process efficiency of incinerator.'''
        return self._process_efficiency

    @process_efficiency.setter
    def process_efficiency(self, process_efficiency):
        if process_efficiency is not None:
            if process_efficiency>=1 or process_efficiency<=0:
                raise ValueError(f'should be between 0 and 1 not {process_efficiency}')
            self._process_efficiency = process_efficiency
        else: 
            raise ValueError('Process efficiency of incinerator expected from user')
            
    @property
    def calorific_value_sludge(self):
        '''[float or NoneType] Calorific value of the dry sludge in KJ/kg, None to calculate from component properties.'''
        return self._calorific_value_sludge

    @calorific_value_sludge.setter
    def calorific_value_sludge(self, calorific_value_sludge):
        self._calorific_value_sludge = calorific_value_sludge
            
    @property
    def calorific_value_fuel(self):
        '''[float or NoneType] Calorific value of the fuel in KJ/kg, None to calculate from component properties.'''
        return self._calorific_value_fuel

    @calorific_value_fuel.setter
    def calorific_value_fuel(self, calorific_value_fuel):
        self._calorific_value_fuel = calorific_value_fuel

    @property
    def heating_value_basis(self):
        '''[str] Basis of the heats of combustion, 'LHV' (lower heating value) or 'HHV' (higher heating value).'''
        return self._heating_value_basis
    @heating_value_basis.setter
    def heating_value_basis(self, basis):
        if basis not in ('LHV', 'HHV'):
            raise ValueError(f"heating_value_basis must be 'LHV' or 'HHV', not {basis!r}")
        self._heating_value_basis = basis
        
    _atomic_mass = {el: periodic_table[el].MW for el in 'CHON'} # g/mol
    heat_per_O2 = 13.1e6 # J/kg O2 consumed, net heat of combustion of organics (Thornton's rule [3], [4])
    cp_air = 1.005 # kJ/kg/K, the heat capacities of dissolved gas components are not for gases
    T_ref = 298.15 # K, reference temperature for the sensible heat of air

    @staticmethod
    def _get_real_mass_per_mol(cmp):
        # Components with a formula have their molecular weight, others (e.g., ASM
        # state variables) have MW as the measured-as mass and i_mass as the real mass
        return cmp.chem_MW if cmp.formula else cmp.i_mass * cmp.MW

    def _get_composition_reaction(self, cmp, use_formula=True):
        """
        Return a :class:`thermosteam.reaction.Reaction` for the complete combustion
        of a component from its composition, for components whose
        combustion is not defined by BioSTEAM (no formula, or products that are
        not components, e.g., P4O10 and Ash for P and metals).

        C becomes CO2, H becomes H2O, N becomes N2, and the O already bound in the
        component reduces the O2 needed. Mineral solids and elements other than
        C, H, O, and N go to the ash.

        - Components with a formula are split into elements by their formula
          (unless `use_formula` is False).
        - Non-organic components without a formula (e.g., X_Ig_ISS) go to the ash.
        - Organic components without a formula use their C and N contents, and H
          and O are solved from the real mass and COD (COD = 2.664 C + 7.937 H - O).
        """
        cmps = self.components
        M = self._atomic_mass
        real = self._get_real_mass_per_mol(cmp) # real mass per mol
        measured = real / cmp.i_mass # measured-as mass per mol (MW if no formula)
        atoms = cmp.atoms if (cmp.formula and use_formula) else {}
        O2_per_C = 2*M['O']/M['C'] # g O2 per g C (to CO2)
        O2_per_H = molecular_weight({'O': 2})/(4*M['H']) # g O2 per g H (to H2O)
        if atoms:
            C, H, O, N = (atoms.get(el, 0) * M[el] for el in 'CHON') # g per mol
            ash = max(cmp.chem_MW - C - H - O - N, 0.)
        elif not cmp.organic:
            C = H = O = N = 0.
            ash = real
        else:
            C, N = cmp.i_C * measured, cmp.i_N * measured
            iss = real*(1-cmp.f_Vmass_Totmass) if cmp.particle_size == 'Particulate' else 0.
            metals = (cmp.i_P + cmp.i_K + cmp.i_Mg + cmp.i_Ca) * measured
            ash = min(max(real - C - N, 0.), iss + metals)
            volatile = max(real - C - N - ash, 0.)
            H = min(max((cmp.i_COD*measured + volatile - O2_per_C*C)/(O2_per_H+1), 0.), volatile)
            O = volatile - H
        ash_cmp = cmps[self.ash_component_ID]
        stoichiometry = {cmp.ID: -1.}
        for ID, mol in ((self.carbon_dioxide_ID, C/M['C']),
                        (self.water_ID, H/(2*M['H'])),
                        (self.nitrogen_ID, N/(2*M['N'])),
                        (self.oxygen_ID, -(C/M['C'] + H/(4*M['H']) - O/(2*M['O']))),
                        (self.ash_component_ID, ash/self._get_real_mass_per_mol(ash_cmp))):
            stoichiometry[ID] = stoichiometry.get(ID, 0.) + mol
        return tmo.reaction.Reaction(stoichiometry, cmp.ID, 1., cmps)

    def _compile_combustion(self):
        """
        Compile the combustion reactions (`combustion_reactions`) of all components,
        their linear operator on the mole flows, and the heating values [J/mol].

        Water, N2, CO2, O2, and the ash component pass through. BioSTEAM's
        combustion reaction of the component is used when available, otherwise
        the reaction is from the component's composition.
        """
        cmps = self.components
        n = len(cmps)
        pass_through = {cmps.index(i) for i in (self.water_ID, self.nitrogen_ID,
                                                self.carbon_dioxide_ID, self.oxygen_ID,
                                                self.ash_component_ID)}
        i_o2 = cmps.index(self.oxygen_ID)
        M_O2 = molecular_weight({'O': 2})
        reactions = []
        for j, cmp in enumerate(cmps):
            if j in pass_through: continue
            rxn = None
            use_formula = True
            if cmp.formula:
                try: rxn = cmp.get_combustion_reaction(cmps)
                except Exception: rxn = None # products that are not components
            if rxn is not None and cmp.organic and cmp.i_COD > 0:
                # The COD of organic components must equal their O2 demand; BioSTEAM's
                # reaction can differ for ions (e.g., acetate), whose formula has no counter ion
                mol = np.zeros(n)
                mol[j] = 1.
                rxn.force_reaction(mol)
                cod = cmp.i_COD * self._get_real_mass_per_mol(cmp) / cmp.i_mass # g COD per mol
                if abs(-mol[i_o2]*M_O2 - cod) > 0.01*cod: rxn, use_formula = None, False
            if rxn is None: rxn = self._get_composition_reaction(cmp, use_formula)
            reactions.append(rxn)
        self.combustion_reactions = rxns = tmo.reaction.ParallelReaction(reactions)
        # All conversions are 100% and each component is the reactant of at most one
        # reaction, so the reactions are a linear operator on the mole flows
        operator = np.zeros((n, n))
        for j in range(n):
            mol = np.zeros(n)
            mol[j] = 1.
            rxns.force_reaction(mol)
            operator[:, j] = mol

        # Heating values [J/mol]: component LHV/HHV for neutral components that have it,
        # otherwise Thornton's rule (net heat per mass of O2 consumed) for the combustion reaction
        i_o2, i_h2o = cmps.index(self.oxygen_ID), cmps.index(self.water_ID)
        Hvap_water = cmps[self.water_ID].Hvap(self.T_ref)
        heat_per_mol_O2 = self.heat_per_O2 * molecular_weight({'O': 2}) / 1000 # J/mol O2
        hhv, lhv = np.zeros(n), np.zeros(n)
        for j, cmp in enumerate(cmps):
            if j in pass_through: continue
            if cmp.formula and not cmp.i_charge and (cmp.HHV or 0) > 0 and (cmp.LHV or 0) > 0:
                hhv[j], lhv[j] = cmp.HHV, cmp.LHV
            else:
                lhv[j] = heat_per_mol_O2 * max(-operator[i_o2, j], 0.)
                hhv[j] = lhv[j] + Hvap_water * max(operator[i_h2o, j], 0.)
        # Flows are in the `measured_as` unit of each component, so the moles are from the
        # real mass (i_mass), which is correct whether or not the component's MW was
        # adjusted to the measured-as unit when compiling the components
        real_mass_per_mol = np.array([self._get_real_mass_per_mol(i) for i in cmps])
        i_mass = np.asarray(cmps.i_mass, dtype=float)
        self._combustion_data = dict(
            operator=operator, hv=dict(HHV=hhv, LHV=lhv),
            i_mass=i_mass, real_mass_per_mol=real_mass_per_mol,
            to_mol=i_mass/real_mass_per_mol, to_mass=real_mass_per_mol/i_mass,
            indices=tuple(cmps.index(i) for i in (self.ash_component_ID, self.carbon_dioxide_ID,
                                                  self.nitrogen_ID, self.water_ID, self.oxygen_ID)),
            )

    def _get_balance(self, mass):
        """
        Return the (ash, CO2, N2, H2O, O2) masses after complete combustion,
        in the same mass unit as the input, from the total mass of each
        component in all influents (array).
        """
        if self._combustion_data is None: self._compile_combustion()
        data = self._combustion_data
        out = (data['operator'] @ (mass*data['to_mol'])) * data['to_mass']
        return tuple(out[i] for i in data['indices'])

    def _run(self):
        sludge, air, fuel = self.ins
        flue_gas, ash = self.outs
        flue_gas.phase = 'g'
        ash.phase = 's'

        if sludge.phase != 'l':
            raise ValueError(f'The phase of incoming sludge is expected to be liquid not {sludge.phase}')
        if air.phase != 'g':
            raise ValueError(f'The phase of air is expected to be gas not {air.phase}')
        if fuel.phase != 'g':
            raise ValueError(f'The phase of fuel is expected to be gas not {fuel.phase}')

        mass_ash, co2, n2, h2o, o2 = self._get_balance(np.asarray(sludge.mass + air.mass + fuel.mass))
        if o2 < 0:
            warn(f'{self.ID} does not receive enough oxygen for complete combustion '
                 f'(short by {-o2:.4g} kg/hr); the flue gas oxygen is set to zero.')
            o2 = 0.
        flue_gas.empty()
        ash.empty()
        flue_gas.set_flow([n2, h2o, co2, o2], 'kg/hr',
                          (self.nitrogen_ID, self.water_ID, self.carbon_dioxide_ID, self.oxygen_ID))
        ash.set_flow(mass_ash, 'kg/hr', (self.ash_component_ID,))

    def _design(self):
        # Energy balance [kJ/hr], following BioSTEAM's BoilerTurbogenerator: the heat of
        # combustion is the sum of the heating values of the feeds, and the efficiency
        # defines the heat content of the flue gas and the heat loss
        if self._combustion_data is None: self._compile_combustion()
        data = self._combustion_data
        hv = data['hv'][self.heating_value_basis] # J/mol = kJ/kmol
        to_mol, i_mass = data['to_mol'], data['i_mass']
        cmps = self.components
        sludge, air, fuel = self.ins
        i_h2o = cmps.index(self.water_ID)

        if self.calorific_value_sludge is None:
            self.heat_sludge = float((sludge.mass*to_mol) @ hv)
        else:
            sludge_mass = sludge.mass * i_mass # real mass
            self.heat_sludge = (sludge_mass.sum()-sludge_mass[i_h2o]) * self.calorific_value_sludge
        if self.calorific_value_fuel is None:
            self.heat_fuel = float((fuel.mass*to_mol) @ hv)
        else:
            self.heat_fuel = float((fuel.mass*i_mass).sum()) * self.calorific_value_fuel
        self.heat_combustion = self.heat_sludge + self.heat_fuel
        self.heat_air = float((air.mass*i_mass).sum()) * self.cp_air * (air.T-self.T_ref)

        # Evaporation of the water in the sludge (heating to the boiling point and vaporization)
        water = cmps[self.water_ID]
        T, Tb = sludge.T, water.Tb
        Cn = water.Cn.l if hasattr(water.Cn, 'l') else water.Cn # liquid heat capacity [J/mol/K]
        heat_up = Cn((T+Tb)/2) * (Tb-T) if T < Tb else 0.
        self.heat_evaporation = float(sludge.mass[i_h2o]*to_mol[i_h2o]) * (water.Hvap(Tb)+heat_up)

        heat_in = self.heat_combustion + self.heat_air
        self.heat_flue_gas = self.process_efficiency * heat_in
        self.heat_loss = heat_in - self.heat_flue_gas

    def _init_state(self):
        self._state = np.zeros(5) # ash, CO2, N2, H2O, O2 [mass/d]
        self._dstate = self._state * 0.
        self._cached_state = self._state.copy()
        self._cached_t = 0

    def _update_state(self):
        cmps = self.components
        flue_gas, ash = self.outs
        idx_ash = cmps.index(self.ash_component_ID)
        idx_gases = cmps.indices([self.carbon_dioxide_ID, self.nitrogen_ID,
                                  self.water_ID, self.oxygen_ID])

        # Float arrays, an integer array would truncate the flows
        if flue_gas.state is None: flue_gas.state = np.append(np.zeros(len(cmps)), 1.)
        if ash.state is None: ash.state = np.append(np.zeros(len(cmps)), 1.)

        flue_gas.state[idx_gases] = self._state[1:]
        ash.state[idx_ash] = self._state[0]

    def _update_dstate(self):
        cmps = self.components
        flue_gas, ash = self.outs
        idx_ash = cmps.index(self.ash_component_ID)
        idx_gases = cmps.indices([self.carbon_dioxide_ID, self.nitrogen_ID,
                                  self.water_ID, self.oxygen_ID])

        if flue_gas.dstate is None: flue_gas.dstate = np.zeros(len(cmps)+1)
        if ash.dstate is None: ash.dstate = np.zeros(len(cmps)+1)

        flue_gas.dstate[idx_gases] = self._dstate[1:]
        ash.dstate[idx_ash] = self._dstate[0]


    @property
    def AE(self):
        if self._AE is None:
            self._compile_AE()
        return self._AE

    def _compile_AE(self):
        _state = self._state
        _dstate = self._dstate
        _update_state = self._update_state
        _update_dstate = self._update_dstate
        _cached_state = self._cached_state
        _get_balance = self._get_balance

        def yt(t, QC_ins, dQC_ins):
            Mass_ins = np.diag(QC_ins[:,-1]) @ QC_ins[:,:-1] # one row per influent
            _state[:] = _get_balance(Mass_ins.sum(axis=0))
            _state[-1] = max(_state[-1], 0.) # no negative oxygen if the air is not enough

            if t > self._cached_t:
                _dstate[:] = (_state - _cached_state)/(t - self._cached_t)
            _cached_state[:] = _state
            self._cached_t = t
            _update_state()
            _update_dstate()
        self._AE = yt
