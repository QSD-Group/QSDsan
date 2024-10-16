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

from .. import SanUnit, WasteStream
import numpy as np
from ..sanunits import WWTpump
from warnings import warn
from ..sanunits._pumping import (
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
    >>> from qsdsan.sanunits import Thickener
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
         Alkalinity : 2.5 mg/L
         COD        : 23873.0 mg/L
         BOD        : 14963.2 mg/L
         TC         : 8298.3 mg/L
         TOC        : 8298.3 mg/L
         TN         : 20363.2 mg/L
         TP         : 367.6 mg/L
         TK         : 68.3 mg/L
         TSS        : 11124.4 mg/L
    outs...
    [0] sludge
    phase: 'l', T: 298.15 K, P: 101325 Pa
    flow (g/hr): S_F    1.56e+03
                    S_NH4  3.11e+03
                    X_OHO  1.47e+04
                    H2O    1.56e+05
        WasteStream-specific properties:
         pH         : 7.0
         COD        : 95050.4 mg/L
         BOD        : 55228.4 mg/L
         TC         : 34369.6 mg/L
         TOC        : 34369.6 mg/L
         TN         : 24354.4 mg/L
         TP         : 1724.0 mg/L
         TK         : 409.8 mg/L
         TSS        : 66748.0 mg/L
    [1] effluent
    phase: 'l', T: 298.15 K, P: 101325 Pa
    flow (g/hr): S_F    8.44e+03
                    S_NH4  1.69e+04
                    X_OHO  300
                    H2O    8.44e+05
        WasteStream-specific properties:
         pH         : 7.0
         COD        : 9978.2 mg/L
         BOD        : 7102.9 mg/L
         TC         : 3208.8 mg/L
         TOC        : 3208.8 mg/L
         TN         : 19584.1 mg/L
         TP         : 102.9 mg/L
         TK         : 1.6 mg/L
         TSS        : 265.9 mg/L

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
    calorific_value_sludge : float 
        The calorific value of influent sludge in KJ/kg. The default value used is 12000 KJ/kg.
    calorific_value_fuel : float 
        The calorific value of fuel employed for combustion in KJ/kg. 
        The default fuel is natural gas with calorific value of 50000 KJ/kg.
        
    Examples
    --------    
    >>> import qsdsan as qs
    >>> cmps = qs.Components.load_default()
    >>> CO2 = qs.Component.from_chemical('S_CO2', search_ID='CO2', 
    ...                                  particle_size='Soluble', 
    ...                                  degradability='Undegradable', 
    ...                                  organic=False)
    >>> cmps_test = qs.Components([cmps.S_F, cmps.S_NH4, cmps.X_OHO, cmps.H2O, 
    ...                            cmps.S_CH4, cmps.S_O2, cmps.S_N2, cmps.S_H2, 
    ...                            cmps.X_Ig_ISS, CO2])
    >>> cmps_test.default_compile()
    >>> qs.set_thermo(cmps_test)
    >>> ws = qs.WasteStream('ws', S_F=10, S_NH4=20, X_OHO=15, H2O=1000)
    >>> natural_gas = qs.WasteStream('nat_gas', phase='g', S_CH4=1000)
    >>> air = qs.WasteStream('air', phase='g', S_O2=210, S_N2=780, S_H2=10)
    >>> from qsdsan.sanunits import Incinerator
    >>> Inc = Incinerator(ID='Inc', ins= (ws, air, natural_gas), 
    ...                   outs=('flu_gas', 'ash'), 
    ...                   isdynamic=False)
    >>> Inc.simulate()
    >>> Inc.show()
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
         Alkalinity : 2.5 mg/L
         COD        : 23873.0 mg/L
         BOD        : 14963.2 mg/L
         TC         : 8298.3 mg/L
         TOC        : 8298.3 mg/L
         TN         : 20363.2 mg/L
         TP         : 367.6 mg/L
         TK         : 68.3 mg/L
         TSS        : 11124.4 mg/L
    [1] air
    phase: 'g', T: 298.15 K, P: 101325 Pa
    flow (g/hr): S_O2  2.1e+05
                    S_N2  7.8e+05
                    S_H2  1e+04
        WasteStream-specific properties: None for non-liquid waste streams
    [2] nat_gas
    phase: 'g', T: 298.15 K, P: 101325 Pa
    flow (g/hr): S_CH4  1e+06
        WasteStream-specific properties: None for non-liquid waste streams
    outs...
    [0] flu_gas
    phase: 'g', T: 298.15 K, P: 101325 Pa
    flow (g/hr): H2O    1e+06
                    S_N2   7.8e+05
                    S_CO2  2.67e+05
        WasteStream-specific properties: None for non-liquid waste streams
    [1] ash
    phase: 's', T: 298.15 K, P: 101325 Pa
    flow (g/hr): X_Ig_ISS  2.37e+05
        WasteStream-specific properties: None for non-liquid waste streams
    
    References
    ----------
    [1] Khuriati, A., P. Purwanto, H. S. Huboyo, Suryono Sumariyah, S. Suryono, and A. B. Putranto. 
    "Numerical calculation based on mass and energy balance of waste incineration in the fixed bed reactor." 
    In Journal of Physics: Conference Series, vol. 1524, no. 1, p. 012002. IOP Publishing, 2020.
    
    [2] Omari, Arthur, Karoli N. Njau, Geoffrey R. John, Joseph H. Kihedu, and Peter L. Mtui. 
    "Mass And Energy Balance For Fixed Bed Incinerators A case of a locally designed incinerator in Tanzania."
    """

    #These are class attributes
    _N_ins = 3
    _N_outs = 2   
    Cp_air = 1 #(Cp = 1 kJ/kg for air)
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, isdynamic=False, 
                  init_with='WasteStream', F_BM_default=None, process_efficiency=0.90, 
                  calorific_value_sludge=12000, calorific_value_fuel=50000, 
                  ash_component_ID='X_Ig_ISS', nitrogen_ID='S_N2', water_ID='H2O',
                  carbon_dioxide_ID='S_CO2', **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo, isdynamic=isdynamic, 
                         init_with=init_with, F_BM_default=F_BM_default)
        
        self.calorific_value_sludge = calorific_value_sludge #in KJ/kg
        self.calorific_value_fuel  = calorific_value_fuel #in KJ/kg (here the considered fuel is natural gas) 
        self.process_efficiency = process_efficiency
        self.ash_component_ID = ash_component_ID
        self.nitrogen_ID = nitrogen_ID
        self.water_ID = water_ID
        self.carbon_dioxide_ID = carbon_dioxide_ID
        self.Heat_air = None 
        self.Heat_fuel = None
        self.Heat_sludge = None
        self.Heat_flue_gas = None
        self.Heat_loss = None
    
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
        '''Calorific value of sludge in KJ/kg.'''
        return self._calorific_value_sludge

    @calorific_value_sludge.setter
    def calorific_value_sludge(self, calorific_value_sludge):
        if calorific_value_sludge is not None: 
            self._calorific_value_sludge = calorific_value_sludge
        else: 
            raise ValueError('Calorific value of sludge expected from user')
            
    @property
    def calorific_value_fuel(self):
        '''Calorific value of fuel in KJ/kg.'''
        return self._calorific_value_fuel

    @calorific_value_fuel.setter
    def calorific_value_fuel(self, calorific_value_fuel):
        if calorific_value_fuel is not None: 
            self._calorific_value_fuel = calorific_value_fuel
        else: 
            raise ValueError('Calorific value of fuel expected from user')
        
    def _run(self):
        
        sludge, air, fuel = self.ins
        flue_gas, ash = self.outs
        flue_gas.phase = 'g'
        ash.phase = 's'
        cmps = self.components
        nitrogen_ID = self.nitrogen_ID
        water_ID = self.water_ID
        carbon_dioxide_ID = self.carbon_dioxide_ID
        ash_cmp_ID = self.ash_component_ID
        
        if sludge.phase != 'l':
            raise ValueError(f'The phase of incoming sludge is expected to be liquid not {sludge.phase}')
        if air.phase != 'g':
            raise ValueError(f'The phase of air is expected to be gas not {air.phase}')
        if fuel.phase != 'g':
            raise ValueError(f'The phase of fuel is expected to be gas not {fuel.phase}')
        
        inf = np.asarray(sludge.mass + air.mass + fuel.mass)
        idx_n2 = cmps.index(nitrogen_ID)
        idx_h2o = cmps.index(water_ID)
        idx_co2 = cmps.index(carbon_dioxide_ID)
        idx_ash = cmps.index(ash_cmp_ID)
        i_mass = cmps.i_mass
        i_iss = i_mass*(1-cmps.f_Vmass_Totmass)
        
        n2 = inf[idx_n2]
        h2o = inf[idx_h2o]
        
        mass_ash = np.sum(inf*i_iss) - n2*i_mass[idx_n2] - h2o*i_mass[idx_h2o]

        # Conservation of mass 
        mass_flue_gas = np.sum(inf*i_mass) - mass_ash
        mass_co2 = mass_flue_gas - n2*i_mass[idx_n2] - h2o*i_mass[idx_h2o]
        
        flue_gas.set_flow([n2, h2o, (mass_co2/i_mass[idx_co2])], 
                          'kg/hr', (nitrogen_ID, water_ID, carbon_dioxide_ID))
        ash.set_flow(mass_ash/i_mass[idx_ash],'kg/hr', (ash_cmp_ID,))
        
        # Energy balance 
        # self.Heat_sludge = sludge.dry_mass*sludge.F_vol*self.calorific_value_sludge/1000 #in KJ/hr (mg/L)*(m3/hr)*(KJ/kg)=KJ/hr*(1/1000)
        # self.Heat_air = np.sum(air.mass*cmps.i_mass)*self.Cp_air #in KJ/hr 
        # self.Heat_fuel = np.sum(fuel.mass*cmps.i_mass)*self.calorific_value_fuel #in KJ/hr 
        # self.Heat_flue_gas = self.process_efficiency*(self.Heat_sludge + self.Heat_air + self.Heat_fuel)
        
        # # Conservation of energy
        # self.Heat_loss = self.Heat_sludge + self.Heat_air + self.Heat_fuel - self.Heat_flue_gas
        
    def _init_state(self):
        self._state = np.zeros(4)
        self._dstate = self._state * 0.
        self._cached_state = self._state.copy()
        self._cached_t = 0

    def _update_state(self):
        cmps = self.components
        flue_gas, ash = self.outs
        idx_ash = cmps.index(self.ash_component_ID)
        idx_gases = cmps.indices([self.carbon_dioxide_ID, self.nitrogen_ID, self.water_ID])

        if flue_gas.state is None: flue_gas.state = np.array([0]*len(cmps)+[1])
        if ash.state is None: ash.state = np.array([0]*len(cmps)+[1])
        
        flue_gas.state[idx_gases] = self._state[1:]
        ash.state[idx_ash] = self._state[0]
        
    def _update_dstate(self):
        cmps = self.components
        flue_gas, ash = self.outs
        idx_ash = cmps.index(self.ash_component_ID)
        idx_gases = cmps.indices([self.carbon_dioxide_ID, self.nitrogen_ID, self.water_ID])

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
        cmps = self.components
        idx_h2o = cmps.index(self.water_ID)
        idx_n2 = cmps.index(self.nitrogen_ID)
        idx_co2 = cmps.index(self.carbon_dioxide_ID)
        idx_ash = cmps.index(self.ash_component_ID)
        i_mass = cmps.i_mass
        i_iss = i_mass*(1-cmps.f_Vmass_Totmass)
        
        def yt(t, QC_ins, dQC_ins):            
            Mass_ins = np.diag(QC_ins[:,-1]) @ QC_ins[:,:-1]
            n2 = Mass_ins[idx_n2]
            h2o = Mass_ins[idx_h2o]
            ash = np.sum(Mass_ins*i_iss) - n2*i_mass[idx_n2] - h2o*i_mass[idx_h2o]
            co2 = np.sum(Mass_ins*i_mass) - ash - n2*i_mass[idx_n2] - h2o*i_mass[idx_h2o]
            
            _state[:] = [ash/i_mass[idx_ash], co2/i_mass[idx_co2], n2, h2o]
            
            if t > self._cached_t:
                _dstate[:] = (_state - _cached_state)/(t - self._cached_t)
            _cached_state[:] = _state
            self._cached_t = t
            _update_state()
            _update_dstate()
        self._AE = yt