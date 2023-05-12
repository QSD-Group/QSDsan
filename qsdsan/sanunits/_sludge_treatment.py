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


__all__ = ('Thickener', 'DewateringUnit', 'Incinerator')

class Thickener(SanUnit):
    
    """
    Thickener based on BSM2 Layout. [1]

    ----------
    ID : str
        ID for the Thickener. The default is ''.
    ins : class:`WasteStream`
        Influent to the clarifier. Expected number of influent is 1. 
    outs : class:`WasteStream`
        Treated effluent and sludge.
    thickener_perc : float
        The percentage of solids in the underflow of the thickener.[1]
    TSS_removal_perc : float
        The percentage of suspended solids removed in the thickener.[1]
    solids_loading_rate : float
        Solid loading rate in the thickener.[2]
    h_cylinderical = float
        Height of cylinder forming the thickener.[2]
        
        
    Examples
    --------
    
    >>> from qsdsan import set_thermo, Components, WasteStream
    >>> cmps = Components.load_default()
    >>> cmps_test = cmps.subgroup(['S_F', 'S_NH4', 'X_OHO', 'H2O'])
    >>> set_thermo(cmps_test)
    >>> ws = WasteStream('ws', S_F = 10, S_NH4 = 20, X_OHO = 15, H2O=1000)
    >>> from qsdsan.sanunits import Thickener
    >>> ps = Thickener(ID='TC', ins= (ws), outs=('Sludge', 'Effluent'))
    >>> ps._run()
    >>> uf, of = ps.outs
    >>> uf.imass['X_OHO']/ws.imass['X_OHO'] # doctest: +ELLIPSIS
    0.98
    >>> ps
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
         COD        : 23643.1 mg/L
         BOD        : 14819.1 mg/L
         TC         : 8218.3 mg/L
         TOC        : 8218.3 mg/L
         TN         : 20167.1 mg/L
         TP         : 364.1 mg/L
         TK         : 67.6 mg/L
    outs...
    [0] Sludge
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (g/hr): S_F    8.46e+03
                     S_NH4  1.69e+04
                     X_OHO  1.47e+04
                     H2O    8.46e+05
        WasteStream-specific properties:
         pH         : 7.0
         COD        : 25857.3 mg/L
         BOD        : 16071.7 mg/L
         TC         : 9029.4 mg/L
         TOC        : 9029.4 mg/L
         TN         : 20291.5 mg/L
         TP         : 406.3 mg/L
         TK         : 78.3 mg/L
    [1] Effluent
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (g/hr): S_F    1.54e+03
                     S_NH4  3.08e+03
                     X_OHO  300
                     H2O    1.54e+05
        WasteStream-specific properties:
         pH         : 7.0
         COD        : 11387.0 mg/L
         BOD        : 7885.7 mg/L
         TC         : 3729.1 mg/L
         TOC        : 3729.1 mg/L
         TN         : 19478.3 mg/L
         TP         : 130.6 mg/L
         TK         : 8.8 mg/L

    References
    ----------
    .. [1] Gernaey, Krist V., Ulf Jeppsson, Peter A. Vanrolleghem, and John B. Copp.
    Benchmarking of control strategies for wastewater treatment plants. IWA publishing, 2014.
    [2] Metcalf, Leonard, Harrison P. Eddy, and Georg Tchobanoglous. Wastewater 
    engineering: treatment, disposal, and reuse. Vol. 4. New York: McGraw-Hill, 1991.
    """
    
    _N_ins = 1
    _N_outs = 2
    _ins_size_is_fixed = False
    _outs_size_is_fixed = False
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, isdynamic=False, 
                  init_with='WasteStream', F_BM_default=None, thickener_perc=7, 
                  TSS_removal_perc=98, solids_loading_rate = 50, h_cylinderical=2, **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo, isdynamic=isdynamic, 
                         init_with=init_with, F_BM_default=F_BM_default)
        
        self.thickener_perc = thickener_perc 
        self.TSS_removal_perc = TSS_removal_perc
        self.solids_loading_rate = solids_loading_rate 
        self.h_cylinderical = h_cylinderical
        self.mixed = WasteStream(thermo=thermo)
        
    @property
    def thickener_perc(self):
        '''tp is the percentage of Suspended Sludge in the underflow of the thickener'''
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
        self.mixed.mix_from(self.ins)
        inf = self.mixed
        _cal_thickener_factor = self._cal_thickener_factor
        if not self.ins: return
        elif inf.isempty(): return
        else: 
            TSS_in = inf.get_TSS()
            thickener_factor = _cal_thickener_factor(TSS_in)
        return thickener_factor
    
    @property
    def thinning_factor(self):
        self.mixed.mix_from(self.ins)
        inf = self.mixed
        TSS_in = inf.get_TSS()
        _cal_thickener_factor = self._cal_thickener_factor
        thickener_factor = _cal_thickener_factor(TSS_in)
        _cal_parameters = self._cal_parameters
        Qu_factor, thinning_factor = _cal_parameters(thickener_factor)
        return thinning_factor
    
    def _cal_thickener_factor(self, TSS_in):
        if TSS_in > 0:
            thickener_factor = self._tp*10000/TSS_in
            if thickener_factor<1:
                thickener_factor=1
            return thickener_factor
        else: return None
            
    def _cal_parameters(self, thickener_factor):
        if thickener_factor<1:
            Qu_factor = 1
            thinning_factor=0
        else:
            Qu_factor = self._TSS_rmv/(100*thickener_factor)
            thinning_factor = (1 - (self._TSS_rmv/100))/(1 - Qu_factor)
        return Qu_factor, thinning_factor
    
    def _update_parameters(self):
        
        # Thickener_factor, Thinning_factor, and Qu_factor need to be 
        # updated again and again. while dynamic simulations 
        
        cmps = self.components 
    
        TSS_in = np.sum(self._state[:-1]*cmps.i_mass*cmps.x)
        _cal_thickener_factor = self._cal_thickener_factor
        self.updated_thickener_factor = _cal_thickener_factor(TSS_in)
        _cal_parameters = self._cal_parameters
        
        updated_thickener_factor = self.updated_thickener_factor
        self.updated_Qu_factor, self.updated_thinning_factor = _cal_parameters(updated_thickener_factor)
        
    def _run(self):
        
        self.mixed.mix_from(self.ins)
        inf = self.mixed
        uf, of = self.outs
        cmps = self.components
        
        TSS_rmv = self._TSS_rmv
        thinning_factor = self.thinning_factor
        thickener_factor = self.thickener_factor
        
        # The following are splits by mass of particulates and solubles 
        
        # Note: (1 - thinning_factor)/(thickener_factor - thinning_factor) = Qu_factor
        Zs = (1 - thinning_factor)/(thickener_factor - thinning_factor)*inf.mass*cmps.s
        Ze = (thickener_factor - 1)/(thickener_factor - thinning_factor)*inf.mass*cmps.s
        
        Xe = (1 - TSS_rmv/100)*inf.mass*cmps.x
        Xs = (TSS_rmv/100)*inf.mass*cmps.x
        
        # e stands for effluent, s stands for sludge 
        Ce = Ze + Xe 
        Cs = Zs + Xs
    
        of.set_flow(Ce,'kg/hr')
        uf.set_flow(Cs,'kg/hr')
        
    def _design(self):
        
        design = self.design_results
        slr = self._slr
                
        design['Area'] = ((self.ins[0].get_TSS()/1000)*self.ins[0].F_vol*24)/slr # in m2
        design['Hydraulic_Loading'] = (self.ins[0].F_vol*24)/design['Area'] #in m3/(m2*day)
        design['Diameter'] = np.sqrt(4*design['Area']/np.pi) #in m
        design['Volume'] = np.pi*np.square(design['Diameter']/2)*self.h_cylinderical #in m3
        design['Curved Surface Area'] = np.pi*design['Diameter']*self.h_cylinderical #in m2
       
    def _init_state(self):
       
        # This function is run only once during dynamic simulations 
    
        # Since there could be multiple influents, the state of the unit is 
        # obtained assuming perfect mixing 
        Qs = self._ins_QC[:,-1]
        Cs = self._ins_QC[:,:-1]
        self._state = np.append(Qs @ Cs / Qs.sum(), Qs.sum())
        self._dstate = self._state * 0.
        
        # To initialize the updated_thickener_factor, updated_thinning_factor
        # and updated_Qu_factor for dynamic simulation 
        self._update_parameters()
        
    def _update_state(self):
        '''updates conditions of output stream based on conditions of the Thickener''' 
        
        # This function is run multiple times during dynamic simulation 
        
        # Remember that here we are updating the state array of size n, which is made up 
        # of component concentrations in the first (n-1) cells and the last cell is flowrate. 
        
        # So, while in the run function the effluent and sludge are split by mass, 
        # here they are split by concentration. Therefore, the split factors are different. 
        
        # Updated intrinsic modelling parameters are used for dynamic simulation 
        thickener_factor = self.updated_thickener_factor
        thinning_factor = self.updated_thinning_factor
        Qu_factor = self.updated_Qu_factor
        cmps = self.components
        
        # For sludge, the particulate concentrations are multipled by thickener factor, and
        # flowrate is multiplied by Qu_factor. The soluble concentrations remains same. 
        uf, of = self.outs
        if uf.state is None: uf.state = np.zeros(len(cmps)+1)
        uf.state[:-1] = self._state[:-1]*cmps.s*1 + self._state[:-1]*cmps.x*thickener_factor
        uf.state[-1] = self._state[-1]*Qu_factor
        
        # For effluent, the particulate concentrations are multipled by thinning factor, and
        # flowrate is multiplied by Qu_factor. The soluble concentrations remains same. 
        if of.state is None: of.state = np.zeros(len(cmps)+1)
        of.state[:-1] = self._state[:-1]*cmps.s*1 + self._state[:-1]*cmps.x*thinning_factor
        of.state[-1] = self._state[-1]*(1 - Qu_factor)

    def _update_dstate(self):
        '''updates rates of change of output stream from rates of change of the Thickener'''
        
        # This function is run multiple times during dynamic simulation 
        
        # Remember that here we are updating the state array of size n, which is made up 
        # of component concentrations in the first (n-1) cells and the last cell is flowrate. 
        
        # So, while in the run function the effluent and sludge are split by mass, 
        # here they are split by concentration. Therefore, the split factors are different. 
        
        # Updated intrinsic modelling parameters are used for dynamic simulation
        thickener_factor = self.updated_thickener_factor
        thinning_factor = self.updated_thinning_factor
        Qu_factor = self.updated_Qu_factor
        cmps = self.components
        
        # For sludge, the particulate concentrations are multipled by thickener factor, and
        # flowrate is multiplied by Qu_factor. The soluble concentrations remains same. 
        uf, of = self.outs
        if uf.dstate is None: uf.dstate = np.zeros(len(cmps)+1)
        uf.dstate[:-1] = self._dstate[:-1]*cmps.s*1 + self._dstate[:-1]*cmps.x*thickener_factor
        uf.dstate[-1] = self._dstate[-1]*Qu_factor
        
        # For effluent, the particulate concentrations are multipled by thinning factor, and
        # flowrate is multiplied by Qu_factor. The soluble concentrations remains same.
        if of.dstate is None: of.dstate = np.zeros(len(cmps)+1)
        of.dstate[:-1] = self._dstate[:-1]*cmps.s*1 + self._dstate[:-1]*cmps.x*thinning_factor
        of.dstate[-1] = self._dstate[-1]*(1 - Qu_factor)
     
    @property
    def AE(self):
        if self._AE is None:
            self._compile_AE()
        return self._AE

    def _compile_AE(self):
        
        # This function is run multiple times during dynamic simulation 
        
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

class DewateringUnit(Thickener):
    
    """
    Dewatering Unit based on BSM2 Layout. [1]
    
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
    number_of_centrifuges : float
        Number of centrifuges in the dewatering unit.[2,3]
    specific_gravity_sludge: float
        Specific gravity of influent sludge from secondary clarifier.[2,3]
    cake density: float
        Density of effleunt dewatered sludge.[2,3]
    centrifugal_force : float
        Centrifugal force in the centrifuge.[2,3]
    rotational_speed : float
        rotational speed of the centrifuge.[2,3]
    polymer_dosage_per_kg_of_sludge : float
        mass of polymer utilised per kg of influent sludge.[2,3]
    h_cylinderical: float
        length of cylinderical portion of dewatering unit.[2,3]
    h_conical: float
        length of conical portion of dewatering unit.[2,3]
    
    References
    ----------
    .. [1] Gernaey, Krist V., Ulf Jeppsson, Peter A. Vanrolleghem, and John B. Copp.
    Benchmarking of control strategies for wastewater treatment plants. IWA publishing, 2014.
    [2] Metcalf, Leonard, Harrison P. Eddy, and Georg Tchobanoglous. Wastewater 
    engineering: treatment, disposal, and reuse. Vol. 4. New York: McGraw-Hill, 1991.
    [3]Design of Municipal Wastewater Treatment Plants: WEF Manual of Practice 
    No. 8 ASCE Manuals and Reports on Engineering Practice No. 76, Fifth Edition. 
    """
    
    _N_ins = 1
    _N_outs = 2
    _ins_size_is_fixed = False
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, isdynamic=False, 
                  init_with='WasteStream', F_BM_default=None, thickener_perc=28, TSS_removal_perc=98, 
                  number_of_centrifuges=1, specific_gravity_sludge=1.03, cake_density=965, 
                  g_factor=2500, rotational_speed = 40, polymer_dosage_per_kg_of_sludge = 0.0075, 
                  h_cylinderical=2, h_conical=1, **kwargs):
        Thickener.__init__(self, ID=ID, ins=ins, outs=outs, thermo=thermo, isdynamic=isdynamic,
                      init_with=init_with, F_BM_default=F_BM_default, thickener_perc=thickener_perc, 
                      TSS_removal_perc=TSS_removal_perc, **kwargs)
        self.number_of_centrifuges=number_of_centrifuges
        self.specific_gravity_sludge=specific_gravity_sludge
        self.cake_density=cake_density #in kg/m3
        self.g_factor = g_factor #unitless, centrifugal acceleration = g_factor*9.81
        self.rotational_speed = rotational_speed #in revolution/sec 
        self.polymer_dosage_per_kg_of_sludge = polymer_dosage_per_kg_of_sludge #in (kg,polymer/kg,sludge) unitless
        self.h_cylinderical = h_cylinderical
        self.h_conical = h_conical
        
    def _design(self):
        sludge_feed_rate = ((self.ins[0].get_TSS()*self.ins[0].F_vol)/1000)/self.number_of_centrifuges #in kg/hr
        
        #TSS_rmv = self._TSS_rmv
        #recovery = 1 - TSS_rmv/100
        #cake_mass_discharge_rate = sludge_feed_rate*recovery #in kg/hr
        #wetcake_mass_discharge_rate = cake_mass_discharge_rate/(self.thickener_perc/100) #in kg/hr
        #cake_density = self.cake_density 
        #wetcake_flowrate = wetcake_mass_discharge_rate/cake_density #in m3/hr
        #volume_reduction_perc= (1 - wetcake_flowrate/(self.ins[0].F_mass/(1000*self.specific_gravity_sludge*self.number_of_centrifuges)))*100
        
        design = self.design_results 
        design['Diameter'] = 2*(self.g_factor*9.81/np.square(2*np.pi*self.rotational_speed)) #in m
        design['Polymer feed rate'] = (self.polymer_dosage_per_kg_of_sludge*sludge_feed_rate) # in kg/hr
        design['Projected Area at Inlet'] = np.pi*np.square(design['Diameter']/2) #in m2
        design['Hydraulic_Loading'] = (self.ins[0].F_vol*24)/design['Projected Area at Inlet'] #in m3/(m2*day)
        design['Volume'] = np.pi*np.square(design['Diameter']/2)*(self.h_cylinderical + (self.h_conical/3)) #in m3
        design['Curved Surface Area'] = np.pi*design['Diameter']/2*(2*self.h_cylinderical + np.sqrt(np.square(design['Diameter']/2) + np.square(self.h_conical))) #in m2
        
class Incinerator(SanUnit):
    
    """
    Fluidized bed incinerator unit for metroWWTP.
    
    Parameters
    ----------
    ID : str
        ID for the Incinerator Unit. The default is ''.
    ins : class:`WasteStream`
        Influent to the Incinerator Unit. Expected number of influent streams are 3. 
        Please remember the order of influents as {wastestream, air, fuel} 
    outs : class:`WasteStream`
        Flue gas and ash. 
    thickener_perc : float
        The percentage of Suspended Sludge in the underflow of the dewatering unit.
    process_efficiency : float
        The process efficiency of the incinerator unit. Expected value between 0 and 1. 
    calorific_value_sludge : float 
        The calorific value of influent sludge in KJ/kg. The default value used is 12000 KJ/kg.
    calorific_value_fuel : float 
        The calorific value of fuel employed for combustion in KJ/kg. 
        The default fuel is natural gas with calofific value of 50000 KJ/kg.
        
    Examples
    --------
    
    >>> import qsdsan as qs
    >>> cmps = qs.Components.load_default()
    >>> CO2 = qs.Component.from_chemical('S_CO2', search_ID='CO2', particle_size='Soluble', degradability='Undegradable', organic=False)
    >>> cmps_test = qs.Components([cmps.S_F, cmps.S_NH4, cmps.X_OHO, cmps.H2O, cmps.S_CH4, cmps.S_O2, cmps.S_N2, cmps.S_H2, cmps.X_Ig_ISS, CO2])
    >>> cmps_test.default_compile()
    >>> qs.set_thermo(cmps_test)
    >>> ws = qs.WasteStream('ws', S_F=10, S_NH4=20, X_OHO=15, H2O=1000)
    >>> natural_gas = qs.WasteStream('nat_gas', phase='g', S_CH4=1000)
    >>> air = qs.WasteStream('air', phase='g', S_O2=210, S_N2=780, S_H2=10)
    >>> from qsdsan.sanunits import Incinerator
    >>> ps = Incinerator(ID='PC', ins= (ws, air, natural_gas), outs=('flu_gas', 'ash'), 
                     isdynamic=True)
    >>> ps._run()
    >>> ps
       
    Incinerator: PC
    ins...
    [0] ws
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (g/hr): S_F    1e+04
                     S_NH4  2e+04
                     X_OHO  1.5e+04
                     H2O    1e+06
        WasteStream-specific properties:
         pH         : 7.0
         COD        : 23643.1 mg/L
         BOD        : 14819.1 mg/L
         TC         : 8218.3 mg/L
         TOC        : 8218.3 mg/L
         TN         : 20167.1 mg/L
         TP         : 364.1 mg/L
         TK         : 67.6 mg/L
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
                     S_CO2  4.77e+05
        WasteStream-specific properties: None for non-liquid waste streams
    [1] ash
        phase: 's', T: 298.15 K, P: 101325 Pa
        flow (g/hr): X_Ig_ISS  2.58e+04
        WasteStream-specific properties: None for non-liquid waste streams
    
    References: 
    ----------
    .. [1] Khuriati, A., P. Purwanto, H. S. Huboyo, Suryono Sumariyah, S. Suryono, and A. B. Putranto. 
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
                  calorific_value_sludge= 12000, calorific_value_fuel=50000, 
                  ash_component_ID = 'X_Ig_ISS', nitrogen_ID = 'S_N2', water_ID = 'H2O',
                  carbon_di_oxide_ID = 'S_CO2', **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo, isdynamic=isdynamic, 
                         init_with=init_with, F_BM_default=F_BM_default)
        
        self.calorific_value_sludge = calorific_value_sludge #in KJ/kg
        self.calorific_value_fuel  = calorific_value_fuel #in KJ/kg (here the considered fuel is natural gas) 
        self.process_efficiency = process_efficiency
        self.ash_component_ID = ash_component_ID
        self.nitrogen_ID = nitrogen_ID
        self.water_ID = water_ID
        self.carbon_di_oxide_ID = carbon_di_oxide_ID
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
        carbon_di_oxide_ID = self.carbon_di_oxide_ID
        
        if sludge.phase != 'l':
            raise ValueError(f'The phase of incoming sludge is expected to be liquid not {sludge.phase}')
        if air.phase != 'g':
            raise ValueError(f'The phase of air is expected to be gas not {air.phase}')
        if fuel.phase != 'g':
            raise ValueError(f'The phase of fuel is expected to be gas not {fuel.phase}')
        
        inf = sludge.mass + air.mass + fuel.mass
        idx_n2 = cmps.index(nitrogen_ID)
        idx_h2o = cmps.index(water_ID)
        
        n2 = inf[idx_n2]
        h2o = inf[idx_h2o]
        
        mass_ash = sum(inf*cmps.i_mass*(1-cmps.f_Vmass_Totmass)) \
               - h2o*cmps.H2O.i_mass*(1-cmps.H2O.f_Vmass_Totmass) \
                   - n2*cmps.N2.i_mass*(1-cmps.N2.f_Vmass_Totmass)

        # Conservation of mass 
        mass_flue_gas = np.sum(inf*cmps.i_mass) - mass_ash
        mass_co2 = mass_flue_gas - n2*cmps.N2.i_mass - h2o*cmps.H2O.i_mass
        flue_gas.set_flow([n2, h2o, (mass_co2/cmps.S_CO2.i_mass)], 
                          'kg/hr', (nitrogen_ID, water_ID, carbon_di_oxide_ID))
        ash_cmp_ID = self.ash_component_ID
        ash_idx = cmps.index(ash_cmp_ID)
        ash.set_flow([mass_ash/cmps.i_mass[ash_idx]/(1-cmps.f_Vmass_Totmass[ash_idx])], 
                     'kg/hr', (ash_cmp_ID))
        
        # Energy balance 
        self.Heat_sludge = sludge.dry_mass*sludge.F_vol*self.calorific_value_sludge/1000 #in KJ/hr (mg/L)*(m3/hr)*(KJ/kg)=KJ/hr*(1/1000)
        self.Heat_air = np.sum(air.mass*cmps.i_mass)*self.Cp_air #in KJ/hr 
        self.Heat_fuel = np.sum(fuel.mass*cmps.i_mass)*self.calorific_value_fuel #in KJ/hr 
        self.Heat_flue_gas = self.process_efficiency*(self.Heat_sludge + self.Heat_air + self.Heat_fuel)
        
        # Conservation of energy
        self.Heat_loss = self.Heat_sludge + self.Heat_air + self.Heat_fuel - self.Heat_flue_gas
        
    def _init_state(self):
        
        sludge, air, fuel = self.ins
        inf = sludge.mass + air.mass + fuel.mass
        self._state = (24*inf)/1000
        self._dstate = self._state * 0.
        self._cached_state = self._state.copy()
        self._cached_t = 0
        
    def _update_state(self):
        cmps = self.components
        for ws in self.outs:
            if ws.state is None:
                ws.state = np.zeros(len(self._state)+1)
        
        nitrogen_ID = self.nitrogen_ID
        water_ID = self.water_ID
        carbon_di_oxide_ID = self.carbon_di_oxide_ID
        
        idx_h2o = cmps.index(water_ID)
        idx_n2 = cmps.index(nitrogen_ID)
        idx_co2 = cmps.index(carbon_di_oxide_ID)
        ash_idx = cmps.index(self.ash_component_ID)
        cmps_i_mass = cmps.i_mass
        cmps_v2tmass = cmps.f_Vmass_Totmass    
        
        inf = self._state
        mass_in_tot = np.sum(inf*cmps_i_mass)       
        self._outs[0].state[idx_h2o] = h2o = inf[idx_h2o]
        self._outs[0].state[idx_n2] = n2 = inf[idx_n2]
        
        mass_ash = np.sum(inf*cmps_i_mass*(1-cmps_v2tmass)) \
                   - h2o*cmps.H2O.i_mass*(1-cmps_v2tmass[idx_h2o]) - n2*cmps.N2.i_mass*(1-cmps_v2tmass[idx_n2])
        mass_flue_gas = mass_in_tot - mass_ash
        mass_co2 = mass_flue_gas - n2 - h2o
        
        self._outs[0].state[idx_co2] = mass_co2/cmps_i_mass[idx_co2]
        self._outs[1].state[ash_idx] = mass_ash/cmps_i_mass[ash_idx]/(1-cmps_v2tmass[ash_idx])

        self._outs[0].state[-1] = 1
        self._outs[1].state[-1] = 1
        
    def _update_dstate(self):
        
        cmps = self.components
        nitrogen_ID = self.nitrogen_ID
        water_ID = self.water_ID
        carbon_di_oxide_ID = self.carbon_di_oxide_ID
        
        idx_h2o = cmps.index(water_ID)
        idx_n2 = cmps.index(nitrogen_ID)
        idx_co2 = cmps.index(carbon_di_oxide_ID)
        ash_idx = cmps.index(self.ash_component_ID)
        d_state = self._dstate
        cmps_i_mass = cmps.i_mass
        cmps_v2tmass = cmps.f_Vmass_Totmass  
        d_n2 = d_state[idx_n2]
        d_h2o = d_state[idx_h2o]
        
        for ws in self.outs:
            if ws.dstate is None:
                ws.dstate = np.zeros(len(self._dstate)+1)
    
        self._outs[0].dstate[idx_n2] = d_n2
        self._outs[0].dstate[idx_h2o] = d_h2o
        
        d_mass_in_tot = np.sum(d_state*cmps_i_mass)
        
        d_mass_ash = np.sum(d_state*cmps_i_mass*(1-cmps_v2tmass)) \
            - d_h2o*cmps.H2O.i_mass*(1-cmps_v2tmass[idx_h2o]) - d_n2*cmps.N2.i_mass*(1-cmps_v2tmass[idx_n2])
        d_mass_flue_gas = d_mass_in_tot - d_mass_ash
        d_mass_co2 = d_mass_flue_gas - d_n2  - d_h2o
        
        self._outs[0].dstate[idx_co2] = d_mass_co2/cmps_i_mass[idx_co2]
        self._outs[1].dstate[ash_idx] = d_mass_ash/cmps_i_mass[ash_idx]/(1-cmps_v2tmass[ash_idx])
        
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

        def yt(t, QC_ins, dQC_ins):            
            # Mass_in is basically the mass flowrate array where each row 
            # corresponds to the flowrates of individual components (in columns)
            # This strcuture is achieved by multiplying the first (n-1) rows of 
            # Q_ins (which corresponds to concentration) to the nth row (which 
            # is the volumetric flowrate)
            Mass_ins = np.diag(QC_ins[:,-1]) @ QC_ins[:,:-1]
            # the _state array is formed by adding each column of the Mass_in
            # array, thus providing the total massflowrate of each component 
            _state[:] = np.sum(Mass_ins, axis=0)
            
            if t > self._cached_t:
                _dstate[:] = (_state - _cached_state)/(t - self._cached_t)
            _cached_state[:] = _state
            self._cached_t = t
            _update_state()
            _update_dstate()
        self._AE = yt