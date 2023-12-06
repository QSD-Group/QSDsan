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
         COD        : 23873.0 mg/L
         BOD        : 14963.2 mg/L
         TC         : 8298.3 mg/L
         TOC        : 8298.3 mg/L
         TN         : 20363.2 mg/L
         TP         : 367.6 mg/L
         TK         : 68.3 mg/L
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

    References
    ----------
    .. [1] Gernaey, Krist V., Ulf Jeppsson, Peter A. Vanrolleghem, and John B. Copp.
    Benchmarking of control strategies for wastewater treatment plants. IWA publishing, 2014.
    .. [2] Chapter-21: Solids Thicknening (Table 21.3). WEF Manual of Practice No. 8. 
    6th Edition. Virginia: McGraw-Hill, 2018. 
    .. [3] Introduction to Wastewater Clarifier Design by Nikolay Voutchkov, PE, BCEE.
    .. [4] Metcalf, Leonard, Harrison P. Eddy, and Georg Tchobanoglous. Wastewater 
    engineering: treatment, disposal, and reuse. Vol. 4. New York: McGraw-Hill, 1991.
    """
    
    _N_ins = 1
    _N_outs = 2
    _ins_size_is_fixed = False
    _outs_size_is_fixed = False
    
    # Costs
    wall_concrete_unit_cost = 1081.73 # $/m3 (Hydromantis. CapdetWorks 4.0. https://www.hydromantis.com/CapdetWorks.html)
    slab_concrete_unit_cost = 582.48 # $/m3 (Hydromantis. CapdetWorks 4.0. https://www.hydromantis.com/CapdetWorks.html)
    stainless_steel_unit_cost=1.8 # Alibaba. Brushed Stainless Steel Plate 304. https://www.alibaba.com/product-detail/brushed-stainless-steel-plate-304l-stainless_1600391656401.html?spm=a2700.details.0.0.230e67e6IKwwFd
    
    pumps = ('sludge',)
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, isdynamic=False, 
                  init_with='WasteStream', F_BM_default=default_F_BM, thickener_perc=7, 
                  TSS_removal_perc=98, solids_loading_rate =4, h_thickener=4, 
                  downward_flow_velocity= 36, F_BM=default_F_BM, **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo, isdynamic=isdynamic, 
                         init_with=init_with)
        self.thickener_perc = thickener_perc 
        self.TSS_removal_perc = TSS_removal_perc
        self.solids_loading_rate = solids_loading_rate 
        self.h_thickener = h_thickener
        self.downward_flow_velocity = downward_flow_velocity
        self.F_BM.update(F_BM)
        self._mixed = WasteStream(f'{ID}_mixed')        
        self._sludge = self.outs[0].copy(f'{ID}_sludge')
        self._thickener_factor = None
        self._thinning_factor = None
        self._Qu_factor = None
        
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
        if self._thickener_factor is None:
            self._mixed.mix_from(self.ins)
            inf = self._mixed
            _cal_thickener_factor = self._cal_thickener_factor
            if not self.ins: return
            elif inf.isempty(): return
            else: 
                TSS_in = inf.get_TSS()
                self._thickener_factor = _cal_thickener_factor(TSS_in)
        return self._thickener_factor
    
    @property
    def thinning_factor(self):
        if self._thinning_factor is None:
            self._Qu_factor, self._thinning_factor = self._cal_Qu_fthin(self.thickener_factor)
        return self._thinning_factor
    
    def _cal_thickener_factor(self, TSS_in):
        if TSS_in > 0:
            thickener_factor = self._tp*10000/TSS_in
            if thickener_factor < 1: thickener_factor = 1
            return thickener_factor
        else: 
            raise ValueError(f'Influent TSS is not valid: ({TSS_in:.2f} mg/L).')
            
    def _cal_Qu_fthin(self, thickener_factor):
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
        self._thickener_factor = f_thick = self._cal_thickener_factor(TSS_in)
        self._Qu_factor, self._thinning_factor = self._cal_Qu_fthin(f_thick)
        
    def _run(self):
        self._mixed.mix_from(self.ins)
        inf = self._mixed
        sludge, eff = self.outs
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
    
        eff.set_flow(Ce,'kg/hr')
        sludge.set_flow(Cs,'kg/hr')
       
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
        thickener_factor = self.thickener_factor
        thinning_factor = self.thinning_factor
        Qu_factor = self._Qu_factor
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
        thickener_factor = self.thickener_factor
        thinning_factor = self.thinning_factor
        Qu_factor = self._Qu_factor
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
         
    def _design_pump(self):
        ID, pumps = self.ID, self.pumps
        self._sludge.copy_like(self.outs[0])
        sludge = self._sludge
        
        ins_dct = {
            'sludge': sludge,
            }
        
        type_dct = dict.fromkeys(pumps, 'sludge')
        inputs_dct = dict.fromkeys(pumps, (1,))
        
        D = self.design_results
        influent_Q = sludge.get_total_flow('m3/hr')/D['Number of thickeners']
        influent_Q_mgd = influent_Q*0.00634 # m3/hr to MGD
       
        for i in pumps:
            if hasattr(self, f'{i}_pump'):
                p = getattr(self, f'{i}_pump')
                setattr(p, 'add_inputs', inputs_dct[i])
            else:
                ID = f'{ID}_{i}'
                capacity_factor=1
                # No. of pumps = No. of influents
                pump = WWTpump(
                    ID=ID, ins= ins_dct[i], pump_type=type_dct[i],
                    Q_mgd=influent_Q_mgd, add_inputs=inputs_dct[i],
                    capacity_factor=capacity_factor,
                    include_pump_cost=True,
                    include_building_cost=False,
                    include_OM_cost=True,
                    )
                setattr(self, f'{i}_pump', pump)

        pipe_ss, pump_ss = 0., 0.
        for i in pumps:
            p = getattr(self, f'{i}_pump')
            p.simulate()
            p_design = p.design_results
            pipe_ss += p_design['Pump pipe stainless steel']
            pump_ss += p_design['Pump stainless steel']
        return pipe_ss, pump_ss
    
    _units = {
        'Design solids loading rate': 'kg/m2/hr',
        'Total mass of solids handled': 'kg',
        'Surface area': 'm2',
        'Thickener diameter': 'm',
        'Number of thickeners' : 'Unitless',
        'Thickener depth': 'm',
        'Conical depth': 'm',
        'Cylindrical depth': 'm',
        'Cylindrical volume': 'm3',
        'Conical volume': 'm3',
        'Thickener volume': 'm3',
        
        'Center feed depth': 'm',
        'Downward flow velocity': 'm/hr',
        'Volumetric flow': 'm3/hr', 
        'Center feed diameter': 'm',
        'Thickener depth': 'm',
        'Volume of concrete wall': 'm3',
        'Volume of concrete slab': 'm3',
        'Stainless steel': 'kg',
        'Pump pipe stainless steel' : 'kg',
        'Pump stainless steel': 'kg',
        'Number of pumps': 'Unitless'
    }
    
    def _design(self):
        
        self._mixed.mix_from(self.ins)
        mixed = self._mixed
        D = self.design_results
        
        # D['Number of thickeners'] = np.ceil(self._mixed.get_total_flow('m3/hr')/self.design_flow)
        D['Design solids loading rate'] = self.solids_loading_rate # in (kg/hr)/m2
        D['Total mass of solids handled'] = (mixed.get_TSS()/1000)*mixed.get_total_flow('m3/hr') # (mg/L)*[1/1000(kg*L)/(mg*m3)](m3/hr) = (kg/hr)
        
        # Common gravity thickener configurations have tanks with diameter between 21-24m (MOP 8)
        diameter_thickener = 24
        number_of_thickeners = 0
        while diameter_thickener >= 22:
            number_of_thickeners += 1
            total_surface_area =  D['Total mass of solids handled']/D['Design solids loading rate'] #m2
            surface_area_thickener = total_surface_area/number_of_thickeners
            diameter_thickener = np.sqrt(4*surface_area_thickener/np.pi)
            
        D['Surface area'] = surface_area_thickener #in m2
        D['Thickener diameter'] = diameter_thickener #in m
        D['Number of thickeners'] = number_of_thickeners
    
        # Common gravity thickener configurations have sidewater depth between 3-4m (MOP 8)
        D['Thickener depth'] = self.h_thickener #in m 
        # The thickener tank floor generally has slope between 2:12 and 3:12 (MOP 8)
        D['Conical depth'] = (2/12)*(D['Thickener diameter']/2)
        D['Cylindrical depth'] = D['Thickener depth'] - D['Conical depth']
        
        # Checks on depth
        if D['Cylindrical depth'] < 0:
            cyl_depth = D['Cylindrical depth']
            RuntimeError(f'Cylindrical depth (= {cyl_depth} ) is negative')
            
        if D['Cylindrical depth'] < D['Conical depth']:
            cyl_depth = D['Cylindrical depth']
            con_depth = D['Conical depth'] 
            RuntimeError(f'Cylindrical depth (= {cyl_depth} ) is lower than Conical depth (= {con_depth})')
             
        D['Cylindrical volume'] = np.pi*np.square(D['Thickener diameter']/2)*D['Cylindrical depth'] #in m3
        D['Conical volume'] = (np.pi/3)*(D['Thickener diameter']/2)**2*D['Conical depth']
        D['Thickener volume'] = D['Cylindrical volume'] + D['Conical volume']
        
        #Check on SOR is pending
        
        # The design here is for center feed of thickener.
        # Depth of the center feed lies between 30-75% of sidewater depth. [2]
        D['Center feed depth'] = 0.5*D['Cylindrical depth']
        # Typical conventional feed wells are designed for an average downflow velocity
        # of 10-13 mm/s and maximum velocity of 25-30 mm/s. [4]
        peak_flow_safety_factor = 2.5 # assumed based on average and maximum velocities
        D['Downward flow velocity'] = self.downward_flow_velocity*peak_flow_safety_factor # in m/hr
        
        D['Volumetric flow'] =  mixed.get_total_flow('m3/hr')/D['Number of thickeners'] # m3/hr
        Center_feed_area = D['Volumetric flow']/D['Downward flow velocity'] # in m2
        D['Center feed diameter'] = np.sqrt(4*Center_feed_area/np.pi) # in m

        #Diameter of the center feed does not exceed 40% of tank diameter [2]
        if D['Center feed diameter'] > 0.40*D['Thickener diameter']:
            cf_dia = D['Center feed diameter'] 
            tank_dia = D['Thickener diameter']
            warn(f'Diameter of the center feed exceeds 40% of tank diameter. It is {cf_dia*100/tank_dia}% of tank diameter')

        # Amount of concrete required
        D_tank = D['Thickener depth']*39.37 # m to inches 
        # Thickness of the wall concrete [m]. Default to be minimum of 1 feet with 1 inch added for every feet of depth over 12 feet.
        thickness_concrete_wall = (1 + max(D_tank-12, 0)/12)*0.3048 # from feet to m
        inner_diameter = D['Thickener diameter']
        outer_diameter = inner_diameter + 2*thickness_concrete_wall
        volume_cylindercal_wall = (np.pi*D['Cylindrical depth']/4)*(outer_diameter**2 - inner_diameter**2)
        D['Volume of concrete wall'] = volume_cylindercal_wall # in m3
        
        # Concrete slab thickness, [ft], default to be 2 in thicker than the wall thickness. (Brian's code)
        thickness_concrete_slab = thickness_concrete_wall + (2/12)*0.3048 # from inch to m
        outer_diameter_cone = inner_diameter + 2*(thickness_concrete_wall + thickness_concrete_slab)
        volume_conical_wall = (np.pi/(3*4))*(((D['Conical depth'] + thickness_concrete_wall + thickness_concrete_slab)*(outer_diameter_cone**2)) - (D['Conical depth']*(inner_diameter)**2))
        D['Volume of concrete slab'] = volume_conical_wall
        
        # Amount of metal required for center feed
        thickness_metal_wall = 0.3048 # equal to 1 feet, in m (!! NEED A RELIABLE SOURCE !!)
        inner_diameter_center_feed = D['Center feed diameter']
        outer_diameter_center_feed = inner_diameter_center_feed + 2*thickness_metal_wall
        volume_center_feed = (3.14*D['Center feed depth']/4)*(outer_diameter_center_feed**2 - inner_diameter_center_feed**2)
        density_ss = 7930 # kg/m3, 18/8 Chromium
        D['Stainless steel'] = volume_center_feed*density_ss # in kg
       
        # Pumps
        pipe, pumps = self._design_pump()
        D['Pump pipe stainless steel'] = pipe
        D['Pump stainless steel'] = pumps
        
        #For thickener 
        D['Number of pumps'] = D['Number of thickeners']
        
    def _cost(self):
        
        self._mixed.mix_from(self.ins)
        mixed = self._mixed
       
        D = self.design_results
        C = self.baseline_purchase_costs
       
        # Construction of concrete and stainless steel walls
        C['Wall concrete'] = D['Number of thickeners']*D['Volume of concrete wall']*self.wall_concrete_unit_cost
        C['Slab concrete'] = D['Number of thickeners']*D['Volume of concrete slab']*self.slab_concrete_unit_cost
        C['Wall stainless steel'] = D['Number of thickeners']*D['Stainless steel']*self.stainless_steel_unit_cost
       
        # Cost of equipment 
        # Source of scaling exponents: Process Design and Economics for Biochemical Conversion of Lignocellulosic Biomass to Ethanol by NREL.
        
        # Scraper 
        # Source: https://www.alibaba.com/product-detail/Peripheral-driving-clarifier-mud-scraper-waste_1600891102019.html?spm=a2700.details.0.0.47ab45a4TP0DLb
        # base_cost_scraper = 2500
        # base_flow_scraper = 1 # in m3/hr (!!! Need to know whether this is for solids or influent !!!)
        thickener_flow = mixed.get_total_flow('m3/hr')/D['Number of thickeners']
        # C['Scraper'] = D['Number of thickeners']*base_cost_scraper*(thickener_flow/base_flow_scraper)**0.6
        # base_power_scraper = 2.75 # in kW
        # THE EQUATION BELOW IS NOT CORRECT TO SCALE SCRAPER POWER REQUIREMENTS 
        # scraper_power = D['Number of thickeners']*base_power_scraper*(thickener_flow/base_flow_scraper)**0.6
        
        # v notch weir
        # Source: https://www.alibaba.com/product-detail/50mm-Tube-Settler-Media-Modules-Inclined_1600835845218.html?spm=a2700.galleryofferlist.normal_offer.d_title.69135ff6o4kFPb
        base_cost_v_notch_weir = 6888
        base_flow_v_notch_weir = 10 # in m3/hr
        C['v notch weir'] = D['Number of thickeners']*base_cost_v_notch_weir*(thickener_flow/base_flow_v_notch_weir)**0.6
       
        # Pump (construction and maintainance)
        pumps = self.pumps
        add_OPEX = self.add_OPEX
        pump_cost = 0.
        building_cost = 0.
        opex_o = 0.
        opex_m = 0.
       
        for i in pumps:
            p = getattr(self, f'{i}_pump')
            p_cost = p.baseline_purchase_costs
            p_add_opex = p.add_OPEX
            pump_cost += p_cost['Pump']
            building_cost += p_cost['Pump building']
            opex_o += p_add_opex['Pump operating']
            opex_m += p_add_opex['Pump maintenance']

        C['Pumps'] = pump_cost*D['Number of thickeners']
        C['Pump building'] = building_cost*D['Number of thickeners']
        add_OPEX['Pump operating'] = opex_o*D['Number of thickeners']
        add_OPEX['Pump maintenance'] = opex_m*D['Number of thickeners']
       
        # Power
        pumping = 0.
        for ID in self.pumps:
            p = getattr(self, f'{ID}_pump')
            if p is None:
                continue
            pumping += p.power_utility.rate
        
        pumping = pumping*D['Number of thickeners']
        
        self.power_utility.rate += pumping
        # self.power_utility.rate += scraper_power

#%% Centrifuge

# Asign a bare module of 1 to all

# default_F_BM = {
#     'Bowl stainless steel': 1, 
#     'Conveyor': 1, 
#     'Pumps': 1,
#         }
# default_F_BM.update(default_WWTpump_F_BM)

class Centrifuge(Thickener):

    """
    Centrifuge based on BSM2 Layout. [1]
    
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
    .. [1] Gernaey, Krist V., Ulf Jeppsson, Peter A. Vanrolleghem, and John B. Copp.
    Benchmarking of control strategies for wastewater treatment plants. IWA publishing, 2014.
    [2] Metcalf, Leonard, Harrison P. Eddy, and Georg Tchobanoglous. Wastewater 
    engineering: treatment, disposal, and reuse. Vol. 4. New York: McGraw-Hill, 1991.
    [3]Design of Municipal Wastewater Treatment Plants: WEF Manual of Practice 
    No. 8 ASCE Manuals and Reports on Engineering Practice No. 76, Fifth Edition. 
    [4] https://www.alibaba.com/product-detail/Multifunctional-Sludge-Dewatering-Decanter-Centrifuge_1600285055254.html?spm=a2700.galleryofferlist.normal_offer.d_title.1cd75229sPf1UW&s=p
    [5] United States Environmental Protection Agency (EPA) 'Biosolids Technology Fact Sheet Centrifuge Thickening and Dewatering'  
    [6] San Diego (.gov) Chapter - 3 'Solids Treatment Facility' 
    (https://www.sandiego.gov/sites/default/files/legacy/mwwd/pdf/mbc/chapterthree.pdf)
    """
    
    _N_ins = 1
    _N_outs = 2
    _ins_size_is_fixed = False
    
    pumps = ('sludge',)
    
    # Costs
    stainless_steel_unit_cost=1.8 # $/Kg (Taken from Joy's METAB code) https://www.alibaba.com/product-detail/brushed-stainless-steel-plate-304l-stainless_1600391656401.html?spm=a2700.details.0.0.230e67e6IKwwFd
    polymer_cost_by_weight = 2.2 # $/Kg (Source: https://www.alibaba.com/product-detail/dewatering-pool-chemicals-cationic-polyacrylamide-cas_1600194474507.html?spm=a2700.galleryofferlist.topad_classic.i5.5de8615c4zGAhg)
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, isdynamic=False, 
                  init_with='WasteStream', F_BM_default=default_F_BM, thickener_perc=28, TSS_removal_perc=98, 
                  solids_feed_rate = 70, g_factor=2500, rotational_speed = 2500, LtoD = 4, F_BM = default_F_BM,
                  polymer_dosage = 20, h_cylindrical=2, h_conical=1, **kwargs):
        
        Thickener.__init__(self, ID=ID, ins=ins, outs=outs, thermo=thermo, isdynamic=isdynamic,
                      init_with=init_with, F_BM_default=1, thickener_perc=thickener_perc, 
                      TSS_removal_perc=TSS_removal_perc, **kwargs)
        
        self.solids_feed_rate = solids_feed_rate
        self.g_factor = g_factor #unitless, centrifugal acceleration = g_factor*9.81
        self.rotational_speed = rotational_speed #in revolution/min
        self.LtoD = LtoD 
        self.polymer_dosage = polymer_dosage #in (lbs,polymer/tonne,solids)
        self.h_cylindrical = h_cylindrical
        self.h_conical = h_conical
        
    _units = {
        'Number of centrifuges': 'ea',
        'Diameter of bowl': 'm',
        'Total length of bowl': 'm',
        'Length of cylindrical portion': 'm',
        'Length of conical portion': 'm',
        'Volume of bowl': 'm3',
        'Stainless steel for bowl': 'kg',
        'Polymer feed rate': 'kg/hr',
        'Pump pipe stainless steel' : 'kg',
        'Pump stainless steel': 'kg',
        'Number of pumps': 'ea'
    }        
    
    def _design_pump(self):
        ID, pumps = self.ID, self.pumps
        self._sludge.copy_like(self.outs[0])
        sludge = self._sludge
        ins_dct = {
            'sludge': sludge,
            }
        type_dct = dict.fromkeys(pumps, 'sludge')
        inputs_dct = dict.fromkeys(pumps, (1,))
        
        D = self.design_results
        influent_Q = sludge.get_total_flow('m3/hr')/D['Number of centrifuges']
        influent_Q_mgd = influent_Q*0.00634 # m3/hr to MGD
       
        for i in pumps:
            if hasattr(self, f'{i}_pump'):
                p = getattr(self, f'{i}_pump')
                setattr(p, 'add_inputs', inputs_dct[i])
            else:
                ID = f'{ID}_{i}'
                capacity_factor=1
                # No. of pumps = No. of influents
                pump = WWTpump(
                    ID=ID, ins= ins_dct[i], pump_type=type_dct[i],
                    Q_mgd=influent_Q_mgd, add_inputs=inputs_dct[i],
                    capacity_factor=capacity_factor,
                    include_pump_cost=True,
                    include_building_cost=False,
                    include_OM_cost=True,
                    )
                setattr(self, f'{i}_pump', pump)

        pipe_ss, pump_ss = 0., 0.
        for i in pumps:
            p = getattr(self, f'{i}_pump')
            p.simulate()
            p_design = p.design_results
            pipe_ss += p_design['Pump pipe stainless steel']
            pump_ss += p_design['Pump stainless steel']
        return pipe_ss, pump_ss
    
    def _design(self):   
        self._mixed.mix_from(self.ins)
        mixed = self._mixed
        
        D = self.design_results 
        TSS_rmv = self._TSS_rmv
        solids_feed_rate = 44.66*self.solids_feed_rate # 44.66 is factor to convert tonne/day to kg/hr
        # Cake's total solids and TSS are essentially the same (pg. 24-6 [3])
        # If TSS_rmv = 98, then total_mass_dry_solids_removed  = (0.98)*(influent TSS mass)
        total_mass_dry_solids_removed = (TSS_rmv/100)*((mixed.get_TSS()*self.ins[0].F_vol)/1000) # in kg/hr
        D['Number of centrifuges'] = np.ceil(total_mass_dry_solids_removed/solids_feed_rate)
        
        
        # HAVE COMMENTED ALL OF THIS SINCE CENTRIFUGE WOULD PROBABLY BE BROUGHT NOT CONSTRUCTED AT THE FACILITY
        
        # k = 0.00000056 # Based on emprical formula (pg. 24-23 of [3])
        # g = 9.81 # m/s2
        # # The inner diameterof the bowl is calculated based on an empirical formula. 1000 is used to convert mm to m.
        # D['Diameter of bowl'] = (self.g_factor*g)/(k*np.square(self.rotational_speed)*1000) # in m
        # D['Total length of bowl'] =  self.LtoD*D['Diameter of bowl'] 
        # # Sanity check: L should be between 1-7 m, diameter should be around 0.25-0.8 m (Source: [4])
        
        # fraction_cylindrical_portion = 0.8
        # fraction_conical_portion = 1 - fraction_cylindrical_portion
        # D['Length of cylindrical portion'] = fraction_cylindrical_portion*D['Total length of bowl']
        # D['Length of conical portion'] =  fraction_conical_portion*D['Total length of bowl']
        # thickness_of_bowl_wall = 0.1 # in m (!!! NEED A RELIABLE SOURCE !!!)
        # inner_diameter = D['Diameter of bowl']
        # outer_diameter = inner_diameter + 2*thickness_of_bowl_wall
        
        # volume_cylindrical_wall = (np.pi*D['Length of cylindrical portion']/4)*(outer_diameter**2 - inner_diameter**2)
        # volume_conical_wall = (np.pi/3)*(D['Length of conical portion']/4)*(outer_diameter**2 - inner_diameter**2)
        # D['Volume of bowl'] = volume_cylindrical_wall + volume_conical_wall # in m3
        
        # density_ss = 7930 # kg/m3, 18/8 Chromium
        # D['Stainless steel for bowl'] = D['Volume of bowl']*density_ss # in kg
        
        polymer_dosage_rate = 0.000453592*self.polymer_dosage # convert from (polymer (lbs)/solids (tonne)) to (polymer (kg)/solids (kg))
        D['Polymer feed rate'] = (polymer_dosage_rate*solids_feed_rate) # in polymer (kg)/hr
        
        # Pumps
        pipe, pumps = self._design_pump()
        D['Pump pipe stainless steel'] = pipe
        D['Pump stainless steel'] = pumps
        
        # For centrifuges 
        D['Number of pumps'] = D['Number of centrifuges']
        
    def _cost(self):
       
        D = self.design_results
        C = self.baseline_purchase_costs
        
        self._mixed.mix_from(self.ins)
        mixed = self._mixed
       
        # HAVE COMMENTED SINCE CENTRIFUGE WOULD PROBABLY BE BROUGHT NOT CONSTRUCTED AT THE FACILITY
        # Construction of concrete and stainless steel walls
        # C['Bowl stainless steel'] = D['Number of centrifuges']*D['Stainless steel for bowl']*self.stainless_steel_unit_cost
        
        # Conveyor 
        # Source: https://www.alibaba.com/product-detail/Sludge-Dewatering-Centrifuge-Decanter-Centrifuge-For_60448094522.html?spm=a2700.galleryofferlist.p_offer.d_title.1c5c5229I5pQeP&s=p
        base_cost_centrifuge = 16000
        base_mass_flow_centrifuge = 80 # in tonne/hr
        thickener_mass_flow = (mixed.get_total_flow('m3/hr')*mixed.get_TSS())/D['Number of centrifuges'] # IN gm/hr
        gm_to_tonne = 0.000001
        thickener_mass_flow = thickener_mass_flow*gm_to_tonne
        C['Centrifuge'] = D['Number of centrifuges']*base_cost_centrifuge*(thickener_mass_flow/base_mass_flow_centrifuge)**0.6
        
        # base_power_conveyor = 55 # in kW
        # THIS IS NOT THE CORRECT EXPRESSION TO SCALE UP POWER OF CENTRIFUGE
        # conveyor_power = D['Number of centrifuges']*base_power_conveyor*(thickener_flow/base_flow_conveyor)**0.6
        
        # Pump (construction and maintainance)
        pumps = self.pumps
        add_OPEX = self.add_OPEX
        pump_cost = 0.
        building_cost = 0.
        opex_o = 0.
        opex_m = 0.
        
        add_OPEX['Polymer'] = D['Number of centrifuges']*D['Polymer feed rate']*self.polymer_cost_by_weight 
       
        for i in pumps:
            p = getattr(self, f'{i}_pump')
            p_cost = p.baseline_purchase_costs
            p_add_opex = p.add_OPEX
            pump_cost += p_cost['Pump']
            building_cost += p_cost['Pump building']
            opex_o += p_add_opex['Pump operating']
            opex_m += p_add_opex['Pump maintenance']

        C['Pumps'] = pump_cost*D['Number of pumps']
        C['Pump building'] = building_cost*D['Number of pumps']
        add_OPEX['Pump operating'] = opex_o*D['Number of pumps']
        add_OPEX['Pump maintenance'] = opex_m*D['Number of pumps']
       
        # Power
        pumping = 0.
        for ID in self.pumps:
            p = getattr(self, f'{ID}_pump')
            if p is None:
                continue
            pumping += p.power_utility.rate
            
        pumping = pumping*D['Number of pumps']
        self.power_utility.rate += pumping
        # self.power_utility.rate += conveyor_power
#%% Incinerator

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
    >>> CO2 = qs.Component.from_chemical('S_CO2', search_ID='CO2', particle_size='Soluble', degradability='Undegradable', organic=False)
    >>> cmps_test = qs.Components([cmps.S_F, cmps.S_NH4, cmps.X_OHO, cmps.H2O, cmps.S_CH4, cmps.S_O2, cmps.S_N2, cmps.S_H2, cmps.X_Ig_ISS, CO2])
    >>> cmps_test.default_compile()
    >>> qs.set_thermo(cmps_test)
    >>> ws = qs.WasteStream('ws', S_F=10, S_NH4=20, X_OHO=15, H2O=1000)
    >>> natural_gas = qs.WasteStream('nat_gas', phase='g', S_CH4=1000)
    >>> air = qs.WasteStream('air', phase='g', S_O2=210, S_N2=780, S_H2=10)
    >>> from qsdsan.sanunits import Incinerator
    >>> Inc = Incinerator(ID='Inc', ins= (ws, air, natural_gas), outs=('flu_gas', 'ash'), 
    ...                 isdynamic=False)
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
         COD        : 23873.0 mg/L
         BOD        : 14963.2 mg/L
         TC         : 8298.3 mg/L
         TOC        : 8298.3 mg/L
         TN         : 20363.2 mg/L
         TP         : 367.6 mg/L
         TK         : 68.3 mg/L
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
        
        inf = np.asarray(sludge.mass + air.mass + fuel.mass)
        idx_n2 = cmps.index(nitrogen_ID)
        idx_h2o = cmps.index(water_ID)
        
        n2 = inf[idx_n2]
        h2o = inf[idx_h2o]
        
        mass_ash = np.sum(inf*cmps.i_mass*(1-cmps.f_Vmass_Totmass)) \
               - h2o*cmps.H2O.i_mass*(1-cmps.H2O.f_Vmass_Totmass) - n2*cmps.N2.i_mass*(1-cmps.N2.f_Vmass_Totmass)

        # Conservation of mass 
        mass_flue_gas = np.sum(inf*cmps.i_mass) - mass_ash
        mass_co2 = mass_flue_gas - n2*cmps.N2.i_mass - h2o*cmps.H2O.i_mass
        
        flue_gas.set_flow([n2, h2o, (mass_co2/cmps.S_CO2.i_mass)], 
                          'kg/hr', (nitrogen_ID, water_ID, carbon_di_oxide_ID))
        ash_cmp_ID = self.ash_component_ID
        ash_idx = cmps.index(ash_cmp_ID)
        ash.set_flow([mass_ash/cmps.i_mass[ash_idx]/(1-cmps.f_Vmass_Totmass[ash_idx])],'kg/hr', (ash_cmp_ID,))
        
        # Energy balance 
        # self.Heat_sludge = sludge.dry_mass*sludge.F_vol*self.calorific_value_sludge/1000 #in KJ/hr (mg/L)*(m3/hr)*(KJ/kg)=KJ/hr*(1/1000)
        # self.Heat_air = np.sum(air.mass*cmps.i_mass)*self.Cp_air #in KJ/hr 
        # self.Heat_fuel = np.sum(fuel.mass*cmps.i_mass)*self.calorific_value_fuel #in KJ/hr 
        # self.Heat_flue_gas = self.process_efficiency*(self.Heat_sludge + self.Heat_air + self.Heat_fuel)
        
        # # Conservation of energy
        # self.Heat_loss = self.Heat_sludge + self.Heat_air + self.Heat_fuel - self.Heat_flue_gas
        
    def _init_state(self):
        
        sludge, air, fuel = self.ins
        inf = np.asarray(sludge.mass + air.mass + fuel.mass)
        self._state = (24*inf)/1000
        self._dstate = self._state * 0.
        self._cached_state = self._state.copy()
        self._cached_t = 0
        
    # def _cost(self):
        
    #     C = self.baseline_purchase_costs
        
    #     sludge, air, fuel = self.ins
        
    #     solids_load_treated = sludge.get_total_flow('m3/hr')*sludge.get_TSS('mg/L')/1000 # in kg/hr
        
    #     C = self.baseline_purchase_costs
        
    #     # Based on regression equations obtained from CapdetWorks
    #     C['Construction and equipment costs'] = 119629*(solids_load_treated)**0.9282 
    #     C['Installed incinerator costs'] = 114834*(solids_load_treated)**0.9284
    #     C['Slab concrete costs'] = 782.28*(solids_load_treated)**0.9111
    #     C['Building costs'] = 4429.8*(solids_load_treated)**0.911
        
    #     # Based on regression equations obtained from CapdetWorks
    #     add_OPEX = self.add_OPEX 
    #     add_OPEX['Material and supply costs'] = 0.0614*(solids_load_treated)**0.9282
    #     add_OPEX['Energy costs'] = 1.3079*(solids_load_treated)**0.9901
        
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
        
        
        