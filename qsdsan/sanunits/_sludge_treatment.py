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

from .. import SanUnit
import numpy as np

__all__ = ('Thickener', 'DewateringUnit',)

class Thickener(SanUnit):
    
    """
    Thickener based on BSM2 Layout. [1]

    Parameters
    ----------
    ID : str
        ID for the Thickener. The default is ''.
    ins : class:`WasteStream`
        Influent to the clarifier. Expected number of influent is 1. 
    outs : class:`WasteStream`
        Treated effluent and sludge.
    thickener_perc : float
        The percentage of Suspended Sludge in the underflow of the thickener.[1]
    TSS_removal_perc : float
        The percentage of suspended solids removed in the thickener.[1]
    solids_loading_rate : float
        Solid loading rate in the thickener.[2]
    h_cylinderical = float
        Height of cylinder forming the thickener.[2]

    References
    ----------
    .. [1] Gernaey, Krist V., Ulf Jeppsson, Peter A. Vanrolleghem, and John B. Copp.
    Benchmarking of control strategies for wastewater treatment plants. IWA publishing, 2014.
    [2] Metcalf, Leonard, Harrison P. Eddy, and Georg Tchobanoglous. Wastewater 
    engineering: treatment, disposal, and reuse. Vol. 4. New York: McGraw-Hill, 1991.
    """
    
    _N_ins = 1
    _N_outs = 2
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, isdynamic=False, 
                  init_with='WasteStream', F_BM_default=None, thickener_perc=7, 
                  TSS_removal_perc=98, solids_loading_rate = 50, h_cylinderical=2, **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo, isdynamic=isdynamic, 
                         init_with=init_with, F_BM_default=F_BM_default)
        
        self.thickener_perc = thickener_perc 
        self.TSS_removal_perc = TSS_removal_perc
        self.solids_loading_rate = solids_loading_rate 
        self.h_cylinderical = h_cylinderical
        
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
        inf, = self.ins
        if not self.ins: return
        elif inf.isempty(): return
        else: 
            TSS_in = inf.get_TSS()
            if TSS_in > 0:
                thickener_factor = self._tp*10000/self.ins[0].get_TSS()
                if thickener_factor<1:
                    thickener_factor=1
                return thickener_factor
            else: return None
    
    @property
    def thinning_factor(self):
        thickener_factor = self.thickener_factor
        if thickener_factor<1:
            thinning_factor=0
        else:
            Qu_factor = self._TSS_rmv/(100*thickener_factor)
            thinning_factor = (1 - (self._TSS_rmv/100))/(1 - Qu_factor)
        return thinning_factor
    
    def _run(self):
        
        inf, = self.ins
        uf, of = self.outs
        cmps = self.components
        
        TSS_rmv = self._TSS_rmv
        thinning_factor = self.thinning_factor
        thickener_factor = self.thickener_factor
        
        Ze = (1 - thinning_factor)/(thickener_factor - thinning_factor)*inf.mass*cmps.s
        Zs = (thickener_factor - 1)/(thickener_factor - thinning_factor)*inf.mass*cmps.s
        
        Xe = (1 - TSS_rmv/100)*inf.mass*cmps.x
        Xs = (TSS_rmv/100)*inf.mass*cmps.x
        
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
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, isdynamic=False, 
                  init_with='WasteStream', F_BM_default=None, thickener_perc=28, TSS_removal_perc=98, 
                  number_of_centrifuges=1, specific_gravity_sludge=1.03, cake_density=965, 
                  centrifugal_force=2500, rotational_speed = 40, polymer_dosage_per_kg_of_sludge = 0.0075, 
                  h_cylinderical=2, h_conical=1, **kwargs):
        Thickener.__init__(self, ID=ID, ins=ins, outs=outs, thermo=thermo, isdynamic=isdynamic,
                      init_with=init_with, F_BM_default=F_BM_default, thickener_perc=thickener_perc, 
                      TSS_removal_perc=TSS_removal_perc, **kwargs)
        self.number_of_centrifuges=number_of_centrifuges
        self.specific_gravity_sludge=specific_gravity_sludge
        self.cake_density=cake_density #in kg/m3
        self.centrifugal_force = centrifugal_force #in Newton
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
        design['Diameter'] = 2*(self.centrifugal_force*9.81/np.square(2*np.pi*self.rotational_speed)) #in m
        design['Polymer feed rate'] = (self.polymer_dosage_per_kg_of_sludge*sludge_feed_rate) # in kg/hr
        design['Projected Area at Inlet'] = np.pi*np.square(design['Diameter']/2) #in m2
        design['Hydraulic_Loading'] = (self.ins[0].F_vol*24)/design['Area'] #in m3/(m2*day)
        design['Volume'] = np.pi*np.square(design['Diameter']/2)*(self.h_cylinderical + (self.h_conical/3)) #in m3
        design['Curved Surface Area'] = np.pi*design['Diameter']/2*(2*self.h_cylinderical + np.sqrt(np.square(design['Diameter']/2) + np.square(self.h_conical))) #in m2