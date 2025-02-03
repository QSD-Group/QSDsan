from .. import SanUnit, WasteStream
from ..sanunits import IdealClarifier
from ..sanunits import WWTpump
from ..sanunits._pumping import default_F_BM as default_WWTpump_F_BM

__all__ = ('AeratedGritChamber',)

F_BM_pump = 1.18*(1 + 0.007/100) # 0.007 is for miscellaneous costs

default_F_BM = {
        'Pumps': F_BM_pump,
        'Pump building': F_BM_pump,
        }

default_equipment_lifetime = {
    'Pumps': 15,
    'Pump pipe stainless steel': 15,
    'Pump stainless steel': 15,
    }

# Assign a bare module of 1 to all
default_F_BM = {
        'Wall concrete': 1.,
        'Slab concrete': 1.,
        'Wall stainless steel': 1.,
        'Scraper': 1,
        'v notch weir': 1,
        'Pumps': 1
        }
default_F_BM.update(default_WWTpump_F_BM)

class AeratedGritChamber(IdealClarifier):
    
    """
    A rectangular aerated grit chamber with a settling model similar to Ideal Clarifier.
    
    Parameters
    ----------
    
    depth : float, optional
        Depth of the grit chamber [m]. Typical depths range from 2 m to 5 m. [1] 
        Default value of 3 m would be used here. 
    freeboard : float, optional
        Freeboard added to the depth of the tank [m]. The default is 0.5 m. [1]
    W_to_D : float, optional
        Design width-to-depth ratio of the grit chamber. Ranges from 1 to 5.
        Default value is 1.5. [1]
    detention_time : float, optional
        The duration for which the wastewater stays in the grit chamber [minutes].
        Typical values range from 2-5 minutes. The default is 3 minutes. [1]
    air_supply : float, optional
        The air supply rate per unit length of the chamber [m3/min/m]. Ranges from 0.2-0.5 m3/min/m.
        The default is 0.48 m3/min/m [1,2] 
    grit_production_rate : float, optional
        The grit production rate for a combined collection system [m3/1000 m3]. 
        Ranges form 0.004-0.20 m3/1000 m3. The default is 0.015 m3/1000 m3. [1]
    Solids removal efficiency : float, optional
        Removal efficiency of total suspended solids (TSS) [%]. 
        Default of 6.68% is used.[3] # Need more reliable and relevant sources.
        
    F_BM : dict
        Equipment bare modules.

    
    References
    ----------
    [1] Metcalf, Leonard, Harrison P. Eddy, and Georg Tchobanoglous. Wastewater 
    engineering: treatment, disposal, and reuse. Vol. 4. New York: McGraw-Hill, 1991.
    [2] The Water Environment Federation (WEF). (2018). Design of Water Resource Recovery Facilities.
    6th ed. McGraw-Hill Education: New York, Chicago, San Francisco, Athens, London, Madrid, Mexico City, 
    Milan, New Delhi, Singapore, Sydney, Toronto.
    https://www.accessengineeringlibrary.com/content/book/9781260031188
    [3] He L, Zhang Y, Song D, Ou Z, Xie Z, Yang S, Guan W, Dong C, Zhang Y. 
    Influence of Pretreatment System on Inorganic Suspended Solids for Influent in Wastewater Treatment Plant. 
    J Environ Public Health. 2022 Sep 28;2022:2768883. doi: 10.1155/2022/2768883. PMID: 36213012; PMCID: PMC9534684.
    [4] Foley, J., De Haas, D., Hartley, K., & Lant, P. (2010). Comprehensive life cycle inventories 
    of alternative wastewater treatment systems. Water Research, 44(5), 1654–1666. 
    https://doi.org/10.1016/j.watres.2009.11.031

    """
    
    _N_ins = 1
    _N_outs = 2
    _ins_size_is_fixed = False
    _outs_size_is_fixed = True
    
    # Costs
    wall_concrete_unit_cost = 1081.73 # $/m3 (Hydromantis. CapdetWorks 4.0. https://www.hydromantis.com/CapdetWorks.html)
    slab_concrete_unit_cost = 582.48 # $/m3 (Hydromantis. CapdetWorks 4.0. https://www.hydromantis.com/CapdetWorks.html)
    stainless_steel_unit_cost=1.8 # Alibaba. Brushed Stainless Steel Plate 304. https://www.alibaba.com/product-detail/brushed-stainless-steel-plate-304l-stainless_1600391656401.html?spm=a2700.details.0.0.230e67e6IKwwFd
    
    pumps = ('sludge',)
    
    def __init__(self, ID='', ins=None, outs=(),
                 depth=3, freeboard=0.5, W_to_D=1.5, detention_time=3, 
                 air_supply=0.48, grit_production_rate=0.015,
                 sludge_flow_rate=200, solids_removal_efficiency=0.0668,
                 thermo=None, isdynamic=False, init_with='WasteStream', 
                 F_BM=default_F_BM, **kwargs):
        super().__init__(ID, ins, outs, thermo,
                         sludge_flow_rate=sludge_flow_rate, 
                         solids_removal_efficiency=solids_removal_efficiency,
                         # thermo=thermo, 
                         isdynamic=isdynamic, 
                         init_with=init_with)
        self.ID = ID
        self._d = depth
        self.freeboard = freeboard
        self.W_to_D = W_to_D
        self._t = detention_time
        self._q_air = air_supply
        self.p_grit = grit_production_rate
        self.solids_removal_efficiency = solids_removal_efficiency
        self.F_BM.update(F_BM)
        self._sludge = WasteStream(f'{ID}_sludge')    

    @property
    def depth(self):
        '''[float] Depth of the grit chamber in m.'''
        return self._d

    @depth.setter
    def depth(self, depth):
        self._d = depth

    @property
    def detention_time(self):
        '''[float] Hydrualic retention time in minutes.'''
        return self._t

    @detention_time.setter
    def detention_time(self, time):
        self._t = time   

    @property
    def air_supply(self):
        '''[float] Air supply rate per unit length of the grit chamber in m3/min/m.'''
        return self._q_air
    
    @air_supply.setter
    def air_supply(self, Q_air):
        self._q_air = Q_air
            
    def _design_pump(self):
        ID, pumps = self.ID, self.pumps
        sludge = self._sludge
        sludge.copy_like(self.outs[1])
        
        ins_dct = {
            'sludge': sludge,
            }
        
        type_dct = dict.fromkeys(pumps, 'sludge')
        inputs_dct = dict.fromkeys(pumps, (1,))
        
        D = self.design_results
        influent_Q = sludge.get_total_flow('m3/hr')
        influent_Q_mgd = influent_Q*0.00634 # m3/hr to MGD
       
        for i in pumps:
            if hasattr(self, f'{i}_pump'):
                p = getattr(self, f'{i}_pump')
                setattr(p, 'add_inputs', inputs_dct[i])
            else:
                ID = f'{ID}_{i}'
                capacity_factor=1
                pump = WWTpump(
                    ID=ID, ins= ins_dct[i], thermo = self.thermo, pump_type=type_dct[i],
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
        'Number of grit chambers': 'ea',
        'Volumetric flow': 'm3/day',
        'Volume of grit chamber': 'm3',
        'Surface area': 'm2',
        'Grit chamber depth': 'm',
        'Grit chamber width': 'm',
        'Grit chamber length': 'm',
        'Hydraulic Retention Time': 'hr', 
        'Air supply': 'm3/hr',
        'Grit production': 'm3/day',
        'Volume of wall concrete': 'm3',
        'Volume of slab concrete': 'm3',
        'Reinforcement steel': 'kg',
        'Stainless steel': 'kg',
        'Pump pipe stainless steel' : 'kg',
        'Pump stainless steel': 'kg',
        'Number of pumps': 'ea'
    }
    
    def _design(self):
        
        D = self.design_results
        # U = self._units
        
        D['Number of grit chambers'] = 2 # 1 in use and 1 redundant.   


        Q_in = self.ins[0].F_vol * 24
        D['Volumetric flow'] = Q_in 
        D['Hydraulic retention time'] = self._t/60
        D['Grit chamber depth'] = self._d

        # Calculating required volume and surface area of grit chamber    
        V_chamber = Q_in * self._t/24
        D['Volume of grit chamber'] = V_chamber
        area = V_chamber/self._d
        D['Surface area'] = area

        # Calculating dimensions of grit chamber
        width = self._d * self.W_to_D
        D['Grit chammber width'] = width
        length = area/width
        D['Grit chamber length'] = length

        # Estimating air supply required
        Q_air = self._q_air * length * 60 # from m3/min to m3/hr
        D['Air supply'] = Q_air

        # Estimating daily grit production
        Q_grit = self.p_grit * Q_in # in m3/day
        D['Grit production'] = Q_grit

        # Estimating concrete and reinforcement steel
        thickness_concrete_wall = (1 + max(self._d/12, 0)) * 0.3048 # from ft to m
        thickness_concrete_slab = thickness_concrete_wall + (2/12)*0.3048 # from inch to m
        VSW = 2 * width * thickness_concrete_wall * (self._d + self.freeboard) # short walls
        VLW = 2 * (length + 2*thickness_concrete_wall) * thickness_concrete_wall * (self._d + self.freeboard) # long walls
        VS = (length + 2*thickness_concrete_wall) * (width + 2*thickness_concrete_wall) * thickness_concrete_slab # slab
        D['Volume of wall concrete'] = VSW + VLW
        D['Volume of slab concrete'] = VS
        D['Reinforcement steel'] = 77.58 * (VSW + VLW + VS) # [4]
       
        # Pumps
        pipe, pumps = self._design_pump()
        D['Pump pipe stainless steel'] = pipe
        D['Pump stainless steel'] = pumps
        
        #For primary clarifier 
        D['Number of pumps'] = D['Number of grit chambers']
   
        # for key in D: print(f'{key}: {D[key]} {U[key]}')
        
    def _cost(self):
        D = self.design_results
        C = self.baseline_purchase_costs
       
        # Construction of concrete and stainless steel walls
        C['Wall concrete'] = D['Number of grit chambers']*D['Volume of wall concrete']*self.wall_concrete_unit_cost
        
        C['Slab concrete'] = D['Number of grit chambers']*D['Volume of slab concrete']*self.slab_concrete_unit_cost
       
        # Cost of equipment 
        
        # Source of scaling exponents: Process Design and Economics for Biochemical Conversion of Lignocellulosic Biomass to Ethanol by NREL.
        
        # Scraper 
        # Source: https://www.alibaba.com/product-detail/Peripheral-driving-clarifier-mud-scraper-waste_1600891102019.html?spm=a2700.details.0.0.47ab45a4TP0DLb
        # base_cost_scraper = 2500
        # base_flow_scraper = 1 # in m3/hr (!!! Need to know whether this is for solids or influent !!!)
        clarifier_flow = D['Volumetric flow']/24
        
        # C['Scraper'] =  D['Number of clarifiers']*base_cost_scraper*(clarifier_flow/base_flow_scraper)**0.6
        
        # base_power_scraper = 2.75 # in kW
        # THE EQUATION BELOW IS NOT CORRECT TO SCALE SCRAPER POWER REQUIREMENTS 
        # scraper_power = D['Number of clarifiers']*base_power_scraper*(clarifier_flow/base_flow_scraper)**0.6
        
        # v notch weir
        # Source: https://www.alibaba.com/product-detail/50mm-Tube-Settler-Media-Modules-Inclined_1600835845218.html?spm=a2700.galleryofferlist.normal_offer.d_title.69135ff6o4kFPb
       
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

        N = D['Number of grit chambers']
        C['Pumps'] = pump_cost*N
        C['Pump building'] = building_cost*N
        add_OPEX['Pump operating'] = opex_o*N
        add_OPEX['Pump maintenance'] = opex_m*N
       
        # Power
        pumping = 0.
        for ID in self.pumps:
            p = getattr(self, f'{ID}_pump')
            if p is None:
                continue
            pumping += p.power_utility.rate
                
        self.power_utility.rate += pumping*N
        # self.power_utility.rate += scraper_power        