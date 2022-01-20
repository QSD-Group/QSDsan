#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems
Copyright (C) 2020, Quantitative Sustainable Design Group

This module is developed by:
    Lewis Rowles <stetsonsc@gmail.com>
    

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/master/LICENSE.txt
for license details.
'''
# %%


import numpy as np
from qsdsan import SanUnit, Construction
from qsdsan.utils.loading import load_data, data_path
import os
__all__ = ('StruvitePrecipitation',)

#path to csv with all the inputs
#data_path = '/Users/lewisrowles/opt/anaconda3/lib/python3.8/site-packages/exposan/biogenic_refinery/_struvite_precipitation.csv' 
data_path += '/sanunit_data/_struvite_precipitation.tsv'

### 
class StruvitePrecipitation(SanUnit):
    '''
    Stuvite Precipitation for P recovery from liquid stream. Solid stuvite is recovered.

    
    Reference documents
    -------------------
    N/A
    
    Parameters
    ----------
    ins : WasteStream (liquid), MagnesiumHydroxide, MagnesiumCarbonate, FilterBag
        
    outs : WasteStream, Struvite
        

        
    References
    ----------
    .. Lohman et al., Advancing Sustainable Sanitation and Agriculture 
    through Investments in Human-Derived Nutrient Systems. 
    Environ. Sci. Technol. 2020, 54, (15), 9217-9227.
    https://dx.doi.org/10.1021/acs.est.0c03764
    
    .. Tarpeh et al., Evaluating ion exchange for nitrogen recovery from 
    source-separated urine in Nairobi, Kenya. Development Engineering. 2018, 
    3, 188–195.
    https://doi.org/10.1016/j.deveng.2018.07.002
    
    '''
    

    def __init__(self, ID='', ins=None, outs=(), **kwargs):
        
        SanUnit.__init__(self, ID, ins, outs)

# load data from csv each name will be self.name    
        data = load_data(path=data_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data
        
        for attr, value in kwargs.items():
            setattr(self, attr, value)




        
# define the number of influent and effluent streams    
    _N_ins = 4
    _N_outs = 2

# in _run: define influent and effluent streams and treatment processes 
    def _run(self):
        waste, magnesium_hydroxide, magnesium_carbonate, bag_filter = self.ins
        treated, struvite = self.outs
        treated.copy_like(self.ins[0])
        magnesium_hydroxide.phase = 's'
        magnesium_carbonate.phase = 's'
        bag_filter = 's'
        struvite.phase = 's'
    
        # recovery of N, K, P
        N_recovered = (waste.imass['P'] * self.P_rec_1 / self.MW_P 
                       * self.N_P_ratio_struvite * self.MW_N)  # (kg N/hr) total N recovered
        P_recovered = waste.imass['P'] * self.P_rec_1  # (kg P/hr) total P recovered
        K_recovered = waste.imass['K'] * self.K_rec_1 # (kg K/hr) total K recovered
        treated.imass['NH3'] =  waste.imass['NH3'] - N_recovered # kg N / hr
        treated.imass['P'] =  waste.imass['P'] - P_recovered # kg P / hr
        treated.imass['K'] =  waste.imass['K'] - K_recovered # kg K / hr
        
        # set values needed for _design and _cost as attributes
        self.volume_treated = waste.F_vol * 1000 * 24 # L liq / d 
        self.quantity_tanks = np.ceil(self.volume_treated / self.cycles_per_day 
                                      / self.reactor_volume) # number of tanks
        
        # !!! add MagnesiumHydroxide, MagnesiumCarbonate, FilterBag as influent streams and Struvite as effluent
        magnesium_hydroxide_demand_time = (waste.imass['P'] / self.MW_P * self.Mg_dose 
                                           * self.Mg_MgOH2_ratio * self.MW_MgOH2) # (kg Mg(OH)2 per hr)
        magnesium_hydroxide.imass['MagnesiumHydroxide'] = magnesium_hydroxide_demand_time
        
        magnesium_carbonate_demand_time = (waste.imass['P'] / self.MW_P * self.Mg_dose 
                                           * self.Mg_MgCO3_ratio * self.MW_MgCO3)  # (kg MgCO3 per hr)
        #assume only magnesium hydroxide is used
        #magnesium_carbonate.imass['MagnesiumCarbonate'] = magnesium_carbonate_demand_time
        magnesium_carbonate.imass['MagnesiumCarbonate'] = 0
        
        struvite_production_time = P_recovered / self.MW_P * self.MW_struvite # kg (NH4)MgPO4•6(H2O) / hr
        struvite.imass['Struvite'] = struvite_production_time

        filter_bag_demand_time = self.quantity_tanks * self.cycles_per_day / self.filter_reuse / 24 # bags/hr
        self.ins[3].imass['FilterBag'] = filter_bag_demand_time # used in place of line below, not sure why its not working
        #bag_filter.imass['FilterBag'] = filter_bag_demand_time          
  
     
    #_design will include all the construction or captial impacts  
    def _design(self):
        design = self.design_results
        # defining the quantities of materials/items
        # note that these items to be to be in the _impacts_items.xlsx
        
        

        design['StainlessSteel'] = SS_quant = self.quantity_tanks * self.reactor_weight # kg SS
        design['PVC'] = PVC_quant = self.quantity_tanks * self.material_P_pipe * self.pvc_mass  # kg PVC
     

        
        self.construction = (
            Construction(item='StainlessSteel', quantity = SS_quant, quantity_unit = 'kg'),
            Construction(item='PVC', quantity = PVC_quant, quantity_unit = 'kg'),
            )
        self.add_construction()
        
    
    #_cost based on amount of steel and stainless plus individual components
    def _cost(self):
        #purchase_costs is used for capital costs
        #can use quantities from above (e.g., self.design_results['StainlessSteel'])
        #can be broken down as specific items within purchase_costs or grouped (e.g., 'Misc. parts')
        self.baseline_purchase_costs['Reactor'] = (self.quantity_tanks * self.cost_P_reactor)
        self.baseline_purchase_costs['Stirrer'] = (self.quantity_tanks * self.cost_P_stirrer)
        self.baseline_purchase_costs['PVC'] = (self.quantity_tanks * self.material_P_pipe * self.cost_P_pipe)
         
        self.F_BM = dict.fromkeys(self.baseline_purchase_costs.keys(), 1)
        
        #certain parts need to be replaced based on an expected lifefime
        #the cost of these parts is considered along with the cost of the labor to replace them
        struvite_replacement_parts_annual_cost = 0 # USD/yr only accounts for time running
        
        struvite_annual_maintenance = 0 #USD/yr only accounts for time running
        

        purchase_costs = 0
        for i in self.baseline_purchase_costs:
            purchase_costs+= self.baseline_purchase_costs[i]
        
        opex_from_cap = 0.35/4 * purchase_costs # $/yr
        
        self.add_OPEX =  (struvite_replacement_parts_annual_cost + struvite_annual_maintenance + opex_from_cap) / (365 * 24) # USD/hr (all items are per hour)
        
        # costs associated with full time opperators can be added in the TEA as staff
      


