
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems
Copyright (C) 2020, Quantitative Sustainable Design Group

This module is developed by:
    Lewis Rowles <stetsonsc@gmail.com>
and adapted by:
    Victoria Morgan <tvlmorgan@gmail.com>
    

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/master/LICENSE.txt
for license details.
'''
# %%
from qsdsan import SanUnit, Construction
from chaospy import distributions as shape
from qsdsan.utils.loading import load_data, data_path
__all__ = ('Ultrafiltration',)
import numpy as np

import os
data_path = os.path.join(data_path, 'sanunit_data/_ultrafiltration_reclaimer.csv')


#To scale the system from 0 - 120 users based on units corresponding to 30 users: 
ppl = 120

# #power demand decreases with number of ultrafiltration units added (calculated by Duke Team)
# power_demand_1 = shape.Triangle(lower=675, midpoint = 750, upper=825)
# power_demand_2 = shape.Triangle(lower=405, midpoint = 540, upper=675)
# power_demand_3 = shape.Triangle(lower=382.5, midpoint = 510, upper=637.5)
# power_demand_4 = shape.Triangle(lower=360, midpoint = 480, upper=600)



# if (ppl <=30): power_demand = power_demand_1
# elif (ppl >30 and ppl <=60): power_demand = power_demand_2
# elif (ppl >60 and ppl <=90): power_demand = power_demand_3
# elif (ppl >90 and ppl <=120): power_demand = power_demand_4
    
### 
class Ultrafiltration(SanUnit):
    '''
    Ultrafiltration for removing suspended solids with automted backwash 
    to prolong filter
    
    '''


    def __init__(self, ID='', ins=None, outs=(), if_gridtied=True, ppl = 1, **kwargs):
        
        SanUnit.__init__(self, ID, ins, outs, F_BM_default=1)
        self.if_gridtied = if_gridtied
        self.ppl = ppl
        self.price_ratio = 1  
        data = load_data(path=data_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data
        
        for attr, value in kwargs.items():
            setattr(self, attr, value)
        
# define the number of influent and effluent streams    
    _N_ins = 1
    _N_outs = 2
    
# #Number of users to account for reclaimer units
#     R = np.ceil(self.ppl / 30)

# in _run: define influent and effluent streams and treatment processes 
    def _run(self):
        waste = self.ins[0]
        treated, retentate = self.outs
        treated.copy_like(self.ins[0])   
        
        self.retentate_prcd = (waste.imass['OtherSS']* self.TSS_removal) #mg/L
        retentate.imass['OtherSS'] = self.retentate_prcd
    
    def _design(self):
        design = self.design_results
        
        
        design['Plastic'] = P_quant = self.Plastic_weight 
        design['Steel'] = S_quant = self.Steel_weight 
        
        
        self.construction = (
              Construction(item='Plastic', quantity = P_quant, quantity_unit = 'kg'),
            Construction(item='Steel', quantity = S_quant, quantity_unit = 'kg')
            )
        self.add_construction(add_cost=False)        
 
    
    def _cost(self):
        C =self.baseline_purchase_costs
        C['Pipes'] = (self.one_in_pipe_SCH40 + self.onehalf_in_pipe_SCH40 + self.three_in_pipe_SCH80)
        C['fittings'] = (self.one_in_elbow_SCH80 + self.one_in_tee_SCH80 + self.one_in_SCH80
            + self.one_onehalf_in_SCH80 + self.onehalf_in_SCH80 + self.three_in_SCH80_endcap + self.one_one_NB_MTA
            + self.one_onehalf_NB_MTA + self.foot_valve + self.one_onehalf_in_SCH80_threadedtee + self.three_in_pipe_clamp
            + self.one_in_pipe_clamp + self.onehalf_in_pipe_clamp + self.two_way_valve + self.UF_brush) 
        C['UF_unit'] = (self.UF_unit) * R               
          
        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio
        
        
        if self.if_gridtied:
            self.power_demand = self.power_demand # kWh/day
        else:
            self.power_demand = 0
            
        self.power_utility(self.power_demand)
    
    def _calc_replacement_cost(self):
        ultrafiltration_replacement_cost = (self.replacement_costs / 20) #USD/yr
        return ultrafiltration_replacement_cost/ (365 * 24) * self.price_ratio # USD/hr (all items are per hour)
       


        

