
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


import numpy as np
from qsdsan import SanUnit, Construction
from qsdsan.utils.loading import load_data, data_path
import os
__all__ = ('Ultrafiltration',)

#path to csv with all the inputs

data_path += 'sanunit_data/_ultrafiltration.tsv'

### 
class Ultrafiltration(SanUnit):
    '''
    Ultrafiltration for removing suspended solids with automted backwash 
    to prolong filter
    
    '''
    

    def __init__(self, ID='', ins=None, outs=(), if_gridtied=True, **kwargs):
        
        SanUnit.__init__(self, ID, ins, outs)
        self.if_gridtied = if_gridtied

# load data from csv each name will be self.name    
        data = load_data(path=data_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data
        
        for attr, value in kwargs.items():
            setattr(self, attr, value)




        
# define the number of influent and effluent streams    
    _N_ins = 1
    _N_outs = 4

# in _run: define influent and effluent streams and treatment processes 
    def _run(self):
        waste = self.ins[0]
        treated, sludge = self.outs
        treated.copy_like(self.ins[0])   
        
        self.sludge_prcd = (waste.imass['OtherSS']* self.TSS_removal) #mg/L
        sludge.imass['OtherSS'] = self.sludge_prcd
    
    def _design(self):
        design = self.design_results
        
        design['Plastic'] = P_quant = self.plastic_weight
        design['PVC'] = PVC_quant = self.PVC_weight
        design['Steel'] = S_quant = self.steel_weight
        
        self.construction = (
            Construction(item='Plastic', quantity = P_quant, quantity_unit = 'kg'),
            Construction(item='PVC', quantity = PVC_quant, quantity_unit = 'kg'),
            Construction(item='Steel', quantity = S_quant, quantity_unit = 'kg'),
            )
        self.add_construction(add_cost=False)        
 
        
     #_cost based on amount of steel and stainless plus individual components
    def _cost(self):

        self.purchase_costs['Pipes'] = (self.one_in_pipe_SCH40 + self.onehalf_in_pipe_SCH40 + self.three_in_pipe_SCH80)
        self.purchase_costs['fittings'] = (self.one_in_elbow_SCH80 + self.one_in_tee_SCH80 + self.one_in_SCH80
            + self.one_onehalf_in_SCH80 + self.onehalf_in_SCH80 + self.three_in_SCH80_endcap + self.one_one_NB_MTA
            + self.one_onehalf_NB_MTA + self.foot_valve + self.one_onehalf_in_SCH80_threadedtee + self.three_in_pipe_clamp
            + self.one_in_pipe_clamp + self.onehalf_in_pipe_clamp + self.two_way_valve + self.UF_brush)   
        self.purchase_costs['UF_unit'] = (self.UF_unit)                                 
            
        self._BM = dict.fromkeys(self.purchase_costs.keys(), 1)
              
        self.add_OPEX =   + self._calc_maintenance_labor_cost()
    
        
    def _calc_replacement_cost(self):
        ion_exchange_replacement_cost = (self.zeolite_bag_mesh * self.zeolite_lifetime) #USD/yr
        return ion_exchange_replacement_cost/ (365 * 24) # USD/hr (all items are per hour)
                  
    def _calc_maintenance_labor_cost(self):
        ion_exchange_maintenance_labor = ((self.labor_maintenance_GAC_replacement * self.wages))
        return ion_exchange_maintenance_labor/ (365 * 24) # USD/hr (all items are per hour)

