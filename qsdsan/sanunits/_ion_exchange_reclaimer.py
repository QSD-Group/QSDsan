
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
__all__ = ('IonExchangeReclaimer',)

#path to csv with all the inputs

data_path += 'sanunit_data/_ion_exchange_reclaimer.csv'

### 
class IonExchangeReclaimer(SanUnit):
    '''
    Ion Exchange for N recovery from liquid stream. Concentrated NH3 is recovered.

    
    Reference documents
    -------------------
    Duke data
    
    Parameters
    ----------
    ins : WasteStream (liquid), Zeolite, GAC, KCl
        
    outs : WasteStream, SpentZeolite, SpentGAC, Concentrated NH3
        

        
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
    _N_ins = 4
    _N_outs = 4

# in _run: define influent and effluent streams and treatment processes 
    def _run(self):
        waste, zeolite_in, gac_in, KCl = self.ins
        treated, zeolite_out, gac_out, conc_NH3 = self.outs
        treated.copy_like(self.ins[0])
        zeolite_in.phase = 's'
        zeolite_out.phase = 's'
        gac_in.phase = 's'
        gac_out.phase = 's'
        conc_NH3.phase = 'l'
        KCl.phase = 's'
        
        
        #from Lena's info 
        # Zeolite
        zeolite_demand_time = self.zeolite_weight / (self.zeolite_lifetime*365*24) # kg zeolite / hr
        #zeolite_demand_influent = zeolite_demand_time / waste.F_vol # kg zeolite / m3 treated
        #zeolite_cost_day = zeolite_demand_time * 24 * (self.zeolite/self.zeolite_weight_per_cost) # $ zeolite / d
        zeolite_in.imass['Zeolite'] = zeolite_demand_time
        zeolite_out.imass['Zeolite'] = zeolite_demand_time
        
        # GAC
        gac_demand_time = self.gac_weight / (self.gac_lifetime*365*24) # kg gac / hr
        #gac_demand_influent = gac_demand_time / waste.F_vol # kg gac / m3 treated
        #gac_cost_day = gac_demand_time * 24 * (self.gac / self.gac_weight) # $ gac / d
        gac_in.imass['GAC'] = gac_demand_time
        gac_out.imass['GAC'] = gac_demand_time
        
        # KCl
        self.KCl_demand_time = (self.KCl_weight * self.regen_freq_per_yr) / (365*24) # kg KCl / hr
        #KCl_cost_day = KCl_demand_time * 24 * (self.KCl / self.KCl_weight_per_cost) # $ KCl / d
        KCl.imass['PotassiumChloride'] = self.KCl_demand_time


        #!!! During storage most N as urea goes to NH3, should that 
        # conversion be added or just use total N here? 
        self.N_removed = waste.imass['NH3'] * (self.TN_removal)
        self.N_recovered = self.N_removed * (self.desorption_recovery_efficiency)  # kg N / hr
        treated.imass['NH3'] =  waste.imass['NH3'] - self.N_removed # kg N / hr
        conc_NH3.imass['NH3'] = self.N_recovered # kg N / hr
        conc_NH3.imass['PotassiumChloride'] = self.KCl_demand_time
        zeolite_out.imass['NH3'] = self.N_removed - self.N_recovered
        
        
    def _design(self):
        design = self.design_results
        
        design['Plastic'] = P_quant = self.Plastic_weight
        design['PVC'] = PVC_quant = self.PVC_weight
        design['Steel'] = S_quant = self.Steel_weight
        
        self.construction = (
            Construction(item='Plastic', quantity = P_quant, quantity_unit = 'kg'),
            Construction(item='PVC', quantity = PVC_quant, quantity_unit = 'kg'),
            Construction(item='Steel', quantity = S_quant, quantity_unit = 'kg'),
            )
        self.add_construction(add_cost=False)        
 
        
     #_cost based on amount of steel and stainless plus individual components
    def _cost(self):

        self.purchase_costs['Pipes'] = (self.four_in_pipe_SCH40 + self.four_in_pipe_SCH80)
        self.purchase_costs['fittings'] = (self.four_in_pipe_SCH80_endcap + self.NRV + self.connector + 
            self.ball_valve + self.three_eight_elbow + self.ten_ten_mm_tee + self.OD_tube + self.four_in_pipe_clamp)
                                           
        self.purchase_costs['GAC_Zeolite'] = (self.GAC_zeolite_mesh)  
            
        self._BM = dict.fromkeys(self.purchase_costs.keys(), 1)
              
        # self.add_OPEX =  self._calc_replacement_cost() + self._calc_maintenance_labor_cost()
    
        
    def _calc_replacement_cost(self):
        ion_exchange_replacement_cost = (self.GAC_zeolite_mesh * self.zeolite_lifetime) #USD/yr
        return ion_exchange_replacement_cost/ (365 * 24) # USD/hr (all items are per hour)
                  
    # def _calc_maintenance_labor_cost(self):
    #     ion_exchange_maintenance_labor = ((self.labor_maintenance_GAC_replacement * self.wages))
    #     return ion_exchange_maintenance_labor/ (365 * 24) # USD/hr (all items are per hour)




       
