
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
__all__ = ('IonExchangeNH3',)

#path to csv with all the inputs
#data_path = '/Users/lewisrowles/opt/anaconda3/lib/python3.8/site-packages/exposan/biogenic_refinery/_ion_exchange_NH3.csv'
#data_path = os.path.abspath(os.path.dirname('_ion_exchange_NH3.csv'))
data_path += 'sanunit_data/_ion_exchange_NH3.tsv'
# data_path += 'sanunit_data/_carbonizer_base.csv'

### 
class IonExchangeNH3(SanUnit):
    '''
    Ion Exchange for N recovery from liquid stream. Concentrated NH3 is recovered.

    
    Reference documents
    -------------------
    N/A
    
    Parameters
    ----------
    ins : WasteStream (liquid), FreshResin, H2SO4
        
    outs : WasteStream, SpentResin, ConcentratedNH3
        

        
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
    _N_ins = 3
    _N_outs = 3

# in _run: define influent and effluent streams and treatment processes 
    def _run(self):
        waste, resin_in, H2SO4 = self.ins
        treated, resin_out, conc_NH3 = self.outs
        treated.copy_like(self.ins[0])
        resin_in.phase = 'l'
        resin_out.phase = 'l'
        conc_NH3.phase = 'l'
        
        
        #!!! During storage most N as urea goes to NH3, should that 
        # conversion be added or just use total N here? 
        N_recovered = waste.imass['NH3'] * (self.N_rec_2) # kg N / hr
        treated.imass['NH3'] =  waste.imass['NH3'] - N_recovered # kg N / hr
        conc_NH3.imass['NH3'] = N_recovered # kg N / hr
        
        # following Terpeh et al. 2018 estimates for regenerating resin
        # assume volume of eluent = volume of urine
        # !!! cost associated with this influent stream neg compared to cost of resin and acid?
        #conc_NH3.F_vol = waste.F_vol
        # 0.1 M H2SO4 (0.65%) used to regenerate resin
        
        # !!! need to add SO4 as a component? 
        # conc_NH3.imass['SO4'] = 0.1 * 98 * conc_NH3.F_vol # kg SO4 / hr
        
        # !!! add resin and H2SO4 as influent streams and spent resin as effluent
        resin_demand_influent = (waste.TN / self.resin_lifetime / self.ad_density / 14) # kg resin / m3 treated
        resin_demand_time = resin_demand_influent * waste.F_vol # kg resin / hr
        resin_cost_day = resin_demand_time * 24 * self.cost_resin # $ resin / d
        resin_in.imass['Polystyrene'] = resin_demand_time
        resin_out.imass['Polystyrene'] = resin_demand_time

        acid_demand_influent = waste.TN * self.vol_H2SO4 / self.ad_density / 14 # L acid / L treated
        acid_demand_time = acid_demand_influent * waste.F_vol * 1000 * 1.83 # kg acid / hr
        acid_cost_day = acid_demand_time * 24 / 1000 * self.cost_H2SO4 # $ acid / d
        H2SO4.imass['H2SO4'] = acid_demand_time

        # set values needed for _design and _cost as attributes
        self.volume_treated = waste.F_vol * 1000 * 24 # L liq / d 
        

        
     
    #_design will include all the construction or captial impacts  
    def _design(self):
        design = self.design_results
        # defining the quantities of materials/items
        # note that these items to be to be in the _impacts_items.xlsx
        
        self.quantity_columns = np.ceil(self.volume_treated / 
                                        self.column_daily_loading_rate) # number of 0.4 m columns
        design['PVC'] = PVC_quant = (self.column_length * self.quantity_columns 
                                     * self.pvc_mass) # kg PVC
        Tubing_quant = (self.tubing_length * self.quantity_columns 
                                           * self.tubing_mass) # kg PE
        Tank_quant = (self.quantity_columns * self.tank_mass / 3) # number of tanks with one tank for three columns
        design['PE'] = PE_quant = Tubing_quant + Tank_quant

        
        self.construction = (
            Construction(item='PVC', quantity = PVC_quant, quantity_unit = 'kg'),
            Construction(item='PE', quantity = PE_quant, quantity_unit = 'kg'),
            )
        self.add_construction()
        
    
    #_cost based on amount of steel and stainless plus individual components
    def _cost(self):
        #purchase_costs is used for capital costs
        #can use quantities from above (e.g., self.design_results['StainlessSteel'])
        #can be broken down as specific items within purchase_costs or grouped (e.g., 'Misc. parts')
        self.baseline_purchase_costs['PVC'] = (self.cost_PVC_column * self.column_length 
                                      * self.quantity_columns)
        self.baseline_purchase_costs['Tubing'] = (self.cost_tubing * self.tubing_length 
                                         * self.quantity_columns)
        self.baseline_purchase_costs['Tank'] = (self.quantity_columns * self.tank_cost / 3) # one tank for three columns
        self.F_BM = dict.fromkeys(self.baseline_purchase_costs.keys(), 1)
        
        #certain parts need to be replaced based on an expected lifefime
        #the cost of these parts is considered along with the cost of the labor to replace them
        ix_replacement_parts_annual_cost = 0 # USD/yr only accounts for time running
        
        ix_annual_maintenance = 0 #USD/yr only accounts for time running
        
        self.add_OPEX =  (ix_replacement_parts_annual_cost + ix_annual_maintenance) / (365 * 24) # USD/hr (all items are per hour)
        
        # costs associated with full time opperators can be added in the TEA as staff
      




       
