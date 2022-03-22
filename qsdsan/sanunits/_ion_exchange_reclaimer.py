
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
from qsdsan.utils.loading import load_data, data_path
__all__ = ('IonExchangeReclaimer',)

data_path += 'sanunit_data/_ion_exchange_reclaimer.csv'

#To scale the system from 0 - 120 users based on units corresponding to 30 users: 
ppl = 120

if (ppl <=30): R = 1 
elif (ppl >30 and ppl <=60): R = 2
elif (ppl >60 and ppl <=90): R = 3
elif (ppl >90 and ppl <=120): R = 4
else: R = 5

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
        
        SanUnit.__init__(self, ID, ins, outs, F_BM_default=1)
        self.if_gridtied = if_gridtied
        self.price_ratio = 1
  
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
        zeolite_demand_time = (self.zeolite_weight / (self.zeolite_lifetime*365*24)) * R # kg zeolite / hr
        zeolite_in.imass['Zeolite'] = zeolite_demand_time
        zeolite_out.imass['Zeolite'] = zeolite_demand_time
        
        # GAC
        gac_demand_time = ((self.gac_weight * self.gac_annual_replacement) / (365*24)) * R # kg gac / hr
        gac_in.imass['GAC'] = gac_demand_time
        gac_out.imass['GAC'] = gac_demand_time
        
        # KCl
        self.KCl_demand_time = ((self.KCl_weight * self.regen_freq_per_yr) / (365*24)) * R # kg KCl / hr
        KCl.imass['PotassiumChloride'] = self.KCl_demand_time

        self.N_removed = waste.imass['NH3'] * (self.TN_removal)
        self.N_recovered = self.N_removed * (self.desorption_recovery_efficiency)  # kg N / hr
        treated.imass['NH3'] =  waste.imass['NH3'] - self.N_removed # kg N / hr
        conc_NH3.imass['NH3'] = self.N_recovered # kg N / hr
        conc_NH3.imass['PotassiumChloride'] = self.KCl_demand_time
        zeolite_out.imass['NH3'] = self.N_removed - self.N_recovered
        
        
    def _design(self):
        design = self.design_results
        
        design['Plastic'] = P_quant = self.Plastic_weight * R
        design['PVC'] = PVC_quant = self.PVC_weight * R
        design['Steel'] = S_quant = self.Steel_weight * R
        
        self.construction = (
            Construction(item='Plastic', quantity = P_quant, quantity_unit = 'kg'),
            Construction(item='PVC', quantity = PVC_quant, quantity_unit = 'kg'),
            Construction(item='Steel', quantity = S_quant, quantity_unit = 'kg'),
            )
        self.add_construction(add_cost=False)        
 
     
    def _calc_labor_cost(self):
        labor_cost = (self.wages * self.labor_maintenance_zeolite_regeneration) * R
        return labor_cost / (365 * 24) # USD/hr (all items are per hour)
    
    def _calc_replacement_cost(self):
        ion_exchange_replacement_cost = (((self.Zeolite_cost * self.gac_annual_replacement)
                                         + (self.KCl_cost * self.KCl_weight + self.regen_freq_per_yr) 
                                         + (self.GAC_cost * self.gac_annual_replacement)
                                         + (self.ion_exchange_replacement_other_parts /10))) * R #USD/yr
        return ion_exchange_replacement_cost/ (365 * 24)  * self.price_ratio # USD/hr (all items are per hour)
                 
    
    def _cost(self):
        
        C = self.baseline_purchase_costs
        C['Pipes'] = (self.four_in_pipe_SCH40 + self.four_in_pipe_SCH80) * R
        C['fittings'] = (self.four_in_pipe_SCH80_endcap + self.NRV + self.connector + 
        self.ball_valve + self.three_eight_elbow + self.ten_ten_mm_tee + self.OD_tube + self.four_in_pipe_clamp) * R
                                           
        C['GAC_Zeolite'] = (self.GAC_zeolite_mesh + self.GAC_cost + self.Zeolite_cost) * R
        C['Regeneration Solution'] = (self.KCl_cost * self.KCl_weight) * R
            
        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio
     
        add_OPEX = self._calc_replacement_cost() 
        self._add_OPEX = {'Additional OPEX': add_OPEX}
    
    




       
