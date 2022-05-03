# !/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Tori Morgan <tvlmorgan@gmail.com>
    Hannah Lohman <hlohman94@gmail.com>
    Lewis Rowles <stetsonsc@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''


import numpy as np
from qsdsan import SanUnit, Construction
from ..utils import ospath, load_data, data_path

__all__ = ('IonExchangeReclaimer',)

ion_exchange_data_path = ospath.join(data_path, 'sanunit_data/_ion_exchange_reclaimer.csv')


class IonExchangeReclaimer(SanUnit):
    '''
    Ion Exchange for N recovery from liquid stream. Concentrated NH3 is recovered.

    The following impact items should be pre-constructed for life cycle assessment:
    Plastic, PVC, Steel.

    Parameters
    ----------
    ppl: int
        Total number of users for scaling of costs.

    ins:
        waste: liquid waste stream to be treated by ion exchange unit
        zeolite_in: zeolite input
        gac_in: GAC input
        KCl: KCl input

    outs:
        treated: treated liquid leaving ion exchange unit
        zeolite_out: spent zeolite
        gac_out: spent GAC
        conc_NH3: concentrated NH3

    References
    ----------
    [1] Lohman et al., Advancing Sustainable Sanitation and Agriculture
    through Investments in Human-Derived Nutrient Systems. 
    Environ. Sci. Technol. 2020, 54, (15), 9217-9227.
    https://dx.doi.org/10.1021/acs.est.0c03764
    
    [2] Tarpeh et al., Evaluating ion exchange for nitrogen recovery from
    source-separated urine in Nairobi, Kenya. Development Engineering. 2018, 
    3, 188–195.
    https://doi.org/10.1016/j.deveng.2018.07.002

    [3] Duke Reclaimer team data
    
    '''

    # Constants
    baseline_ppl = 30  # baseline population served by Reclaimer
    exponent_scale = 0.6  # exponential scaling constant

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream', ppl=1, **kwargs):
        
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with, F_BM_default=1)

        self.ppl = ppl

        self.qty_reclaimers = np.ceil(self.ppl / self.baseline_ppl)  # number of reclaimer units required
        self.qty_reclaimers = self.qty_reclaimers.astype(int)  # convert from float to integer
  
        data = load_data(path=ion_exchange_data_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)

    _N_ins = 4
    _N_outs = 4

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
        
        # Zeolite
        zeolite_demand_time = (self.zeolite_weight / (self.zeolite_lifetime * 365 * 24)) * self.qty_reclaimers  # kg zeolite/hr
        zeolite_in.imass['Zeolite'] = zeolite_demand_time
        zeolite_out.imass['Zeolite'] = zeolite_demand_time
        
        # GAC
        gac_demand_time = (self.gac_weight / (self.gac_lifetime * 365 * 24)) * self.qty_reclaimers  # kg GAC/hr
        gac_in.imass['GAC'] = gac_demand_time
        gac_out.imass['GAC'] = gac_demand_time
        
        # KCl
        self.KCl_demand_time = (self.KCl_weight / (self.KCl_regeneration_freq * 365 * 24)) * self.qty_reclaimers  # kg KCl/hr
        KCl.imass['PotassiumChloride'] = self.KCl_demand_time

        self.N_removed = waste.imass['NH3'] * self.TN_removal
        self.N_recovered = self.N_removed * self.desorption_recovery_efficiency  # kg N / hr
        treated.imass['NH3'] = waste.imass['NH3'] - self.N_removed  # kg N / hr
        conc_NH3.imass['NH3'] = self.N_recovered  # kg N / hr
        conc_NH3.imass['PotassiumChloride'] = self.KCl_demand_time
        zeolite_out.imass['NH3'] = self.N_removed - self.N_recovered

    def _design(self):
        design = self.design_results
        
        design['Plastic'] = P_quant = self.Plastic_weight * (self.ppl / self.baseline_ppl)  # linear scale
        design['PVC'] = PVC_quant = self.PVC_weight * (self.ppl / self.baseline_ppl)  # linear scale
        design['Steel'] = S_quant = self.Steel_weight * (self.ppl / self.baseline_ppl)  # linear scale
        
        self.construction = (
            Construction(item='Plastic', quantity=P_quant, quantity_unit='kg'),
            Construction(item='PVC', quantity=PVC_quant, quantity_unit='kg'),
            Construction(item='Steel', quantity=S_quant, quantity_unit='kg'),
            )
        self.add_construction(add_cost=False)        

    def _cost(self):
        
        C = self.baseline_purchase_costs
        C['Pipes'] = (self.four_in_pipe_SCH40 + self.four_in_pipe_SCH80)
        C['Fittings'] = (self.four_in_pipe_SCH80_endcap + self.NRV + self.connector +
                         self.ball_valve + self.three_eight_elbow + self.ten_ten_mm_tee + self.OD_tube +
                         self.four_in_pipe_clamp)
                                           
        C['GAC_Zeolite'] = self.GAC_zeolite_mesh + (self.GAC_cost * self.gac_weight) + (self.Zeolite_cost * self.zeolite_weight)
        C['Regeneration Solution'] = (self.KCl_cost * self.KCl_weight)

        # Exponentially scale capital cost with number of users
        scale = (self.ppl / self.baseline_ppl) ** self.exponent_scale
        for equipment, cost in C.items():
            C[equipment] = cost * scale

        self.add_OPEX = self._calc_replacement_cost() + self._calc_maintenance_labor_cost()

    def _calc_replacement_cost(self):
        scale = (self.ppl / self.baseline_ppl) ** self.exponent_scale
        zeolite_replacement_cost = self.Zeolite_cost * self.zeolite_weight / self.zeolite_lifetime  # USD/year
        gac_replacement_cost = self.GAC_cost * self.gac_weight / self.gac_lifetime  # USD/year
        kcl_replacement_cost = self.KCl_cost * self.KCl_weight / self.KCl_regeneration_freq  # USD/year
        # replacement parts for everything else, assume 2-6% of capital as annual cost
        other_replacement_cost = self.om_capital_ratio * (self.four_in_pipe_SCH40 + self.four_in_pipe_SCH80 +
                                                          self.four_in_pipe_SCH80_endcap + self.NRV + self.connector +
                                                          self.ball_valve + self.three_eight_elbow +
                                                          self.ten_ten_mm_tee + self.OD_tube + self.four_in_pipe_clamp +
                                                          self.GAC_zeolite_mesh)  # USD/year
        ion_exchange_replacement_cost = scale * (zeolite_replacement_cost + gac_replacement_cost + kcl_replacement_cost
                                                 + other_replacement_cost) / (365 * 24)  # USD/hr
        return ion_exchange_replacement_cost

    def _calc_maintenance_labor_cost(self):
        scale = (self.ppl / self.baseline_ppl) ** self.exponent_scale
        labor_cost = (self.wages * self.labor_maintenance_zeolite_regeneration) * scale  # USD/year
        return labor_cost / (365 * 24)  # USD/hr
