#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 20 10:52:16 2021

@author: lewisrowles
"""

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
from . import SludgeSeparator
from qsdsan.utils.loading import load_data, data_path
import os
__all__ = ('ScrewPress',)

su_data_path = os.path.join(data_path,'sanunit_data/_screw_press.tsv')


### 
class ScrewPress(SludgeSeparator):
    '''
    Screw Press is used for dewatering where sludge, conditioned with cationic 
    polymer, is fed into the unit. Sludge is continuously dewatered as it travels
    along the screw.
    
    Reference documents
    -------------------
    Tchobanoglous, G.; Stensel, H. D.; Tsuchihashi, R.; Burton, F.; Abu-Orf, 
    M.; Bowden, G.; Pfrang, W. Wastewater Engineering: Treatment and Resource 
    Recovery, 5th ed.; Metcalf & Eddy, Inc., AECOM, McGraw-Hill: New York, 2014.
    
    Parameters
    ----------
    ins : WasteStream
        Waste for treatment, e.g., wastewater or latrine sludge
        
        SanStream
        ***should polymer addition also be included as a stream?
        make polymer a sanstream 
        
        
    outs : WasteStream
        Liquids, solids.


        
    References
    ----------
    .. N/A
    
    '''
    

    def __init__(self, ID='', ins=None, outs=(),thermo=None, init_with='WasteStream',
                 split=None, settled_frac=None,
                 if_N2O_emission=False, **kwargs):
        
        SludgeSeparator.__init__(self, ID, ins, outs, thermo, init_with, 
                                 split, settled_frac, F_BM_default=1)
        self.price_ratio = 1 
        
# load data from csv each name will be self.name    
        data = load_data(path=su_data_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data
        
        for attr, value in kwargs.items():
            setattr(self, attr, value)


        #**** how can I use this distribution for solids capture? 
        #change name in csv to settled_frac
        self.dewatering_solids_capture = self.settled_frac

        
# define the number of influent and effluent streams    
    _N_ins = 2
    _N_outs = 2
    
# in _run: define influent and effluent streams and treatment processes 
    def _run(self):
        waste, polymer = self.ins
        liq, cake_sol = self.outs
        SludgeSeparator._run(self)
        
        #*** is this needed in both places???? Will it be pulled into the function? Isnt cake solids TS needed too?
        self.dewatering_solids_capture = self.settled_frac
        liq, cake_sol = self._adjust_solid_water(waste, liq, cake_sol)
        
        # need the flowrate of the liquid stream or use equations below with dewatering_solids_flowrate
        self.dewatering_solids_flowrate = liq.F_vol
        self.waste_sol_flow = cake_sol.F_mass
        
        #breakpoint()
        # note to self change these once confirm outputs from above function...
        
        # Does the funtion above for self._adjust_solid_water take care of these equations? 
        # self.waste_sol_flow = waste.F_mass # kg TS/hr
        # dewatering_solids_conc = (self.waste_sol_flow) * self.dewatering_solids_capture    # kg TS/hr
        # dewatering_solids_flowrate = dewatering_solids_conc / (1.06 * self.cake_solids_TS * 1000) # m3 / hr
        # liq.F_vol = self.dewatering_centrate_flowrate = (waste.F_vol - dewatering_solids_flowrate)   # m3 / hr
        # dewatering_centrate_conc = (self.waste_sol_flow) * (1 - self.dewatering_solids_capture) # kg / d
        # dewatering_total_solids = dewatering_solids_conc / dewatering_solids_flowrate # kg / m3
        
        waste_TS = waste.F_mass * (1 - (waste.imass['H2O'] / waste.F_mass))
        self.polymer_demand = self.dewatering_polymer_dose * waste_TS / 1000 # kg polymer / hr
        polymer.imass['Polyacrylamide'] = self.polymer_demand 
        
        
        
        
        
    #_design will include all the construction or captial impacts  
    def _design(self):
        design = self.design_results
        
        # defining the quantities of materials/items
        # note that these items to be to be in the _impacts_items.xlsx
        # add function to calculate the number of screw presses required based on 
        # influent flowrate 
        design['Steel'] = S_quant = (self.dewatering_screw_press_steel)

        
        self.construction = (
            Construction(item='Steel', quantity = S_quant, quantity_unit = 'kg'),
            )
        self.add_construction()
        

    #_cost based on amount of steel and stainless plus individual components
    def _cost(self):
        
        #purchase_costs is used for capital costs
        #can use quantities from above (e.g., self.design_results['StainlessSteel'])
        #can be broken down as specific items within purchase_costs or grouped (e.g., 'Misc. parts')
        self.baseline_purchase_costs['Screw Press'] = (self.dewatering_screw_press_cost) * self.price_ratio
        
        
        #dewatered_water_treatment_cost = self.dewatering_solids_flowrate * self.wwt_cost
        #OM_costs = (dewatered_water_treatment_cost) # $/hr


        #need to be a cost per hour
        #self.add_OPEX =  0 # USD/hr (all items are per hour)

       
        power_demand = (self.dewatering_energy_demand * (self.waste_sol_flow)) #kW (all items are per hour)

        self.power_utility(power_demand)
        #breakpoint()
        



    


