#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 21 19:55:59 2021

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
from qsdsan.utils.loading import load_data, data_path
import os

__all__ = ('DryerFromHHX',)

#path to csv with all the inputs

su_data_path = os.path.join(data_path, 'sanunit_data/_dryer_from_hhx.tsv')
### 
class DryerFromHHX(SanUnit):
    '''
    Dryer is used in combination with hydronic heat exchanger.

    
    Reference documents
    -------------------
    N/A
    
    Parameters
    ----------
    ins : WasteStream
        Dewatered solids, heat.
    outs : WasteStream 
        Dried solids, fugitive N2O, fugitive CH4.
            set both as WasteStream 
            others could be:
            Gases consist of SO2_emissions, NOx_emissions, CO_emissions, 
            Hg_emissions, Cd_emissions, As_emissions, Dioxin_Furans_emissions.

        
    References
    ----------
    .. N/A
    
    '''
    

    def __init__(self, ID='', ins=None, outs=(), moisture_content_out=0.35, **kwargs):
        
        SanUnit.__init__(self, ID, ins, outs, F_BM_default=1)
        self.moisture_content_out = moisture_content_out
        self.price_ratio = 1 

# load data from csv each name will be self.name    
        data = load_data(path=su_data_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data
        
        for attr, value in kwargs.items():
            setattr(self, attr, value)




        
# define the number of influent and effluent streams    
    _N_ins = 2
    _N_outs = 3

# in _run: define influent and effluent streams and treatment processes 
    def _run(self):
        waste_in, heat_in = self.ins
        waste_out, N2O, CH4  = self.outs
        waste_out.copy_like(self.ins[0])
        heat_in.phase = N2O.phase = CH4.phase = 'g'
        
        # calculate heat needed to dry to 35%
        moisture_content_in = (waste_in.imass['H2O'] / waste_in.F_mass) # fraction
        TS_in = waste_in.F_mass - waste_in.imass['H2O'] # kg TS dry/hr

        # set necessary moisture content of effluent as 35%
        waste_out.imass['H2O'] = ((self.moisture_content_out)/(moisture_content_in)) * TS_in # fraction
        water_to_dry = waste_in.imass['H2O'] - waste_out.imass['H2O'] # kg water/hr
        heat_needed_to_dry_35 = (water_to_dry) * self.energy_required_to_dry_sludge # MJ/hr
        
        #!!! heat_supplied isn't calculated from heat_in and heat loss?
        heat_supplied = self.dryer_heat_transfer_coeff * self.area_surface * ((self.water_out_temp + 273) - (self.feedstock_temp+ 273))

        # if heat_needed_to_dry_35 > heat_supplied:
        #     breakpoint()
        
        # set to use COD or C for carbon based on influent composition
        # if waste_in.COD == 0:
        #     C_in = waste_in.imass['C']
            
        # else:
        #     C_in = self.carbon_COD_ratio * waste_in.COD
        
        # emissions
        drying_CO2_to_air = (self.drying_CO2_emissions * self.carbon_COD_ratio
                             * waste_in.COD * waste_in.F_vol / 1000) # kg CO2 /hr

        # add conversion factor for COD to TC?
        drying_CH4_to_air = (self.drying_CH4_emissions * self.carbon_COD_ratio
                             * waste_in.COD * waste_in.F_vol / 1000) # kg CH4 /hr

        drying_NH3_to_air = self.drying_NH3_emissions * waste_in.TN * waste_in.F_vol / 1000 # kg NH3 /hr
        drying_N2O_to_air = drying_NH3_to_air * self.NH3_to_N2O # kg N2O /hr
        CH4.imass['CH4'] = drying_CH4_to_air
        N2O.imass['N2O'] = drying_N2O_to_air
        # reduce COD and TN in waste_out based on emissions
        
        waste_out._COD = (waste_in.COD * (waste_in.F_vol/waste_out.F_vol)) - ((drying_CO2_to_air + drying_CH4_to_air) / self.carbon_COD_ratio)
 

        # waste_out.imass['C'] -= ((drying_CO2_to_air + drying_CH4_to_air))

        waste_out.imass['N'] -= drying_NH3_to_air


        # jacket loss data and calculate losses due to convection and radiation
        self.jacket_heat_loss_conv = self.heat_transfer_coeff * (self.water_air_hx_temp1 - self.ambient_temp) * self.water_air_hx_area1 / 1000 / 0.2778 # MJ/hr
        self.jacket_heat_loss_radiation = self.radiative_emissivity * 5.67e-8 * ((self.water_air_hx_temp1 + 273)**4 - (self.ambient_temp + 273)**4) / 1000 / 0.2778 # MJ/hr
        self.jacket_heat_loss_sum = self.jacket_heat_loss_conv + self.jacket_heat_loss_radiation # MJ/hr
        #these can be calculated for pyrolysis and catalytic converter also 

        

    #_cost based on amount of steel and stainless plus individual components
    def _cost(self):
        
        
        #power_demand = (self.drying_energy_demand) #kW
        power_demand = 0
        self.power_utility(power_demand) #kW







