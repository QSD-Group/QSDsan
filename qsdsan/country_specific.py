#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import models as m
import pandas as pd
from systems import price_dct
from chaospy import distributions as shape
import biosteam as bst
from biosteam.evaluation import Model
from exposan import biogenic_refinery as br
from qsdsan.utils import (
    load_data, data_path,
    AttrSetter, AttrFuncSetter, DictAttrSetter,
    FuncGetter,
    time_printer
    )
from qsdsan import currency, ImpactItem

systems = br.systems
sys_dct = systems.sys_dct
price_dct = systems.price_dct
GWP_dct = systems.GWP_dct
GWP = systems.GWP

import os
su_data_path = os.path.join(data_path, 'sanunit_data/')

# add energy_price as a parameter to the new_model

input_dct = {
    'China': {'energy_GWP': ('uniform', 0.745, 0.712133646, 0.848557635), 
              'energy_price': 0.084,
#             'operator': 50.08,
              'wages': ('uniform', 6.26, 4.695, 7.825), 
              'construction': 31.12,
              'e_cal': 3191,
              'p_anim': 40,
              'p_veg': 60.63,
              'price_ratio': 0.610,
              'household_size': 3,
              # 'N_fertilizer_price': 0.94,
              # 'P_fertilizer_price': 1.74,
              # 'K_fertilizer_price': 1.12
               },
 #    'India': {'energy_GWP': ('uniform', 0.852, 0.805105241, 0.958663233), 
 #              'energy_price': 0.081,
 #              'operator': 29.12,
    #           'wages': ('uniform', 3.64, 3.094, 4.186), 
 #              'construction': 16.85,
 #              'e_cal': 2533,
 #              'p_anim': 15,
 #              'p_veg': 48.35,
 #              'price_ratio': 0.300,
 #              # 'household_size': 5,
 #              # 'N_fertilizer_price': 0.16,
 #              # 'P_fertilizer_price': 0.57,
 #              # 'K_fertilizer_price': 0.44
 #              },
 #    'South Africa': {'energy_GWP': ('uniform', 0.955, 0.909088322, 1.086711278),
 #              'energy_price': 0.14,
 #              'operator': 23.6,
    #           'wages': ('uniform', 2.95, 2.2125, 3.6875), 
 #              'construction': 11.76,
 #              'e_cal': 2899,
 #              'p_anim': 36.03,
 #              'p_veg': 48.33,
 #              'price_ratio': 0.460,
 #              # 'household_size': 3,
 #              # 'N_fertilizer_price': 0.81,
 #              # 'P_fertilizer_price': 0.78, #source: https://www.namc.co.za/wp-content/uploads/2021/07/Trends-in-selected-Agricultural-input-prices-June-2021.pdf
 #              # 'K_fertilizer_price': 0.87
 #              },
 #    'Senegal': {'energy_GWP': ('uniform', 0.939, 0.850772468, 1.016179461),  
 #              'energy_price': 0.186,
 #              'operator': 29.12,
    #           'wages': ('uniform', 3.64, 3.094, 4.186), 
 #              'construction': 16.85,
 #              'e_cal': 2545,
 #              'p_anim': 13.69,
 #              'p_veg': 48.67,
 #              'price_ratio': 0.408,
 # #              'household_size': 9,
 # #              'N_fertilizer_price': 1.40,
 # #              'P_fertilizer_price': 1.41, #source: https://doi.org/10.1371/journal.pone.0227764
 # # !             'K_fertilizer_price': 
 #              },
 #    'Uganda': {'energy_GWP': ('uniform', 0.159, 0.1113, 0.19875), 
 #              'energy_price': 0.184,
 #              'operator': 10.64,
 #              'wages': ('uniform', 1.33, 0.9975, 1.6625), 
 #              'construction': 4.80,
 #              'e_cal': 1981,
 #              'p_anim': 12.25,
 #              'p_veg': 34.69,
 #              'price_ratio': 0.348,
 #              # 'household_size': 5,
 #              # 'N_fertilizer_price': 1.79,
 #              # 'P_fertilizer_price': 3.97,
 #              # 'K_fertilizer_price': 1.33
 #              },
    }

b = price_dct['Electricity']

def run_country(dct):
    results = {}
    

    
    # all_paramsA = modelA.get_parameters()

    
    for country, country_dct in dct.items():
        sysA = systems.sysA
        sysA.simulate()
        streams = sys_dct['stream_dct'][sysA.ID]
        # operator
        #sysA._TEA.annual_labor = country_dct['operator']* 3*365
        
        modelA = Model(sysA, m.add_metrics(sysA))
        # paramA = modelA.parameter
        
        # Shared parameters
        modelA = m.add_shared_parameters(sysA, modelA, country_specific=True)
        param = modelA.parameter
        
        
        #price ratio
        i = country_dct['price_ratio']
        # i=1
        price_dct = systems.price_dct
        old_price_ratio = systems.price_ratio
        for stream in ('Concrete', 'Steel', 'Polymer', 'Resin', 'FilterBag', 'MgOH2', 'MgCO3', 'H2SO4', 'biochar'):
            # new_price = price_dct[stream] = price_dct[stream] * i/old_price_ratio 
            # streams[stream].price = new_price
            old_price = price_dct[stream]
            new_price = old_price * i/old_price_ratio
            if stream=='Concrete': print(f'\n old price is {old_price}, new price is {new_price}')
            streams[stream].price = new_price
        systems.price_ratio = i
        for u in sysA.units:
            if hasattr(u, 'price_ratio'):
                u.price_ratio = i
        
        # wages
        kind, low_val, peak_val, max_val = country_dct['wages']
        b=peak_val
        if kind == 'triangle':
            D = shape.Triangle(lower=low_val, midpoint=peak_val, upper=max_val)
        else:
            D = shape.Uniform(lower=low_val,upper=max_val)
        @param(name='Labor wages', element=systems.A1, kind='coupled', units='USD/h',
               baseline=b, distribution=D)
        def set_labor_wages(i):
            labor_cost = 0 
            sysA = systems.sysA
            for u in sysA.units:
                if hasattr(u, '_calc_maintenance_labor_cost'):
                    u.wages = i
                    labor_cost += u._calc_maintenance_labor_cost()
            sysA.TEA.annual_labor = labor_cost
        
        #energy_GWP
        kind, low_val, peak_val, max_val = country_dct['energy_GWP']
        b=peak_val
        if kind == 'triangle':
            D = shape.Triangle(lower=low_val, midpoint=peak_val, upper=max_val)
        else:
            D = shape.Uniform(lower=low_val,upper=max_val)
            
        @param(name='Electricity CF', element='LCA', kind='isolated',
                units='kg CO2-eq/kWh', baseline=b, distribution=D)
        def set_electricity_CF(i):
            GWP_dct['Electricity'] = ImpactItem.get_item('e_item').CFs['GlobalWarming'] = i
            
        # energy_price
        bst.PowerUtility.price = country_dct['energy_price']
        
        # household size
        unit = sysA.path[0]
        b = country_dct['household_size']
        D = shape.Normal(mu=b, sigma=1.8)
        @param(name='Household size', element=unit, kind='coupled', units='cap/household',
               baseline=b, distribution=D)
        def set_household_size(i):
            systems.household_size = max(1, i)

        
        # Diet and excretion
        A1 = systems.A1
        # p_anim
        A1.p_anim = country_dct['p_anim']
        
        # p_veg
        A1.p_veg = country_dct['p_veg']
        
        # e_cal
        A1.e_cal = country_dct['e_cal']
        
        path = su_data_path + '_excretion.tsv'
        data = load_data(path)
        m.batch_setting_unit_params(data, modelA, A1, exclude=('e_cal','p_anim','p_veg')) #adds distributions to tsv file
        
        # Pit latrine and conveyance
        # modelA = m.add_pit_latrine_parameters(sysA, modelA)
        
        # Industrial control panel
        A4 = systems.A4
        path = su_data_path + '_industrial_control_panel.tsv'
        data = load_data(path)
        m.batch_setting_unit_params(data, modelA, A4)
        
        # Housing biogenic refinery
        A5 = systems.A5
        
        # construction
        A5.const_wage = country_dct['construction']
        
        path = su_data_path + '_housing_biogenic_refinery.tsv'
        data = load_data(path)
        m.batch_setting_unit_params(data, modelA, A5)
        
        # Screw press
        A6 = systems.A6
        path = su_data_path + '_screw_press.tsv'
        data = load_data(path)
        m.batch_setting_unit_params(data, modelA, A6)
        
        # Liquid treatment bed
        A7 = systems.A7
        path = su_data_path + '_liquid_treatment_bed.tsv'
        data = load_data(path)
        m.batch_setting_unit_params(data, modelA, A7)
        
        # Carbonizer base
        A8 = systems.A8
        path = su_data_path + '_carbonizer_base.tsv'
        data = load_data(path)
        m.batch_setting_unit_params(data, modelA, A8)
        
        # Pollution control device
        A9 = systems.A9
        path = su_data_path + '_pollution_control_device.tsv'
        data = load_data(path)
        m.batch_setting_unit_params(data, modelA, A9)
        
        # Oil heat exchanger
        A10 = systems.A10
        path = su_data_path + '_oil_heat_exchanger.tsv'
        data = load_data(path)
        m.batch_setting_unit_params(data, modelA, A10)
        
        # Hydronic heat exchanger
        A11 = systems.A11
        path = su_data_path + '_hydronic_heat_exchanger.tsv'
        data = load_data(path)
        m.batch_setting_unit_params(data, modelA, A11)
        
        # Dryer from HHx
        A12 = systems.A12
        path = su_data_path + '_dryer_from_hhx.tsv'
        data = load_data(path)
        m.batch_setting_unit_params(data, modelA, A12)

        

        results[country] = m.run_uncertainty(model=modelA, seed=5, N=100) #!!! switch back to 1000 for final 
        del sysA, modelA
        
    return results


results = run_country(dct=input_dct)


#need results folder in the file where this is located
def save_uncertainty_results(results):
    import os
    path = os.path.dirname(os.path.realpath(__file__))
    path = '/opt/anaconda3/envs/op-model/lib/python3.8/site-packages/exposan/biogenic_refinery'
    path += '/results'
    if not os.path.isdir(path):
         os.mkdir(path)
    del os


    for country, dct in results.items():
        file_name = path+'/'+country+'.xlsx'
        if dct['parameters'] is None:
            raise ValueError('No cached result, run model first.')
        with pd.ExcelWriter(file_name) as writer:
            dct['parameters'].to_excel(writer, sheet_name='Parameters')
            dct['data'].to_excel(writer, sheet_name='Uncertainty results')
            if 'percentiles' in dct.keys():
                dct['percentiles'].to_excel(writer, sheet_name='Percentiles')
            dct['spearman'].to_excel(writer, sheet_name='Spearman')
            # model.table.to_excel(writer, sheet_name='Raw data')

save_uncertainty_results(results)

 # b = systems.get_operator_daily_wage()
 #    D = shape.Triangle(lower=14.55, midpoint=b, upper=43.68)
 #    @param(name='Operator daily wages', element='TEA', kind='cost', units='USD/d',
 #          baseline=b, distribution=D)
 #    def set_operator_daily_wage(i):
 #        sys._TEA.annual_labor = i* 3*365 
        
        
 #        b = price_dct['Electricity']
 #    D = shape.Triangle(lower=0.04, midpoint=b, upper=0.1)
 #    @param(name='Electricity price', element='TEA', kind='isolated',
 #           units='$/kWh', baseline=b, distribution=D)
 #    def set_electricity_price(i):
 #        PowerUtility.price = i
        
        
 #        b = GWP_dct['Electricity']
 #    D = shape.Uniform(lower=0.106, upper=0.121)
 #    @param(name='Electricity CF', element='LCA', kind='isolated',
 #               units='kg CO2-eq/kWh', baseline=b, distribution=D)
 #    def set_electricity_CF(i):
 #        GWP_dct['Electricity'] = ImpactItem.get_item('E_item').CFs['GlobalWarming'] = i
        
        
