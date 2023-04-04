#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

This is a pH solver, currently under development.

@author: Jianan Feng <jiananf2@illinois.edu>
"""

import sympy as sym, numpy as np
from math import log, sqrt, lcm
from qsdsan.utils.ph_chemical_inventory import chemical_inventory, precipitation_inventory

__all__ = ('pH_solver',)

def pH_solver(kw=10**-14,
              activity=True,
              check_precipitation = True,
              weak_chemicals={},
              other_chemicals={}):
    '''
    weak_chemicals={'H3PO4':0.5,
                      'NH4':0.5}

    other_chemicals = {
     (('Mg', 2),): (0.5,)
      }
    
    e.g.,
    from qsdsan.utils.ph import pH_solver
    pH_solver(precise=6, activity=True,
              weak_chemicals={'PO4':0.05, 'NH3':0.05},
              other_chemicals={(('Cl', -1),): (0.05,), (('Mg', 2),): (0.05,)})
    
    '''
    
    if weak_chemicals == {} and other_chemicals == {}:
        raise ValueError('no chemical is given')
    
    chemicals = {}
    
    for chemical in weak_chemicals.items():
        try:
            chemicals[list(chemical_inventory[chemical[0]].keys())[0]] = list(chemical_inventory[chemical[0]].values())[0]+(chemical[1],)
        except KeyError:
            raise KeyError(f'{chemical[0]} is not in the chemical inventory, check the inventory or manually import the chemical')

    for chemical in other_chemicals.items():
        if chemical[0][0][0] == 'SO4':
            raise ValueError('SO4 is treated as HSO4, which is a weak acid, please add in weak_chemicals')
        chemicals[chemical[0]] = chemical[1]

    H = sym.symbols('H')
    OH = sym.symbols('OH')

    eqn_list=[]
    eqn_dict={}
    eqn_dict['H'] = H
    eqn_dict['OH'] = OH
    eqn_list.append(sym.Eq(H*OH, kw))

    for chemical in list(chemicals):
        for i in range(len(chemical)):
            eqn_dict[chemical[i][0]] = sym.symbols(f'{chemical[i][0]}')
            
    # mass balance
            try:
                MB+=eqn_dict[chemical[i][0]]
            except NameError:
                MB=eqn_dict[chemical[i][0]]
        eqn_list.append(sym.Eq(MB, chemicals[chemical][-1]))
        MB-=MB

    # equilibrium
        if len(chemical) > 1:
            for i in range(len(chemical)-1):
                eqn_list.append(sym.Eq(eqn_dict[chemical[i+1][0]]*H - chemicals[chemical][i]*eqn_dict[chemical[i][0]],0))
                
    # charge balance
    chemical_ion = {}
    for chemical in chemicals.keys():
        for i in chemical:
            chemical_ion[i[0]] = i[1]
    chemical_ion['H'] = 1
    chemical_ion['OH'] = -1
    tot_charge = 0
    for chemical in chemical_ion.keys():
        tot_charge += chemical_ion[chemical]*eqn_dict[chemical]
    eqn_list.append(sym.Eq(tot_charge, 0))

    # for i in [10**-2, 10**-1, 1, 10]: # TODO what if i=9 and leave the loop? check?
    #     try:
    #         ans = sym.nsolve(eqn_list, tuple(eqn_dict.values()), [i]*len(eqn_dict), dict=True, maxsteps=100)
    #     except Exception:
    #         pass
    #     else:
    #         if all(i>=0 for i in list(ans[0].values())) == True:
    #             break
    # ans = sym.nsolve(eqn_list, tuple(eqn_dict.values()), list(ans[0].values()), dict=True, maxsteps=100) # repeat once to increase precise
    ans = sym.nsolve(eqn_list, tuple(eqn_dict.values()), [1]*len(eqn_dict), dict=True, maxsteps=100)
    
    if not check_precipitation:
    
        pH = -log(ans[0][sym.symbols('H')], 10)
        
        print(f'\npH: {pH:.2f}\n\nions in the system: {ans}')
    
    else:   
        # get all ions   
        ions = {}
        for item in ans[0].items():
            ions[str(item[0])] =item[1]
        ions['OH'] = kw/ions['H']
    
        all_precipitates = []
        ions, precipitate, chemicals, check = precipitation_iterator(ions, chemicals, chemical_ion, all_precipitates, kw)
        while check == 1:
            all_precipitates.append(precipitate)
            ions, precipitate, chemicals, check = precipitation_iterator(ions, chemicals, chemical_ion, all_precipitates, kw) # !!! keep being the first precipitate, start from here (for 4/2/2023)
        
        pH = -log(ions['H'], 10)
    
        ions_round = {}
        for item in ions.items():
            ions_round[item[0]] = precision_round(float(item[1]), digits=2)
    
        print(f'\npH: {pH:.2f}\n\nions in the system: {ions_round}\n\nprecipitates: {all_precipitates}')
    
        # return pH, ions, all_precipitates

def precipitation_iterator(ions, chemicals, chemical_ion, existed_precipitate, kw):
    
        # separate ions to cations and anions and also remove neutral ones
        cations = {}
        anions = {}
        for ion in ions.items():
            if chemical_ion[ion[0]] > 0:
                cations[ion[0]] = ion[1]
            else:
                anions[ion[0]] = ion[1]

        eqn_list=[sym.Eq(sym.symbols('OH')*sym.symbols('H'), kw),]
        eqn_dict={
            'H': sym.symbols('H'),
            'OH': sym.symbols('OH'),
                  }
        
        precipitation_SI = {}
        
        # identify possible precipitations:
        for cation in cations.items():
            if cation[0] in list(precipitation_inventory.keys()):
                for anion in anions.items():
                    if anion[0] in precipitation_inventory[cation[0]]:
                        index = precipitation_inventory[cation[0]].index(anion[0])
                        ksp = precipitation_inventory[cation[0]][index+1]
                        
                        least_common = lcm(abs(chemical_ion[cation[0]]), abs(chemical_ion[anion[0]]))
                        cation_coef = int(least_common/abs(chemical_ion[cation[0]]))
                        anion_coef = int(least_common/abs(chemical_ion[anion[0]]))
                        
                        if (ions[cation[0]]**cation_coef)*(ions[anion[0]]**anion_coef) > 1.01*ksp: # allow 1% calculation error, consistent with below sym.nsolve until relative error is smaller than 0.01
                            # storage saturation index (SI) = log(K/Ksp), where K is the concentratilon products (not consider activity for now)
                            SI = log((ions[cation[0]]**cation_coef)*(ions[anion[0]]**anion_coef)/ksp, 10)
                            if f'({cation[0]}){cation_coef}({anion[0]}){anion_coef}' not in existed_precipitate:
                                precipitation_SI[sym.symbols(f'({cation[0]}){cation_coef}({anion[0]}){anion_coef}')] = [(cation[0], cation_coef), (anion[0], anion_coef), ksp, SI]
        
        # breakpoint()
        for ion in ions.items():
            if ion[1] < 0:
                ions[ion[0]] = 0
        # breakpoint()
        if len(precipitation_SI) == 0:
            return ions, '', {}, 0
        else:
            
          
            
            
            # substract from the MB, but just for the ones with the largest SI
            SI_list = [precipitation[1][-1] for precipitation in precipitation_SI.items()]
            first_precipitation_index = SI_list.index(max(SI_list))
            first_precipitation = list(precipitation_SI.values())[first_precipitation_index]
            precipitate_name = f'({first_precipitation[0][0]}){first_precipitation[0][1]}({first_precipitation[1][0]}){first_precipitation[1][1]}'

            eqn_list.append(sym.Eq((sym.symbols(first_precipitation[0][0])**first_precipitation[0][1])*(sym.symbols(first_precipitation[1][0])**first_precipitation[1][1]), first_precipitation[-2]))
            eqn_dict[f'({first_precipitation[0][0]}){first_precipitation[0][1]}({first_precipitation[1][0]}){first_precipitation[1][1]}'] = sym.symbols(f'({first_precipitation[0][0]}){first_precipitation[0][1]}({first_precipitation[1][0]}){first_precipitation[1][1]}')
            
            for chemical in chemicals.items():
                for ion in first_precipitation[:-2]:
                    if ion[0] in list(dict(chemical[0]).keys()):
                        new_tot = chemicals[chemical[0]][-1] - ion[1]*sym.symbols(precipitate_name)
                        tempo_list = list(chemicals[chemical[0]])
                        tempo_list[-1] = new_tot
                        tempo_dict = {chemical[0]: tuple(tempo_list)}
                        chemicals.update(tempo_dict)
                                
            # TODO forget activitiy for precipitation for now
    
            for chemical in list(chemicals):
                for i in range(len(chemical)):
                    eqn_dict[chemical[i][0]] = sym.symbols(f'{chemical[i][0]}')
                    
            # mass balance
                    try:
                        MB+=eqn_dict[chemical[i][0]]
                    except NameError:
                        MB=eqn_dict[chemical[i][0]]
                eqn_list.append(sym.simplify(sym.Eq(MB, chemicals[chemical][-1])))
                MB-=MB
    
            # equilibrium
                if len(chemical) > 1:
                    for i in range(len(chemical)-1):
                        eqn_list.append(sym.Eq(eqn_dict[chemical[i+1][0]]*sym.symbols('H') - chemicals[chemical][i]*eqn_dict[chemical[i][0]],0))
    
            chemical_ion = {}
            for chemical in chemicals.keys():
                for i in chemical:
                    chemical_ion[i[0]] = i[1]
            chemical_ion['H'] = 1
            chemical_ion['OH'] = -1
            tot_charge = 0
            for chemical in chemical_ion.keys():
                tot_charge += chemical_ion[chemical]*eqn_dict[chemical]
            eqn_list.append(sym.Eq(tot_charge, 0))
            
            dict_i_relative_error = {}
            
            # how about except precipitate, other initial guesses follow the previous results?
            
            # for i in [10**-10, 10**-9, 10**-8, 10**-7, 10**-6, 10**-5, 10**-4, 10**-3, 10**-2, 10**-1, 1, 10]: # TODO what if i=9 and leave the loop? check?
            for i in [10**-5, 10**-2.5, 1]: # TODO what if i=9 and leave the loop? check?
                try:
                    ans = sym.nsolve(eqn_list, tuple(eqn_dict.values()), [i]*len(eqn_dict), dict=True, maxsteps=100)
                except Exception:
                    pass
                else:
                    # H and OH should be strictly larger than 0  
                    if all(j >= -0.05 for j in list(ans[0].values())) == True: # !!! is it reasonable to set -0.05 here as tolerance? we may also need to change this values to 0 in the final answer?

                    # add priority: what if some have negative values but less error?
                    
                        if ans[0][sym.symbols('H')] > 0:
                            if ans[0][sym.symbols('OH')] > 0:
                                # store relative error and select the minimum one to iterate below
                                dict_i_relative_error[i] = abs(((ans[0][sym.symbols(first_precipitation[0][0])]**first_precipitation[0][1]*ans[0][sym.symbols(first_precipitation[1][0])]**first_precipitation[1][1])-first_precipitation[-2])/first_precipitation[-2])
                                if abs(((ans[0][sym.symbols(first_precipitation[0][0])]**first_precipitation[0][1]*ans[0][sym.symbols(first_precipitation[1][0])]**first_precipitation[1][1])-first_precipitation[-2])/first_precipitation[-2]) <= 0.01:
                                    break
                        
                        
            min_error_index = list(dict_i_relative_error.values()).index(min(list(dict_i_relative_error.values())))
            
            i = list(dict_i_relative_error.keys())[min_error_index]
            
            ans = sym.nsolve(eqn_list, tuple(eqn_dict.values()), [i]*len(eqn_dict), dict=True, maxsteps=100)

            # iteration_times = 0
            # while abs(((ans[0][sym.symbols(first_precipitation[0][0])]**first_precipitation[0][1]*ans[0][sym.symbols(first_precipitation[1][0])]**first_precipitation[1][1])-first_precipitation[-2])/first_precipitation[-2]) > 0.01:
            ans = sym.nsolve(eqn_list, tuple(eqn_dict.values()), list(ans[0].values()), dict=True, maxsteps=100)
                # iteration_times += 1
                # if iteration_times == 5:
                #     break

            # update chemicals
            for chemical in chemicals.items():
                for ion in first_precipitation[:-2]:
                    if ion[0] in list(dict(chemical[0]).keys()):
                        new_tot = max(0, chemicals[chemical[0]][-1] + ion[1]*sym.symbols(precipitate_name) - ion[1]*ans[0][sym.symbols(precipitate_name)])
                        tempo_list = list(chemicals[chemical[0]])
                        tempo_list[-1] = new_tot
                        tempo_dict = {chemical[0]: tuple(tempo_list)}
                        chemicals.update(tempo_dict)
            
            ans[0].pop(sym.symbols(precipitate_name))

            ions = {}
            for item in ans[0].items():
                ions[str(item[0])] = item[1]
            
            return ions, precipitate_name, chemicals, 1  # TODO when to stop use this iterator and return final pH?
        
def precision_round(number, digits=2):
    power = '{:e}'.format(number).split('e')[1]
    return round(number, -(int(power) - digits))