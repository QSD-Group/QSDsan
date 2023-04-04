#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

This is a pH solver, currently under development.

@author: Jianan Feng <jiananf2@illinois.edu>
"""

import sympy as sym, numpy as np
from math import log, sqrt, lcm
from qsdsan.utils.ph_chemical_inventory import chemical_inventory, precipitation_inventory
import time

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

    ans = sym.nsolve(eqn_list, tuple(eqn_dict.values()), [1]*len(eqn_dict), dict=True, maxsteps=100, verify=False)

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
            ions, precipitate, chemicals, check = precipitation_iterator(ions, chemicals, chemical_ion, all_precipitates, kw)
        
        pH = -log(ions['H'], 10)
    
        ions_round = {}
        for item in ions.items():
            ions_round[item[0]] = precision_round(float(item[1]), digits=2)
    
        print(f'\npH: {pH:.2f}\n\nions in the system: {ions_round}\n\nprecipitates: {all_precipitates}')
    
        return pH, ions, all_precipitates

def precipitation_iterator(ions, chemicals, chemical_ion, existed_precipitate, kw):
    
    guess = max(chemical[-1] for chemical in list(chemicals.values()))
    # separate ions to cations and anions and also remove neutral ones
    cations = {}
    anions = {}
    for ion in ions.items():
        if chemical_ion[ion[0]] > 0:
            cations[ion[0]] = ion[1]
        if chemical_ion[ion[0]] < 0:
            anions[ion[0]] = ion[1]

    eqn_list=[sym.Eq(sym.symbols('OH')*sym.symbols('H'), kw),]
    eqn_dict={
        'H': sym.symbols('H'),
        'OH': sym.symbols('OH'),
              }

    precipitation_SI = {}
    
    SI_struvite = 0
    if '(Mg)1(NH4)1(PO4)1' not in existed_precipitate:
        try:
            if cations['Mg']*cations['NH4']*anions['PO4'] > 2.5*10**-13:
                    SI_struvite = log(cations['Mg']*cations['NH4']*anions['PO4']/(2.5*10**-13), 10)
        except KeyError:
            pass

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
                    
                    if (ions[cation[0]]**cation_coef)*(ions[anion[0]]**anion_coef) > ksp:
                        # storage saturation index (SI) = log(K/Ksp), where K is the concentratilon products (not consider activity for now)
                        SI = log((ions[cation[0]]**cation_coef)*(ions[anion[0]]**anion_coef)/ksp, 10)
                        if f'({cation[0]}){cation_coef}({anion[0]}){anion_coef}' not in existed_precipitate:
                            precipitation_SI[sym.symbols(f'({cation[0]}){cation_coef}({anion[0]}){anion_coef}')] = [(cation[0], cation_coef), (anion[0], anion_coef), ksp, SI]

    for ion in ions.items():
        if ion[1] < 0:
            ions[ion[0]] = 0

    if len(precipitation_SI) == 0 and SI_struvite == 0:
        return ions, '', {}, 0
    else:

        # substract from the MB, but just for the ones with the largest SI
        SI_list = [precipitation[1][-1] for precipitation in precipitation_SI.items()]
        
        if SI_list == []:
            first_precipitation = [('Mg', 1),('NH4', 1),('PO4', 1), 2.5*10**-13, SI_struvite]
            precipitate_name = '(Mg)1(NH4)1(PO4)1'
            eqn_list.append(sym.Eq((sym.symbols(first_precipitation[0][0])**first_precipitation[0][1])*(sym.symbols(first_precipitation[1][0])**first_precipitation[1][1])*(sym.symbols(first_precipitation[2][0])**first_precipitation[2][1]), first_precipitation[-2]))
            eqn_dict[precipitate_name] = sym.symbols(precipitate_name)
            
        elif max(SI_list) < SI_struvite:
            first_precipitation = [('Mg', 1),('NH4', 1),('PO4', 1), 2.5*10**-13, SI_struvite]
            precipitate_name = '(Mg)1(NH4)1(PO4)1'
            eqn_list.append(sym.Eq((sym.symbols(first_precipitation[0][0])**first_precipitation[0][1])*(sym.symbols(first_precipitation[1][0])**first_precipitation[1][1])*(sym.symbols(first_precipitation[2][0])**first_precipitation[2][1]), first_precipitation[-2]))
            eqn_dict[precipitate_name] = sym.symbols(precipitate_name)

        else:
            first_precipitation_index = SI_list.index(max(SI_list))
            first_precipitation = list(precipitation_SI.values())[first_precipitation_index]
            precipitate_name = f'({first_precipitation[0][0]}){first_precipitation[0][1]}({first_precipitation[1][0]}){first_precipitation[1][1]}'

            eqn_list.append(sym.Eq((sym.symbols(first_precipitation[0][0])**first_precipitation[0][1])*(sym.symbols(first_precipitation[1][0])**first_precipitation[1][1]), first_precipitation[-2]))
            eqn_dict[precipitate_name] = sym.symbols(precipitate_name)
        
        for chemical in chemicals.items():
            for ion in first_precipitation[:-2]:
                if ion[0] in list(dict(chemical[0]).keys()):
                    new_tot = chemicals[chemical[0]][-1] - ion[1]*sym.symbols(precipitate_name)
                    tempo_list = list(chemicals[chemical[0]])
                    tempo_list[-1] = new_tot
                    tempo_dict = {chemical[0]: tuple(tempo_list)}
                    chemicals.update(tempo_dict)

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

        ans = sym.nsolve(eqn_list, tuple(eqn_dict.values()), [guess]*len(eqn_dict), dict=True, maxsteps=100, verify=False)      
        
        while ans[0][sym.symbols('H')] <= 0 or ans[0][sym.symbols('OH')] <= 0:
            guess *= 0.1
            ans = sym.nsolve(eqn_list, tuple(eqn_dict.values()), [guess]*len(eqn_dict), dict=True, maxsteps=100, verify=False)      

        ans = sym.nsolve(eqn_list, tuple(eqn_dict.values()), list(ans[0].values()), dict=True, maxsteps=100, verify=False)

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
        return ions, precipitate_name, chemicals, 1
        
def precision_round(number, digits=2):
    power = '{:e}'.format(number).split('e')[1]
    return round(number, -(int(power) - digits))