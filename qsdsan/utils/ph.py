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
              precise = 2,
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
        raise RuntimeError('no chemical is given')
    
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

    eqn_list=[]
    eqn_dict={}

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

    pH = iterator(eqn_list, eqn_dict, chemicals, kw=kw, precipitation=check_precipitation)
    for i in range(precise+1):
        pH = iterator(eqn_list, eqn_dict, chemicals,*pH[:2], kw=kw, precipitation=check_precipitation)
    
    while abs((pH[2]-pH[3])/pH[2]) > 0.01 or abs((pH[2]-pH[3])/pH[3]) > 0.01:
        pH = iterator(eqn_list, eqn_dict, chemicals,*pH[:2], kw=kw, precipitation=check_precipitation)

    if activity == True:
        ions_sym = pH[4]
        I = (pH[2]+pH[3])/2
        if I>0.5:
            Warning('I is larger than 0.5, Davies equation may lead to wrong answer')
        delta_pH = 0.51*(sqrt(I)/(1+sqrt(I))-0.3*I) # Davies equation requires I<0.5
        pH = round((min(round(pH[0]+delta_pH, precise), round(pH[1]+delta_pH, precise))), precise)
    else:
        ions_sym = pH[4]
        pH = min(round(pH[0], precise), round(pH[1], precise))
    
    ions = {}
    for ion_sym in ions_sym.items():
        ions[str(ion_sym[0])] = ion_sym[1]
    ions['OH'] = kw/ions['H']
    
    
    
    
    return pH, ions
    
def iterator(eqn_list=None, eqn_dict=None, chemicals={}, minimum=-10, maximum=24, kw=None, precipitation=None):
    
    step = max(maximum-minimum+1, 11)
    # step = 11
    H = sym.symbols('H')
    charge_dict={}
    for potential_pH in np.linspace(minimum, maximum, step):
        eqn_list.append(sym.Eq(H,10**(-potential_pH)))
        
        eqn_dict['H'] = H

        ans = sym.solve(eqn_list,tuple(eqn_dict.values()), dict=True)
        
        ans=ans[0]
        
        charge = 10**(-potential_pH)-10**(potential_pH+log(kw, 10))
        I = 0.5*(10**(-potential_pH)+10**(potential_pH+log(kw, 10)))
        for chemical in chemicals:
            for i in range(len(chemical)):
                for j in range (len(ans)):
                    if chemical[i][0] == str(list(ans.keys())[j]):
                        charge += chemical[i][1]*float(ans[list(ans.keys())[j]])
                        I += 0.5*float(ans[list(ans.keys())[j]])*(chemical[i][1]**2)
        del eqn_list[-1]
        charge_dict[potential_pH] = (charge, I, ans)
        
        if len(charge_dict) > 1:
            charge_values = list(charge_dict.values())
            if charge_values[-1][0]*charge_values[-2][0] <= 0:
                break
            
    if precipitation:   
        # get all ions   
        ions = {}
        for item in ans.items():
            ions[str(item[0])] =item[1]
        ions['OH'] = kw/ions['H']
        
        # get ion:charge dict
        chemical_ion = {}
        for chemical in chemicals.keys():
            for i in chemical:
                chemical_ion[i[0]] = i[1]
        chemical_ion['H'] = 1
        chemical_ion['OH'] = -1
    
        # separate ions to cations and anions and also remove neutral ones
        cations = {}
        anions = {}
        for ion in ions.items():
            if chemical_ion[ion[0]] > 0:
                cations[ion[0]] = ion[1]
            else:
                anions[ion[0]] = ion[1]
    
    
    
        eqn_list=[]
        eqn_dict={}
        
        # identify possible precipitations:
        for cation in cations.items():
            if cation[0] in list(precipitation_inventory.keys()):
                for anion in anions.items():
                    if anion[0] in precipitation_inventory[cation[0]]:
                        index = precipitation_inventory[cation[0]].index(anion[0])
                        ksp = precipitation_inventory[cation[0]][index+1]
                        
                        least_common = lcm(abs(chemical_ion[cation[0]]), abs(chemical_ion[anion[0]]))
                        cation_cof = int(least_common/abs(chemical_ion[cation[0]]))
                        anion_cof = int(least_common/abs(chemical_ion[anion[0]]))
                            

                        
                        if (ions[cation[0]]**cation_cof)*(ions[anion[0]]**anion_cof) > ksp:
                            # add new equations here?
                            eqn_list.append(sym.Eq((sym.symbols(cation[0])**cation_cof)*(sym.symbols(anion[0])**anion_cof),ksp))
                            
                            # substract from the MB:
                            for chemical in chemicals.items():
                                if cation[0] in list(dict(chemical[0]).keys()):
                                    new_tot = chemicals[chemical[0]][-1] - cation_cof*sym.symbols(f'({cation[0]}){cation_cof}({anion[0]}){anion_cof}')
                                    tempo_list = list(chemicals[chemical[0]])
                                    tempo_list[-1] = new_tot
                                    tempo_dict = {chemical[0]: tuple(tempo_list)}
                                    chemicals.update(tempo_dict)
                            
                                if anion[0] in list(dict(chemical[0]).keys()):
                                    new_tot = chemicals[chemical[0]][-1] - anion_cof*sym.symbols(f'({cation[0]}){cation_cof}({anion[0]}){anion_cof}')
                                    tempo_list = list(chemicals[chemical[0]])
                                    tempo_list[-1] = new_tot
                                    tempo_dict = {chemical[0]: tuple(tempo_list)}
                                    chemicals.update(tempo_dict)
                            
                            
                            
# forget activitiy for precipitation for now

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

                            
                            
                            
                            
                        
                        print(cation[0]+anion[0]+str(precipitation_inventory[cation[0]][index+1]))
                    
            
            
            
            
            
            
            
            
            
            
            
    return list(charge_dict.keys())[-2], list(charge_dict.keys())[-1], list(charge_dict.values())[-2][1], list(charge_dict.values())[-1][1], list(charge_dict.values())[-2][2]



