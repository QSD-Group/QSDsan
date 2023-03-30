#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

This is a pH solver, currently under development.

@author: Jianan Feng <jiananf2@illinois.edu>
"""

import sympy as sym, numpy as np
from math import log, sqrt

__all__ = ('pH_solver',)

def pH_solver(kw=10**-14,
              activity=True,
              check_precipitation = True,
              precise = 2,
              chemicals={}):
    '''
    The chemical input should be a dict, e.g.,
    chemicals = {
     (('H3PO4', 0), ('H2PO4', -1), ('HPO4', -2), ('PO4', -3)): (1.1*10**-2, 2.0*10**-7, 3.6*10**-13, 0.5),
     (('NH4', 1), ('NH3', 0)): (1.76*10**-5, -0.5),
     (('Mg', 2),): (0.5,)
      }
    '''
    
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

    pH = iterator(eqn_list, eqn_dict, chemicals,0,14, kw=kw)
    for i in range(precise+1):
        pH = iterator(eqn_list, eqn_dict, chemicals,*pH[:2], kw=kw)
    
    while abs((pH[2]-pH[3])/pH[2]) > 0.01 or abs((pH[2]-pH[3])/pH[3]) > 0.01:
        pH = iterator(eqn_list, eqn_dict, chemicals,*pH[:2], kw=kw)

    if activity == True:
        I = (pH[2]+pH[3])/2 
        delta_pH = 0.51*(sqrt(I)/(1+sqrt(I))-0.3*I) # Davies equation requires I<0.5
        pH = round((min(round(pH[0], precise), round(pH[1], precise)) + delta_pH), precise)
    else:
        pH = min(round(pH[0], precise), round(pH[1], precise))

    return pH
    
def iterator(eqn_list=None, eqn_dict=None, chemicals={}, minimum=0, maximum=14, kw=None):
    
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
        charge_dict[potential_pH] = (charge, I)
        
        if len(charge_dict) > 1:
            charge_values = list(charge_dict.values())
            if charge_values[-1][0]*charge_values[-2][0] <= 0:
                break
            
    return list(charge_dict.keys())[-2],list(charge_dict.keys())[-1], list(charge_dict.values())[-2][1], list(charge_dict.values())[-1][1]







