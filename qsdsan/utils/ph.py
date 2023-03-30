#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

This is a pH solver, currently under development.

@author: Jianan Feng <jiananf2@illinois.edu>
"""

import sympy as sym

def pH_solver(temperature=298.15,
              activity=True,
              check_precipitation = True,
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

    # check charge balance under different potential pH values
    for potential_pH in (1,2,3,4,5,6,7,8,9,10,11,12,13,14):
        eqn_list.append(sym.Eq(H,10**(-potential_pH)))
        
        eqn_dict['H'] = H
        
        ans = sym.solve(eqn_list,tuple(eqn_dict.values()))
        
        charge = 0
        for chemical in chemicals:
            for i in range(len(chemical)):
                for j in range (len(ans)):
                    if chemical[i][0] == str(list(ans.keys())[j]):
                        charge += chemical[i][1]*float(ans[list(ans.keys())[j]])
        
        print(f'When pH={potential_pH}, the net charge is {charge}')
        del eqn_list[-1]
    
    return pH