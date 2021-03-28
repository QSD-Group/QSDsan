#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems
Copyright (C) 2020, Quantitative Sustainable Design Group

This module is developed by:
    Joy Cheung <joycheung1994@gmail.com>
    Yalin Li

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''


from sympy import symbols, sympify, simplify, Matrix, solve
from sympy.parsing.sympy_parser import parse_expr
import numpy as np

__all__ = ('get_stoichiometric_coeff', )

#%%

def split_coefficient(nID, sign):
    if len(nID.split(']')) > 1:
        coeff, ID = nID.split(']')[0].split('[')[1], nID.split(']')[1]
        try: return sign * float(coeff), ID
        except ValueError: 
            if sign == 1: return coeff, ID
            elif sign == -1: return '-('+coeff+')', ID
    else:
        for i, letter in enumerate(nID):
            if letter != 'e' and letter.isalpha(): break
        if i: 
            ID = nID[i:]
            n = sign * float(nID[:i])
        else: 
            ID = nID
            n = sign
    return n, ID

def extract_coefficients(nIDs, dct, sign):
    for nID in nIDs:
        n, ID = split_coefficient(nID, sign)
        if ID in dct:
            raise ValueError('Components can only appear once in a reaction; '
                            f'multiple instances of {repr(ID)} found')
        dct[ID] = n        

def str2dct(reaction) -> dict:
    reaction = reaction.replace(' ', '')
    left, right = reaction.split('->')
    reactants = left.split('+')
    products = right.split('+')
    dct = {}
    extract_coefficients(reactants, dct, -1)
    extract_coefficients(products, dct, 1)
    return dct

def dct2arr(dct, components):
    arr = np.zeros(components.size)
    cmp_index = components._index
    for ID, coefficient in dct.items():
        arr[cmp_index[ID]] = coefficient
    return arr

def dct2list(dct, components):
    lst = [0] * components.size
    cmp_index = components._index
    for ID, coefficient in dct.items():
        lst[cmp_index[ID]] = coefficient
    return lst    

def get_ic(cmps, conservation_for):        
    '''
    return conversion factors as a sympy matrix to solve for unknown stoichiometric coefficients.
    '''
    if conservation_for:
        arr = getattr(cmps, 'i_'+conservation_for[0])
        for c in conservation_for[1:]:
            arr = np.vstack((arr, getattr(cmps, 'i_'+c)))
        return Matrix(arr.tolist())
    else: return None

def symbolize(coeff_dct, components, conserved_for, parameters):
    n = sum([v in ('?', '-(?)') for v in coeff_dct.values()])
    if n > 0:
        unknowns = symbols('unknown0:%s' % n)
        i = 0
        v_arr = []
        IDs = sorted(coeff_dct)
        for cmp in IDs:
            if coeff_dct[cmp] in ('?', '-(?)'):
                v_arr.append('unknown%s' % i)
                i += 1
            else: v_arr.append(coeff_dct[cmp])
        v = Matrix(sympify(v_arr, parameters))
        ic = get_ic(components.subgroup(IDs), conserved_for)
        sol = solve(ic * v, unknowns)
        coeff_dct = dict(zip(IDs, simplify(v.subs(sol))))
        del unknowns                    
    else:
        coeff_dct = {k: simplify(parse_expr(v, local_dict=parameters)) for k, v in coeff_dct.items()}
    return coeff_dct
    
def get_stoichiometric_coeff(reaction, ref_component, components, conserved_for, parameters):
    '''
    parse input reaction to get array of symbolic expressions (or a function to return numpy array with parameters as kwargs) 
    or numpy array of input values for stoichiometric coefficients.
    '''
    isa = isinstance
    coeff_dct = None    

    if isa(reaction, np.ndarray):
        if len(reaction) == components.size: coeff = reaction
        else: 
            raise ValueError("The array of stoichiometric coefficients must have "
                             "the same length as Components {repr(components))}.")
    elif isa(reaction, dict):
        if all([isa(v, (float, int)) for v in reaction.values()]):
            coeff = dct2arr(reaction, components)
        else: 
            coeff_dct = reaction
    elif isa(reaction, str):
        coeff_dct = str2dct(reaction)
    else:
        raise ValueError("reaction must be either a str or a dict or an array; "
                         f"not a '{type(reaction).__name__}' object")        
    if coeff_dct:
        coeff_dct = symbolize(coeff_dct, components, conserved_for, parameters)
        coeff = dct2list(coeff_dct, components)
    if ref_component:
        normalize_factor = abs(coeff[components._index[ref_component]])
        if isa(coeff, np.ndarray): coeff /= normalize_factor
        else: coeff[:] = [v/normalize_factor for v in coeff]
    return coeff

#%%
# from qsdsan import Components

# # hydrolysis
# rxn = 'XB_Subst -> [1-fsi]SF + [fsi]SU_E + [?]SNH4 + [?]SPO4 + [?]SCO3'
# ref = 'XB_Subst'
# cmps = Components.load_default()
# conserved_for = ('COD', 'N', 'P', 'charge')
# parameters = {'fsi': symbols('fsi')}

# # anoxic growth on SF
# rxn = '[1/yh]SF + [(1-yh)/yh/2.86]SNO3 + [?]SPO4 + [?]SNH4 + [?]SCO3 -> [(1-yh)/yh/2.86]SN2 + XOHO'
# ref = 'XOHO'
# # cmps = Components.load_default()
# # conserved_for = ('COD', 'N', 'P', 'charge')
# parameters = {'yh': symbols('yh'),}

# stoichiometry = get_stoichiometric_coeff(rxn, ref, cmps, conserved_for, parameters)
# dict(zip(cmps.IDs, stoichiometry))

#%%
# nID = '[1-fsi]SF'
# coeffient = '1-fsi'

# split_coefficient(nID, 1)
