# -*- coding: utf-8 -*-
'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems
Copyright (C) 2020, Quantitative Sustainable Design Group

This module is developed by:
    Joy Cheung <joycheung1994@gmail.com>

This module is under the UIUC open-source license. Please refer to 
https://github.com/QSD-Group/QSDsan/blob/master/LICENSE.txt
for license details.
'''

import thermosteam as tmo
from ._parse import get_stoichiometric_coeff
from thermosteam.utils import chemicals_user
from sympy import symbols, Matrix
import numpy as np
    
@chemicals_user        
class Process():
    
    def __init__(self, reaction, ref_component, rate_equation=None, components=None, 
                 conserved_for=('COD', 'N', 'P', 'charge'), parameters=None):
        
        self._stoichiometry = []
        self.ref_component = ref_component
        self._components = self._load_chemicals(components)
        self._conserved_for = conserved_for
        self._parameters = {p: symbols(p) for p in parameters}
        
        self._stoichiometry = get_stoichiometric_coeff(reaction, self.ref_component, self._components, self._conserved_for, self._parameters)
        self._parse_rate_eq(rate_equation)
                
    def get_conversion_factors(self, as_matrix=False):        
        '''
        return conversion factors as a numpy ndarray to check conservation
        or return them as a sympy matrix to solve for unknown stoichiometric coefficients.
        '''
        if self._conservation_for:
            cmps = self._components
            arr = getattr(cmps, 'i_'+self._conservation_for[0])
            for c in self._conservation_for[1:]:
                arr = np.vstack((arr, getattr(cmps, 'i_'+c)))
            if as_matrix: return Matrix(arr.tolist())
            return arr
        else: return None
                    
    def check_conservation(self, tol=1e-8):
        '''check conservation for given tuple of materials subject to conservation. '''
        isa = isinstance
        if isa(self._stoichiometry, np.ndarray):
            ic = self.get_conversion_factors()
            v = self._stoichiometry
            ic_dot_v = ic @ v
            conserved_arr = np.isclose(ic_dot_v, np.zeros(ic_dot_v.shape), atol=tol)
            if not conserved_arr.all(): 
                materials = self._conserved_for
                unconserved = [(materials[i], ic_dot_v[i]) for i, conserved in enumerate(conserved_arr) if not conserved]
                raise RuntimeError("The following materials are unconserved by the "
                                   "stoichiometric coefficients. A positive value "
                                   "means the material is created, a negative value "
                                   "means the material is destroyed:\n "
                                   + "\n ".join([f"{material}: {value:.2f}" for material, value in unconserved]))
        else: 
            raise RuntimeError("Can only check conservations with numerical "
                               "stoichiometric coefficients.")
    
    #TODO: what happens to rate eq?
    def reverse(self):
        '''reverse the process as to flip the signs of all components.'''
        pass
        
    @property
    def ref_component(self):
        return getattr(self._components, self._ref_component)    
    @ref_component.setter
    def ref_component(self, ref_cmp):
        if ref_cmp: 
            self._ref_component = ref_cmp
            self._normalize_stoichiometry(ref_cmp)
            self._normalize_rate_eq(ref_cmp)

    @property
    def conserved_for(self):
        return self._conserved_for
    @conserved_for.setter
    def conserved_for(self, materials):
        self._conserved_for = materials
    
    @property
    def parameters(self):
        return sorted(self._parameters)
    
    #TODO: append new symbols to the dictionary
    def append_parameters(self):
        pass
    
    #TODO: set parameter values (and evaluate coefficients and rate??)
    def set_parameters(self):
        pass
    
    @property
    def stoichiometry(self):
        return dict(zip(self._components, self._stoichiometry))
        
    @property
    def rate_equation(self):
        return self._rate_equation
    
    #TODO: parse rate equations into symbols
    def _parse_rate_eq(self, eq):
        pass    
    
    def _normalize_stoichiometry(self, new_ref):
        isa = isinstance
        factor = abs(self._stoichiometry[self._components._index[new_ref]])
        if isa(self._stoichiometry, np.ndarray):
            self._stoichiometry /= factor
        elif isa(self._stoichiometry, list):
            self._stoichiometry = [v/factor for v in self._stoichiometry]
    
    def _normalize_rate_eq(self, new_ref):
        factor = self._stoichiometry[self._components._index[new_ref]]
        self._rate_equation *= factor

# #%%
# chemicals = tmo.Chemicals(['H2O', 'H2', 'O2'], cache=True)
# tmo.settings.set_thermo(chemicals)
# reaction = tmo.Reaction('H2O -> H2 + 0.5O2', reactant='H2O', X=0.7)

# # from sympy.solvers.solveset import linsolve
# from sympy import symbols, Matrix, solve
# # unknowns = symbols('v1:7')
# f, v1, v2, v3 = symbols('fsi v11 v12 v13')
# M = Matrix([[1,0,0,1,0,1],[.03, 1, 0, .01, 0, .04], [.01, 0, 1, 0, 0, .01], [0, .07, -.04839, 0, -1, 0]])
# v = Matrix([1-f, v1, v2, f, v3, -1])
# sol = solve(M * v, v1, v2, v3)

# # import numpy as np
# import os
# os.chdir("C:/Users/joy_c/Dropbox/PhD/Research/QSD/codes_developing/QSD-for-WaSH/sanitation")
# from qsdsan import Component, Components, WasteStream
# cmps = Components.load_default()
# tmo.settings.set_thermo(cmps)
# # react = cmps.subgroup(('SF', 'SNH4', 'SPO4', 'SU_E', 'SCO3', 'XB_Subst'))
# # M = Matrix(np.vstack((react.i_COD, react.i_N, react.i_P, react.i_charge)).tolist())
# p1 = Process('')
# p1.get_conversion_factors()
# p1._stoichiometry = np.ones(len(cmps.IDs))
# p1.check_conservation()


# dct = {'XS':-1, 
#        'SF':'1-fsi',
#        'SU_E':'fsi',
#        'SNH4':'?',
#        'SPO4':'?',
#        'SCO3':'?'}
# n_v = sum(v == '?' for v in dct.values())
# sp.var('v0:%s' % n_v)
