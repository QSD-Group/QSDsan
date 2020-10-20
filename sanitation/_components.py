#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 21:12:51 2020

@author: yalinli_cabbi
"""

import thermosteam as tmo

__all__ = ('Components', 'CompiledComponents')



# %%

# =============================================================================
# Define the Components class
# =============================================================================

# Component should not exist in the chemical database, otherwise should just use
# that chemical
class Components(tmo.Chemicals):
    '''
    A subclass of the Chemicals object in the thermosteam package, contains Component objects as attributes, is compatible with Chemical objects
    '''

    def __setattr__(self, ID, component):
        raise TypeError("can't set attribute; use <Components>.append instead")
    
    def __setitem__(self, ID, component):
        raise TypeError("can't set item; use <Components>.append instead")
        



# %%

# =============================================================================
# Define the CompiledComponents class
# =============================================================================
class CompiledComponents(tmo.CompiledChemicals):
    '''
    A subclass of the CompiledChemicals object in the thermosteam package, contains Component objects as attributes, is compatible with Chemical objects
    '''
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    