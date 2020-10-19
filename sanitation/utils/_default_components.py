#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 08:33:29 2020

@author: yalinli_cabbi
"""

import thermosteam as tmo
from .._component import Component

# A Chemicals object containing waste treatment-relevant characteristics,
# note that this component has not been compiled to allow for addition of
# new Chemical objects
__all__ = ('load_default_components',)



# %%

def load_default_components():
    
    # =============================================================================
    # Initialize the Chemicals object
    # =============================================================================
    components = tmo.Chemicals(())

    def component_database(ID, search_ID, **properties):
        component = Component(ID, search_ID=search_ID, **properties)
        components.append(component)
        return component
    
    def component_copied(ID, reference, **properties):
        component = reference.copy(ID)
        for i, j in properties.items():
            component.i = j
        components.append(component)
        return component
    
    def component_new(ID, **properties):
        component = Component(ID, **properties)
        components.append(component)
        return component
    
    # =============================================================================
    # Element-related
    # =============================================================================
    
    # Nitrogen-related
    N_tot = component_database('N_tot', search_ID='N', phase='l',
                               description='Total nitrogen')
    N_TKN = component_copied('N_TKN', N_tot, description='Total Kjeldahl nitrogen')
    N_NH3 = component_copied('N_NH3', N_tot, description='Total ammonia as nitrogen')
    
    # Phosphorus-related
    P_tot = component_database('P_tot', search_ID='P', phase='l',
                               description='Total phosphorus')
    
    # Potassium-related
    K_tot = component_database('K_tot', search_ID='K', phase='l',
                               description='Total potassium')
    
    # Magnesium-related
    Mg_tot = component_database('Mg_tot', search_ID='Mg', phase='l',
                                description='Total magnesium')
    
    # Calcium-related
    Ca_tot = component_database('Ca_tot', search_ID='Ca', phase='l',
                                description='Total calcium')
    
    # =============================================================================
    # Bulk Properties
    # =============================================================================
    
    COD = component_new('COD', phase='l', MW=32,
                        description='Chemical oxygen demand as O2',
                        degradability='Chemical', organic=True)
    BOD5 = component_copied('BOD5', COD,
                            description='Five-day biological oxygen demand as O2',
                        degradability='Biological', organic=True)
    BODu = component_copied('BODu', COD,
                            description='Ultimate biological oxygen demand as O2',
                            degradability='Biological', organic=True)
    
    H2O = tmo.Chemical('H2O')
    components.append(H2O)

    return components










