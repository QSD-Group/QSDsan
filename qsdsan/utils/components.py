#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>

Part of the code is based on modules in BioSTEAM,
https://github.com/BioSTEAMDevelopmentGroup/biosteam


This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

from thermosteam.functional import rho_to_V
from . import auom


__all__ = ('add_V_from_rho', 'correct_param_with_T', 'default_component_dict', )


def add_V_from_rho(component, rho, rho_unit='kg/m3'):
    '''
    Add a constant molar volume model to the component with the given rho.

    Parameters
    ----------
    component : obj
        The component for which the molar volume model will be added.
    rho : float
        The density of the component.
    rho_unit : str
        Unit of the density rho.

    Examples
    --------
    >>> from qsdsan import Component
    >>> from qsdsan.utils import add_V_from_rho
    >>> P4O10 = Component('P4O10', phase='s', organic=False,
    ...                   particle_size='Particulate', degradability='Undegradable')
    >>> add_V_from_rho(P4O10, 2.39, rho_unit='g/mL') # http://www.chemspider.com/Chemical-Structure.14128.html
    >>> P4O10.V
    VolumeSolid(CASRN="16752-60-6 (P4O10)", MW=283.889048, extrapolation="linear", method="USER_METHOD")
    >>> P4O10.V(330) # doctest: +ELLIPSIS
    0.0001187...
    '''
    rho = auom(rho_unit).convert(rho, 'kg/m3')
    V_model = rho_to_V(rho, component.MW)
    try: component.V.add_model(V_model)
    except AttributeError:
        handle = getattr(component.V, component.locked_state)
        handle.add_model(V_model)


def correct_param_with_T(k0, T0, T, theta):
    r'''
    Correct the coefficient value considering the temperature effect by:

    .. math::

        k = k_0 * \theta^{T-T_0}

    Parameters
    ----------
    k0 : float
        Coefficient value at a certain temperature.
    T0 : float
        Temperature for the known coefficient value.
    T : float
        New temperature will the coefficient value will be calculated for.
    theta :
        Temperature coefficient.

    Examples
    --------
    >>> # Correct the maximum specific growth rate for heterotrophs from 20 to 30°C
    >>> from qsdsan.utils import correct_param_with_T
    >>> mu_hat_30 = correct_param_with_T(k0=6, T0=20+273.15, T=30+273.15, theta=1.08)
    >>> mu_hat_30 # doctest: +ELLIPSIS
    12.95...
    '''
    return k0*theta**(T-T0)



def default_component_dict(dct={}, *, cmps, gas=None, soluble=None, solid=None):
    '''
    Fill out a dict with default values of components based on its property
    (values already in the dict will NOT be updated).

    Parameters
    ----------
    dct : dict
        Dict to be updated with the default values.
    g : value
        Default value for gas components (components in `cmps.gases`).
    s : value
        Default value for soluble components (components NOT in `cmps.gases` nor `cmps.solids`).
    x : value
        Default value for solid components (components in `cmps.solids`).

    Examples
    --------
    >>> from qsdsan import Components
    >>> from qsdsan.utils import default_component_dict
    >>> cmps = Components.load_default()
    >>> split_dct = default_component_dict(dct={}, cmps=cmps, gas=1, soluble=0.95, solid=0)
    >>> split_dct
    {'S_H2': 1,
     'S_CH4': 1,
     'S_CH3OH': 0.95,
     'S_Ac': 0.95,
     'S_Prop': 0.95,
     'S_F': 0.95,
     'S_U_Inf': 0.95,
     'S_U_E': 0.95,
     'C_B_Subst': 0.95,
     'C_B_BAP': 0.95,
     'C_B_UAP': 0.95,
     'C_U_Inf': 0.95,
     'X_B_Subst': 0,
     'X_OHO_PHA': 0,
     'X_GAO_PHA': 0,
     'X_PAO_PHA': 0,
     'X_GAO_Gly': 0,
     'X_PAO_Gly': 0,
     'X_OHO': 0,
     'X_AOO': 0,
     'X_NOO': 0,
     'X_AMO': 0,
     'X_PAO': 0,
     'X_MEOLO': 0,
     'X_FO': 0,
     'X_ACO': 0,
     'X_HMO': 0,
     'X_PRO': 0,
     'X_U_Inf': 0,
     'X_U_OHO_E': 0,
     'X_U_PAO_E': 0,
     'X_Ig_ISS': 0,
     'X_MgCO3': 0,
     'X_CaCO3': 0,
     'X_MAP': 0,
     'X_HAP': 0,
     'X_HDP': 0,
     'X_FePO4': 0,
     'X_AlPO4': 0,
     'X_AlOH': 0,
     'X_FeOH': 0,
     'X_PAO_PP_Lo': 0,
     'X_PAO_PP_Hi': 0,
     'S_NH4': 0.95,
     'S_NO2': 0.95,
     'S_NO3': 0.95,
     'S_PO4': 0.95,
     'S_K': 0.95,
     'S_Ca': 0.95,
     'S_Mg': 0.95,
     'S_CO3': 0.95,
     'S_N2': 1,
     'S_O2': 1,
     'S_CAT': 0.95,
     'S_AN': 0.95,
     'H2O': 0.95}
    '''
    gases = cmps.gases
    solids = cmps.solids
    for cmp in cmps:
        ID = cmp.ID
        if ID not in dct:
            if cmp in gases:
                dct[ID] = gas
            elif cmp in solids:
                dct[ID] = solid
            else:
                dct[ID] = soluble
    return dct