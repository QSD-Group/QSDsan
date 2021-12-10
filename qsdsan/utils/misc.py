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
    '''
    gases = cmps.gases
    solids = cmps.solids
    for i in cmps:
        ID = i.ID
        if ID not in dct:
            if ID in gases:
                dct[ID] = gas
            elif ID in solids:
                dct[ID] = solid
            else:
                dct[ID] = soluble
    return dct