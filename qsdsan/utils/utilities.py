#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

from . import auom

__all__ = ('sum_system_utility',)


def sum_system_utility(system, operating_hours=None, exclude_units=(),
                       utility='power', result_unit=None,
                       calculate_net_utility=False):
    '''
    Sum up the select utility of a system
    (power in kWh/yr, heating/cooling duty in GJ/yr,
    note that the duty will be negative for cooling utilities).

    Parameters
    ----------
    system : obj
        The system whose utility usage will be calculated.
    operating_hours : float
        Operating hours for the utility, will use `system.operating_hours`
        if not provided here and `operating_hours` has been set for the system,
        will be set to 1 if neither information is available.
    exclude_units : Iterable(obj)
        Units within the system that will be excluded from the calculation.
    utility : str
        Utility of interest, can be either "power", "heating", or "cooling".
    result_unit : str
        If provided, the results will be converted to this unit.
    calculate_net_utility : bool
        Whether to calculate the net utility usage
        (e.g., subtract the amount produced by some units such as CHP).

    Examples
    --------
    >>> from qsdsan.utils import load_example_cmps, load_example_sys, sum_system_utility
    >>> sys = load_example_sys(load_example_cmps())
    >>> sys.simulate()
    >>> sum_system_utility(sys, utility='heating', result_unit='kJ/yr') # doctest: +NUMBER
    463359.8752661339
    >>> sum_system_utility(sys, utility='cooling', result_unit='GJ/yr') # doctest: +NUMBER
    0.0
    >>> # Exclude a certain unit
    >>> sum_system_utility(sys, utility='power') # doctest: +NUMBER
    8026.734351097134
    >>> sum_system_utility(sys, utility='power', exclude_units=(sys.units[0],)) # doctest: +NUMBER
    5795.725960018872
    '''

    hrs = operating_hours or getattr(system, 'operating_hours') or 1.
    utility = utility.lower()
    try: iter(exclude_units)
    except: exclude_units = (exclude_units,)
    units = [i for i in system.units if i not in exclude_units]
    get = getattr
    if utility in ('power', 'electricity'):
        unit = 'kWh/yr' if not result_unit else result_unit
        attr = 'consumption' if not calculate_net_utility else 'rate'
        tot = sum([get(i.power_utility, attr) for i in units if i.power_utility])*hrs
        return auom('kWh/hr').convert(tot, unit)
    elif utility == 'heating':
        unit = 'GJ/yr' if not result_unit else result_unit
        hu = sum([i.heat_utilities for i in units], ())
        if not calculate_net_utility:
            tot = sum([i.duty for i in hu if i.flow*i.duty>0])/1e6*hrs
        else:
            tot = sum([i.duty for i in hu if i.duty>0])/1e6*hrs
        return auom('GJ/yr').convert(tot, unit)
    elif utility == 'cooling':
        unit = 'GJ/yr' if not result_unit else result_unit
        cu = sum([i.heat_utilities for i in units], ())
        if not calculate_net_utility:
            tot = sum([i.duty for i in cu if i.flow*i.duty<0])/1e6*hrs
        else:
            tot = sum([i.duty for i in cu if i.duty<0])/1e6*hrs
        return auom('GJ/yr').convert(tot, unit)
    else:
        raise ValueError('`utility` can only be "power", "heating", or "cooling", '
                         f'not "{utility}".')