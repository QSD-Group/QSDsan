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

'''Select pipe size from volumetric flow rate and velocity.'''

import numpy as np
from math import pi


__all__ = (
    'calculate_excavation_volume',
    'calculate_pipe_material',
    'select_pipe',
    'cost_pump',
           )


# %%

# =============================================================================
# Excavation
# =============================================================================

def calculate_excavation_volume(L, W, D, excav_slope, constr_acess):
    '''
    Calculate the excavation volume needed as a frustum.

    Note that the units of the parameters do no matter as along as they are consistent.

    Parameters
    ----------
    L : float
        Length at the bottom.
    W : float
        Width at the bottom.
    D : float
        Depth.
    excav_slope : float
        Slope (horizontal/vertical).
    constr_acess : float
        Construction access on top of the length/width.
    '''
    L_bottom = L + 2*constr_acess
    W_bottom = W + 2*constr_acess
    diff = D * excav_slope
    return 0.5*(L_bottom*W_bottom+(L_bottom+diff)*(W_bottom+diff))


# %%

# =============================================================================
# Piping
# =============================================================================

def calculate_pipe_material(OD, t, ID, L, density=None):
    '''
    Calculate the total material needed of pipes.

    Parameters
    ----------
    OD : float
        Outer diameter of the pipe, [in].
    t : float
        Thickness of the pipe, [in].
    ID : float
        Inner diameter of the pipe, [in].
    L : float
        Length of the pipe, [ft].
    density : float
        Density of the material, [kg/ft3].

    Returns
    -------
    Quantity of the material in ft3 (if no density provided) or kg (if density provided).
    '''
    V = pi/4 * ((OD/12)**2-(ID/12)**2) * L
    quantity = V if not density else V*density
    return quantity


# Based on ANSI (American National Standards Institute) pipe chart
# the original code has a bug (no data for 22) and has been fixed here
boundaries = np.concatenate([
    np.arange(1/8, 0.5, 1/8),
    np.arange(0.5, 1.5, 1/4),
    np.arange(1.5, 5,   0.5),
    np.arange(5,   12,  1),
    np.arange(12,  36,  2),
    np.arange(36,  54,  6)
    ])

size = boundaries.shape[0]

pipe_dct = {
    1/8 : (0.405,  0.049), # OD (outer diameter), t (wall thickness)
    1/4 : (0.540,  0.065),
    3/8 : (0.675,  0.065),
    1/2 : (0.840,  0.083),
    3/4 : (1.050,  0.083),
    1   : (1.315,  0.109),
    1.25: (1.660,  0.109),
    1.5 : (1.900,  0.109),
    2   : (2.375,  0.109),
    2.5 : (2.875,  0.120),
    3   : (3.500,  0.120),
    3.5 : (4.000,  0.120),
    4   : (4.500,  0.120),
    4.5 : (5.000,  0.120),
    5   : (5.563,  0.134),
    6   : (6.625,  0.134),
    7   : (7.625,  0.134),
    8   : (8.625,  0.148),
    9   : (9.625,  0.148),
    10  : (10.750, 0.165),
    11  : (11.750, 0.165),
    12  : (12.750, 0.180),
    14  : (14.000, 0.188),
    16  : (16.000, 0.199),
    18  : (18.000, 0.188),
    20  : (20.000, 0.218),
    22  : (22.000, 0.250),
    24  : (24.000, 0.250),
    26  : (26.000, 0.250),
    28  : (28.000, 0.250),
    30  : (30.000, 0.312),
    32  : (32.000, 0.312),
    34  : (34.000, 0.312),
    36  : (36.000, 0.312),
    42  : (42.000, 0.312),
    48  : (48.000, 0.312)
    }


def select_pipe(Q, v):
    '''
    Select pipe based on Q (flow in ft3/s) and velocity (ft/s).

    Parameters
    ----------
    Q : float
        Flow rate of the fluid, [ft3/s] (cfs).
    v : float
        Velocity of the fluid, [ft/s].

    Returns
    -------
    Outer diameter, thickness, and inner diameter of the pipe (three floats), all in ft.
    '''
    A = Q / v # cross-section area
    d = (4*A/np.pi) ** 0.5 # minimum inner diameter, [ft]
    d *= 12 # minimum inner diameter, [in]
    d_index = np.searchsorted(boundaries, d, side='left') # a[i-1] < v <= a[i]
    d_index = d_index-1 if d_index==size else d_index # if beyond the largest size
    OD, t = pipe_dct[boundaries[d_index]]
    ID = OD - 2*t # inner diameter, [in]
    return OD, t, ID


# %%

# =============================================================================
# Pumping
# =============================================================================

#!!! Maybe move this to the Pump unit!
def cost_pump(unit=None, Q_mgd=None, recir_ratio=None,
              building_unit_cost=90):
    '''
    Calculate the cost of the pump and pump building for a unit
    based on its `Q_mgd` (hydraulic flow in million gallons per day),
    `recir_ratio` (recirculation ratio) attributes.

    Optionally, can left `unit` as None and directly provide
    `Q_mgd` and `recir_ratio` values.

    Parameters
    ----------
    unit : obj
        :class:`~.SanUnit` with the attribute `Q_mgd` and `recir_ratio`
        (if `Q_mgd` and/or `recir_ratio` are not provided).
    Q_mgd : float
        Hydraulic flow [mgd] (million gallons per day).
    recir_ratio : float
        Recirculation ratio.
    building_unit_cost : float
        Unit cost of the pump building, [USD/ft2].

    Returns
    -------
    Costs of the pumps and pump building (two floats) in USD.
    '''
    Q_mgd = Q_mgd or unit.Q_mgd
    recir_ratio = recir_ratio or unit.recir_ratio

    # Installed pump cost, this is a fitted curve
    pumps = 2.065e5 + 7.721*1e4*Q_mgd

    # Design capacity of intermediate pumps, gpm,
    # 2 is the excess capacity factor to handle peak flows
    GPMI = 2 * Q_mgd * 1e6 / 24 / 60

    # Design capacity of recirculation pumps, gpm
    GPMR = recir_ratio * Q_mgd * 1e6 / 24 / 60

    building = 0.
    for GPM in (GPMI, GPMR):
        if GPM == 0:
            N = 0
        else:
            N = 1 # number of buildings
            GPMi = GPM
            while GPMi > 80000:
                N += 1
                GPMi = GPM / N

        PBA = N * (0.0284*GPM+640) # pump building area, [ft]
        building += PBA * building_unit_cost

    return pumps, building