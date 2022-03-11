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

import numpy as np
from math import pi

__all__ = (
    'calculate_concrete_volume',
    'concrete',
    'calculate_excavation_volume',
    'excavation',
    'calculate_pipe_material',
    'select_pipe',
    )


# %%

# =============================================================================
# Concrete
# =============================================================================

def calculate_concrete_volume(L, W, D, t_wall, t_slab, include_cover):
    '''
    Calculate needed volume of concrete.

    Note that the units of the parameters do no matter as along as they are consistent.

    Parameters
    ----------
    L : float
        Length of the concrete wall (without considering the wall thickness), [ft].
    W : float
        Width of the concrete wall (without considering the wall thickness), [ft].
    D : float
        Depth of the concrete wall, [ft].
    t_wall : float
        Thickness of the concrete wall, [ft].
    t_slab : float
        Thickness of the concrete slab, [ft].
    include_cover : bool
        Whether to include a cover on the top (t will be the same as the floor).
    '''
    # Wall concrete
    L_w_t = L + t_wall
    W_w_t = W + t_wall
    L_wall = 2 * (L_w_t+W_w_t)
    V_wall = L_wall * D * t_wall
    # Slab concrete
    A_slab = L_w_t * W_w_t if not include_cover else 2 * L_w_t * W_w_t
    V_slab = A_slab * (t_wall+t_slab)
    return V_wall, V_slab


def concrete(ID, L_concrete=0, W_concrete=0, D_concrete=0,
             t_wall=None, t_slab=None, include_cover=False,
             wall_concrete_unit_cost=24, # about $650/yd3, 650/27
             slab_concrete_unit_cost=13, # about $350/yd3, 350/27
             F_BM=1, lifetime=None):
    '''
    Handy decorator to add concrete usage-related properties and design/cost functions
    to a :class:`qs.SanUnit`.

    Parameters
    ----------
    ID : str
        ID for which concrete will be used.
    L_concrete : float
        Length of the concrete wall (without considering the wall thickness), [ft].
    W_concrete : float
        Width of the concrete wall (without considering the wall thickness), [ft].
    D_concrete : float
        Depth of the concrete wall, [ft].
    t_wall : float
        Thickness of the concrete wall, [ft],
        if not provided, will be calculated based on `D_concrete` (1 ft per 12 ft depth, minimum 1 ft).
    t_slab : float
        Thickness of the concrete slab, [ft],
        if not provided, will be calculated based on `t_wall` (t_wall+2/12 ft).
    wall_concrete_unit_cost : float
        Unit cost of the wall concrete, [$/ft3].
    slab_concrete_unit_cost : float
        Unit cost of the slab concrete, [$/ft3].
    include_cover : bool
        Whether to include a cover on the top (t will be the same as the floor).
    F_BM : float
        Bare module factor of the concrete.
    lifetime : int
        Lifetime of the concrete.

    Examples
    --------
    >>> from qsdsan import SanUnit, set_thermo
    >>> from qsdsan.utils import load_example_cmps, concrete
    >>> # Create the class
    >>> @concrete('Pump building', L_concrete=20, W_concrete=10, D_concrete=5)
    ... @concrete('Blower building', L_concrete=20, W_concrete=10, D_concrete=5,
    ...           t_wall=2, t_slab=3, include_cover=True)
    ... class FakeUnit(SanUnit):
    ...     def _design(self):
    ...         self._add_concrete_design()
    ...     def _cost(self):
    ...         self._add_concrete_cost()
    >>> # Use the class
    >>> cmps = load_example_cmps()
    >>> set_thermo(cmps)
    >>> U1 = FakeUnit('U1', ins='ws1', outs='ws2')
    >>> U1.simulate()
    >>> U1.results() # doctest: +SKIP

    References
    ----------
    [1] Shoener et al., Design of Anaerobic Membrane Bioreactors for the
    Valorization of Dilute Organic Carbon Waste Streams.
    Energy Environ. Sci. 2016, 9 (3), 1102–1112.
    https://doi.org/10.1039/C5EE03715H.
    '''
    return lambda cls: add_concrete(cls, ID, L_concrete, W_concrete, D_concrete,
                                    t_wall, t_slab, include_cover,
                                    wall_concrete_unit_cost, slab_concrete_unit_cost,
                                    F_BM, lifetime)


def add_concrete(cls, ID, L_concrete, W_concrete, D_concrete,
                 t_wall, t_slab, include_cover,
                 wall_concrete_unit_cost, slab_concrete_unit_cost,
                 F_BM, lifetime):
    def set_L_concrete(cls, i): cls._L_concrete = i
    def set_W_concrete(cls, i): cls._W_concrete = i
    def set_D_concrete(cls, i): cls._D_concrete = i
    def set_t_wall(cls, i): cls._t_wall = i
    def set_t_slab(cls, i): cls._t_slab = i
    def set_include_cover(cls, i): cls._include_cover = i
    def set_wall_concrete_unit_cost(cls, i): cls._wall_concrete_unit_cost = i
    def set_slab_concrete_unit_cost(cls, i): cls._slab_concrete_unit_cost = i
    prop_dct = {
        'L_concrete':
            (L_concrete, # value of the property
             lambda cls: cls._L_concrete, # getter
             set_L_concrete, # setter
             '[dict(float)] Length of the concrete wall (without considering the wall thickness), [ft].'), # documentation
        'W_concrete':
            (W_concrete,
             lambda cls: cls._W_concrete,
             set_W_concrete,
             '[dict(float)] Length of the concrete wall (without considering the wall thickness), [ft].'),
        'D_concrete':
            (D_concrete,
             lambda cls: cls._D_concrete,
             set_D_concrete,
             '[dict(float)] Depth of the concrete wall, [ft].'),
        't_wall':
            (t_wall,
             lambda cls: cls._t_wall,
             set_t_wall,
             '[dict(float)] Thickness of the concrete wall, [ft], if not provided, will be calculated based on `D_concrete` (1 ft per 12 ft depth, minimum 1 ft).'),
        't_slab':
            (t_slab,
             lambda cls: cls._t_slab,
             set_t_wall,
             '[dict(float)] Thickness of the concrete slab, [ft], if not provided, will be calculated based on `t_wall` (t_wall+2/12 ft).'),
        'include_cover':
            (include_cover,
             lambda cls: cls._include_cover,
             set_include_cover,
             '[dict(bool)] Whether to include a cover on the top (t will be the same as the floor).'),
        'wall_concrete_unit_cost':
            (wall_concrete_unit_cost,
             lambda cls: cls._wall_concrete_unit_cost,
             set_wall_concrete_unit_cost,
             '[dict(float)] Unit cost of the wall concrete, [$/ft3].'),
        'slab_concrete_unit_cost':
            (slab_concrete_unit_cost,
             lambda cls: cls._slab_concrete_unit_cost,
             set_slab_concrete_unit_cost,
             '[dict(float)] Unit cost of the slab concrete, [$/ft3].'),
        }
    # Add attributes and properties
    get_cls = getattr
    set_cls = setattr
    try:
        items = get_cls(cls, '_concrete')
        items.append(ID)
    except AttributeError:
        items = cls._concrete = [ID]
    for prop_name, val in prop_dct.items():
        # Firstly set the attributes with a `_` prefix.
        attr_name = f'_{prop_name}'
        try:
            existing = get_cls(cls, attr_name)
            existing[ID] = val[0]
        except AttributeError:
            set_cls(cls, attr_name, {ID: val[0]})
            # Then add the corresponding properties with documentation.
            prop = property(fget=val[1], fset=val[2], doc=val[3])
            set_cls(cls, prop_name, prop)

    wall_n_slab = ('wall', 'slab')
    for wall_or_slab in wall_n_slab:
        key = f'{ID} {wall_or_slab} concrete'
        cls._F_BM_default[key] = F_BM
        cls._default_equipment_lifetime[key] = lifetime

    # Add convenient functions for instances of the class
    try: get_cls(cls, '_add_concrete_design')
    except AttributeError:
        def _add_concrete_design(self):
            D = self.design_results
            units = self._units
            L_dct, W_dct, D_dct = get_cls(self, 'L_concrete'), \
                get_cls(self, 'W_concrete'), get_cls(self, 'D_concrete')
            t_wall_dct, t_slab_dct,cover_dct = get_cls(self, 't_wall'), \
                get_cls(self, 't_slab'), get_cls(self, 'include_cover')
            for item in items:
                D_concrete = D_dct[item]
                t_wall = t_wall_dct[item] = t_wall_dct[item] or (1+max(D_concrete-12, 0)/12)
                t_slab = t_slab_dct[item] = t_slab_dct[item] or t_wall+2/12
                Vs = calculate_concrete_volume(
                    L_dct[item], W_dct[item], D_concrete, t_wall, t_slab, cover_dct[item])
                for wall_or_slab, V in zip(wall_n_slab, Vs):
                    key = f'{item} {wall_or_slab} concrete'
                    units[key] = 'ft3'
                    D[key] = V
        cls._add_concrete_design = _add_concrete_design

        def _add_concrete_cost(self):
            D = self.design_results
            C = self.baseline_purchase_costs
            cost_dcts = get_cls(self, 'wall_concrete_unit_cost'), \
                get_cls(self, 'slab_concrete_unit_cost')
            for item in items:
                for wall_or_slab, cost_dct in zip(wall_n_slab, cost_dcts):
                    key = f'{item} {wall_or_slab} concrete'
                    C[key] = D[key] * cost_dct[item]
        cls._add_concrete_cost = _add_concrete_cost

    return cls


# %%

# =============================================================================
# Excavation
# =============================================================================

def calculate_excavation_volume(L, W, D, excav_slope, constr_access):
    '''
    Calculate the excavation volume needed as a frustum.

    Note that the units of the parameters do no matter as along as they are consistent.

    Parameters
    ----------
    L : float
        Excavation length at the bottom.
    W : float
        Excavation width at the bottom.
    D : float
        Excavation depth.
    excav_slope : float
        Excavation slope (horizontal/vertical).
    constr_acess : float
        Extra (i.e., on top of the length/width) room for construction access.
    '''
    L_bottom = L + 2*constr_access
    W_bottom = W + 2*constr_access
    diff = D * excav_slope
    return 0.5*(L_bottom*W_bottom+(L_bottom+diff)*(W_bottom+diff))*D


def excavation(ID, L_excav=0, W_excav=0, D_excav=0,
               excav_slope=1.5, constr_access=3, excav_unit_cost=0.011, # about $0.3/yd3, 0.3/27
               F_BM=1, lifetime=None):
    '''
    Handy decorator to add excavation-related properties and design/cost functions
    to a :class:`qs.SanUnit`.

    The excavation volume is calculated as a frustum.

    Parameters
    ----------
    ID : str
        ID of this excavation activity.
    L_excav : float
        Excavation length at the bottom, [ft].
    W_excav : float
        Excavation width at the bottom, [ft].
    D_excav : float
        Excavation depth. [ft].
    excav_slope : float
        Excavation slope (horizontal/vertical).
    constr_access : float
        Extra (i.e., on top of the length/width) room for construction access, [ft].
    excav_unit_cost : float
        Unit cost of the excavation activity, [$/ft3].
    F_BM : float
        Bare module factor of this excavation activity.
    lifetime : int
        Lifetime of this excavation activity.

    Examples
    --------
    >>> from qsdsan import SanUnit, set_thermo
    >>> from qsdsan.utils import load_example_cmps, excavation
    >>> # Create the class
    >>> @excavation('Reactor building', L_excav=100, W_excav=150, D_excav=20,
    ...             excav_unit_cost=8/27)
    ... @excavation('Pump and blower building', L_excav=30, W_excav=15, D_excav=10)
    ... class FakeUnit(SanUnit):
    ...     def _design(self):
    ...         self._add_excavation_design()
    ...     def _cost(self):
    ...         self._add_excavation_cost()
    >>> # Use the class
    >>> cmps = load_example_cmps()
    >>> set_thermo(cmps)
    >>> U1 = FakeUnit('U1', ins='ws1', outs='ws2')
    >>> U1.simulate()
    >>> U1.results() # doctest: +SKIP

    References
    ----------
    [1] Shoener et al., Design of Anaerobic Membrane Bioreactors for the
    Valorization of Dilute Organic Carbon Waste Streams.
    Energy Environ. Sci. 2016, 9 (3), 1102–1112.
    https://doi.org/10.1039/C5EE03715H.
    '''
    return lambda cls: add_excavation(cls, ID, L_excav, W_excav, D_excav,
                                      excav_slope, constr_access, excav_unit_cost,
                                      F_BM, lifetime)

def add_excavation(cls, ID, L_excav, W_excav, D_excav,
                   excav_slope, constr_access, excav_unit_cost,
                   F_BM, lifetime):
    def set_L_excav(cls, i): cls._L_excav = i
    def set_W_excav(cls, i): cls._W_excav = i
    def set_D_excav(cls, i): cls._D_excav = i
    def set_excav_slope(cls, i): cls._excav_slope = i
    def set_constr_access(cls, i): cls._constr_access = i
    def set_excav_unit_cost(cls, i): cls._excav_unit_cost = i
    prop_dct = {
        'L_excav':
            (L_excav, # value of the property
             lambda cls: cls._L_excav, # getter
             set_L_excav, # setter
             '[dict(float)] Excavation length at the bottom, [ft].'), # documentation
        'W_excav':
            (W_excav,
             lambda cls: cls._W_excav,
             set_W_excav,
             '[dict(float)] Excavation width at the bottom, [ft].'),
        'D_excav':
            (D_excav,
             lambda cls: cls._D_excav,
             set_D_excav,
             '[dict(float)] Excavation depth at the bottom, [ft].'),
        'excav_slope':
            (excav_slope,
             lambda cls: cls._excav_slope,
             set_excav_slope,
             '[dict(float)] Slope for excavation (horizontal/vertical).'),
        'constr_access':
            (constr_access,
             lambda cls: cls._constr_access,
             set_constr_access,
             '[dict(float)] Extra (i.e., on top of the length/width) room for construction access, [ft].'),
        'excav_unit_cost':
            (excav_unit_cost,
             lambda cls: cls._excav_unit_cost,
             set_excav_unit_cost,
             '[dict(float)] Unit cost of the excavation activity, [$/ft3].'),
        }
    # Add attributes and properties
    get_cls = getattr
    set_cls = setattr
    try:
        items = get_cls(cls, '_excavation')
        items.append(ID)
    except AttributeError:
        items = cls._excavation = [ID]
    for prop_name, val in prop_dct.items():
        # Firstly set the attributes with a `_` prefix.
        attr_name = f'_{prop_name}'
        try:
            existing = get_cls(cls, attr_name)
            existing[ID] = val[0]
        except AttributeError:
            set_cls(cls, attr_name, {ID: val[0]})
            # Then add the corresponding properties with documentation.
            prop = property(fget=val[1], fset=val[2], doc=val[3])
            set_cls(cls, prop_name, prop)

    key = f'{ID} excavation'
    cls._F_BM_default[key] = F_BM
    cls._default_equipment_lifetime[key] = lifetime

    # Add convenient functions for instances of the class
    try: get_cls(cls, '_add_excavation_design')
    except AttributeError:
        def _add_excavation_design(self):
            D = self.design_results
            units = self._units
            L_dct, W_dct, D_dct = \
                get_cls(self, 'L_excav'), get_cls(self, 'W_excav'), get_cls(self, 'D_excav')
            excav_slope_dct, constr_access_dct = \
                get_cls(self, 'excav_slope'), get_cls(self, 'constr_access')
            for item in items:
                key = f'{item} excavation'
                units[key] = 'ft3'
                D[key] = calculate_excavation_volume(
                    L_dct[item], W_dct[item], D_dct[item],
                    excav_slope_dct[item], constr_access_dct[item])
        cls._add_excavation_design = _add_excavation_design

        def _add_excavation_cost(self):
            D = self.design_results
            C = self.baseline_purchase_costs
            excav_unit_cost_dct = get_cls(self, 'excav_unit_cost')
            for item in items:
                key = f'{item} excavation'
                C[key] = D[key] * excav_unit_cost_dct[item]
        cls._add_excavation_cost = _add_excavation_cost
    return cls


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