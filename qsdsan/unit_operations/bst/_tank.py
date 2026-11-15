#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:

    Yalin Li <mailto.yalin.li@gmail.com>

Part of this module is based on the biosteam package:
https://github.com/BioSTEAMDevelopmentGroup/biosteam

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

from warnings import warn
from math import pi
from inspect import signature, Parameter
from biosteam.units import Tank as BSTTank, MixTank as BSTMixTank, StorageTank as BSTStorageTank
from biosteam.exceptions import bounds_warning, DesignWarning
from biosteam.units.design_tools import flash_vessel_design
from biosteam.units.design_tools.specification_factors import material_densities_lb_per_ft3
from ... import SanUnit, Construction
from ...utils import auom


__all__ = ('Tank', 'MixTank', 'StorageTank', )

_lb_to_kg = auom('lb').conversion_factor('kg')
_m_to_ft = auom('m').conversion_factor('ft')
_Pa_to_psi = auom('Pa').conversion_factor('psi')

# Compatibility with BioSTEAM versions before/after the `design_input` feature:
# the latter drops the `kind` argument of `bounds_warning` and the `_default_*`
# class attributes of `Tank`, taking the defaults from design inputs instead.
_bounds_warning_has_kind = 'kind' in signature(bounds_warning).parameters

def _bounds_warning(source, key, value, units, bounds, kind='cost'):
    if _bounds_warning_has_kind:
        bounds_warning(source, key, value, units, bounds, kind)
    else:
        bounds_warning(source, key, value, units, bounds)


def _get_tank_default(cls, name):
    '''Default value of a tank design parameter (e.g., ``tau``) for ``cls``.'''
    legacy = getattr(cls, f'_default_{name}', Parameter.empty)
    if legacy is not Parameter.empty: return legacy
    for klass in cls.__mro__: # design inputs may be overridden by plain properties
        default = getattr(klass.__dict__.get(name), 'default', Parameter.empty)
        if default is not Parameter.empty: return default
    raise AttributeError(f"{cls.__name__} has no default for '{name}'")


class Tank(SanUnit, BSTTank, isabstract=True):
    '''
    Similar to the :class:`biosteam.units.Tank`,
    but can be initialized with :class:`qsdsan.SanStream` and :class:`qsdsan.WasteStream`.

    See Also
    --------
    `biosteam.units.Tank <https://biosteam.readthedocs.io/en/latest/API/units/tank.html>`_
    '''

    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                 vessel_type=None, tau=None, V_wf=None,
                 vessel_material=None, kW_per_m3=0.,
                 init_with='WasteStream', F_BM_default=None,
                 include_construction=True,):

        SanUnit.__init__(self, ID, ins, outs, thermo,
                         init_with=init_with, F_BM_default=F_BM_default,
                         include_construction=include_construction,)

        cls = type(self)
        self.vessel_type = vessel_type or _get_tank_default(cls, 'vessel_type')
        self.tau = tau or _get_tank_default(cls, 'tau')
        self.V_wf = V_wf or _get_tank_default(cls, 'V_wf')
        self.vessel_material = vessel_material or _get_tank_default(cls, 'vessel_material')
        self.kW_per_m3 = kW_per_m3 or _get_tank_default(cls, 'kW_per_m3')


class MixTank(Tank, BSTMixTank):
    '''
    Similar to the :class:`biosteam.units.MixTank`,
    but can be initialized with :class:`qsdsan.SanStream` and :class:`qsdsan.WasteStream`.

    .. note::

        For dynamic simulation, CSTR should be used.

    See Also
    --------
    `biosteam.units.MixTank <https://biosteam.readthedocs.io/en/latest/API/units/tank.html>`_
    '''


class StorageTank(Tank, BSTStorageTank):
    '''
    Similar to the :class:`biosteam.units.MixTank`, but can calculate material usage.

    See Also
    --------
    `biosteam.units.StorageTank <https://biosteam.readthedocs.io/en/latest/API/units/tank.html>`_
    '''

    _units = {'Diameter': 'ft',
              'Length': 'ft',
              'Wall thickness': 'in',
              'Weight': 'lb'}
    _bounds = {'Vertical vessel weight': (4200, 1e6),
               'Horizontal vessel weight': (1e3, 9.2e5),
               'Horizontal vessel diameter': (3, 21),
               'Vertical vessel length': (12, 40)}
    _vessel_material = 'Stainless steel'

    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 vessel_type=None, tau=None, V_wf=None,
                 vessel_material=None, kW_per_m3=0.,
                 init_with='WasteStream', F_BM_default=None,
                 include_construction=True, length_to_diameter=2):
        Tank.__init__(self, ID=ID, ins=ins, outs=outs, thermo=thermo,
                      init_with=init_with, F_BM_default=F_BM_default,
                      include_construction=include_construction,
                      vessel_type=vessel_type, tau=tau, V_wf=V_wf,
                      vessel_material=vessel_material, kW_per_m3=kW_per_m3,)
        self.length_to_diameter = length_to_diameter

    def _init_lca(self):
        item_name = self.vessel_material.replace(' ', '_')
        self.construction = [
            Construction(item_name.lower(), linked_unit=self, item=item_name, quantity_unit='kg'),
            ]


    def _design(self):
        BSTStorageTank._design(self)
        D = self.design_results

        Diameter = (4*D['Total volume']/pi/self.length_to_diameter)**(1/3)
        Diameter *= _m_to_ft # convert from m to ft
        L = Diameter * self.length_to_diameter # ft
        D.update(self._horizontal_vessel_design(self.ins[0].P*_Pa_to_psi, Diameter, L))
        D['Material'] = self.vessel_material
        if self.include_construction: self.construction[0].quantity = D['Weight']*_lb_to_kg


    def _horizontal_vessel_design(self, pressure, diameter, length) -> dict:
        pressure = pressure
        diameter = diameter
        length = length
        # Calculate vessel weight and wall thickness
        if self.vessel_material == 'Carbon steel':
            rho_M = material_densities_lb_per_ft3[self.vessel_material]
        else:
            rho_M = material_densities_lb_per_ft3['Stainless steel 304']
        if pressure < 14.68:
            warn('vacuum pressure vessel ASME codes not implemented yet; '
                 'wall thickness may be inaccurate and stiffening rings may be '
                 'required', category=DesignWarning)
        VW, VWT = flash_vessel_design.compute_vessel_weight_and_wall_thickness(
            pressure, diameter, length, rho_M)
        _bounds_warning(self, 'Horizontal vessel weight', VW, 'lb',
                        self._bounds['Horizontal vessel weight'], 'cost')
        _bounds_warning(self, 'Horizontal vessel diameter', diameter, 'ft',
                        self._bounds['Horizontal vessel diameter'], 'cost')
        Design = {}
        Design['Vessel type'] = 'Horizontal'
        Design['Length'] = length # ft
        Design['Diameter'] = diameter # ft
        Design['Weight'] = VW # lb
        Design['Wall thickness'] = VWT # in
        return Design

    @property
    def vessel_material(self):
        return self._vessel_material
    @vessel_material.setter
    def vessel_material(self, i):
        exist_material = getattr(self, '_vessel_material', None)
        BSTStorageTank.vessel_material.fset(self, i)
        if i and exist_material == i: return # type doesn't change, no need to reload construction items
        self._init_lca()