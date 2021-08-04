#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Smiti Mittal <smitimittal@gmail.com>
    Yalin Li <zoe.yalin.li@gmail.com>
    Anna Kogler

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''


from .. import Equipment

__all__ = ('Electrode',)


class Electrode(Equipment):
    '''
    Electrodes to be used in an electrochemical cell.
    Refer to the example in :class:`ElectroChemCell` for how to use this class.

    Parameters
    ----------
    N : int
        Number of units of the given electrode.
    electrode_type : str
        Type of the electrode, can only be "anode" or "cathode".
    material: str
        Material of the electrode.
    unit_cost: float
        Unit cost of the electrode, will use default cost (if available)
        if not provided.
    surface_area : float
        Surface area of the electrode in m2.

    See Also
    --------
    :class:`ElectroChemCell`

    '''

    # Include all attributes (no properties) in addition to the ones in the
    # parent `Equipment` class
    # Using __slots__ can improve computational efficiency when the class does not
    # have many attributes
    __slots__ = ('_type', '_material', '_N', '_unit_cost', 'surface_area')

    _default_unit_cost = {'graphite': 50}


    def __init__(self, name=None, # when left as None, will be the same as the class name
                 design_units={}, # if no value, then should be an empty dict
                 F_BM=1.,
                 lifetime=10000, lifetime_unit='hr', N=0,
                 electrode_type='anode',
                 material='graphite', unit_cost=0.1, surface_area=1):
        Equipment.__init__(self=self, name=name, design_units=design_units, F_BM=F_BM, lifetime=lifetime, lifetime_unit=lifetime_unit)
        self.N = N
        self.electrode_type = electrode_type
        self.unit_cost = unit_cost
        self.material = material
        self.surface_area = surface_area


    # All subclasses of `Equipment` must have a `_design` and a `_cost` method
    def _design(self):
        design = {
            'Type of electrode' : self.electrode_type,
            f'Number of {self.electrode_type}': self.N,
            f'Material of {self.electrode_type}': self.material,
            f'Surface area of {self.electrode_type}': self.surface_area
            }
        self.design_units = {f'Surface area of {self.electrode_type}': 'm2'}
        return design

    # All subclasses of `Equipment` must have a `_cost` method, which returns the
    # purchase cost of this equipment
    def _cost(self):
        return self.unit_cost*self.N


    # You can use property to add checks
    @property
    def N(self):
        '''[str] Number of units of the electrode.'''
        return self._N
    @N.setter
    def N(self, i):
        self._N = int(i)

    @property
    def electrode_type(self):
        '''[str] Type of the electrode, either "anode" or "cathode".'''
        return self._type
    @electrode_type.setter
    def electrode_type(self, i):
        if i.lower() in ('anode', 'cathode'):
            self._type = i
        else:
            raise ValueError(f'Electrode can only be "anode" or "cathode", not "{i}".')

    @property
    def unit_cost(self):
        '''[float] Cost of one electrode.'''
        if self._unit_cost:
            return self._unit_cost
        cost = self._default_unit_cost.get(self.material)
        return cost or 0.
    @unit_cost.setter
    def unit_cost(self, i):
        self._unit_cost = i

    @property
    def material(self):
        '''[str] Material of the electrode.'''
        return self._material
    @material.setter
    def material(self, i):
        self._material = i.lower()