#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Smiti Mittal <smitimittal@gmail.com>
    Yalin Li <zoe.yalin.li@gmail.com>
    Anna Kogler
    
This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/master/LICENSE.txt
for license details.
'''

# %%

import math
from .. import Equipment, SanUnit, Component, WasteStream

__all__ = ('Electrode',)

#%%

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
    __slots__ = ('_type', '_material', '_N', 'unit_cost', 'surface_area')

    def __init__(self, name=None, # when left as None, will be the same as the class name
                 design_units={}, # if no value, then should be an empty dict
                 F_BM=1.,
                 lifetime=10000, lifetime_unit='hr', N=0,
                 electrode_type='anode', # note that I prefer not to use 'type' because it's a builtin function
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
            f'Type of electrode' : self.electrode_type,
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
        try:
            self._N = int(i)
        except:
            raise ValueError(f'N must be an integer')

    @property
    def electrode_type(self):
        '''[str] Type of the electrode, either "anode" or "cathode".'''
        return self._type
    @electrode_type.setter
    def electrode_type(self, i):
        if i.lower() in ('anode', 'cathode'):
            self._type = i
        else:
            raise ValueError(f'Electrode can only be "anode" or "cathode", not {i}.')

    @property
    def material(self):
        '''[str] Material of the electrode.'''
        return self._material
    @material.setter
    def material(self, i):
        material = i.lower()
        if material == 'graphite':
            # You can have some default unit cost based on the material,
            # I'm just making up numbers
            # But be careful that by doing this, you might overwriter users' input
            if not self.unit_cost:
                self.unit_cost = 50
        self._material = material
