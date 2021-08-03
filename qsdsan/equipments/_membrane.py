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

__all__ = ('Membrane',)


class Membrane(Equipment):
    '''
    Membranes to be used in an electrochemical cell.
    Refer to the example in :class:`ElectroChemCell` for how to use this class.

    Parameters
    ----------
    N : int
        Number of units of the given membrane.
    material: str
        Material of the membrane.
    unit_cost: float
        Unit cost of the membrane per m2, will use default cost (if available)
        if not provided.
    surface_area : float
        Surface area of the membrane in m2.

    See Also
    --------
    :class:`ElectroChemCell`

    '''
    __slots__ = ('_N', 'name', 'unit_cost', 'material', 'surface_area')

    def __init__(self, name=None, # when left as None, will be the same as the class name
                 design_units={},
                 F_BM=1., lifetime=10000, lifetime_unit='hr', N=0,
                 material='polypropylene', unit_cost=0.1, surface_area=1):
        Equipment.__init__(self=self, name=name, design_units=design_units,
                           F_BM=F_BM, lifetime=lifetime, lifetime_unit=lifetime_unit)
        self.name = name
        self.N = N
        self.unit_cost = unit_cost
        self.material = material
        self.surface_area = surface_area


    def _design(self):
        design = {
            f'Number of {self.name}': self.N,
            f'Material of {self.name}': self.material,
            f'Surface area of {self.name}': self.surface_area
            }
        self.design_units = {f'Surface area of {self.name}': 'm2'}
        return design


    # All subclasses of `Membrane` must have a `_cost` method, which returns the
    # purchase cost of this equipment
    def _cost(self):
        return self.unit_cost*self.N*self.surface_area


    @property
    def N(self):
        '''[int] Number of units of the electrode.'''
        return self._N
    @N.setter
    def N(self, i):
        self._N = int(i)