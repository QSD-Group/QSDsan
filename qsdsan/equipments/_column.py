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

__all__ = ('Column',)


class Column(Equipment):
    '''
    Columns to be used in an electrochemical cell.
    Refer to the example in :class:`~.sanunits.ElectroChemCell` for how to use this class.

    Parameters
    ----------
    N : int
        Number of units of the given column.
    material: str
        Material of the column.
    unit_cost: float
        Unit cost of the column per m2, will use default cost (if available)
        if not provided.
    surface_area : float
        Surface area of the column in m2.

    See Also
    --------
    :class:`~.sanunits.ElectroChemCell`

    '''

    def __init__(self, ID='', linked_unit=None,
                 units={
                     'Number of columns': '',
                     'Material of the column': '',
                     'Surface area of columns': 'm2',
                     },
                 F_BM=1., lifetime=10000, lifetime_unit='hr', N=0,
                 material='resin', unit_cost=0.1, surface_area=1):
        Equipment.__init__(self=self, ID=ID, linked_unit=linked_unit, units=units,
                           F_BM=F_BM, lifetime=lifetime, lifetime_unit=lifetime_unit)
        self.N = N
        self.unit_cost = unit_cost
        self.material = material
        self.surface_area = surface_area


    def _design(self):
        design = {
            'Number of columns': self.N,
            'Material of the column': self.material,
            'Surface area of columns': self.surface_area
            }
        return design


    def _cost(self):
        return self.unit_cost*self.N*self.surface_area


    @property
    def N(self):
        '''[str] Number of units of the columns.'''
        return self._N
    @N.setter
    def N(self, i):
        self._N = int(i)