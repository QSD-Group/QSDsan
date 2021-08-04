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

__all__ = ('Machine',)


class Machine(Equipment):
    '''
    Supplementary machines to be used in an electrochemical process.
    Refer to the example in :class:`ElectroChemCell` for how to use this class.

    Parameters
    ----------
    N : int
        Number of units of the given machine.
    unit_cost: float
        Unit cost of the machine

    See Also
    --------
    :class:`ElectroChemCell`

    '''
    __slots__ = ('_N', 'name', 'unit_cost')

    def __init__(self, name=None, # when left as None, will be the same as the class name
                 design_units={},
                 F_BM=1., lifetime=10000, lifetime_unit='hr', N=0,
                 unit_cost=0.1):
        Equipment.__init__(self=self, name=name, design_units=design_units,
                           F_BM=F_BM, lifetime=lifetime, lifetime_unit=lifetime_unit)
        self.name = name
        self.N = N
        self.unit_cost = unit_cost

    # All subclasses of `Machine` must have a `_design` and a `_cost` method
    def _design(self):
        design = {
            f'Number of {self.name}': self.N,
            }
        return design

    # All subclasses of `Membrane` must have a `_cost` method, which returns the
    # purchase cost of this equipment
    def _cost(self):
        return self.unit_cost*self.N


    @property
    def N(self):
        '''[int] Number of units of the electrode.'''
        return self._N
    @N.setter
    def N(self, i):
        self._N = int(i)