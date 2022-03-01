#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Smiti Mittal <smitimittal@gmail.com>
    Yalin Li <zoe.yalin.li@gmail.com>
    Anna Kogler
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''


from .. import Equipment
# from .utils import auom

__all__ = ('Membrane', 'HollowFiberMembrane')


class Membrane(Equipment):
    '''
    Membranes to be used in an electrochemical cell.
    Refer to the example in :class:`~.sanunits.ElectroChemCell` for how to use this class.

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
    :class:`~.sanunits.ElectroChemCell`

    '''

    def __init__(self, ID='', linked_unit=None,
                 units={
                     'Number of membrane modules': '',
                     'Material of membrane': '',
                     'Surface area of membrane': 'm2'
                     },
                 F_BM=1., lifetime=10000, lifetime_unit='hr', N=0,
                 material='polypropylene', unit_cost=0.1, surface_area=1):
        Equipment.__init__(self=self, ID=ID, linked_unit=linked_unit, units=units,
                           F_BM=F_BM, lifetime=lifetime, lifetime_unit=lifetime_unit)
        self.N = N
        self.unit_cost = unit_cost
        self.material = material
        self.surface_area = surface_area


    def _design(self):
        design = {
            'Number of membrane modules': self.N,
            'Material of membrane': self.material,
            'Surface area of membrane': self.surface_area
            }
        return design


    def _cost(self):
        return self.unit_cost*self.N*self.surface_area


    @property
    def N(self):
        '''[int] Number of units of the membrane.'''
        return self._N
    @N.setter
    def N(self, i):
        self._N = int(i)
        

class HollowFiberMembrane(Equipment):
    
    _PermSelect_pricing = {
        # model: (area [cm^2], price [USD], min Q_gas [scfm], max Q_gas [scfm])
        'PDMSXA-10': (10, 117, 4e-5, 4e-3), 
        'PDMSXA-1000': (1000, 304, 4e-3, 0.4),
        'PDMSXA-7500': (7500, 541, 2e-2, 1.0),
        'PDMSXA-2.1': (2.1e4, 1045, 4e-2, 2.0)
        }
    
    scfm_to_kmolphr = 0.002641 * 453.59237 * 60 / 1000
    
    def __init__(self, ID='', linked_unit=None,
                 units=dict(), F_BM=1.15, lifetime=5, lifetime_unit='yr',
                 material='polydimethylsiloxane'):
        Equipment.__init__(self=self, ID=ID, linked_unit=linked_unit, units=units,
                           F_BM=F_BM, lifetime=lifetime, lifetime_unit=lifetime_unit)
        self.material = material
        
    def _design(self):
        for ws in self.linked_unit.outs:
            if ws.phase == 'g':
                Q_mol = ws.F_mol # kmol/hr
                break
        key = None
        N = 0
        while not key: 
            N += 1
            Q = Q_mol / N
            for md, data in self._PermSelect_pricing.items():
                Q_max = data[-1] * self.scfm_to_kmolphr
                if Q < Q_max:
                    key = md
        return {'Material': self.material,
                'Number of membrane units': N,
                'Membrane unit model':key}
                
    def _cost(self):
        design = self._design()
        N = design['Number of membrane units']
        key = design['Membrane unit model']
        p = self._PermSelect_pricing[key][1]
        return N * p
        