#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/master/LICENSE.txt
for license details.
'''


# %%

import pandas as pd
from thermosteam.utils import copy_maybe
from . import currency, ImpactIndicator, ImpactItem
from ._units_of_measure import auom
from .utils.formatting import format_number as f_num

indicators = ImpactIndicator._indicators

__all__ = ('Construction',)


class Construction:
    '''Construction cost and environmental impacts.'''

    __slots__ = ('_item', '_quantity')
    
    def __init__(self, item=None, quantity=0., unit=''):
        self.item = item
        self._update_quantity(quantity, unit)

    def _update_quantity(self, quantity=0., unit=''):
        if not unit or unit == self.item.functional_unit:
            self._quantity = float(quantity)
        else:
            converted = auom(unit).convert(float(quantity), self.item.functional_unit)
            self._quantity = converted
           
    def __repr__(self):
        return f'<Construction: {self.item.ID}>'
    
    def show(self):
        '''Show basic information about this :class:`Construction` object.'''
        item = self.item
        impacts = self.impacts
        info = f'Construction : {item.ID}'
        info += f'\nQuantity     : {f_num(self.quantity)} {self.item.functional_unit}'
        info += f'\nTotal cost   : {f_num(self.cost)} {currency}'
        info += '\nTotal impacts:'
        print(info)
        if len(impacts) == 0:
            print(' None')
        else:
            index = pd.Index((i.ID+' ('+i.unit+')' for i in self.indicators))
            df = pd.DataFrame({
                'Impacts': tuple(self.impacts.values())
                },
                index=index)
            # print(' '*15+df.to_string().replace('\n', '\n'+' '*15))
            print(df.to_string())
        
    _ipython_display_ = show
    
    def copy(self):
        new = Construction.__new__(Construction)
        for slot in Construction.__slots__:
            value = getattr(self, slot)
            #!!! Not sure if this will cause problem because two objects pointing to the same one
            setattr(new, slot, copy_maybe(value))
        return new
    __copy__ = copy
    
    @property
    def item(self):
        '''[ImpactItem] Item associated with this construction activity.'''
        return self._item
    @item.setter
    def item(self, i):
        if isinstance(i, str):
            i = ImpactItem._items[i]
        elif i is not ImpactItem:
            raise TypeError('Only <ImpactItem> or  <ImpactItem>.ID can be set, '
                            f'not {type(i).__name__}.')
        self._item = i

    @property
    def indicators(self):
        ''' [tuple] ImpactIndicators associated with the construction item.'''
        return self.item.indicators

    @property
    def quantity(self):
        '''[float] Quantity of this construction item.'''
        return self._quantity
    @quantity.setter
    def quantity(self, quantity, unit=''):
        self._update_quantity(quantity, unit)

    @property
    def price(self):
        '''[float] Unit price of the item.'''
        return self.item.price

    @property
    def cost(self):
        '''[float] Total cost of this construction item.'''
        return self.quantity*self.price

    @property
    def impacts(self):
        '''[dict] Total impacts of this construction item.'''
        impacts = {}
        for indicator, CF in self.item.CFs.items():
            impacts[indicator] = self.quantity*CF
        return impacts


















    
    
    
    