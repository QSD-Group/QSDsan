#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Sanitation Explorer: Sustainable design of non-sewered sanitation technologies
Copyright (C) 2020, Sanitation Explorer Development Group

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the UIUC open-source license. Please refer to 
https://github.com/QSD-for-WaSH/sanitation/blob/master/LICENSE.txt
for license details.

'''


# %%

from . import currency, ImpactIndicator, ImpactItem
from ._units_of_measure import auom
from .utils.formatting import format_number as f_num

indicators = ImpactIndicator._indicators

__all__ = ('Construction',)


class Construction:
    '''Construction cost and environmental impacts'''

    __slots__ = ('_item', '_quantity', '_price')
    
    def __init__(self, item=None, quantity=0., quantity_unit='', price=0., price_unit=currency):
        self.item = item
        self._update_quantity(quantity, quantity_unit)
        self._update_price(price, price_unit)

    def _update_quantity(self, quantity=0., unit=''):
        if not unit or unit == self.item._functional_unit:
            self._quantity = float(quantity)
        else:
            converted = auom(unit).convert(float(quantity), self.item.functional_unit)
            self._quantity = converted

    def _update_price(self, price=0., unit=''):
        if not unit or unit == currency:
            self._price = float(price)
        else:
            converted = auom(unit).convert(float(price), currency)
            self._price = converted
            
    def __repr__(self):
        item = self.item
        impacts = self.impacts
        info = f'Construction: {item.ID}'
        info += f'\n Quantity     : {f_num(self.quantity)} {self.item.functional_unit}'
        info += f'\n Total cost   : {f_num(self.cost)} {currency}'
        info += '\n Total impacts:'
        if len(impacts) == 0:
            info += ' None'        
        else:
            for indicator in impacts.keys():
                formated = f_num(impacts[indicator])
                unit = indicators[indicator].unit
                info += f'\n     {indicator}: {formated} {unit}'
        return info
        
    @property
    def item(self):
        '''[ImpactItem] Item associated with this construction.'''
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
    def quantity(self):
        '''[float] Quantity of this construction item.'''
        return self._quantity
    @quantity.setter
    def quantity(self, quantity, unit=''):
        self._update_quantity(quantity, unit)

    @property
    def price(self):
        '''[float] Unit price of the item.'''
        return self._price
    @price.setter
    def price(self, price, unit=''):
        self._update_price(price, unit)

    @property
    def cost(self):
        '''[float] Total cost of this construction item.'''
        return self.quantity*self.price

    @property
    def impacts(self):
        '''[dict] Total impacts of this construction item.'''
        impacts = {}
        for indicator, CF in self.item.CFs.items():
            impacts[indicator.ID] = self.quantity*CF
        return impacts


















    
    
    
    