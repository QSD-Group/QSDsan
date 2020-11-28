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

Ref:
    [1] Trimmer et al., Navigating Multidimensional Social–Ecological System
        Trade-Offs across Sanitation Alternatives in an Urban Informal Settlement.
        Environ. Sci. Technol. 2020, 54 (19), 12641–12653.
        https://doi.org/10.1021/acs.est.0c03296.

'''



# %%

import pandas as pd
from . import ImpactIndicator
from ._units_of_measure import ureg, RelativeUnitsOfMeasure
from .utils.loading import data_path
from .utils.formatting import format_number as f_num

indicators = ImpactIndicator._indicators
_parse_unit = ImpactIndicator._parse_unit
data_path += 'lca_data/_construction_item.xlsx'

__all__ = ('ConstructionItem', 'Construction')


class ConstructionItem:
    '''For construction material and activity items.'''
    
    _items = {}
    _default_data = None
    
    __slots__ = ('_ID', '_kind', '_functional_unit', '_price', '_price_unit', '_CFs')
 
    def __new__(cls, ID, kind='material', functional_unit='kg', price=0.,
                price_unit='USD', **indicator_CFs):
        if ID in cls._items.keys():
            raise ValueError(f'The ID {ID} is in use by {cls._items[ID]}')
        self = super().__new__(cls)
        self._ID = str(ID)
        self._kind = str(kind)
        self._functional_unit = getattr(ureg, str(functional_unit))
        self._CFs = {}
        if str(price) == 'nan':
            self._price = 0.            
        else:
            self._price = price
        if str(price_unit) == 'nan':
            self._price_unit = ureg.USD
        else:
            self._price_unit = getattr(ureg, price_unit)
        for CF, value in indicator_CFs.items():
            try:
                CF_value, CF_unit = value # unit provided for CF
                self.add_indicator_CF(CF, CF_value, CF_unit)
            except:
                self.add_indicator_CF(CF, value)
        
        cls._items[ID] = self
        return self
    
    # This makes sure it won't be shown as memory location of the object
    def __repr__(self):
        return f'<ConstructionItem: {self.ID}>'

    
    #TODO: DataFrame can make it look nicer
    def show(self):
        info = f'ConstructionItem: {self.ID} [per {self.functional_unit}]'
        #!!! Make it possible to customerize price unit
        info += f'\n price: {f_num(self.price)} [self.price_unit]'
        info += '\n CFs:'
        CFs = self.CFs
        if len(CFs) == 0:
            info += ' None'
        else:
            for indicator in self.CFs.keys():
                #!!! Currently Eutrofication values are fake
                info += f'\n     {indicator.ID}: {f_num(CFs[indicator])} {indicator.unit}'
        print(info)    
    
    _ipython_display_ = show


    def add_indicator_CF(self, indicator, CF_value, CF_unit=''):
        if type(indicator) is str:
            indicator = indicators[indicator]
        try: CF_unit2 = CF_unit.replace(' eq', '-eq')
        except: pass
        if CF_unit and CF_unit != indicator.unit and CF_unit2 != indicator.unit:
            try:
                CF_value = RelativeUnitsOfMeasure(_parse_unit(CF_unit)[0]). \
                    convert(CF_value, indicator._ureg_unit)
            except:
                raise ValueError(f'Conversion of the given unit {CF_unit} to '
                                 f'the defaut unit {indicator.unit} is not supported.')
        self._CFs[indicator] = CF_value
    
    def update_indicator_CF(self, indicator, CF_value, CF_unit):
        self.add_indicator_CF(indicator, CF_value, CF_unit)

    
    #!!! Are the values GWP100 from ref [1]?
    @classmethod
    def load_default_items(cls, path=data_path):
        if cls._default_data is not None:
            data_file = cls._default_data
        else: data_file = pd.ExcelFile(data_path)
        items = {}
        for sheet in data_file.sheet_names:
            data = data_file.parse(sheet, index_col=0)
            if sheet == 'info':
                for item in data.index:
                    if item in cls._items.keys():
                        items[item] = cls._items[item]
                    else:
                        new = cls.__new__(cls, ID=item,
                                          kind=data.loc[item]['kind'],
                                          functional_unit=data.loc[item]['functional_unit'],
                                          price=data.loc[item]['price'],
                                          price_unit=data.loc[item]['price_unit'])
                        items[item] = new
            else:
                for item in data.index:
                    old = items[item]
                    old.add_indicator_CF(indicator=sheet,
                                         CF_value=float(data.loc[item]['expected']),
                                         CF_unit=data.loc[item]['unit'])
        cls._default_data = data_file
    
    @classmethod
    def get_all_items(cls):
        return set(i for i in cls._items.values())
    
    @property
    def ID(self):
        '''[str] ID of the item.'''
        return self._ID
    
    @property
    def functional_unit(self):
        '''[str] Functional unit of the item.'''
        try:
            return self._functional_unit.format_babel()
        except:
            return self._functional_unit
    @functional_unit.setter
    def functional_unit(self, i):
        self._functional_unit = getattr(ureg, i)
    
    def _update_price(self, price=0., unit=''):
        if not unit or unit == self.price_unit:
            self._price = float(price)
        else:
            converted = RelativeUnitsOfMeasure(getattr(ureg, unit)). \
                    convert(float(price), self.item.price_unit)
            self._price = converted

    @property
    def price(self):
        '''[float] Price of the material or activity per functional unit.'''
        return self._price
    @price.setter
    def price(self, price, unit=''):
        self._update_price(price, unit)
        
    @property
    def price_unit(self):
        '''[float] Unit of the item price.'''
        try:
            return self._price_unit.format_babel()
        except:
            return self._price_unit
    @price.setter
    def price(self, i):
        self._price = float(i)
    
    @property
    def CFs(self):
        '''[dict] Characterization factors of the material or activity for different impact indicators.'''
        return self._CFs













    
    
    
    
class Construction:
    '''
    A class for the calculation of cost and environmental impacts associated with
    construction materials and activity items.
    
    '''

    __slots__ = ('_item', '_quantity',)
    
    def __init__(self, item=None, quantity=0., unit=''):
        self._item = item
        self._quantity = quantity
        self._update_quantity(quantity, unit)

    def _update_quantity(self, quantity=0., unit=''):
        if not unit or unit == self.item._functional_unit:
            self._quantity = float(quantity)
        else:
            converted = RelativeUnitsOfMeasure(getattr(ureg, unit)). \
                    convert(float(quantity), self.item._functional_unit)
            self._quantity = converted
            
    def __repr__(self):
        item = self.item
        impacts = self.impacts
        info = f'Construction: {item.ID}'
        info += f'\n Quantity    : {f_num(self.quantity)} {item.functional_unit}'
        info += f'\n Total cost  : {f_num(self.cost)} {item.price_unit}'
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
    def quantity(self):
        '''[float] Quantity of the construction item.'''
        return self._quantity
    @quantity.setter
    def quantity(self, quantity, unit=''):
        self._update_quantity(quantity, unit)
        
    @property
    def item(self):
        '''[ConstructionItem] Item associated with this construction.'''
        return self._item
    @item.setter
    def item(self, i):
        if i is not ConstructionItem:
            raise TypeError('Only <ConstructionItem> can be set as item, '
                            f'not {type(i).__name__}.')
        self._item = i
  
    @property
    def unit(self):
        '''[str] Unit of the construction, the same as the functional unit of the ConstructionItem.'''
        return self.item.functional_unit  
  
    @property
    def cost(self):
        '''[float] Total cost of the item during the construction.'''
        return self.quantity*self.item.price

    @property
    def impacts(self):
        '''[dict] Unit of the construction, the same as the functional unit of the ConstructionItem.'''
        impacts = {}
        for indicator, CF in self.item.CFs.items():
            impacts[indicator.ID] = self.quantity*CF
        return impacts









ConstructionItem.load_default_items()








    
    
    
    