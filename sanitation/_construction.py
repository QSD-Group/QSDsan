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

indicators = ImpactIndicator.indicators
_parse_unit = ImpactIndicator._parse_unit
data_path += 'lca_data/_construction_item.xlsx'

__all__ = ('ConstructionItem', 'Construction')


class ConstructionItem:
    '''For construction material and activity items.'''
    
    items = {}
    _default_data = None
    
    __slots__ = ('_ID', '_kind', '_functional_unit', '_price', '_CFs')
 
    def __new__(cls, ID, kind='material', functional_unit='kg', price=0.):
        if ID in cls.items.keys():
            raise ValueError(f'The item {ID} exists, choose another ID or use '
                             'the previously created item in <ConstructionItem>.items')
        self = super().__new__(cls)
        self._ID = ID
        self._kind = kind
        self._functional_unit = getattr(ureg, functional_unit)
        self._CFs = {}
        self._price = price
        cls.items[ID] = self
        return self
    
    # This makes sure it won't be shown as memory location of the object
    def __repr__(self):
        return f'<ConstructionItem: {self.ID}>'

    
    #TODO: DataFrame can make it look nicer
    def show(self):
        info = f'ConstructionItem: {self.ID} [per {self.functional_unit}]'
        info += f'\n price: {self.price}'
        info += '\n CFs:'
        CFs = self.CFs
        if len(CFs) == 0:
            info += ' None'
        else:
            for indicator in self.CFs.keys():
                #!!! Need to add units for different categories
                info += f'\n     {indicator.ID}: {CFs[indicator]} [{indicator.unit}]'
        print(info)    
    
    _ipython_display_ = show


    def add_indicator_CF(self, indicator, CF_value, CF_unit):
        if type(indicator) is str:
            indicator = indicators[indicator]
        try: CF_unit2 = CF_unit.replace(' eq', '-eq')
        except: pass
        if CF_unit != indicator.unit and CF_unit2 != indicator.unit:
            try:
                CF_value = RelativeUnitsOfMeasure(_parse_unit(CF_unit)[0]). \
                    convert(CF_value, indicator._ureg_unit)
            except:
                raise ValueError(f'Conversion of the given unit {CF_unit} to '
                                 f'the defaut unit {indicator.unit} is not supported.')
        self._CFs[indicator] = CF_value

    
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
                    if item in cls.items.keys():
                        items[item] = cls.items[item]
                    else:
                        new = cls.__new__(cls, ID=item,
                                          kind=data.loc[item]['kind'],
                                          functional_unit=data.loc[item]['functional_unit'])
                        items[item] = new
            else:
                for item in data.index:
                    old = items[item]
                    old.add_indicator_CF(indicator=sheet,
                                         CF_value=float(data.loc[item]['expected']),
                                         CF_unit=data.loc[item]['unit'])
        cls._default_data = data_file
    
    
    @property
    def ID(self):
        '''ID of the item.'''
        return self._ID
    
    @property
    def functional_unit(self):
        '''Functional unit of the item.'''
        return self._functional_unit
    @functional_unit.setter
    def functional_unit(self, i):
        self._functional_unit = getattr(ureg, i)
    
    @property
    def price(self):
        '''[float] Price of the material or activity per functional unit.'''
        return self._price
    @price.setter
    def price(self, i):
        self._price = float(i)
    
    @property
    def CFs(self):
        '''[dict] Characterization factors of the material or activity for different impact indicators.'''
        return self._CFs













    
    
    
    
class Construction:
    '''
    A class to calculate the cost and environmental impacts associated with
    construction materials and activity items.
    
    '''

    __slots__ = ('quantity',)
    
    def __init__(self):
        self.quantity = 0.



ConstructionItem.load_default_items()








    
    
    
    