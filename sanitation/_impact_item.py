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

import pandas as pd
from . import ImpactIndicator
from ._units_of_measure import auom, parse_unit
from .utils.loading import data_path
from .utils.formatting import format_number as f_num

indicators = ImpactIndicator._indicators
data_path += '_impact_item.xlsx'

__all__ = ('ImpactItem',)


class ImpactItem:
    '''A class for calculation of environmental impacts.'''
    
    _items = {}
    _default_data = None
    
    __slots__ = ('_ID', '_kind', '_functional_unit', '_CFs', '_CF_values')
 
    def __new__(cls, ID, kind='material', functional_unit='kg', **indicator_CFs):
        if ID in cls._items.keys():
            raise ValueError(f'The ID {ID} is in use by {cls._items[ID]}')
        self = super().__new__(cls)
        self._ID = str(ID)
        self._kind = str(kind)
        self._functional_unit = auom(functional_unit)
        self._CFs = {}
        self._CF_values = {}
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
        return f'<ImpactItem: {self.ID}>'

    
    #TODO (maybe): DataFrame can make it look nicer
    def show(self):
        info = f'ImpactItem: {self.ID} [per {self.functional_unit}]'
        info += '\n Characterization factors:'
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
                CF_value = auom(parse_unit(CF_unit)[0]). \
                    convert(CF_value, indicator._ureg_unit.units)
            except:
                raise ValueError(f'Conversion of the given unit {CF_unit} to '
                                 f'the defaut unit {indicator.unit} is not supported.')
        self._CFs[indicator] = CF_value
        self._CF_values[indicator.ID] = CF_value

    
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
                                          functional_unit=data.loc[item]['functional_unit'])
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
        return self._functional_unit.units
    @functional_unit.setter
    def functional_unit(self, i):
        self._functional_unit = auom(i)
    
    @property
    def CFs(self):
        '''
        [dict] Characterization factors of the item for different impact indicators,
        dict keys are ImpactIndicator objects.
        '''
        return self._CFs
    @CFs.setter
    def CFs(self, indicator, CF_value, CF_unit=''):
        self.add_indicator_CF(indicator, CF_value, CF_unit)

    @property
    def CF_values(self):
        '''
        [dict] Characterization factors of the item for different impact indicators,
        dict keys are the IDs of ImpactIndicator objects.
        '''
        return self._CF_values
    @CF_values.setter
    def CF_values(self, indicator, CF_value, CF_unit=''):
        self.add_indicator_CF(indicator, CF_value, CF_unit)


ImpactItem.load_default_items()











