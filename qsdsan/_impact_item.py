#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems
Copyright (C) 2020, Quantitative Sustainable Design Group

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the UIUC open-source license. Please refer to 
https://github.com/QSD-Group/QSDsan/blob/master/LICENSE.txt
for license details.

Ref (default item CFs):
    [1] Trimmer et al., Navigating Multidimensional Social–Ecological System
        Trade-Offs across Sanitation Alternatives in an Urban Informal Settlement.
        Environ. Sci. Technol. 2020, 54 (19), 12641–12653.
        https://doi.org/10.1021/acs.est.0c03296.
'''


# %%

import pandas as pd
from warnings import warn
from thermosteam import Stream
from thermosteam.utils import copy_maybe
from . import currency, WasteStream, ImpactIndicator
from ._units_of_measure import auom, parse_unit
from .utils.loading import data_path

indicators = ImpactIndicator._indicators
data_path += '_impact_item.xlsx'

__all__ = ('ImpactItem', 'StreamImpactItem')


class ImpactItem:
    '''
    A class for calculation of environmental impacts.
    
    Parameters
    ----------
    ID : str
        ID of the ImpactItem. 
    functional_unit : str
        Functional unit of the ImpactItem.
    price : float
        Price of the item per functional unit.
    price_unit : str
        Unit of the price.
    **indicator_CFs : kwargs, ImpactIndicator or str = float or (float, unit)
        ImpactIndicators and their characteriziation factors.
    
    '''
    
    _items = {}
    _default_data = None
    
    __slots__ = ('_ID', '_functional_unit', '_price', '_CFs')
    
    def __init__(self, ID, functional_unit='kg', price=0., price_unit='', **indicator_CFs):
        
        if ID in ImpactItem._items.keys() and ImpactItem._items[ID] is not self:
            raise ValueError(f'The ID {ID} is in use by {ImpactItem._items[ID]}')
        self._ID = ID
        self._functional_unit = auom(functional_unit)
        self._update_price(price, price_unit)
        self._CFs = {}
        for CF, value in indicator_CFs.items():
            try:
                CF_value, CF_unit = value # unit provided for CF
                self.add_indicator_CF(CF, CF_value, CF_unit)
            except:
                self.add_indicator_CF(CF, value)
        ImpactItem._items[ID] = self
    
    
    # This makes sure it won't be shown as memory location of the object
    def __repr__(self):
        return f'<ImpactItem: {self.ID}>'

    def show(self):
        info = f'ImpactItem      : {self.ID} [per {self.functional_unit}]'
        info += f'\nPrice           : {self.price} {currency}'
        info += '\nImpactIndicators:'
        print(info)
        if len(self.CFs) == 0:
            print(' None')
        else:
            index = pd.Index((i.ID+' ('+i.unit+')' for i in self.indicators))
            df = pd.DataFrame({
                'Characterization factors': tuple(self.CFs.values())
                },
                index=index)
            print(df.to_string())
        
    _ipython_display_ = show


    def _update_price(self, price=0., unit=''):
        if not unit or unit == currency:
            self._price = float(price)
        else:
            converted = auom(unit).convert(float(price), currency)
            self._price = converted

    def add_indicator_CF(self, indicator, CF_value, CF_unit=''):
        if isinstance(indicator, str):
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
        self._CFs[indicator.ID] = CF_value

    def copy(self):
        new = ImpactItem.__new__(ImpactItem)
        for slot in ImpactItem.__slots__:
            value = getattr(self, slot)
            #!!! Not sure if this will cause problem because two objects pointing to the same one
            setattr(new, slot, copy_maybe(value))
        return new
    __copy__ = copy
    
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
                        new = cls.__new__(cls)
                        new.__init__(ID=item,
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
        return tuple(set(i for i in cls._items.values()))
    
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
    def indicators(self):
        ''' [tuple] ImpactIndicators associated with the item.'''
        return tuple(indicators[i] for i in self.CFs.keys())
    
    @property
    def price(self):
        '''Price of the item per functional unit.'''
        return self._price
    @price.setter
    def price(self, price, unit=''):
        self._update_price(price, unit)
    
    @property
    def CFs(self):
        '''[dict] Characterization factors of the item for different impact indicators.'''
        return self._CFs
    @CFs.setter
    def CFs(self, indicator, CF_value, CF_unit=''):
        self.add_indicator_CF(indicator, CF_value, CF_unit)


# %%

class StreamImpactItem(ImpactItem):
    '''
    A class for calculation of environmental impacts associated with WasteStreams.
    
    Parameters
    ----------
    linked_ws : WasteStream or str
        The associated WasteStream for environmental impact calculation.
    **indicator_CFs : kwargs
        ImpactIndicators and their characteriziation factors.
    
    '''

    __slots__ = ('_ID', '_linked_ws', '_functional_unit', '_CFs')

    def __init__(self, linked_ws, **indicator_CFs):
        
        self._linked_ws = None
        self.linked_ws = linked_ws
        ID = self.linked_ws.ID + '_item'
        self._ID = ID
        self._functional_unit = auom('kg')
        self._CFs = {}
        for CF, value in indicator_CFs.items():
            try:
                CF_value, CF_unit = value # unit provided for CF
                self.add_indicator_CF(CF, CF_value, CF_unit)
            except:
                self.add_indicator_CF(CF, value)
        
        ImpactItem._items[ID] = self


    def __repr__(self):
        return f'<StreamImpactItem: WasteStream {self.linked_ws}>'


    def show(self):
        info = f'StreamImpactItem: [per {self.functional_unit}]'        
        info += f'\nLinked to     : {self.linked_ws}'
        info += f'\nPrice           : {self.price} {currency}'
        info += '\nImpactIndicators:'
        print(info)
        if len(self.CFs) == 0:
            print(' None')
        else:
            index = pd.Index((i.ID+' ('+i.unit+')' for i in self.indicators))
            df = pd.DataFrame({
                'Characterization factors': tuple(self.CFs.values())
                },
                index=index)
            # print(' '*18+df.to_string().replace('\n', '\n'+' '*18))
            print(df.to_string())

    
    _ipython_display_ = show

    def copy(self, new_ws):
        new = StreamImpactItem.__new__(StreamImpactItem)
        for slot in StreamImpactItem.__slots__:
            if slot == '_linked_ws': continue
            value = getattr(self, slot)
            #!!! Not sure if this will cause problem because two objects pointing to the same one
            setattr(new, slot, copy_maybe(value))
        new.linked_ws = new_ws
        return new
    __copy__ = copy


    @property
    def linked_ws(self):
        '''[WasteStream] or [str] The associated WasteStream for environmental impact calculation.'''
        return self._linked_ws
    @linked_ws.setter
    def linked_ws(self, i):
        if self._linked_ws:
            self._linked_ws._impact_item = None
        if not isinstance(i, WasteStream) and not isinstance(i, Stream) and i is not None:
            if isinstance(i, str):
                try:
                    i = getattr(WasteStream.registry, i)
                except:
                    raise ValueError(f'The WasteStream ID {i} not '
                                     'found in <WasteStream>.registry')
            else:
                raise TypeError('linked_ws must be a WasteStream or '
                                f'the ID of WasteStream, not {type(i).__name__}.')
        if i is not None:
            if i.impact_item and i.impact_item is not i:
                msg = f'The origin StreamImpactItem linked to WasteStream {i} ' \
                    'is replaced with the current one.'
                warn(message=msg, stacklevel=3)
                i.impact_item.linked_ws = None
            i._impact_item = self
        self._linked_ws = i


    @property
    def functional_unit(self):
        '''[str] Functional unit of the item, set to 'kg'.'''
        return self._functional_unit.units
    
    @property
    def price(self):
        '''[float] Price of the linked WasteStream.'''
        if self.linked_ws:
            return self.linked_ws.price
        else: return 0.



ImpactItem.load_default_items()








