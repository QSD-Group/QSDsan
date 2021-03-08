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
from warnings import warn
from thermosteam.utils import copy_maybe
from . import currency, WasteStream, ImpactIndicator
from ._units_of_measure import auom, parse_unit
from .utils.loading import data_path
from .utils.formatting import format_number as f_num

indicators = ImpactIndicator._indicators
data_path += '_impact_item.xlsx'
isinstance = isinstance
getattr = getattr

__all__ = ('ImpactItem', 'StreamImpactItem')

def check_source(item):
    if item.source:
        raise ValueError(f'This ImpactItem is copied from {item.source.ID}, '
                         'value cannot be set.')

class ImpactItem:
    '''
    A class for calculation of environmental impacts.
    
    Parameters
    ----------
    ID : str
        ID of the impact item. If no ID is provided, this item will not be
        saved in the ImpactItem dict.
    functional_unit : str
        Functional unit of the impact item.
    price : float
        Price of the item per functional unit.
    price_unit : str
        Unit of the price.
    source : ImpactItem
        If provided, all attributions and properties of this impact item will
        be copied from the provided source.
    **indicator_CFs : kwargs, :class:`ImpactIndicator` or str = float or (float, unit)
        Impact indicators and their characteriziation factors.
    
    '''
    
    _items = {}
    _default_data = None
    
    __slots__ = ('_ID', '_functional_unit', '_price', '_CFs', '_source')
    
    def __init__(self, ID=None, functional_unit='kg', price=0., price_unit='',
                 source=None, **indicator_CFs):
        
        self._ID = ID
        if source:
            self.source = source
        else:
            self._source = None
            self._functional_unit = auom(functional_unit)
            self._update_price(price, price_unit)
            self._CFs = {}
            for CF, value in indicator_CFs.items():
                try:
                    CF_value, CF_unit = value # unit provided for CF
                    self.add_indicator_CF(CF, CF_value, CF_unit)
                except:
                    self.add_indicator_CF(CF, value)
            if ID:
                if ID in ImpactItem._items.keys():
                    old = ImpactItem._items[ID]
                    for i in old.__slots__:
                        if not getattr(old, i) == getattr(self, i):
                            raise ValueError(f'The ID {ID} is in use by {ImpactItem._items[ID]}, '\
                                             'use another ID instead.')
                else:
                    ImpactItem._items[ID] = self
    
    # This makes sure it won't be shown as memory location of the object
    def __repr__(self):
        return f'<ImpactItem: {self.ID}>'

    def show(self):
        '''Show basic information of this ``ImpactItem`` object'''
        info = f'ImpactItem      : {self.ID} [per {self.functional_unit}]'
        if self.source:
            info += f'\nSource          : {self.source.ID}'
        info += f'\nPrice           : {f_num(self.price)} {currency}'
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
        check_source(self)
        if not unit or unit == currency:
            self._price = float(price)
        else:
            converted = auom(unit).convert(float(price), currency)
            self._price = converted

    def add_indicator_CF(self, indicator, CF_value, CF_unit=''):
        '''Add an indicator charactorization factor for this :class:`ImpactItem` object.'''
        check_source(self)
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

    def copy(self, new_ID=None, set_as_source=False):
        '''
        Return a new :class:`ImpactItem` object with the same settings.
        Set the original :class:`ImpactItem` as the source for the new one if
        `set_as_source` is True.
        '''        
        new = ImpactItem.__new__(ImpactItem)
        new.ID = new_ID
        if set_as_source:
            new.source = self
        else:
            for slot in ImpactItem.__slots__:
                value = getattr(self, slot)
                setattr(new, slot, copy_maybe(value))
        return new
    __copy__ = copy

    
    
    #!!! Are the values GWP100 from ref [1]?
    @classmethod
    def load_default_items(cls, path=data_path):
        '''
        Load all default indicators as in /data/_impact_item.xlsx from Trimmer et al. [1]_
        
        References
        ----------
        .. [1] Trimmer et al., Navigating Multidimensional Social–Ecological System
            Trade-Offs across Sanitation Alternatives in an Urban Informal Settlement.
            Environ. Sci. Technol. 2020, 54 (19), 12641–12653.
            https://doi.org/10.1021/acs.est.0c03296.
        
        '''
        if cls._default_data is not None:
            data_file = cls._default_data
        else: data_file = pd.ExcelFile(data_path, engine='openpyxl')
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
    def get_item(cls, ID):
        '''Get an item by its ID.'''
        return cls._items[ID]
    
    @classmethod
    def get_all_items(cls):
        '''Get a tuple of all impact items'''
        return tuple(set(i for i in cls._items.values()))
    
    @property
    def source(self):
        '''
        [ImpactItem] If provided, all attributions and properties of this
        impact item will be copied from the provided source.
        '''
        return self._source
    @source.setter
    def source(self, i):
        if not isinstance(i, ImpactItem):
            raise ValueError('`source` can only be an `ImpactItem`, '
                             f'not {type(i).__name__}.')
        self._source = i
    
    @property
    def ID(self):
        '''
        [str] ID of the item. If no ID is provided, this item will not be
        saved in the ImpactItem dict.
        '''
        return self._ID
    @ID.setter
    def ID(self, i):
        self._ID = i
    
    @property
    def functional_unit(self):
        '''[str] Functional unit of the item.'''
        if self.source: return self.source._functional_unit.units
        return self._functional_unit.units
    @functional_unit.setter
    def functional_unit(self, i):
        check_source(self)
        self._functional_unit = auom(i)
    
    @property
    def indicators(self):
        ''' [tuple] :class:`ImpactIndicator` objects associated with the item.'''
        return tuple(indicators[i] for i in self.CFs.keys())
    
    @property
    def price(self):
        '''Price of the item per functional unit.'''
        if self.source: return self.source.price
        return self._price
    @price.setter
    def price(self, price, unit=''):
        self._update_price(price, unit)
    
    @property
    def CFs(self):
        '''[dict] Characterization factors of the item for different impact indicators.'''
        if self.source: return self.source._CFs
        return self._CFs
    @CFs.setter
    def CFs(self, indicator, CF_value, CF_unit=''):
        check_source(self)
        self.add_indicator_CF(indicator, CF_value, CF_unit)


# %%

class StreamImpactItem(ImpactItem):
    '''
    A class for calculation of environmental impacts associated with chemical
    inputs and emissions.
    
    Parameters
    ----------
    linked_stream : WasteStream
        The associated :class:`WasteStream` for environmental impact calculation.
    **indicator_CFs : kwargs
        ImpactIndicators and their characteriziation factors.
    
    '''

    __slots__ = ('_ID', '_linked_stream', '_functional_unit', '_CFs', '_source')

    def __init__(self, ID=None, linked_stream=None, source=None,
                 **indicator_CFs):

        self._linked_stream = None
        self.linked_stream = linked_stream
        if not ID and linked_stream:
            ID = self.linked_stream.ID + '_item'
        self._ID = ID
        if source:
            self.source = source
        else:
            self._source = None
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
        if self.linked_stream:
            return f'<StreamImpactItem: WasteStream {self.linked_stream}>'
        else:
            return '<StreamImpactItem: no linked WasteStream>'

    def show(self):
        '''Show basic information about this :class:`StreamImpactItem` object.'''
        info = f'StreamImpactItem: [per {self.functional_unit}]'    
        if self.linked_stream:
            info += f'\nLinked to       : {self.linked_stream}'
        else:
            info += '\nLinked to       : None'
        if self.source:
            info += f'\nSource          : {self.source.ID}'
        info += f'\nPrice           : {f_num(self.price)} {currency}'
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


    def copy(self, new_ID=None, stream=None, set_as_source=False):
        '''
        Return a new :class:`StreamImpactItem` object with the same settings,
        link the new :class:`StreamImpactItem` to the provided stream if given,
        set the original :class:`StreamImpactItem` as the source for the new one if
        `set_as_source` is True.
        '''
        new = StreamImpactItem.__new__(StreamImpactItem)
        new.ID = new_ID
        new._linked_stream = None # initiate this attribute
        if set_as_source:
            new.source = self
        else:
            for slot in StreamImpactItem.__slots__:
                if slot == '_linked_stream': continue
                value = getattr(self, slot)
                setattr(new, slot, copy_maybe(value))
        if stream:
            stream.impact_item = new
            if not new_ID:
                new.ID = f'{stream.ID}_item'
        return new
    __copy__ = copy

    @property
    def source(self):
        '''
        [StreamImpactItem] If provided, all attributions and properties of this
        :class:`StreamImpactItem` will be copied from the provided source.
        
        .. note::
            Since the price is copied from the price of the `linked_stream`, it
            can be different form the source.
        '''
        return self._source
    @source.setter
    def source(self, i):
        if not isinstance(i, StreamImpactItem):
            raise ValueError('source can only be a StreamImpactItem, '
                             f'not a {type(i).__name__}.')
        self._source = i

    @property
    def linked_stream(self):
        '''
        [WasteStream] The associated :class:`WasteStream` for environmental impact calculation,
        can be set by either the :class:`WasteStream` object or its ID
        '''
        try: return self._linked_stream
        except: breakpoint()
    @linked_stream.setter
    def linked_stream(self, new_ws):
        isa = isinstance
        if new_ws and not isa(new_ws, WasteStream):
            if isa(new_ws, str):
                try:
                    new_ws = getattr(WasteStream.registry, new_ws)
                except:
                    raise ValueError(f'The ``WasteStream` ID {new_ws} not '
                                     'found in <WasteStream>.registry')
            else:
                raise TypeError('`linked_stream` must be a `WasteStream` or '
                                f'the ID of a `WasteStream`, not {type(new_ws).__name__}.')
        if self._linked_stream:
            old_ws = self._linked_stream
            self._linked_stream.impact_item = None
            warn(f'`ImpactItem` {self.ID} is unlinked from {old_ws.ID} and ' \
                 f'linked to {new_ws.ID}.', stacklevel=2)
        if new_ws:
            if hasattr(self, '_ID'):
                if new_ws.impact_item and new_ws.impact_item.ID != self.ID:
                    msg = f'The original `StreamImpactItem` linked to `WasteStream` {new_ws} ' \
                        f'is replaced with {self}.'
                    warn(message=msg, stacklevel=2)
            new_ws._impact_item = self
        self._linked_stream = new_ws

    @property
    def ID(self):
        '''
        [str] ID of the item. If no ID is provided but its `linked_stream` is provided,
        the ID will be set as ID of the `linked_stream` with a suffix '_item'.
        '''
        return self._ID
    @ID.setter
    def ID(self, i):
        self._ID = i

    @property
    def functional_unit(self):
        '''[str] Functional unit of the item, set to 'kg'.'''
        return auom('kg')
    
    @property
    def price(self):
        '''[float] Price of the linked `WasteStream`.'''
        if self.linked_stream:
            return self.linked_stream.price
        else: return 0.



ImpactItem.load_default_items()








