#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''


# %%

import pandas as pd
from warnings import warn
from thermosteam.utils import copy_maybe
from . import currency, SanStream, WasteStream, ImpactIndicator
from ._units_of_measure import auom, parse_unit
from .utils.formatting import format_number as f_num

indicators = ImpactIndicator._indicators
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
    indicator_CFs : kwargs, :class:`ImpactIndicator` or str = float or (float, unit)
        Impact indicators and their characteriziation factors.
    
    '''
    
    _items = {}
    
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

    def deregister(self):
        '''Remove this :class:`ImpactItem` from the record.'''
        ID = self.ID
        self._items.pop(ID)
        print(f'The impact item "{ID}" has been removed from the record.')
    
    
    @classmethod
    def load_items_from_excel(cls, path):
        '''
        Load impact items from an Excel file.
        
        This Excel should have multiple sheets:
            
            - The "info" sheet should have two columns: "ID" (e.g., Cement) \
            and "functional_unit" (e.g., kg) of different impact items.
            
            - The remaining sheets should contain characterization factors of \
            impact indicators.
            
                - Name of the sheet should be the ID (e.g., GlobalWarming) or \
                alias (e.g., GWP) of the indicator.
                
                - Each sheet should have at least two columns: "unit" (e.g., kg CO2-eq) \
                and "expected" (values) of the CF.
                
                - You can also have additional columns to be used for other purpose \
                (e.g., uncertainty analysis).
        
        .. note::
            
            This function is just one way to batch-load impact items,
            you can always write your own function that fits your datasheet format,
            as long as it provides all the information to construct new impact items.
        
        
        Parameters
        ----------
        path : str
            Complete path of the Excel file.

        Tip
        ---
        Refer to the `Bwaise system <https://github.com/QSD-Group/EXPOsan/tree/main/exposan/bwaise/data>`_
        in the ``Exposan`` repository for a sample file.
        '''
        if not (path.endswith('.xls') or path.endswith('.xlsx')):
            raise ValueError('Only Excel files ends with ".xlsx" or ".xls" can be interpreted.')
        
        data_file = pd.ExcelFile(path, engine='openpyxl')
        items = cls._items
        for sheet in data_file.sheet_names:
            data = data_file.parse(sheet, index_col=0)

            if sheet == 'info':
                for item in data.index:
                    if item in items.keys():
                        raise ValueError(f'The impact item "{item}" has been added.')
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
    linked_stream : :class:`SanStream`
        The associated :class:`SanStream` for environmental impact calculation.
    source : :class:`StreamImpactItem`
        If provided, all attributions and properties of this
        :class:`StreamImpactItem` will be copied from the provided source.
    indicator_CFs : kwargs
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
            kind = type(self.linked_stream).__name__
            return f'<StreamImpactItem: {kind} {self.linked_stream}>'
        else:
            return '<StreamImpactItem: no linked stream>'

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
            raise ValueError('source can only be a StreamImpactItem, ' \
                             f'not a {type(i).__name__}.')
        self._source = i

    @property
    def linked_stream(self):
        '''
        [SanStream] The associated :class:`SanStream` for environmental impact calculation,
        can be set by either the :class:`SanStream` object or its ID.
        '''
        return self._linked_stream
        
    @linked_stream.setter
    def linked_stream(self, new_s):
        if new_s and not isinstance(new_s, SanStream):
            if isinstance(new_s, str):
                try:
                    new_s = getattr(SanStream.registry, new_s)
                except:
                    try:
                        new_s = getattr(WasteStream.registry, new_s)
                    except:
                        raise ValueError(f'The ID "{new_s}" not found in registry.')
            else:
                raise TypeError('`linked_stream` must be a `SanStream` or '
                                f'the ID of a `SanStream`, not {type(new_s).__name__}.')
        
        if self._linked_stream:
            old_s = self._linked_stream
            self._linked_stream.impact_item = None
            warn(f'`ImpactItem` {self.ID} is unlinked from {old_s.ID} and ' \
                 f'linked to {new_s.ID}.', stacklevel=2)
        if new_s:
            if hasattr(self, '_ID'):
                if new_s.impact_item and new_s.impact_item.ID != self.ID:
                    msg = f'The original `StreamImpactItem` linked to stream {new_s} ' \
                        f'is replaced with {self}.'
                    warn(message=msg, stacklevel=2)
            new_s._impact_item = self
        self._linked_stream = new_s

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
        '''[float] Price of the linked stream.'''
        if self.linked_stream:
            return self.linked_stream.price
        else: return 0.







