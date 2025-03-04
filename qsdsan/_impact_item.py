#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Yalin Li <mailto.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''


# %%

import sys
import pandas as pd
from warnings import warn
from thermosteam.utils import registered
from . import currency, CHECK_IMPACT_ITEM_CONSISTENCY, SanStream, WasteStream, ImpactIndicator
from .utils import (
    auom, parse_unit, copy_attr,
    format_number as f_num
    )


__all__ = ('ImpactItem', 'StreamImpactItem')


def raise_inconsistency_error(ID, attr, exist_val, new_val):
    warn(f'The impact item "{ID}" already exists '
         f'and has a different {attr} ({exist_val}) '
         f'than the provided one ({new_val}).')


def check_source(item, return_item=False, raise_error=True):
    if hasattr(item, '_source'):
        if item.source:
            if raise_error:
                raise ValueError(f'This impact item is copied from {item.source.ID}, '
                                 'value cannot be set.')
            else: item = item.source
    if return_item: return item


@registered(ticket_name='item')
class ImpactItem:
    '''
    A class for calculation of environmental impacts.

    Parameters
    ----------
    ID : str
        ID of the impact item.
    functional_unit : str
        Functional unit of the impact item.
    price : float
        Price of the item per functional unit.
    source : :class:`ImpactItem`
        If provided, all attributions and properties of this impact item will
        be copied from the provided source.
    indicator_CFs : kwargs
        Impact indicators and their characterization factors (CFs),
        can be in the form of str=float or str=(float, unit).

    Tip
    ---
    :class:`ImpactItem` should be used for environmental impacts associated with
    construction and transportation.
    For impacts associated with streams (e.g., chemicals, wastes, emissions),
    use :class:`StreamImpactItem` instead.


    Examples
    --------
    Firstly make impact indicators.

    >>> import qsdsan as qs
    >>> GWP = qs.ImpactIndicator('GlobalWarming', alias='GWP', unit='kg CO2-eq')
    >>> FEC = qs.ImpactIndicator('FossilEnergyConsumption', alias='FEC', unit='MJ')

    We can make impact items in different ways (numbers are made up).

    >>> Steel = qs.ImpactItem('Steel', 'kg', GWP=2.55)
    >>> Steel.show()
    ImpactItem      : Steel [per kg]
    Price           : None USD
    ImpactIndicators:
                               Characterization factors
    GlobalWarming (kg CO2-eq)                      2.55
    >>> # Unit will be automatically converted to match the unit of the impact indicator
    >>> Electricity = qs.ImpactItem('Electricity', functional_unit='kWh',
    ...                             GWP=(480, 'g CO2-eq'), FEC=(5926, 'kJ'))
    >>> Electricity.show()
    ImpactItem      : Electricity [per kWh]
    Price           : None USD
    ImpactIndicators:
                                  Characterization factors
    GlobalWarming (kg CO2-eq)                         0.48
    FossilEnergyConsumption (MJ)                      5.93
    >>> # Note that 5.93 appearing instead of 5.926 is for nicer print
    >>> Electricity.CFs
    {'GlobalWarming': 0.48, 'FossilEnergyConsumption': 5.926}
    >>> # Get all impact items
    >>> qs.ImpactItem.get_all_items()
    {'Steel': <ImpactItem: Steel>, 'Electricity': <ImpactItem: Electricity>}

    You can make copies of impact items and choose to link to the source or not.

    >>> Steel2 = Steel.copy('Steel2', set_as_source=True)
    >>> Steel2.CFs['GlobalWarming']
    2.55
    >>> Steel3 = Steel.copy('Steel3', set_as_source=False)
    >>> Steel3.CFs['GlobalWarming']
    2.55
    >>> # Once linked, CFs of the copy will update with the source
    >>> Steel.CFs['GlobalWarming'] = 2
    >>> Steel.CFs['GlobalWarming']
    2
    >>> Steel2.CFs['GlobalWarming']
    2
    >>> Steel3.CFs['GlobalWarming']
    2.55
    >>> # Update the copy won't update the source
    >>> Steel2.CFs['GlobalWarming'] = 5
    >>> Steel.CFs['GlobalWarming']
    2
    >>> Steel2.CFs['GlobalWarming']
    2

    Manage the registry.

    >>> qs.ImpactItem.get_all_items()
    {'Steel': <ImpactItem: Steel>,
     'Electricity': <ImpactItem: Electricity>,
     'Steel2': <ImpactItem: Steel2>,
     'Steel3': <ImpactItem: Steel3>}
    >>> Steel2.deregister()
    The impact item "Steel2" has been removed from the registry.
    >>> qs.ImpactItem.get_all_items()
    {'Steel': <ImpactItem: Steel>,
     'Electricity': <ImpactItem: Electricity>,
     'Steel3': <ImpactItem: Steel3>}
    >>> Steel2.register()
    The impact item "Steel2" has been added to the registry.
    >>> qs.ImpactItem.get_all_items()
    {'Steel': <ImpactItem: Steel>,
     'Electricity': <ImpactItem: Electricity>,
     'Steel3': <ImpactItem: Steel3>,
     'Steel2': <ImpactItem: Steel2>}
    >>> qs.ImpactItem.clear_registry()
    All impact items have been removed from the registry.
    >>> qs.ImpactItem.get_all_items()
    {}
    >>> # Clear all registries for testing purpose
    >>> from qsdsan.utils import clear_lca_registries
    >>> clear_lca_registries()
    '''

    _items = {}

    __slots__ = ('_ID', '_functional_unit', '_price', '_CFs', '_source')

    def _format_CF_vals(self, **indicator_CFs):
        CF_dct = {}
        for indicator, value in indicator_CFs.items():
            try: # unit provided for CF
                CF_value, CF_unit = value
                CF_dct[indicator] = CF_value, CF_unit
            except Exception as e: # no unit
                if 'unpack' in str(sys.exc_info()[1]):
                    CF_value, CF_unit = value, ''
                    CF_dct[indicator] = CF_value, CF_unit
                else:
                    raise e
        return CF_dct

    def _register_with_consistency_check(self, ID, functional_unit, price, source, CF_dct):
        '''Check consistency in attr values.'''
        exist_item = ImpactItem.get_item(ID)
        if exist_item and CHECK_IMPACT_ITEM_CONSISTENCY:
            self._ID = ID
            val_dct = {
                'functional_unit': functional_unit,
                'price': price,
                'source': source,
                }
            for attr, new_val in val_dct.items():
                exist_val = getattr(exist_item, attr)
                if exist_val != new_val and str(exist_val) != str(new_val): # functional_unit as auom
                    raise_inconsistency_error(ID, attr, exist_val, new_val)

            for indicator, (value, unit) in CF_dct.items():
                indicator_ID = ImpactIndicator.get_indicator(indicator).ID
                exist_CF = exist_item.CFs.get(indicator_ID)
                if exist_CF: # indicator already added
                    if exist_CF != value:
                        raise_inconsistency_error(ID, 'CF value', exist_CF, value)
                else: # add the CF value for the new indicator
                    exist_item.add_indicator(indicator, value, unit)
            for slot in ImpactItem.__slots__:
                if slot[1:] in val_dct: continue
                setattr(self, slot, getattr(exist_item, slot))
            return exist_item
        else:
            self._register(ID)
            return self

    def __init__(self, ID='', functional_unit='', price=0.,
                 source=None, **indicator_CFs):
        CF_dct = self._format_CF_vals(**indicator_CFs)

        returned = self._register_with_consistency_check(
            ID, functional_unit, price, source, CF_dct
            )
        if returned is not self: return

        if source: self.source = source
        else:
            self._functional_unit = auom(functional_unit)
            self._source = None
            self._CFs = {}
            self.price = price
            for indicator, (value, unit) in CF_dct.items():
                self.add_indicator(indicator, value, unit)


    # This makes sure it won't be shown as memory location of the object
    def __repr__(self):
        return f'<ImpactItem: {self.ID}>'


    def show(self):
        '''Show basic information of this impact item.'''
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


    def add_indicator(self, indicator, CF_value, CF_unit=''):
        '''
        Add an indicator with characterization factor values.

        Parameters
        ----------
        indicator : obj or str
            :class:`ImpactIndicator` or its ID.
        CF_value : float
            Characterization factor of the indicator.
        CF_unit : str
            Unit of the characterization factor value.
        '''
        if not hasattr(self, '_source'): self._source = None
        source_item = check_source(self, True)
        if isinstance(indicator, str):
            ind = ImpactIndicator.get_indicator(indicator)
            if ind is None: raise RuntimeError(f'The ImpactIndicator "{indicator}" has not been defined.')
            indicator = ind

        CF_unit2 = CF_unit.replace(' eq', '-eq')

        if CF_unit and CF_unit != indicator.unit and CF_unit2 != indicator.unit:
            try:
                CF_value = auom(parse_unit(CF_unit)[0]). \
                    convert(CF_value, indicator._ureg_unit.units)
            except:
                raise ValueError(f'Conversion of the given unit {CF_unit} to '
                                  f'the default unit {indicator.unit} is not supported.')

        if not hasattr(source_item, '_CFs'): source_item._CFs = {}
        source_item._CFs[indicator.ID] = CF_value

    def remove_indicator(self, indicator):
        '''
        Remove an indicator from this impact item.

        Parameters
        ----------
        indicator : str or :class:`~.ImpactIndicator`
            :class:`~.ImpactIndicator` or its ID.
        '''
        ID = indicator if isinstance(indicator, str) else indicator.ID
        check_source(self, True).CFs.pop(ID)
        print(f'The impact indicator "{ID}" has been removed.')


    def copy(self, new_ID='', set_as_source=False):
        '''
        Return a new :class:`ImpactItem` object with the same settings.

        Parameters
        ----------
        new_ID : str
            ID of the new impact item.
        set_as_source : bool
            Whether to set the original impact item as the source.

        Examples
        --------
        >>> import qsdsan as qs
        >>> GWP = qs.ImpactIndicator('GlobalWarming', alias='GWP', unit='kg CO2-eq')
        >>> Steel = qs.ImpactItem('Steel', 'kg', GWP=2.55)
        >>> Steel.show()
        ImpactItem      : Steel [per kg]
        Price           : None USD
        ImpactIndicators:
                                   Characterization factors
        GlobalWarming (kg CO2-eq)                      2.55
        >>> Steel_cp = Steel.copy(set_as_source=True)
        >>> Steel_cp.show()
        ImpactItem      : item1 [per kg]
        Source          : Steel
        Price           : None USD
        ImpactIndicators:
                                   Characterization factors
        GlobalWarming (kg CO2-eq)                      2.55
        >>> Steel_cp.source is Steel
        True
        '''
        cls = self.__class__
        new = cls.__new__(cls)
        new.__init__(new_ID)

        if set_as_source:
            if self.source:
                new.source = self.source
            else:
                new.source = self
        else:
            new = copy_attr(new, self, skip=('_ID', '_source'))
            new.source = None

        return new

    __copy__ = copy

    def register(self, print_msg=True):
        '''Add this impact item to the registry.'''
        self.registry.register_safely(self.ID, self)
        if print_msg:
            print(f'The impact item "{self.ID}" has been added to the registry.')

    def deregister(self, print_msg=True):
        '''Remove this impact item from the registry.'''
        self.registry.discard(self.ID)
        if print_msg:
            print(f'The impact item "{self.ID}" has been removed from the registry.')

    @classmethod
    def clear_registry(cls, print_msg=True):
        '''Remove all existing impact items from the registry.'''
        cls.registry.clear()
        if print_msg:
            print('All impact items have been removed from the registry.')


    @classmethod
    def get_all_items(cls):
        '''Get all impact items.'''
        return cls.registry.data

    @classmethod
    def get_item(cls, ID):
        '''Get an item by its ID.'''
        return cls.get_all_items().get(ID)

    @classmethod
    def _load_from_df(cls, name, df):
        if name.lower() == 'info':
            if 'kind' not in df.columns: df['kind'] = 'ImpactItem'
            for num in df.index:
                kind = df.iloc[num].kind.lower()
                if kind == 'streamimpactitem':
                    new = StreamImpactItem.__new__(StreamImpactItem)
                else: new = cls.__new__(cls)
                new.__init__(ID=df.iloc[num].ID,
                             functional_unit=df.iloc[num].functional_unit)
        else:
            for num in df.index:
                item = cls.get_item(df.iloc[num].ID)
                item.add_indicator(indicator=name,
                                   CF_value=float(df.iloc[num].expected),
                                   CF_unit=df.iloc[num].unit)

    @classmethod
    def load_items_from_excel(cls, path_or_dict, index_col=None):
        '''Same as :func:`load_from_file`, has been deprecated.'''
        warn('`load_items_from_excel` has been deprecated, '
             'please use `load_from_file` instead.', stacklevel=2)
        cls.load_from_file(path_or_dict, index_col)


    @classmethod
    def load_from_file(cls, path_or_dict, index_col=None):
        '''
        Load impact items from an Excel file or a :class:`pandas.DataFrame`.

        This Excel should have multiple sheets:

            - The "info" sheet should have three columns: "ID" (e.g., Cement) \
             "functional_unit" (e.g., kg), and "kind" ("ImpactItem" or "StreamImpactItem") \
             of different impact items.

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
        path_or_dict : str or dict of :class:`pandas.DataFrame`
            A dict of DataFrame or complete path of the datasheet in xls/xlsx.
        index_col : None or int
            Index column of the :class:`pandas.DataFrame`.

        Tip
        ---
        Refer to the `Bwaise system <https://github.com/QSD-Group/EXPOsan/tree/main/exposan/bwaise/data>`_
        in the `Exposan` repository for a sample file.
        '''
        if isinstance(path_or_dict, str):
            if not (path_or_dict.endswith('.xls') or path_or_dict.endswith('.xlsx')):
                raise ValueError('Only Excel files ends with ".xlsx" or ".xls" can be interpreted.')

            data_file = pd.ExcelFile(path_or_dict, engine='openpyxl')

            for sheet_name in data_file.sheet_names:
                data = data_file.parse(sheet_name, index_col=index_col)
                cls._load_from_df(sheet_name, data)
        else:
            for k, v in path_or_dict.items():
                cls._load_from_df(k, v)

    @property
    def source(self):
        '''
        [:class:`ImpactItem`] If provided, all attributions and properties of this
        impact item will be copied from the provided source.
        ID of the impact item can be provided instead of the object.
        '''
        if hasattr(self, '_source'): return self._source
        else:
            self._source = None
            return self._source
    @source.setter
    def source(self, i):
        if not isinstance(i, ImpactItem):
            if isinstance(i, str): i = self.get_item(i)
            elif not i: i = None
            else:
                raise ValueError('`source` can only be an `ImpactItem` or its ID, '
                                 f'not {type(i).__name__}.')
        self._source = i

    @property
    def ID(self):
        '''[str] ID of this item.'''
        return self._ID
    @ID.setter
    def ID(self, ID):
        self._ID = ID

    @property
    def functional_unit(self):
        '''[str] Functional unit of this item.'''
        source_item = check_source(self, True, False)
        if not hasattr(source_item, '_functional_unit'):
            source_item._functional_unit = auom('kg')
        return source_item._functional_unit.units
    @functional_unit.setter
    def functional_unit(self, i):
        check_source(self, True)._functional_unit = auom(i)

    @property
    def indicators(self):
        ''' [tuple] Impact indicators associated with this item.'''
        return tuple(ImpactIndicator.get_all_indicators(True)[i]
                     for i in self.CFs.keys())

    @property
    def price(self):
        '''Price of this item per functional unit.'''
        source_item = check_source(self, True, False)
        if not hasattr(source_item, '_price'): source_item._price = 0.
        return source_item._price
    @price.setter
    def price(self, price):
        source_item = check_source(self, True)
        source_item._price = float(price)

    @property
    def CFs(self):
        '''[dict] Characterization factors of this item for different impact indicators.'''
        if self.source: return self.source._CFs.copy()
        elif hasattr(self, '_CFs'): return self._CFs
        else:
            self._CFs = None
            return self._CFs
    @CFs.setter
    def CFs(self, indicator, CF_value, CF_unit=''):
        check_source(self, True).add_indicator(indicator, CF_value, CF_unit)

    @property
    def registered(self):
        '''[bool] If this impact item is registered in the record.'''
        data = self.registry.data.get(self.ID)
        return True if data is self else False


# %%

class StreamImpactItem(ImpactItem):
    '''
    A class for calculation of environmental impacts associated with streams
    (e.g., chemical inputs, emissions).

    Parameters
    ----------
    ID : str
        ID of the item. If no ID is provided but its `linked_stream` is provided,
        the ID will be set as ID of the `linked_stream` with a suffix '_item'
    linked_stream : :class:`SanStream`
        The associated :class:`SanStream` for environmental impact calculation.
    functional_unit : str
        Functional unit of the impact item. The default is 'kg'.
    source : :class:`StreamImpactItem`
        If provided, all attributions and properties of this
        :class:`StreamImpactItem` will be copied from the provided source.
    flow_getter : callable|float|int
        If none specified, default to `SanStream.F_mass`.
    indicator_CFs : kwargs
        ImpactIndicators and their characterization factors (CFs).

    Tip
    ---
    For environmental impacts associated with construction and transportation,
    use :class:`ImpactItem` instead.

    Examples
    --------
    Refer to :class:`ImpactItem` for general features.
    Below is about the additional features for :class:`StreamImpactItem`.

    Assume we want to account for the global warming potential for methane:

    >>> # Make impact indicators
    >>> import qsdsan as qs
    >>> GWP = qs.ImpactIndicator('GlobalWarming', alias='GWP', unit='kg CO2-eq')
    >>> FEC = qs.ImpactIndicator('FossilEnergyConsumption', alias='FEC', unit='MJ')
    >>> # Make an stream impact item
    >>> methane_item = qs.StreamImpactItem('methane_item', GWP=28)
    >>> methane_item.show()
    StreamImpactItem: methane_item [per kg]
    Linked to       : None
    Price           : None USD
    ImpactIndicators:
                               Characterization factors
    GlobalWarming (kg CO2-eq)                        28
    >>> # Make a stream and link the stream to the impact item
    >>> cmps = qs.utils.create_example_components()
    >>> qs.set_thermo(cmps)
    >>> methane = qs.SanStream('methane', Methane=1, units='kg/hr',
    ...                        stream_impact_item=methane_item)
    >>> methane_item.show()
    StreamImpactItem: methane_item [per kg]
    Linked to       : methane
    Price           : None USD
    ImpactIndicators:
                               Characterization factors
    GlobalWarming (kg CO2-eq)                        28

    We can make copies of the impact item, and link it to the original one.

    >>> methane2 = methane.copy('methane2')
    >>> methane_item2 = methane_item.copy('methane_item2', stream=methane2,
    ...                                   set_as_source=True)
    >>> methane_item2.CFs['GlobalWarming']
    28
    >>> methane_item.CFs['GlobalWarming'] = 1
    >>> methane_item2.CFs['GlobalWarming']
    1

    We can also add or remove impact indicators.

    >>> methane_item.remove_indicator('GlobalWarming')
    The impact indicator "GlobalWarming" has been removed.
    >>> methane_item2.show()
    StreamImpactItem: methane_item2 [per kg]
    Linked to       : methane2
    Source          : methane_item
    Price           : None USD
    ImpactIndicators:
     None
    >>> methane_item.add_indicator('GlobalWarming', 28)
    >>> methane_item2.show()
    StreamImpactItem: methane_item2 [per kg]
    Linked to       : methane2
    Source          : methane_item
    Price           : None USD
    ImpactIndicators:
                               Characterization factors
    GlobalWarming (kg CO2-eq)                        28
    '''

    __slots__ = ('_ID', '_linked_stream', '_functional_unit', 
                 '_CFs', '_source', '_flow_getter')

    def __init__(self, ID='', linked_stream=None, functional_unit='kg', 
                 source=None, flow_getter=None, **indicator_CFs):
        self._linked_stream = None
        self.linked_stream = linked_stream

        CF_dct = self._format_CF_vals(**indicator_CFs)

        price = linked_stream.price if linked_stream else 0.
        returned = self._register_with_consistency_check(
            ID, 'kg', price, source, CF_dct
            )
        if returned is not self: return
        self._register(ID)

        if not ID and linked_stream:
            ID = self.linked_stream.ID + '_item'
        self._ID = ID

        if source:
            self.source = source
        else:
            self._source = None
            self._functional_unit = auom(functional_unit)
            # self._functional_unit = auom('kg')
            self._CFs = {}
            self.flow_getter = flow_getter
            for indicator, (value, unit) in CF_dct.items():
                self.add_indicator(indicator, value, unit)        


    def __repr__(self):
        return f'<StreamImpactItem: {self.ID}>'


    def show(self):
        '''Show basic information about this :class:`StreamImpactItem` object.'''
        info = f'StreamImpactItem: {self.ID} [per {self.functional_unit}]'
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


    def copy(self, new_ID='', stream=None, set_as_source=False):
        '''
        Return a new :class:`StreamImpactItem` object with the same settings.

        Parameters
        ----------
        new_ID : str
            ID of the new stream impact item.
        stream : :class:`~.SanStream`
            Linked stream to the copy.
        set_as_source : bool
            Whether to set the original impact item as the source.
        '''
        cls = self.__class__
        new = cls.__new__(cls)
        new.ID = new_ID
        new._linked_stream = None # initiate this attribute

        if set_as_source:
            if self.source:
                new.source = self.source
            else:
                new.source = self
        else:
            new = copy_attr(new, self, skip=('_ID', '_linked_stream', '_source'))
            new.source = None

        if stream:
            stream.stream_impact_item = new
            if not new_ID:
                new.ID = f'{stream.ID}_item'

        return new

    __copy__ = copy
    
    @property
    def flow_getter(self):
        '''[callable] A function taking a `SanStream` object and returns the 
        flow [functional unit/hr] for impact calculation. If None specified, default to
        `SanStream.F_mass`.'''
        if not hasattr(self, '_flow_getter'):
            try: return self.source._flow_getter
            except: self.flow_getter = None
        return self._flow_getter
    @flow_getter.setter
    def flow_getter(self, f):
        if f is None:
            self._flow_getter = lambda ws: ws.F_mass
        else:
            if callable(f):
                nargs = f.__code__.co_argcount
                if nargs != 1:
                    raise ValueError('flow_getter must take exactly 1 positional argument,'
                                     'which is expected to be a `SanStream` object')
                self._flow_getter = f
            elif isinstance(f, (float, int)):
                self._flow_getter = lambda ws: f
            else:
                raise TypeError('flow_getter must be a callable, or a number, or None, '
                                f'not {type(f)}')

    @property
    def source(self):
        '''
        [:class:`StreamImpactItem`] If provided, all attributions and properties of this
        :class:`StreamImpactItem` will be copied from the provided source.
        ID of the impact item can be provided instead of the object.

        .. note::

            Since the price is copied from the price of the `linked_stream`, it
            can be different from the source.
        '''
        if hasattr(self, '_source'): return self._source
        else:
            self._source = None
            return self._source
    @source.setter
    def source(self, i):
        if isinstance(i, str): i = self.get_item(i)

        if not isinstance(i, StreamImpactItem):
            if not i: i = None
            else:
                raise ValueError('`source` can only be a StreamImpactItem, ' \
                                 f'not a {type(i).__name__}.')
        self._source = i


    @property
    def linked_stream(self):
        '''
        [:class:`SanStream`] The associated :class:`SanStream` for environmental impact calculation,
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
            self._linked_stream.stream_impact_item = None
            if old_s is not new_s:
                warn(f'`StreamImpactItem` {self.ID} is unlinked from {old_s.ID} and ' \
                     f'linked to {new_s.ID}.')
        if new_s:
            if new_s.stream_impact_item and new_s.stream_impact_item is not self:
                if hasattr(self, '_ID'):
                    warn(f'The original `StreamImpactItem` linked to stream {new_s.ID} '
                         f'is replaced with {self.ID}.')
                else:
                    warn(f'The original `StreamImpactItem` linked to stream {new_s.ID} '
                         f'is replaced upon the creation of a new stream.')
            new_s._stream_impact_item = self
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

    # @property
    # def functional_unit(self):
    #     '''[str] Functional unit of the item, set to 'kg'.'''
    #     return auom('kg')

    @property
    def price(self):
        '''[float] Price of the linked stream.'''
        if self.linked_stream:
            return self.linked_stream.price
        else:
            return 0.
    _price = price