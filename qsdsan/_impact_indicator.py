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

from ._units_of_measure import parse_unit
from .utils.loading import load_data

__all__ = ('ImpactIndicator', )


class ImpactIndicator:
    '''
    To handle different impact indicators in life cycle assessment.
    
    Parameters
    ----------
    ID : str
        ID of the ImpactIndicator.
    alias : str
        Alternative ID of the ImpactIndicator.
        
        .. note::

            "synonym" was used bfore v0.2.2 it is still supported but may be
            deprecated in the future.

    method : str
        Impact assessment method, e.g., 'TRACI'.
    category : str
        Category of the ImpactIndicator, e.g., 'human health'.
    unit : str
        Unit of the ImpactIndicator, e.g., 'kg CO2-eq'.
    description : str
        Supplementary explanation.
    '''
    
    _indicators = {}
    
    __slots__ = ('_ID', '_alias', '_method', '_category', '_unit', '_ureg_unit',
                 '_unit_remaining', '_description')

    def __init__(self, ID, alias='', method='', category='', unit='', description='',
                 **kwargs):
        
        if ID in ImpactIndicator._indicators.keys():
            raise ValueError(f'The ID "{ID}" is currently in use.')
        self._ID = ID
        self._unit = str(unit)
        self._ureg_unit, self._unit_remaining = parse_unit(unit)
        self._method = method
        self._category = category
        self._description = description
        ImpactIndicator._indicators[ID] = self

        if 'synonym' in kwargs.keys():
            synonym = kwargs['synonym']
            if (not alias or str(alias)=='nan'):
                raise DeprecationWarning('`synonym` has been changed to `alias` for qsdsan v0.2.2 and above.')
                alias = synonym
            else:
                raise DeprecationWarning('`synonym` has been changed to `alias` for qsdsan v0.2.2 and above, ' \
                                         f'the given `synonym` "{synonym}" is ignored as `alias` "{alias}" is provided.')
        
        if alias and str(alias) != 'nan':
            self.set_alias(alias)

    def __repr__(self):
        return f'<ImpactIndicator: {self.ID}>'

    def show(self):
        '''Show basic information about this indicator.'''
        if self.unit:
            info = f'ImpactIndicator: {self.ID} as {self.unit}'
        else:
            info = f'ImpactIndicator: {self.ID}'
        line =   '\n Alias(es)   : '
        aliases = self.get_alias()
        if aliases:
            for alias in aliases[:-1]:
                line += alias + '; '
            line += aliases[-1]
            if len(line) > 40: line = line[:40] + '...'
            info += line
        info += f'\n Method     : {self.method or None}'
        info += f'\n Category   : {self.category or None}'
        line =  f'\n Description: {self.description or None}'
        if len(line) > 40: line = line[:40] + '...'
        info += line
        print(info)
    
    _ipython_display_ = show
    
    def set_alias(self, alias):
        '''
        Give the indicator an alias.

        Parameters
        ----------
        alias : str
            New alias of the indicator.

        '''
        dct = ImpactIndicator._indicators
        if alias in dct.keys() and dct[alias] is not self:
            raise ValueError(f'The alias "{alias}" is already in use.')
        else:
            dct[alias] = self
    
    def get_alias(self):
        '''Return all aliases of the indicator as a list.'''
        return tuple(i for i, j in ImpactIndicator._indicators.items()
                     if j==self and i != self.ID)

    def deregister(self):
        '''Remove this indicator from the record.'''
        ID = self.ID
        self._indicators.pop(ID)
        print(f'The indicator "{ID}" has been removed from the record.')


    @classmethod
    def load_indicators_from_file(cls, path):
        '''
        Load indicators from a datasheet.
        
        The first row of this datasheet should have "indicator" 
        (must have value as it is used as the ID, e.g., GlobalWarming),
        "alias" (e.g., GWP), "unit" (e.g., kg CO2-eq), "method" (e.g., TRACI),
        "category" (e.g., environmental impact), and "description".
        Aside from "indicator", other information is optional.
        
        Each row should be a data entry.
        
        .. note::
            
            This function is just one way to batch-load impact indicators,
            you can always write your own function that fits your datasheet format,
            as long as it provides all the information to construct the indicators.
        
        
        Parameters
        ----------
        path : str
            Complete path of the datasheet, currently support tsv, csv, and xls/xlsx.
        
        Tip
        ---
        [1] tsv is preferred as it shows up on GitHub.
        
        [2] Refer to the `Bwaise system <https://github.com/QSD-Group/EXPOsan/tree/main/exposan/bwaise/data>`_
        in the ``Exposan`` repository for a sample file.
        '''
        data = load_data(path=path)
        for indicator in data.index:
            if indicator in cls._indicators.keys():
                raise ValueError(f'The indicator "{indicator}" has been added.')
            else:
                new = cls.__new__(cls)
                dct = {}
                for k in ('alias', 'unit', 'method', 'category', 'description'):
                    try:
                        dct[k] = data.loc[indicator][k]
                    except KeyError:
                        dct[k] = ''

                new.__init__(ID=indicator, **dct)
                cls._indicators[indicator] = new


    @classmethod
    def get_indicator(cls, ID):
        '''Get an indicator by its ID.'''
        return cls._indicators[ID]

    @classmethod
    def get_all_indicators(cls):
        '''Get all defined indicators.'''
        return tuple(i for i in set([i for i in ImpactIndicator._indicators.values()]))

    @property
    def ID(self):
        '''ID of the impact indicator.''' 
        return self._ID

    @property
    def unit(self):
        '''Unit of the impact indicator.'''    
        return self._unit
    @unit.setter
    def unit(self, i):
        self._unit = str(i)
        self._ureg_unit, self._unit_remaining = parse_unit(i)

    @property
    def method(self):
        '''Impact assessment method of the indicator.'''    
        return self._method
    @method.setter
    def method(self, i):
        self._method = i

    @property
    def category(self):
        '''Impact category of the indicator.'''    
        return self._category
    @category.setter
    def category(self, i):
        self._category = i

    @property
    def description(self):
        '''Description of the impact indicator.'''    
        return self._description
    @description.setter
    def description(self, i):
        self._description = i





