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

from ._units_of_measure import parse_unit
from .utils.loading import load_data, data_path
data_path += '_impact_indicator.tsv'

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
        Category of the ImpactIndicator, e.g., 'human healt'.
    unit : str
        Unit of the ImpactIndicator, e.g., 'kg CO2-eq'.
    description : str
        Supplementary explanation.
    '''
    
    _indicators = {}
    _default_data = None
    
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
            raise DeprecationWarning('`synonym` has been renamed as `alias`.')
            if not alias:
                alias = kwargs['synonym']
            else:
                raise ValueError('`synonym` and `alias` are both provided.')
        
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
        line =   '\n alias(es)   : '
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
        ID : str
            Original ID.
        alias : str
            New alias of the indicator.

        '''
        dct = ImpactIndicator._indicators
        if alias in dct.keys() and dct[alias] is not self:
            raise ValueError(f'The alias "{alias}" already in use.')
        else:
            dct[alias] = self
    
    def get_alias(self):
        '''Return all alias(es) of the indicator as a list.'''
        return tuple(i for i, j in ImpactIndicator._indicators.items()
                     if j==self and i != self.ID)


    @classmethod
    def load_default_indicators(cls):
        '''Load all default indicators as in /data/_impact_indicator.xlsx.'''
        if cls._default_data is not None:
            data = cls._default_data
        else: data = load_data(path=data_path)
        for indicator in data.index:
            if indicator in cls._indicators.keys():
                continue
            else:
                new = cls.__new__(cls)
                new.__init__(ID=indicator,
                             alias=data.loc[indicator]['alias'],
                             unit=data.loc[indicator]['unit'],
                             method=data.loc[indicator]['method'],
                             category=data.loc[indicator]['category'],
                             description=data.loc[indicator]['description'])
                cls._indicators[indicator] = new
        cls._default_data = data


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




# ImpactIndicator.load_default_indicators()



