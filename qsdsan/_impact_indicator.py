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
data_path += '_impact_indicator.csv'

__all__ = ('ImpactIndicator', )


class ImpactIndicator:
    '''
    To handle different impact indicators in life cycle assessment.
    
    Parameters
    ----------
    ID : str
        ID of the ImpactIndicator.
    synonym : str
        Alternative ID of the ImpactIndicator.
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
    
    __slots__ = ('_ID', '_synonym', '_method', '_category', '_unit', '_ureg_unit',
                 '_unit_remaining', '_description')

    def __init__(self, ID, synonym='', method='', category='', unit='', description=''):
        
        if ID in ImpactIndicator._indicators.keys():
            raise ValueError(f'The ID "{ID}" is currently in use.')
        self._ID = ID
        self._unit = str(unit)
        self._ureg_unit, self._unit_remaining = parse_unit(unit)
        self._method = method
        self._category = category
        self._description = description
        ImpactIndicator._indicators[ID] = self
        if synonym and str(synonym) != 'nan':
            self.set_synonym(synonym)

    def __repr__(self):
        return f'<ImpactIndicator: {self.ID}>'

    def show(self):
        '''Show basic information about this indicator.'''
        if self.unit:
            info = f'ImpactIndicator: {self.ID} as {self.unit}'
        else:
            info = f'ImpactIndicator: {self.ID}'
        line =   '\n Synonyms   : '
        synonyms = self.get_synonym()
        if synonyms:
            for synonym in synonyms[:-1]:
                line += synonym + '; '
            line += synonyms[-1]
            if len(line) > 40: line = line[:40] + '...'
            info += line
        info += f'\n Method     : {self.method or None}'
        info += f'\n Category   : {self.category or None}'
        line =  f'\n Description: {self.description or None}'
        if len(line) > 40: line = line[:40] + '...'
        info += line
        print(info)
    
    _ipython_display_ = show
    
    def set_synonym(self, synonym):
        '''
        Give the indicator a synonym.

        Parameters
        ----------
        ID : str
            Original ID.
        synonym : str
            New synonym of the indicator.

        '''
        dct = ImpactIndicator._indicators
        if synonym in dct.keys() and dct[synonym] is not self:
            raise ValueError(f'The synonym "{synonym}" already in use.')
        else:
            dct[synonym] = self
    
    def get_synonym(self):
        '''Return all synonyms of the indicator as a list.'''
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
                             synonym=data.loc[indicator]['synonym'],
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




ImpactIndicator.load_default_indicators()



