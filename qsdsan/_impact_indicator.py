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

from warnings import warn
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

            "synonym" was used bfore v0.2.2 it is still supported, but may be
            deprecated in the future.

    method : str
        Impact assessment method, e.g., 'TRACI'.
    category : str
        Category of the ImpactIndicator, e.g., 'human health'.
    unit : str
        Unit of the ImpactIndicator, e.g., 'kg CO2-eq'.
    description : str
        Supplementary explanation.
    register : bool
        Whether to add to the registry.
        
    Examples
    --------
    Make an impact indicator for global warming potential.
    
    >>> import qsdsan as qs
    >>> GWP = qs.ImpactIndicator('GlobalWarming', method='TRACI',
    ...                          category='environmental impact',
    ...                          unit='kg CO2-eq',
    ...                          description='Effect of climate change measured as \
    ...                                      global warming potential.')
    
    See relevant information.
    
    >>> GWP.show()
    ImpactIndicator: GlobalWarming as kg CO2-eq
     Method     : TRACI
     Category   : environmental impact
     Description: Effect of climate change ...
     
    
    Add, retrieve, and remove alias.
    
    >>> GWP.set_alias('GWP')
    >>> GWP.set_alias('GlobalWarmingPotential')
    >>> GWP.get_alias()
    ['GWP', 'GlobalWarmingPotential']
    >>> GWP.remove_alias('GlobalWarmingPotential')
    >>> GWP.get_alias()
    ['GWP']

    Manage the registry.

    >>> GWP.deregister()
    The impact indicator "GlobalWarming" has been removed from the registry.
    >>> qs.ImpactIndicator.get_all_indicators()
    []
    >>> GWP.register()
    The impact indicator "GlobalWarming" has been added to the registry.    
    >>> qs.ImpactIndicator.get_all_indicators()
    [<ImpactIndicator: GlobalWarming>]
    '''
    
    _indicators = {}
    
    __slots__ = ('_ID', '_alias', '_method', '_category', '_unit', '_ureg_unit',
                 '_unit_remaining', '_description', '_registered')

    def __init__(self, ID, alias='', method='', category='', unit='', description='',
                 register=True, **kwargs):
        
        if register:
            if ID in self._indicators.keys():
                warn(f'The impact indicator "{ID}" is replaced in the registry.')
        
        self._ID = ID
        self._unit = str(unit)
        self._ureg_unit, self._unit_remaining = parse_unit(unit)
        self._method = method
        self._category = category
        self._description = description
        self._indicators[ID] = self
        self._registered = True

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
        '''Show basic information about this impact indicator.'''
        if self.unit:
            info = f'ImpactIndicator: {self.ID} as {self.unit}'
        else:
            info = f'ImpactIndicator: {self.ID}'
        line =   '\n Alias(es)  : '
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
        Give the impact indicator an alias.

        Parameters
        ----------
        alias : str
            New alias of the impact indicator.
        '''
        dct = self._indicators
        if alias in dct.keys() and dct[alias] is not self:
            warn(f'The alias "{alias}" is already in use.')
            return

        else:
            dct[alias] = self
    
    def get_alias(self):
        '''Return all aliases of the impact indicator as a list.'''

        return sorted([i for i, j in self._indicators.items()
                       if j==self and i != self.ID])

    def remove_alias(self, alias):
        '''
        Remove a certain alias.
        
        Parameters
        ----------
        alias : str
            The alias of the impact indicator to be removed.
        '''
        if self._indicators[alias] is not self:
            warn('"{alias}" is not the alias of {self.ID}.')
            return
        else:
            self._indicators.pop(alias)

    def register(self):
        '''Add this impact indicator to the registry.'''
        
        ID = self.ID
        if self._registered:
            warn(f'The impact indicator "{ID}" is already in registry.')
            return
        else:
            IDs = (ID, *self.get_alias())
            for i in IDs:
                self._indicators[ID] = self
        
        print(f'The impact indicator "{ID}" has been added to the registry.')

    def deregister(self):
        '''Remove this impact indicator from the registry.'''
        ID = self.ID
        IDs = (ID, *self.get_alias())
        for i in IDs:
            self._indicators.pop(i)
        self._registered = False
            
        print(f'The impact indicator "{ID}" has been removed from the registry.')


    @classmethod
    def load_indicators_from_file(cls, path):
        '''
        Load impact indicator from a datasheet.
        
        The first row of this datasheet should have "indicator" 
        (must have value as it is used as the ID, e.g., GlobalWarming),
        "alias" (e.g., GWP), "unit" (e.g., kg CO2-eq), "method" (e.g., TRACI),
        "category" (e.g., environmental impact), and "description".
        Aside from "indicator", other information is optional.
        
        Each row should be a data entry.
        
        .. note::
            
            This function is just one way to batch-load impact indicators,
            you can always write your own function that fits your datasheet format,
            as long as it provides all the information to construct the impact indicator.
        
        
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
                raise ValueError(f'The impact indicator "{indicator}" has been added.')
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
        '''Get an impact indicator by its ID.'''
        return cls._indicators[ID]

    @classmethod
    def get_all_indicators(cls, as_dict=False):
        '''
        Get all defined impact indicator.
        
        Parameters
        ----------
        as_dict : bool
            False returns a list and True returns a dict.
        '''
        if as_dict:
            return cls._indicators
        
        return sorted(set([i for i in cls._indicators.values()]), key=lambda i: i.ID)

    @property
    def ID(self):
        '''ID of the impact indicator.''' 
        return self._ID
    @ID.setter
    def ID(self, ID):
        old_ID = self.ID
        self._indicators[ID] = self._indicators.pop(old_ID)
        self._ID = ID

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
        '''Impact assessment method of the impact indicator.'''    
        return self._method
    @method.setter
    def method(self, i):
        self._method = i

    @property
    def category(self):
        '''Impact category of the impact indicator.'''    
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

    @property
    def registered(self):
        '''[bool] If this impact indicator is registered in the records.'''
        return self._registered



