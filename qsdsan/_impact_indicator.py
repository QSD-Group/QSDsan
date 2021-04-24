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
from thermosteam.utils import registered
from ._units_of_measure import parse_unit
from .utils.loading import load_data

__all__ = ('ImpactIndicator', )


@registered(ticket_name='ind')
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
     Alias      : None
     Method     : TRACI
     Category   : environmental impact
     Description: Effect of climate change ...
    >>> # Add an alias
    >>> GWP.alias = 'GWP'
    >>> GWP.show()
    ImpactIndicator: GlobalWarming as kg CO2-eq
     Alias      : GWP
     Method     : TRACI
     Category   : environmental impact
     Description: Effect of climate change ...
    >>> # Add another impact indicator
    >>> FEC = qs.ImpactIndicator('FossilEnergyConsumption', alias='FEC', unit='MJ')
    >>> # Get all impact indicators
    >>> qs.ImpactIndicator.get_all_indicators()
    {'GlobalWarming': <ImpactIndicator: GlobalWarming>,
     'FossilEnergyConsumption': <ImpactIndicator: FossilEnergyConsumption>}

    Manage the registry.

    >>> GWP.deregister()
    The impact indicator "GlobalWarming" has been removed from the registry.
    >>> qs.ImpactIndicator.get_all_indicators()
    {'FossilEnergyConsumption': <ImpactIndicator: FossilEnergyConsumption>}
    >>> GWP.register()
    The impact indicator "GlobalWarming" has been added to the registry.    
    >>> qs.ImpactIndicator.get_all_indicators()
    {'FossilEnergyConsumption': <ImpactIndicator: FossilEnergyConsumption>,
     'GlobalWarming': <ImpactIndicator: GlobalWarming>}
    >>> qs.ImpactIndicator.clear_registry()
    All impact indicators have been removed from registry.
    >>> qs.ImpactIndicator.get_all_indicators()
    {}
    '''
    
    # _indicators = ImpactIndicator.registry.data
    
    __slots__ = ('_ID', '_alias', '_method', '_category', '_unit', '_ureg_unit',
                 '_unit_remaining', '_description')

    def __init__(self, ID, alias='', method='', category='', unit='', description='',
                 # register=True,
                 **kwargs):
        
        # if register:
        #     if ID in self._indicators.keys():
        #         warn(f'The impact indicator "{ID}" is replaced in the registry.')
        
        # self._ID = ID
        self._register(ID)
        self.alias = alias
        
        self._unit = str(unit)
        self._ureg_unit, self._unit_remaining = parse_unit(unit)
        self._method = method
        self._category = category
        self._description = description
        # self._indicators[ID] = self
        # self._registered = True

        if 'synonym' in kwargs.keys():
            synonym = kwargs['synonym']
            if (not alias or str(alias)=='nan'):
                raise DeprecationWarning('`synonym` has been changed to `alias` for qsdsan v0.2.2 and above.')
                alias = synonym
            else:
                raise DeprecationWarning('`synonym` has been changed to `alias` for qsdsan v0.2.2 and above, ' \
                                         f'the given `synonym` "{synonym}" is ignored as `alias` "{alias}" is provided.')
        

            
            # self.set_alias(alias)
        
        #!!! temp
        

    def __repr__(self):
        return f'<ImpactIndicator: {self.ID}>'

    #!!! Temporary, remove after making sure registry works
    # @classmethod
    # def _indicators(cls):
    #     return {i.ID:i for i in cls.registry}

    def show(self):
        '''Show basic information about this impact indicator.'''
        if self.unit:
            info = f'ImpactIndicator: {self.ID} as {self.unit}'
        else:
            info = f'ImpactIndicator: {self.ID}'
        alias = self.alias if self.alias else 'None'
        line =   f'\n Alias      : {alias}'
        
        # aliases = self.get_alias()
        # if aliases:
        #     for alias in aliases[:-1]:
        #         line += alias + '; '
        #     line += aliases[-1]
        if len(line) > 40: line = line[:40] + '...'
        info += line
        info += f'\n Method     : {self.method or None}'
        info += f'\n Category   : {self.category or None}'
        line =  f'\n Description: {self.description or None}'
        if len(line) > 40: line = line[:40] + '...'
        info += line
        print(info)
    
    _ipython_display_ = show
    
    


    def register(self):
        '''Add this impact indicator to the registry.'''
        self.registry.register_safely(self.ID, self)
        
        
        # ID = self.ID
        # if self._registered:
        #     warn(f'The impact indicator "{ID}" is already in registry.')
        #     return
        # else:
        #     IDs = (ID, *self.get_alias())
        #     for i in IDs:
        #         self._indicators[ID] = self
        
        print(f'The impact indicator "{self.ID}" has been added to the registry.')

    def deregister(self):
        '''Remove this impact indicator from the registry.'''
        self.registry.discard(self.ID)
        
        # ID = self.ID
        # IDs = (ID, *self.get_alias())
        # for i in IDs:
        #     self._indicators.pop(i)
        # self._registered = False
            
        print(f'The impact indicator "{self.ID}" has been removed from the registry.')


    @classmethod
    def clear_registry(cls):
        '''Remove all existing impact indicators from the registry.'''
        cls.registry.clear()
        # for i in cls._indicators.values():
        #     i._registered = False
        # cls._indicators = {}
        print('All impact indicators have been removed from registry.')


    #!!! PAUSED, TO REMOVE _indicators

    @classmethod
    def get_all_indicators(cls, include_alias=False):
        '''
        Get all defined impact indicator as a dict.
        
        Parameters
        ----------
        include_alias : bool
            If True, aliases will be included as keys in the dict as well.
        '''

        if not include_alias:
            return cls.registry.data

        else:
            dct = cls.registry.data.copy()
            dct.update(cls._get_alias_dct(cls))
            return dct

        # if not include_alias:
        #     return cls.registry.data
        
        # else:
        #     dct = 
        #     return cls.registry.data
        
        # else:
        #     # all_indicators = set(cls.registry.data.items())
        #     breakpoint()
        #     return {i.ID:i for i in sorted(set(cls.registry.data.items()),
        #                                    key=lambda j: j.ID)}
        # if as_dict:
        #     cls.registry.data
        #     return {i.ID:i for i in cls.registry}
        
        # return sorted(set([i for i in cls.registry.values()]), key=lambda i: i.ID)

    @classmethod
    def get_indicator(cls, ID_or_alias):
        '''Get an impact indicator by its ID or alias.'''
        dct = cls.get_all_indicators(True)
        return dct.get(ID_or_alias)

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
        in the `Exposan` repository for a sample file.
        '''
        data = load_data(path=path)
        for indicator in data.index:
            # if indicator in cls.get_all_indicators(True).keys():
            #     raise ValueError(f'The impact indicator "{indicator}" has been added.')
            # else:
            new = cls.__new__(cls)
            kwargs = {}
            for k in ('alias', 'unit', 'method', 'category', 'description'):
                try:
                    kwargs[k] = data.loc[indicator][k]
                except KeyError:
                    kwargs[k] = ''

            new.__init__(ID=indicator, **kwargs)

            # cls._indicators[indicator] = new

    def _get_alias_dct(cls):
        dct = {}
        for i in cls.registry.data.values():
            if i.alias:
                dct[i.alias] = i
        return dct

    @property
    def ID(self):
        '''[str] ID of this impact indicator.'''
        return self._ID
    @ID.setter
    def ID(self, ID):
        self._ID = ID

    @property
    def alias(self):
        '''[str] Alias of this impact indicator.'''
        if not hasattr(self, '_alias'): # for initiation
            self._alias = None
        return self._alias
    @alias.setter
    def alias(self, alias):
        alias = None if str(alias) == 'nan' else alias
        alias_dct = self._get_alias_dct()
        
        if alias:
            if not isinstance(alias, str):
                raise TypeError(f'`alias` can only be a str, not {type(alias).__name__}.')

            if alias in alias_dct.keys():
                old_ind = alias_dct[alias]
                if old_ind.ID != self.ID:
                    warn(f'The alias "{alias}" is now being used for "{self.ID}", ' \
                         f'instead of {old_ind.ID}.')
                    old_ind._alias = None
            self._alias = alias
        
        else:
            self._alias = None
    
    
    # def set_alias(self, alias):
    #     '''
    #     Give the impact indicator an alias.

    #     Parameters
    #     ----------
    #     alias : str
    #         Alias of the impact indicator.
    #     '''
    #     if alias and str(alias) != 'nan':
    #         self._alias = alias
    #         self.registry.register_safely(alias, self)
        
        
    #     self.register
        
    #     dct = self.registry.data
    #     if alias in dct.keys() and dct[alias] is not self:
    #         warn(f'The alias "{alias}" is already in use.')
    #         return

    #     else:
    #         dct[alias] = self
    
    # def get_alias(self):
    #     '''Return all aliases of the impact indicator as a list.'''

    #     return sorted([i for i, j in self.registry.data.items()
    #                    if j==self and i != self.ID])

    # def remove_alias(self, alias):
    #     '''
    #     Remove a certain alias.
        
    #     Parameters
    #     ----------
    #     alias : str
    #         The alias of the impact indicator to be removed.
    #     '''
    #     if self.registry.data[alias] is not self:
    #         warn('"{alias}" is not the alias of {self.ID}.')
    #         return
    #     else:
    #         self._indicators.pop(alias)


    @property
    def unit(self):
        '''[str] Unit of this impact indicator.'''    
        return self._unit
    @unit.setter
    def unit(self, i):
        self._unit = str(i)
        self._ureg_unit, self._unit_remaining = parse_unit(i)

    @property
    def method(self):
        '''[str] Impact assessment method of this impact indicator.'''    
        return self._method
    @method.setter
    def method(self, i):
        self._method = i

    @property
    def category(self):
        '''[str] Impact category of this impact indicator.'''    
        return self._category
    @category.setter
    def category(self, i):
        self._category = i

    @property
    def description(self):
        '''[str] Description of this impact indicator.'''    
        return self._description
    @description.setter
    def description(self, i):
        self._description = i

    @property
    def registered(self):
        '''[bool] If this impact indicator is registered in the record.'''
        data = self.registry.data.get(self.ID)
        return True if data else False

