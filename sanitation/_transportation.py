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

from . import currency, ImpactIndicator, ImpactItem
from ._units_of_measure import auom
from .utils.formatting import format_number as f_num

indicators = ImpactIndicator._indicators

__all__ = ('Transportation',)


class Transportation:
    '''Transportation cost and environmental impacts.'''

    __slots__ = ('_item', '_load_type', '_load', '_distance', '_interval',
                 '_fee', '_default_units', '_simulated_in')
    
    def __init__(self, item=None,
                 load_type='mass', load=1., load_unit='kg',
                 distance=1., distance_unit='km',
                 interval=1., interval_unit='day',
                 fee=0., fee_unit=currency):
        self.item = item
        self._default_units = {
            'distance': 'km',
            'interval': 'day',
            'fee': currency
            }
        if load_type == 'mass':
            self._default_units['load'] = 'kg'
            self._default_units['quantity'] = 'kg*km'
        elif load_type == 'volume':
            self._default_units['load'] = 'm3'
            self._default_units['quantity'] = 'm3*km'
        else:
            raise ValueError("load_type can only be 'mass' or 'volume', "
                             f'not {load_type}.')
        self._update_value('load', load, load_unit)
        self._update_value('distance', distance, distance_unit)
        self._update_value('interval', interval, interval_unit)
        self._update_value('fee', fee, fee_unit)

        try:
            auom(str(load_unit)+'*'+str(distance_unit)).convert(1, self.item.functional_unit)
        except:
            raise ValueError(f'Units of load {load_unit} and distance {distance_unit} '
                             f'do not match the item functional_unit {self.item.functional_unit}.')
        
    
    def _update_value(self, var, value, unit=''):
        default_unit = self._default_units[var]
        if not unit or unit == default_unit:
            setattr(self, '_'+var, value)
        else:
            converted = auom(unit).convert(float(value), unit)
            setattr(self, '_'+var, converted)

    def __repr__(self):
        return f'<Transportation: {self.item.ID}>'

    def show(self):
        item = self.item
        impacts = self.impacts
        du = self._default_units
        info = f'Transportation: {item.ID}'
        info += f"\n Load         : {f_num(self.load)} {du['load']}"
        info += f"\n Distance     : {f_num(self.distance)} {du['distance']}"
        info += f"\n Interval     : {f_num(self.interval)} {du['interval']}"
        info += f'\n Total cost   : {f_num(self.cost)} {currency}'
        info += '\n Total impacts:'
        if len(impacts) == 0:
            info += ' None'        
        else:
            for indicator in impacts.keys():
                formated = f_num(impacts[indicator])
                unit = indicators[indicator].unit
                info += f'\n     {indicator}: {formated} {unit}'
        print(info)

    _ipython_display_ = show

    @property
    def item(self):
        '''[ImpactItem] Item associated with this transportation.'''
        return self._item
    @item.setter
    def item(self, i):
        if isinstance(i, str):
            i = ImpactItem._items[i]
        elif i is not ImpactItem:
            raise TypeError('Only <ImpactItem> or  <ImpactItem>.ID can be set, '
                            f'not {type(i).__name__}.')
        self._item = i

    @property
    def load_type(self):
        '''[str] Either 'mass' or 'm3'.'''
        return self._load_type
    @load_type.setter
    def load_type(self, i):
        if i == 'mass':
            self._default_units['load'] = 'kg'
            self._default_units['quantity'] = 'kg*km'
        elif i == 'volume':
            self._default_units['load'] = 'm3'
            self._default_units['quantity'] = 'm3*km'
        else:
            raise ValueError("load_type can only be 'mass' or 'volume', "
                             f'not {i}.')
        self._load_type = i

    @property
    def load(self):
        '''
        [float] Transportation load each trip.

        Note
        ----
            Set this to 1 and let the load_unit match the functional unit of the item
            if load does not affect fee and impacts.

        '''
        return self._load
    @load.setter
    def load(self, load, unit=''):
        self._update_value('load', load, unit)

    @property
    def distance(self):
        '''
        [float] Transportation distance each trip.

        Note
        ----
            Set this to 1 and let the distance_unit match the functional unit of the item
            if distance does not affect fee and impacts.

        '''
        return self._distance
    @distance.setter
    def distance(self, distance, unit=''):
        self._update_value('distance', distance, unit)

    @property
    def quantity(self):
        '''[float] Quantity of item functional unit.'''
        quantity = auom(self._default_units['quantity']). \
            convert(self.load*self.distance, self.item.functional_unit)
        return quantity

    @property
    def interval(self):
        '''[float] Time between trips.'''
        return self._interval
    @interval.setter
    def interval(self, interval, unit=''):
        self._update_value('interval', interval, unit)

    @property
    def fee(self):
        '''[float] Fee per functional unit of the item.'''
        return self._fee
    @fee.setter
    def fee(self, fee, unit=''):
        self._update_value('fee', fee, unit)

    @property
    def cost(self):
        '''[float] Total cost per trip.'''
        return self.fee*self.quantity

    @property
    def impacts(self):
        '''[dict] Total impacts of this construction item.'''
        impacts = {}
        for indicator, CF in self.item.CFs.items():
            impacts[indicator.ID] = self.quantity*CF
        return impacts

    # @property
    # def simulated_in(self):
    #     '''[str] ID of the SanUnit associated with this Transportation object.'''
    #     return self._simulated_in
    # @simulated_in.setter
    # def simulated_in(self, i):
    #     self._simulated_in = i






