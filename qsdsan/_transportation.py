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
from thermosteam.utils import registered
from . import currency, ImpactItem
from .utils import (
    auom, copy_attr,
    format_number as f_num
    )

__all__ = ('Transportation',)


@registered(ticket_name='Trans')
class Transportation:
    '''
    Transportation activity for cost and environmental impact calculations.

    Parameters
    ----------
    ID : str
        ID of this transportation activity.
    item : :class:`ImpactItem`
        Impact item associated with this transportation activity.
    load_type : str
        Can be either 'mass' or 'volume'.
    load : float
        Quantity of the load per trip.
    load_unit : str
        Unit of the load.
    distance : float
        Distance per trip.
    distance_unit : str
        Unit of the distance.
    interval : float
        Distance per trip.
    interval_unit : str
        Unit of the transportation interval.

    Examples
    --------
    >>> import qsdsan as qs
    >>> # Make impact indicator
    >>> GWP = qs.ImpactIndicator('GlobalWarming', alias='GWP', unit='kg CO2-eq')
    >>> FEC = qs.ImpactIndicator('FossilEnergyConsumption', alias='FEC', unit='MJ')
    >>> # Assuming transporting 1 kg of goods for 1 km emits 10 kg CO2-eq
    >>> Trucking = qs.ImpactItem('Trucking', 'kg*km', GWP=10, FEC=5)
    >>> # Make a transportation activity for transporting 1000 kg goods for 1 mile every day
    >>> shipping = qs.Transportation('shipping', item=Trucking,
    ...                              load_type='mass', load=1, load_unit='tonne',
    ...                              distance='1', distance_unit='mile',
    ...                              interval='1', interval_unit='day')
    >>> shipping.show()
    Transportation: shipping
    Impact item   : Trucking [per trip]
    Load          : 1000 kg
    Distance      : 1.61 km
    Interval      : 24 hr
    Total cost    : None USD
    Total impacts :
                                  Impacts
    GlobalWarming (kg CO2-eq)    1.61e+04
    FossilEnergyConsumption (MJ) 8.05e+03
    '''

    __slots__ = ('_ID', '_item', '_load_type', '_load', '_distance', '_interval',
                 'default_units')

    def __init__(self, ID='', item=None,
                 load_type='mass', load=1., load_unit='kg',
                 distance=1., distance_unit='km',
                 interval=1., interval_unit='hr'):
        if ID == '': # this is only to auto-generate ID for ones that don't have
            self._register(ID)
        else:
            self._ID = ID
        self.item = item
        self.default_units = {
            'distance': 'km',
            'interval': 'hr',
            }

        self._load_type = load_type
        if load_type == 'mass':
            self.default_units['load'] = 'kg'
            self.default_units['quantity'] = 'kg*km'
        elif load_type == 'volume':
            self.default_units['load'] = 'm3'
            self.default_units['quantity'] = 'm3*km'
        else:
            raise ValueError("load_type can only be 'mass' or 'volume', "
                             f'not {load_type}.')

        self._update_value('load', load, load_unit)
        self._update_value('distance', distance, distance_unit)
        self._update_value('interval', interval, interval_unit)

        if item:
            try:
                auom(str(load_unit)+'*'+str(distance_unit)).convert(1, self.item.functional_unit)
            except:
                raise ValueError(f'Units of `load` {load_unit} and `distance` {distance_unit} '
                                 f'do not match the item `functional_unit` {self.item.functional_unit}.')


    def _update_value(self, var, value, unit=''):
        default_unit = self.default_units[var]
        if not unit or unit == default_unit:
            setattr(self, '_'+var, value)
        else:
            converted = auom(unit).convert(float(value), default_unit)
            setattr(self, '_'+var, converted)

    def __repr__(self):
        return f'<Transportation: {self.ID}>'

    def show(self):
        item = self.item
        impacts = self.impacts
        du = self.default_units
        info = f'Transportation: {self.ID}'
        info += f'\nImpact item   : {item.ID} [per trip]'
        info += f"\nLoad          : {f_num(self.load)} {du['load']}"
        info += f"\nDistance      : {f_num(self.distance)} {du['distance']}"
        info += f"\nInterval      : {f_num(self.interval)} {du['interval']}"
        info += f'\nTotal cost    : {f_num(self.cost)} {currency}'
        info += '\nTotal impacts :'
        print(info)
        if len(impacts) == 0:
            print(' None')
        else:
            index = pd.Index((i.ID+' ('+i.unit+')' for i in self.indicators))
            df = pd.DataFrame({
                'Impacts': tuple(self.impacts.values())
                },
                index=index)
            # print(' '*16+df.to_string().replace('\n', '\n'+' '*16))
            print(df.to_string())

    _ipython_display_ = show # funny that _ipython_display_ and _ipython_display behave differently

    def copy(self, new_ID=''):
        new = Transportation.__new__(Transportation)
        new.__init__(new_ID)
        new = copy_attr(new, self, skip=('_ID',))
        return new

    __copy__ = copy


    @property
    def item(self):
        '''[:class:`ImpactItem`] Item associated with this transportation activity.'''
        return self._item
    @item.setter
    def item(self, i):
        if not i:
            i = None
        elif isinstance(i, str):
            i = ImpactItem.get_item(i)
        elif not isinstance(i, ImpactItem):
            raise TypeError('Only `ImpactItem` or the ID of `ImpactItem` can be set, '
                            f'not {type(i).__name__}.')
        self._item = i

    @property
    def load_type(self):
        '''[str] Either "mass" or "volume".'''
        return self._load_type
    @load_type.setter
    def load_type(self, i):
        if i == 'mass':
            self.default_units['load'] = 'kg'
            self.default_units['quantity'] = 'kg*km'
        elif i == 'volume':
            self.default_units['load'] = 'm3'
            self.default_units['quantity'] = 'm3*km'
        else:
            raise ValueError('load_type can only be "mass" or "volume", '
                             f'not {i}.')
        self._load_type = i

    @property
    def indicators(self):
        ''' [tuple] Impact indicators associated with the transportation item.'''
        return self.item.indicators

    @property
    def load(self):
        '''
        [float] Transportation load each trip.

        .. note::

            Set this to 1 and let the load_unit match the functional unit of the item
            if load does not affect price and impacts.

        '''
        return self._load
    @load.setter
    def load(self, load, unit=''):
        self._update_value('load', load, unit)

    @property
    def distance(self):
        '''
        [float] Transportation distance each trip.

        .. note::

            Set this to 1 and let the `distance_unit` match the functional unit of the item
            if distance does not affect price and impacts.

        '''
        return self._distance
    @distance.setter
    def distance(self, distance, unit=''):
        self._update_value('distance', distance, unit)

    @property
    def quantity(self):
        '''[float] Quantity of item functional unit.'''
        quantity = auom(self.default_units['quantity']). \
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
    def price(self):
        '''[float] Unit price of the item.'''
        return self.item.price

    @property
    def cost(self):
        '''[float] Total cost per trip.'''
        return self.price*self.quantity

    @property
    def impacts(self):
        '''[dict] Total impacts of this transportation item.'''
        impacts = {}
        for indicator, CF in self.item.CFs.items():
            impacts[indicator] = self.quantity*CF
        return impacts