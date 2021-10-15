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

__all__ = ('Construction',)


@registered(ticket_name='Constr')
class Construction:
    '''
    Construction activity for cost and environmental impact calculations.

    Parameters
    ----------
    ID : str
        ID of this construction activity.
    item : :class:`ImpactItem`
        Impact item associated with this consturction activity.
    quantity : float
        Quantity of the impact item involved in this construction activity.
    lifetime : float
        Lifetime of the constructed item.
    lifetime_unit : str
        Unit of the lifetime.

    Examples
    --------
    >>> import qsdsan as qs
    >>> # Make impact indicators
    >>> GWP = qs.ImpactIndicator('GlobalWarming', alias='GWP', unit='kg CO2-eq')
    >>> FEC = qs.ImpactIndicator('FossilEnergyConsumption', alias='FEC', unit='MJ')
    >>> # Make impact item
    >>> Steel = qs.ImpactItem('Steel', 'kg', GWP=2.55, FEC=0.5)
    >>> # Make a construction activity that uses 100 g of steel
    >>> steel_100_g = qs.Construction('steel_100_g', item=Steel, quantity=100,
    ...                               quantity_unit='g')
    >>> steel_100_g.show()
    Construction : steel_100_g
    Impact item  : Steel
    Lifetime     : None yr
    Quantity     : 0.1 kg
    Total cost   : None USD
    Total impacts:
                                  Impacts
    GlobalWarming (kg CO2-eq)       0.255
    FossilEnergyConsumption (MJ)     0.05
    '''

    __slots__ = ('_ID', '_item', '_quantity', '_lifetime')

    def __init__(self, ID='', item=None, quantity=0., quantity_unit='',
                 lifetime=None, lifetime_unit='yr'):
        if ID == '': # this is only to auto-generate ID for ones that don't have
            self._register(ID)
        else:
            self._ID = ID
        self.item = item
        self._update_quantity(quantity, quantity_unit)
        self._lifetime = None
        if lifetime:
            self._lifetime = auom(lifetime_unit).convert(lifetime, 'yr')

    def _update_quantity(self, quantity=0., quantity_unit=''):
        if not quantity_unit or quantity_unit == self.item.functional_unit:
            self._quantity = float(quantity)
        else:
            converted = auom(quantity_unit).convert(float(quantity), self.item.functional_unit)
            self._quantity = converted

    def __repr__(self):
        return f'<Construction: {self.ID}>'

    def show(self):
        '''Show basic information about this :class:`Construction` object.'''
        item = self.item
        impacts = self.impacts
        info = f'Construction : {self.ID}'
        info += f'\nImpact item  : {item.ID}'
        info += f'\nLifetime     : {f_num(self.lifetime)} yr'
        info += f'\nQuantity     : {f_num(self.quantity)} {item.functional_unit}'
        info += f'\nTotal cost   : {f_num(self.cost)} {currency}'
        info += '\nTotal impacts:'
        print(info)
        if len(impacts) == 0:
            print(' None')
        else:
            index = pd.Index((i.ID+' ('+i.unit+')' for i in self.indicators))
            df = pd.DataFrame({
                'Impacts': tuple(self.impacts.values())
                },
                index=index)
            # print(' '*15+df.to_string().replace('\n', '\n'+' '*15))
            print(df.to_string())

    _ipython_display_ = show

    def copy(self, new_ID=''):
        new = Construction.__new__(Construction)
        new.__init__(new_ID)
        new = copy_attr(new, self, skip=('_ID',))
        return new

    __copy__ = copy

    @property
    def lifetime(self):
        '''[float] Lifetime of this construction activity.'''
        return self._lifetime
    @lifetime.setter
    def lifetime(self, lifetime, unit='yr'):
        if lifetime is None:
            self.lifetime = lifetime
        else:
            self._lifetime = auom(unit).convert(lifetime, 'yr')

    @property
    def item(self):
        '''[:class:`ImpactItem`] The impact item associated with this construction activity.'''
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
    def indicators(self):
        ''' [tuple] Impact indicators associated with the construction item.'''
        return self.item.indicators

    @property
    def quantity(self):
        '''[float] Quantity of this construction item.'''
        return self._quantity
    @quantity.setter
    def quantity(self, quantity, unit=''):
        self._update_quantity(quantity, unit)

    @property
    def price(self):
        '''[float] Unit price of the item.'''
        return self.item.price

    @property
    def cost(self):
        '''[float] Total cost of this construction item.'''
        return self.quantity*self.price

    @property
    def impacts(self):
        '''[dict] Total impacts of this construction activity over its lifetime.'''
        impacts = {}
        for indicator, CF in self.item.CFs.items():
            impacts[indicator] = self.quantity*CF
        return impacts