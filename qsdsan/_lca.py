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
import numpy as np
import pandas as pd
from math import ceil
from collections.abc import Iterable
from warnings import warn
from . import ImpactIndicator, ImpactItem, Stream, SanStream, SanUnit
from .utils import (
    auom,
    clear_lca_registries,
    format_number as f_num
    )

__all__ = ('LCA',)


class LCA:
    '''
    For life cycle assessment (LCA) of a System.

    Parameters
    ----------
    system : :class:`biosteam.System`
        System for which this LCA is conducted for.
    lifetime : int
        Lifetime of the LCA in years.
    indicators : Iterable(obj)
        `ImpactIndicator` objects or their IDs/aliases.
    uptime_ratio : float
        Fraction of time that the system is operating.
    annualize_construction : bool
        Used in the case that the lifetime of this LCA (e.g., 10 years)
        is not divisible by the lifetime of certain equipment (e.g., 8 years).
        If True, then the impacts from construction will be annualized using
        the lifetime of the equipment;
        if False, then the total number of the equipment needed throughout this
        LCA will be calculated using `ceil(LCA lifetime/equipment lifetime)`.
    simulate_system : bool
        Whether to simulate the system before creating the LCA object.
    simulate_kwargs : dict
        Keyword arguments for system simulation (used when `simulate_system` is True).
    item_quantities : kwargs, :class:`ImpactItem` or str = float/callable or (float/callable, unit)
        Other :class:`ImpactItem` objects (e.g., electricity) and their quantities.
        Note that callable functions are used so that quantity of items can be updated.


    Examples
    --------
    A system should be constructed prior to LCA, here we import a pre-constructed one.

    >>> import qsdsan as qs
    >>> from qsdsan.utils import create_example_system
    >>> sys = create_example_system()
    >>> sys.diagram() # doctest: +SKIP
    >>> sys.simulate()
    >>> sys.show()
    System: sys
    ins...
    [0] salt_water
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): H2O   111
                        NaCl  0.856
    [1] methanol
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Methanol  0.624
    [2] ethanol
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Ethanol  0.217
    outs...
    [0] alcohols
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Methanol  0.624
                        Ethanol   0.217
    [1] waste_brine
        phase: 'l', T: 350 K, P: 101325 Pa
        flow (kmol/hr): H2O   88.8
                        NaCl  0.684

    And we also need to specify the impact indicators that we are interested in.

    >>> GWP = qs.ImpactIndicator('GlobalWarming', alias='GWP', unit='kg CO2-eq')
    >>> FEC = qs.ImpactIndicator('FossilEnergyConsumption', alias='FEC', unit='MJ')

    There are four different types of impacts in `QSDsan`:
    construction, transportation, stream, and others.

    Note that it is best to add the impact items when developing the unit module,
    (i.e., typically in the `_design` function, but can also in `_run`  or `_cost`)
    but for illustrative purpose, we add it after the system is constructed.

    Construction is mainly used for impacts that only occur once per lifetime
    of the equipment or the unit.

    For example, assume we want to consider the amount of stainless steel
    and concrete used in constructing the MixTank M1.

    >>> # Make the impact item, numbers are made up
    >>> SS = qs.ImpactItem('SS', functional_unit='kg', GWP=3, FEC=50)
    >>> Concrete = qs.ImpactItem('Concrete', functional_unit='kg', GWP=4, FEC=30)
    >>> # Specify the amount of stainless steel and concrete used in the unit
    >>> SS_constr_M1 = qs.Construction(item=SS, quantity=100)
    >>> Concrete_constr_M1 = qs.Construction(item=Concrete, quantity=50)
    >>> # Retrieve the unit from the registry
    >>> flowsheet = qs.Flowsheet.flowsheet.default
    >>> M1 = flowsheet.unit.M1
    >>> # Add the construction activity
    >>> M1.construction = (SS_constr_M1, Concrete_constr_M1)

    Transportation activity can be added in a similar manner, assuming that
    stainless steel and concrete are delivered by truck from 500 km away.

    The interval set below is assuming a system lifetime of 10 year
    and this delivery is only needed once for the entire lifetime.

    >>> lifetime = 10
    >>> Trucking = qs.ImpactItem('Trucking', functional_unit='kg*km',
    ...                          GWP=0.5, FEC=1.5)
    >>> total_quantity = SS_constr_M1.quantity + Concrete_constr_M1.quantity
    >>> Trans_M1 = qs.Transportation(item=Trucking, load_type='mass',
    ...                              load=total_quantity, distance=500,
    ...                              interval=lifetime, interval_unit='yr')
    >>> M1.transportation = Trans_M1

    We can als consider the impacts associated with chemicals and emissions.
    For example, assume the acquisition of methanol, ethanol and disposal of
    the waste brine all have impacts, but the generated alcohols can be treated
    as a product therefore have credits with

    >>> # Retrieve streams
    >>> methanol = flowsheet.stream.methanol
    >>> ethanol = flowsheet.stream.ethanol
    >>> alcohols = flowsheet.stream.alcohols
    >>> waste_brine = flowsheet.stream.waste_brine
    >>> # Create `StreamImpactItem` and link to the streams
    >>> methanol_item = qs.StreamImpactItem(linked_stream=methanol, GWP=2, FEC=13)
    >>> ethanol_item = qs.StreamImpactItem(linked_stream=ethanol, GWP=2.1, FEC=25)
    >>> alcohols_item = qs.StreamImpactItem(linked_stream=alcohols, GWP=-0.2, FEC=-5)
    >>> brine_item = qs.StreamImpactItem(linked_stream=waste_brine, GWP=2, FEC=3)

    Finally, there might be other impacts we want to include in the LCA,
    for example, the electricity needed to operate the system.

    We can use add those additional items when creating the `LCA` object.

    >>> # Get the electricity usage of the system throughout the lifetime,
    >>> # note that the default power utility unit is hr
    >>> total_power = sys.power_utility.rate*24*365*lifetime
    >>> # Create an impact item for the electricity
    >>> e_item = qs.ImpactItem('e_item', 'kWh', GWP=1.1, FEC=24)
    >>> # Create the LCA object
    >>> lca = qs.LCA(system=sys, lifetime=10, e_item=total_power)

    Now we can look at the total impacts associate with this system.

    >>> lca.show() # doctest: +ELLIPSIS
    LCA: sys (lifetime 10 yr)
    ...
    >>> # Retrieve impacts associated with a specific indicator
    >>> lca.get_total_impacts()[GWP.ID] # doctest: +ELLIPSIS
    349737809...
    >>> # Or breakdowns of the different category
    >>> lca.get_impact_table('Construction') # doctest: +SKIP
    >>> # Below is for testing purpose, you do not need it
    >>> lca.get_impact_table('Construction').to_dict() # doctest: +ELLIPSIS
    {'Quantity': ...
    >>> lca.get_impact_table('Transportation').to_dict() # doctest: +ELLIPSIS
    {'Quantity': ...
    >>> lca.get_impact_table('Stream').to_dict() # doctest: +ELLIPSIS
    {'Mass [kg]': ...
    >>> lca.get_impact_table('Construction').to_dict() # doctest: +ELLIPSIS
    {'Quantity': ...

    You can also allocate the impact based on mass, energy, value, or a ratio you like

    >>> lca.get_allocated_impacts(sys.products, allocate_by='mass')['waste_brine']['FossilEnergyConsumption'] # doctest: +ELLIPSIS
    46018554...
    >>> lca.get_allocated_impacts(sys.products, allocate_by='energy')['alcohols']['GlobalWarming'] # doctest: +ELLIPSIS
    11063012...
    >>> alcohols.price = 5
    >>> waste_brine.price = 1
    >>> GWP_alcohols = lca.get_allocated_impacts(sys.products, allocate_by='value')['alcohols']['GlobalWarming']
    >>> GWP_brine = lca.get_allocated_impacts(sys.products, allocate_by='value')['waste_brine']['GlobalWarming']
    >>> GWP_alcohols + GWP_brine # doctest: +ELLIPSIS
    5469809...
    >>> lca.get_total_impacts(exclude_streams=sys.products)['GlobalWarming'] # doctest: +ELLIPSIS
    5469809...
    >>> # Clear all registries for testing purpose
    >>> from qsdsan.utils import clear_lca_registries
    >>> clear_lca_registries()

    See Also
    --------
    `SanUnit and System <https://qsdsan.readthedocs.io/en/latest/tutorials/SanUnit_and_System.html>`_

    `TEA and LCA <https://qsdsan.readthedocs.io/en/latest/tutorials/TEA_and_LCA.html>`_
    '''

    __slots__ = ('_system',  '_lifetime', '_uptime_ratio',
                 '_construction_units', '_transportation_units',
                 '_lca_streams', '_indicators',
                 '_other_items', '_other_items_f', 'annualize_construction')


    def __init__(self, system, lifetime,
                 indicators=(), uptime_ratio=1, annualize_construction=False,
                 simulate_system=True, simulate_kwargs={},
                 **item_quantities):
        if simulate_system: system.simulate(**simulate_kwargs)
        self._construction_units = set()
        self._transportation_units = set()
        self._lca_streams = set()
        self._update_system(system)
        self.lifetime = lifetime
        self.indicators = indicators
        self.uptime_ratio = uptime_ratio
        self.annualize_construction = annualize_construction
        self._other_items = {}
        self._other_items_f = {}
        for item, val in item_quantities.items():
            if item == 'lifetime_unit': continue # legacy codes
            try:
                f_quantity, unit = val # unit provided for the quantity
            except Exception as e:
                if 'unpack' in str(sys.exc_info()[1]):
                    f_quantity = val
                    unit = ''
                else:
                    raise e
            self.add_other_item(item, f_quantity, unit)

    clear_lca_registries = clear_lca_registries

    def _update_system(self, system):
        for u in system.units:
            if not isinstance (u, SanUnit):
                continue
            if getattr(u, 'construction', []):
                self._construction_units.add(u)
            if getattr(u, 'transportation', []):
                self._transportation_units.add(u)
        self._construction_units = sorted(self._construction_units,
                                          key=lambda u: u.ID)
        self._transportation_units = sorted(self._transportation_units,
                                            key=lambda u: u.ID)
        for s in (i for i in system.feeds+system.products):
            if not hasattr(s, 'stream_impact_item'):
                continue
            if s.stream_impact_item:
                self._lca_streams.add(s)
        self._lca_streams = sorted(self._lca_streams, key=lambda s: s.ID)
        self._system = system
        try: # for older versions of biosteam without the `_LCA` attribute
            system._LCA = self
        except AttributeError:
            pass
        

    def add_other_item(self, item, f_quantity, unit=''):
        '''Add other :class:`ImpactItem` in LCA.'''
        if isinstance(item, str):
            item = ImpactItem.get_item(item)
            if item is None:
                raise ValueError(f'No ImpactItem with the ID {item}.')
        fu = item.functional_unit
        if not callable(f_quantity):
            f = lambda: f_quantity
        else:
            nargs = f_quantity.__code__.co_argcount
            if nargs == 0: f = f_quantity
            else: f = lambda: f_quantity(self)
            # f = f_quantity
        quantity = f()
        if unit and unit != fu:
            try:
                quantity = auom(unit).convert(quantity, fu)
            except:
                raise ValueError(f'Conversion of the given unit {unit} to '
                                 f'item functional unit {fu} is not supported.')
        self._other_items_f[item.ID] = {'item':item, 'f_quantity':f, 'unit':unit}
        self.other_items[item.ID] = {'item':item, 'quantity':quantity}


    def refresh_other_items(self):
        '''Refresh quantities of other items using the given functions.'''
        for item_ID, record in self._other_items_f.items():
            item, f_quantity, unit = record.values()
            self.other_items[item_ID]['quantity'] = f_quantity()


    def __repr__(self):
        return f'<LCA: {self.system}>'

    def show(self):
        '''Show basic information of this :class:`LCA` object.'''
        lifetime = self.lifetime
        info = f'LCA: {self.system} (lifetime {f_num(lifetime)} yr)'
        info += '\nImpacts:'
        print(info)
        if len(self.indicators) == 0:
            print(' None')
        else:
            index = pd.Index((i.ID+' ('+i.unit+')' for i in self.indicators))
            df = pd.DataFrame({
                'Construction': tuple(self.total_construction_impacts.values()),
                'Transportation': tuple(self.total_transportation_impacts.values()),
                'Stream': tuple(self.total_stream_impacts.values()),
                'Others': tuple(self.total_other_impacts.values()),
                'Total': tuple(self.total_impacts.values())
                },
                index=index)
            # print(' '*9+df.to_string().replace('\n', '\n'+' '*9))
            print(df.to_string())

    _ipython_display_ = show


    def get_construction_impacts(self, units=None, annual=False):
        '''
        Return all construction-related impacts for the given units.

        Parameters
        ----------
        units : Iterable(obj)
            Unit operations considered for impacts
            (will default to all unit operations in the system).
        annual : bool
            If True, will return the annual impacts considering `uptime_ratio`
            instead of across the system lifetime.
        '''
        units = self.construction_units if units is None else units
        annualize = self.annualize_construction
        if not isinstance(units, Iterable) or isinstance(units, str):
            units = (units,)
        time = self.lifetime
        impacts = dict.fromkeys((i.ID for i in self.indicators), 0.)
        for i in units:
            if not isinstance(i, SanUnit):
                continue
            for j in i.construction:
                impact = j.impacts
                if j.lifetime is not None: # this equipment has a lifetime
                    constr_lifetime = j.lifetime
                    ratio = ceil(time/constr_lifetime) if not annualize else time/constr_lifetime
                else: # equipment doesn't have a lifetime
                    if i.lifetime and not isinstance(i.lifetime, dict): # unit has a uniform lifetime
                        constr_lifetime = i.lifetime
                        ratio = ceil(time/constr_lifetime) if not annualize else time/constr_lifetime
                    else: # no lifetime, assume just need one
                        ratio = 1.
                for m, n in impact.items():
                    if m not in impacts.keys():
                        continue
                    impacts[m] += n*ratio
        if annual == True:
            lifetime = self.lifetime
            for i, j in impacts.items(): impacts[i] = j/lifetime
        return impacts

    def get_transportation_impacts(self, units=None, annual=False):
        '''
        Return all transportation-related impacts for the given unit.
        
        Parameters
        ----------
        units : Iterable(obj)
            Unit operations considered for impacts
            (will default to all unit operations in the system).
        annual : bool
            If True, will return the annual impacts considering `uptime_ratio`
            instead of across the system lifetime.
        '''
        units = self.transportation_units if units is None else units
        if not isinstance(units, Iterable):
            units = (units,)
        time = self.lifetime_hr
        impacts = dict.fromkeys((i.ID for i in self.indicators), 0.)
        for i in units:
            if not isinstance(i, SanUnit):
                continue
            for j in i.transportation:
                impact = j.impacts
                for m, n in impact.items():
                    if m not in impacts.keys():
                        continue
                    impacts[m] += n*time/j.interval
        if annual == True:
            lifetime = self.lifetime
            for i, j in impacts.items(): impacts[i] = j/lifetime
        return impacts


    def get_stream_impacts(self, stream_items=None, exclude_streams=None, kind='all', annual=False, **kwargs):
        '''
        Return all stream-related impacts for the given streams.
        
        Parameters
        ----------
        stream_items : Iterable(obj)
            Streams considered for impacts
            (will default to all streams in the system with `StreamImpactItem`).
        annual : bool
            If True, will return the annual impacts considering `uptime_ratio`
            instead of across the system lifetime.
        '''
        if 'exclude' in kwargs.keys() and exclude_streams is None: exclude_streams=kwargs['exclude']
        isa = isinstance
        if stream_items == None:
            stream_items = self.stream_inventory
        if not isa(stream_items, Iterable):
            stream_items = (stream_items,)
        if not isa(exclude_streams, Iterable):
            exclude_streams = (exclude_streams,)
        impacts = dict.fromkeys((i.ID for i in self.indicators), 0.)
        time = self.lifetime_hr if annual is False else 365*24*self.uptime_ratio
        for j in stream_items:
            # In case that ws instead of the item is given
            if isa(j, Stream):
                if not isa(j, SanStream):
                    continue
                ws = j
                if j.stream_impact_item:
                    j = ws.stream_impact_item
                else: continue
            else:
                ws = j.linked_stream

            if ws in exclude_streams: continue
            
            F_mass = j.flow_getter(ws)
            for m, n in j.CFs.items():
                if kind in ('all', 'total', 'net'):
                    pass
                elif kind in ('direct', 'direct_emission'):
                    n = max(n, 0)
                elif kind == 'offset':
                    n = min(n, 0)
                else:
                    raise ValueError('kind can only be "all", "direct_emission", or "offset", '
                                     f'not "{kind}".')
                if m not in impacts.keys():
                    continue
                impacts[m] += n*time*F_mass
        return impacts

    def get_other_impacts(self, annual=False):
        '''
        Return all additional impacts from "other" :class:`ImpactItems` objects
        based on defined quantity.
        
        Parameters
        ----------
        annual : bool
            If True, will return the annual impacts considering `uptime_ratio`
            instead of across the system lifetime.
        '''
        self.refresh_other_items()
        impacts = dict.fromkeys((i.ID for i in self.indicators), 0.)
        other_dct = self.other_items
        for i in other_dct.keys():
            item = ImpactItem.get_item(i)
            for m, n in item.CFs.items():
                if m not in impacts.keys():
                    continue
                impacts[m] += n*other_dct[i]['quantity']
        if annual == True:
            lifetime = self.lifetime
            for i, j in impacts.items(): impacts[i] = j/lifetime
        return impacts

    def get_total_impacts(
            self,
            operation_only=False,
            exclude_streams=None,
            annual=False,
            **kwargs
            ):
        '''
        Return total impacts, normalized to a certain time frame.
        
        Parameters
        ----------
        operation_only : bool
            If True, then no construction impacts will be included.
        exclude_streams : Iterable(obj)
            Streams to be excluded from the LCA.
        annual : bool
            If True, will return the annual impacts considering `uptime_ratio`
            instead of across the system lifetime.
        '''
        if 'exclude' in kwargs.keys() and exclude_streams is None: exclude_streams=kwargs['exclude']
        impacts = dict.fromkeys((i.ID for i in self.indicators), 0.) 
        trans = self.get_transportation_impacts(self.transportation_units, annual=annual)
        ws_impacts = self.get_stream_impacts(stream_items=self.stream_inventory,
                                             exclude_streams=exclude_streams,
                                             annual=annual)
        other = self.get_other_impacts(annual=annual)
        if operation_only == False:
            constr = self.get_construction_impacts(self.construction_units, annual=annual)
            categories = (constr, trans, ws_impacts, other)
        else: categories = (trans, ws_impacts, other)
        
        for i in categories:
            for m, n in i.items():
                if m not in impacts.keys():
                    continue
                impacts[m] += n
                
        return impacts

    def get_allocated_impacts(self, streams=(), allocate_by='mass',
                              operation_only=False, annual=False):
        '''
        Allocate total impacts to one or multiple streams.

        Note that original impacts assigned to the streams will be excluded,
        i.e., the total impact for allocation will be calculated using
        `LCA.get_total_impacts(exclude_streams=streams)`.

        Parameters
        ----------
        streams : Iterable(obj)
            One or a Iterable of streams. Note that impacts of these streams will be
            excluded in calculating the total impacts.
        allocate_by : str, Iterable, or function to generate an Iterable
            If provided as a str, can be "mass" (`F_mass`), "energy" (`HHV`),
            or 'value' (`F_mass`*`price`) to allocate the impacts accordingly.
            If provided as an Iterable (no need to normalize so that sum of the Iterable is 1),
            will allocate impacts according to the Iterable.
            If provided as a function,  will call the function to generate an
            Iterable to allocate the impacts accordingly.
        operation_only : bool
            If True, then no construction impacts will be included.
        annual : bool
            If True, will return the annual impacts considering `uptime_ratio`
            instead of across the system lifetime.
            

        .. note::

            Energy of the stream will be calculated as the sum of HHVs of all components
            in the stream.

        '''
        if not isinstance(streams, Iterable):
            streams = (streams,)
        impact_dct = self.get_total_impacts(operation_only=operation_only, exclude_streams=streams, annual=annual)
        impact_vals = np.array([i for i in impact_dct.values()])
        allocated = {}
        if len(streams) == 1:
            if not isinstance(streams[0], SanStream):
                return None
            return impact_dct
        if allocate_by == 'mass':
            ratios = np.array([i.F_mass for i in streams])
        elif allocate_by == 'energy':
            ratios = np.array([i.HHV for i in streams])
        elif allocate_by == 'value':
            ratios = np.array([i.F_mass*i.price for i in streams])
        elif iter(allocate_by):
            ratios = allocate_by
        elif callable(allocate_by):
            ratios = allocate_by()
        else:
            raise ValueError('allocate_by can only be "mass", "energy", "value", '
                             'an Iterable (with the same length as `streams`), '
                             'or a function to generate an Iterable.')
        if ratios.sum() == 0:
            raise ValueError('Calculated allocation ratios are all zero, cannot allocate.')
        ratios = ratios/ratios.sum()
        for n, s in enumerate(streams):
            if not isinstance(s, SanStream):
                continue
            if not s in self.system.streams:
                raise ValueError(f'`WasteStream` {s} not in the system.')
            allocated[s.ID] = dict(zip(impact_dct.keys(), ratios[n]*impact_vals))
        return allocated


    def get_unit_impacts(
            self, units,
            exclude_streams=None,
            operation_only=False,
            annual=False,
            **kwargs
            ):
        '''
        Return total impacts with certain units.
        
        Parameters
        ----------
        units : Iterable(obj)
            Unit operations to be included in the calculation.
        operation_only : bool
            If True, then no construction impacts will be included.
        exclude_streams : Iterable(obj)
            Streams to be excluded from the LCA.
        annual : bool
            If True, will return the annual impacts considering `uptime_ratio`
            instead of across the system lifetime.
        '''
        if 'exclude' in kwargs.keys() and exclude_streams is None: exclude_streams=kwargs['exclude']
        if not isinstance(units, Iterable):
            units = (units,)
        
        trans = self.get_transportation_impacts(units, annual=annual)
        stream_items = set(i for i in
                       sum((tuple(unit.ins+unit.outs) for unit in units), ())
                       if i.stream_impact_item)

        s = self.get_stream_impacts(stream_items=stream_items,
                                    exclude_streams=exclude_streams,
                                    annual=annual)
        other = self.get_other_impacts(annual=annual)
        tot = s.copy()
        for m in tot.keys():
            tot[m] += trans[m] + s[m] + other[m]
        if not operation_only:
            constr = self.get_construction_impacts(units, annual=annual)
            for m in tot.keys(): tot[m] += constr[m]          
            
        return tot

    def _append_cat_sum(self, cat_table, cat, tot):
        num = len(cat_table)
        cat_table.loc[num] = '' # initiate a blank spot for value to be added later

        for i in self.indicators:
            cat_table[f'{i.ID} [{i.unit}]'][num] = tot[i.ID]
            cat_table[f'Category {i.ID} Ratio'][num] = 1

        if cat in ('construction', 'transportation'):
            cat_table.rename(index={num: ('Sum', 'All')}, inplace=True)
            cat_table.index = \
                pd.MultiIndex.from_tuples(cat_table.index,
                                          names=[cat.capitalize(), 'SanUnit'])
        else:
            cat_table.rename(index={num: 'Sum'}, inplace=True)

        return cat_table

    def get_impact_table(self, category, annual=False):
        '''
        Return a :class:`pandas.DataFrame` table for the given impact category.
        
        Parameters
        ----------
        category : str
            Can be 'construction', 'transportation', 'stream', or 'other'.
        annual : bool
            If True, will return the annual impacts considering `uptime_ratio`
            instead of across the system lifetime.
        '''
        time = self.lifetime_hr
        cat = category.lower()
        tot_f = getattr(self, f'get_{cat}_impacts')
        kwargs = {'annual': annual} if cat != 'other' else {}
        tot = tot_f(**kwargs)

        if cat in ('construction', 'transportation'):
            units = sorted(getattr(self, f'_{cat}_units'),
                              key=(lambda su: su.ID))
            items = sorted(set(i.item for i in getattr(self,  f'{cat}_inventory')),
                           key=(lambda item: item.ID))
            if len(items) == 0:
                return f'No {cat}-related impacts.'

            # Note that item_dct = dict.fromkeys([item.ID for item in items], []) won't work
            item_dct = dict.fromkeys([item.ID for item in items])
            for item_ID in item_dct.keys():
                item_dct[item_ID] = dict(SanUnit=[], Quantity=[])
            for su in units:
                if not isinstance(su, SanUnit):
                    continue
                for i in getattr(su, cat):
                    item_dct[i.item.ID]['SanUnit'].append(su.ID)
                    if cat == 'transportation':
                        item_dct[i.item.ID]['Quantity'].append(i.quantity*time/i.interval)
                    else: # construction
                        lifetime = i.lifetime or su.lifetime or self.lifetime
                        if isinstance(lifetime, dict): # in the case the the equipment is not in the unit lifetime dict
                            lifetime = lifetime.get(i.item.ID) or self.lifetime
                        constr_ratio = self.lifetime/lifetime if self.annualize_construction else ceil(self.lifetime/lifetime)
                        item_dct[i.item.ID]['Quantity'].append(i.quantity*constr_ratio)

            dfs = []
            for item in items:
                dct = item_dct[item.ID]
                dct['SanUnit'].append('Total')
                dct['Quantity'] = np.append(dct['Quantity'], sum(dct['Quantity']))
                if dct['Quantity'].sum() == 0.: dct['Item Ratio'] = 0
                else: dct['Item Ratio'] = dct['Quantity']/dct['Quantity'].sum()*2
                for i in self.indicators:
                    if i.ID in item.CFs:
                        dct[f'{i.ID} [{i.unit}]'] = impact = dct['Quantity']*item.CFs[i.ID]
                        dct[f'Category {i.ID} Ratio'] = impact/tot[i.ID]
                    else:
                        dct[f'{i.ID} [{i.unit}]'] = dct[f'Category {i.ID} Ratio'] = 0
                df = pd.DataFrame.from_dict(dct)
                index0 = f'{item.ID} [{item.functional_unit}]'
                df.set_index([pd.MultiIndex.from_arrays(
                    [(index0,)*len(dct['SanUnit'])], names=(category,)),
                    'SanUnit'],
                    inplace=True)
                dfs.append(df)

            table = pd.concat(dfs)
            return self._append_cat_sum(table, cat, tot)

        ind_head = sum(([f'{i.ID} [{i.unit}]',
                         f'Category {i.ID} Ratio'] for i in self.indicators), [])

        if cat in ('stream', 'streams'):
            headings = ['Stream', 'Mass [kg]', *ind_head]
            item_dct = dict.fromkeys(headings)
            for key in item_dct.keys():
                item_dct[key] = []
            for ws_item in self.stream_inventory:
                ws = ws_item.linked_stream
                item_dct['Stream'].append(ws.ID)
                mass = ws_item.flow_getter(ws) * time
                item_dct['Mass [kg]'].append(mass)
                for ind in self.indicators:
                    if ind.ID in ws_item.CFs.keys():
                        impact = ws_item.CFs[ind.ID]*mass
                        item_dct[f'{ind.ID} [{ind.unit}]'].append(impact)
                        item_dct[f'Category {ind.ID} Ratio'].append(impact/tot[ind.ID])
                    else:
                        item_dct[f'{ind.ID} [{ind.unit}]'].append(0)
                        item_dct[f'Category {ind.ID} Ratio'].append(0)
            table = pd.DataFrame.from_dict(item_dct)
            table.set_index(['Stream'], inplace=True)
            return self._append_cat_sum(table, cat, tot)

        elif cat == 'other':
            headings = ['Other', 'Quantity', *ind_head]
            item_dct = dict.fromkeys(headings)
            for key in item_dct.keys():
                item_dct[key] = []
            for other_ID in self.other_items.keys():
                other = self.other_items[other_ID]['item']
                item_dct['Other'].append(f'{other_ID} [{other.functional_unit}]')
                quantity = self.other_items[other_ID]['quantity']
                item_dct['Quantity'].append(quantity)
                for ind in self.indicators:
                    if ind.ID in other.CFs.keys():
                        impact = other.CFs[ind.ID]*quantity
                        item_dct[f'{ind.ID} [{ind.unit}]'].append(impact)
                        item_dct[f'Category {ind.ID} Ratio'].append(impact/tot[ind.ID])
                    else:
                        item_dct[f'{ind.ID} [{ind.unit}]'].append(0)
                        item_dct[f'Category {ind.ID} Ratio'].append(0)

            table = pd.DataFrame.from_dict(item_dct)
            table.set_index(['Other'], inplace=True)
            return self._append_cat_sum(table, cat, tot)

        raise ValueError(
            'category can only be "Construction", "Transportation", "Stream", or "Other", ' \
            f'not "{category}".')


    def save_report(self, file=None, sheet_name='LCA',
                    n_row=0, row_space=2, annual=False):
        '''
        Save all LCA tables as an Excel file.
        
        Parameters
        ----------
        file : str
            Path and name of the excel file,
            will use the ID of the system with an '_lca' suffix, if not provided.
        annual : bool
            If True, will return the annual impacts considering `uptime_ratio`
            instead of across the system lifetime.
        '''
        if not file:
            file = f'{self.system.ID}_lca.xlsx'
        tables = [self.get_impact_table(cat, annual=annual)
                  for cat in ('Construction', 'Transportation',
                              'Stream', 'other')]
        with pd.ExcelWriter(file) as writer:
            for table in tables:
                if isinstance(table, str): continue
                table.to_excel(writer, sheet_name=sheet_name, startrow=n_row)
                n_row += table.shape[0] + row_space + len(table.columns.names) # extra lines for the heading


    @property
    def system(self):
        '''[biosteam.System] The system linked to this LCA.'''
        return self._system
    @system.setter
    def system(self, i):
        self._update_system(i)

    @property
    def lifetime(self):
        '''[int] Lifetime of the system, [yr].'''
        return self._lifetime
    @lifetime.setter
    def lifetime(self, lifetime):
        self._lifetime = int(lifetime)

    @property
    def lifetime_hr(self):
        '''[float] Lifetime of the system in hours, [hr].'''
        return self._lifetime*365*24*self.uptime_ratio

    @property
    def uptime_ratio(self):
        '''[float] Fraction of time that the system is operating.'''
        return self._uptime_ratio
    @uptime_ratio.setter
    def uptime_ratio(self, i):
        if 0 <=i<= 1:
            self._uptime_ratio = float(i)
        else:
            raise ValueError('uptime_ratio must be in [0,1].')

    @property
    def indicators(self):
        '''
        [list] All impact indicators associated with this LCA object.
        If not `ImpactIndicator` has been added, then will be defaulted to
        sum of the `ImpactIndicator` objects added to the system associated
        with this LCA (e.g., associated with construction, streams, etc.

        '''
        if self._indicators:
            return self._indicators

        if not self.construction_inventory:
            constr = set()
        else:
            constr = set(sum((i.indicators for i in self.construction_inventory
                              if i is not None), ()))
        if not self.transportation_inventory:
            trans = set()
        else:
            trans = set(sum((i.indicators for i in self.transportation_inventory
                             if i is not None), ()))
        if not self.stream_inventory:
            ws = set()
        else:
            ws = set(sum((i.indicators for i in self.stream_inventory
                          if i is not None), ()))
        if not self.other_items:
            other = set()
        else:
            other = set(sum((ImpactItem.get_item(i).indicators
                             for i in self.other_items.keys()), ()))
        tot = constr.union(trans, ws, other)
        if len(tot) == 0:
            warn('No `ImpactIndicator` has been added.')
        return list(tot)
    @indicators.setter
    def indicators(self, i):
        if not (isinstance(i, Iterable) and not isinstance(i, str)):
            i = (i,)
        inds = []
        for ind in i:
            if isinstance(ind, str):
                ind = ImpactIndicator.get_indicator(ind)
            if not isinstance(ind, ImpactIndicator):
                raise TypeError(f'{ind} is not an `ImpactIndicator` or ID/alias of an `ImpactIndicator`.')
            inds.append(ind)
        self._indicators = inds

    @property
    def construction_units(self):
        '''[set] All units in the linked system with construction activity.'''
        return self._construction_units

    @property
    def construction_inventory(self):
        '''[tuple] All construction activities.'''
        return sum((i.construction for i in self.construction_units), [])

    @property
    def total_construction_impacts(self):
        '''[dict] Total impacts associated with construction activities.'''
        return self.get_construction_impacts(self.construction_units)

    @property
    def transportation_units(self):
        '''[set] All units in the linked system with transportation activity.'''
        return self._transportation_units

    @property
    def transportation_inventory(self):
        '''[tuple] All transportation activities.'''
        return sum((i.transportation for i in self.transportation_units), [])

    @property
    def total_transportation_impacts(self):
        '''[dict] Total impacts associated with transportation activities.'''
        return self.get_transportation_impacts(self.transportation_units)

    @property
    def lca_streams(self):
        '''[set] All streams in the linked system with impacts.'''
        return self._lca_streams

    @property
    def stream_inventory(self):
        '''[tuple] All chemical inputs, fugitive gases, waste emissions, and products.'''
        return tuple(i.stream_impact_item for i in self.lca_streams)

    @property
    def total_stream_impacts(self):
        '''[dict] Total impacts associated with `WasteStreams` (e.g., chemicals, emissions).'''
        return self.get_stream_impacts(stream_items=self.stream_inventory)

    @property
    def other_items (self):
        '''[dict] Other impact items (e.g., electricity) and their quantities.'''
        return self._other_items
    @other_items.setter
    def other_items(self, item, f_quantity, unit=''):
        self.add_other_item(item, f_quantity, unit)

    @property
    def total_other_impacts(self):
        '''[dict] Total impacts associated with other ImpactItems (e.g., electricity).'''
        return self.get_other_impacts()

    @property
    def total_impacts(self):
        '''[dict] Total impacts of the entire system (construction, transportation, and wastestream).'''
        return self.get_total_impacts()