#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems
Copyright (C) 2020, Quantitative Sustainable Design Group

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the UIUC open-source license. Please refer to 
https://github.com/QSD-Group/QSDsan/blob/master/LICENSE.txt
for license details.

TODO:
    Add a function to save LCA details.
'''


# %%

import pandas as pd
from . import ImpactItem, StreamImpactItem
from ._units_of_measure import auom
from .utils.formatting import format_number as f_num

items = ImpactItem._items

__all__ = ('LCA',)


class LCA:
    '''
    For life cycle assessment (LCA) of a System.
    
    Parameters
    ----------
    system : biosteam.System
        System for which this LCA is conducted for.
    lifetime : float
        Lifetime of the LCA.
    lifetime_unit : str
        Unit of lifetime.
    uptime_ratio : float
        Fraction of time that the plant is operating.
    **item_quantities : kwargs, ImpactItem or str = float or (float, unit)
        Other ImpactItems (e.g., electricity) and their quantities.
    
    '''
    
    __slots__ = ('_system',  '_lifetime', '_uptime_ratio',
                 '_construction_units', '_transportation_units',
                 '_lca_waste_streams',
                 '_impact_indicators', '_other_items')
    
    
    def __init__(self, system, lifetime, lifetime_unit='yr', uptime_ratio=1,
                 **item_quantities):

        self._construction_units = set()
        self._transportation_units = set()
        self._lca_waste_streams = set()
        self._update_system(system)
        self._update_lifetime(lifetime, lifetime_unit)
        self.uptime_ratio = uptime_ratio
        self._other_items = {}
        for item, quantity in item_quantities.items():
            try:
                q_number, q_unit = quantity # unit provided for the quantity
                self.add_other_item(item, q_number, q_unit)
            except:
                self.add_other_item(item, quantity)
                
    
    def _update_system(self, system):
        for unit in system.units:
            if unit.construction:
                self._construction_units.add(unit)
            if unit.transportation:
                self._transportation_units.add(unit)
        self._construction_units = sorted(self._construction_units,
                                          key=lambda i: type(i).__name__)
        self._transportation_units = sorted(self._transportation_units,
                                            key=lambda i: type(i).__name__)
        for ws in system.streams:
            if ws.impact_item:
                self._lca_waste_streams.add(ws)
        self._lca_waste_streams = sorted(self._lca_waste_streams,
                                          key=lambda i: i.ID)
        self._system = system
    
    def _refresh_ws_items(self, system=None):
        if not system:
            system = self.system
        for item in items.values():
            if isinstance(item, StreamImpactItem):
                ws = item.linked_ws
                ws._impact_item = item
                if ws in system.streams and ws not in self._lca_waste_streams:
                    self._lca_waste_streams.add(ws)
        self._lca_waste_streams = sorted(self._lca_waste_streams,
                                          key=lambda i: i.ID)
        
    
    def _update_lifetime(self, lifetime=0., unit='yr'):
        if not unit or unit == 'yr':
            self._lifetime = float(lifetime)
        else:
            converted = auom(unit).convert(float(lifetime), 'yr')
            self._lifetime = converted
    
    def add_other_item(self, item, quantity, unit='', return_dict=None):
        if isinstance(item, str):
            item = items[item]
        fu = item.functional_unit
        if unit and unit != fu:
            try:
                quantity = auom(unit).convert(quantity, fu)
            except:
                raise ValueError(f'Conversion of the given unit {unit} to '
                                 f'item functional unit {fu} is not supported.')
        if not return_dict:
            self.other_items[item.ID] = quantity
        else:
            return_dict[item.ID] = quantity
            return return_dict
      
    def __repr__(self):
        return f'<LCA: {self.system}>'

    def show(self, lifetime_unit='yr'):
        lifetime = auom('yr').convert(self.lifetime, lifetime_unit)
        info = f'LCA: {self.system} (lifetime {f_num(lifetime)} {lifetime_unit})'
        info += '\nImpacts:'
        print(info)
        if len(self.indicators) == 0:
            print(' None')
        else:
            index = pd.Index((i.ID+' ('+i.unit+')' for i in self.indicators))
            df = pd.DataFrame({
                'Construction': tuple(self.construction_impacts.values()),
                'Transportation': tuple(self.transportation_impacts.values()),
                'WasteStream': tuple(self.waste_stream_impacts.values()),
                'Others': tuple(self.other_impacts.values()),
                'Total': tuple(self.total_impacts.values())
                },
                index=index)
            # print(' '*9+df.to_string().replace('\n', '\n'+' '*9))
            print(df.to_string())
    
    _ipython_display_ = show
    
    
    def _get_constr_impacts(self, iterable):
        impacts = dict.fromkeys((i.ID for i in self.indicators), 0.)
        for i in iterable:
            for j in i.construction:
                impact = j.impacts
                for m, n in impact.items():
                    impacts[m] += n
        return impacts
    
    def _get_trans_impacts(self, iterable):
        impacts = dict.fromkeys((i.ID for i in self.indicators), 0.)
        hr = self.lifetime*365*24*self.uptime_ratio
        for i in iterable:
            for j in i.transportation:
                impact = j.impacts
                for m, n in impact.items():
                    impacts[m] += n*hr/j.interval
        return impacts
    
    
    def _get_ws_impacts(self, iterable=None, exclude=None):
        impacts = dict.fromkeys((i.ID for i in self.indicators), 0.)
        hr = self.lifetime*365*24*self.uptime_ratio
        if not iterable:
            iterable = self.waste_stream_inventory
        if None in iterable:
            self._refresh_ws_items()
            iterable = self.waste_stream_inventory
        for j in iterable:
            # breakpoint()
            ws = j.linked_ws
            if ws is exclude: continue
            for m, n in j.CFs.items():
                impacts[m] += n*hr*ws.F_mass
        return impacts
    
    def _get_add_impacts(self, **item_quantities):
        impacts = dict.fromkeys((i.ID for i in self.indicators), 0.)
        for j, k in item_quantities.items():
            if isinstance(j, str):
                item = items[j]
            for m, n in item.CFs.items():
                impacts[m] += n*k
        return impacts
    
    def _get_total_impacts(self, exclude=None):
        impacts = dict.fromkeys((i.ID for i in self.indicators), 0.)
        ws_impacts = self._get_ws_impacts(exclude)
        for i in (self.construction_impacts,
                  self.transportation_impacts,
                  ws_impacts,
                  self.other_impacts):
            for m, n in i.items():
                impacts[m] += n
        return impacts
        
    def get_normalized_impacts(self, waste_stream):
        '''[dict] Normalize all impacts based on the mass flow of a WasteStream.'''
        assert waste_stream in self.system.streams, \
               f'WasteStream {waste_stream} not in the System {self.system}.'
        impacts = self._get_total_impacts(exclude=waste_stream)
        for i in impacts.values():
            i /= waste_stream.F_mass
        return impacts
    
    def get_units_impacts(self, units, exclude=None, **item_quantities):
        try: iter(units)
        except: units = (units,)
        constr = self._get_constr_impacts(units)
        trans = self._get_trans_impacts(units)
        ws_list = set()
        for i in (unit.ins+unit.outs for unit in units):
            if i in self._lca_waste_streams:
                ws_list.add(i)
        ws = self._get_ws_impacts(iterable=ws_list, exclude=exclude)
        add = self._get_add_impacts(**item_quantities)
        for m, n in constr.items():
            n += trans[m] + ws[m] + add[m]
        return constr
    
    @property
    def system(self):
        '''[biosteam.System] The System linked to this LCA.'''
        return self._system
    @system.setter
    def system(self, i):
        self._update_system(i)
    
    @property
    def lifetime(self):
        '''[float] Lifetime of the system, [yr].'''
        return self._lifetime
    @lifetime.setter
    def lifetime(self, lifetime, unit='yr'):
        self._update_lifetime(lifetime, unit)
    
    @property
    def uptime_ratio(self):
        '''[float] Fraction of time that the plant is operating.'''
        return self._uptime_ratio
    @uptime_ratio.setter
    def uptime_ratio(self, i):
        if 0 <= i <= 1:
            self._uptime_ratio = float(i)
        else:
            raise ValueError('uptime_ratio must be in [0,1].')
    
    @property
    def indicators(self):
        '''[set] All ImpactIndicators associated with this LCA.'''
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
        if not self.waste_stream_inventory:
            ws = set()
        else:
            ws = set(sum((i.indicators for i in self.waste_stream_inventory
                          if i is not None), ()))
        if not self.other_items:
            add = set()
        else:
            add = set(sum((items[i].indicators for i in self.other_items.keys()), ()))
        return constr.union(trans, ws, add)
    
    @property
    def construction_units(self):
        '''[set] All units in the linked System with constrution activity.'''
        return self._construction_units
    
    @property
    def construction_inventory(self):
        '''[tuple] All construction activities.'''
        return sum((i.construction for i in self.construction_units), ())
    
    @property
    def construction_impacts(self):
        '''[dict] Total impacts associated with construction activities.'''
        return self._get_constr_impacts(self.construction_units)
    
    @property
    def transportation_units(self):
        '''[set] All units in the linked System with transportation activity.'''
        return self._transportation_units
    
    @property
    def transportation_inventory(self):
        '''[tuple] All transportation activities.'''
        return sum((i.transportation for i in self.transportation_units), ())
    
    @property
    def transportation_impacts(self):
        '''[dict] Total impacts associated with transportation activities.'''
        return self._get_trans_impacts(self.transportation_units)
    
    @property
    def lca_waste_streams(self):
        '''[set] All WasteStreams in the linked System with StreamImpactItems.'''
        return self._lca_waste_streams
    
    @property
    def waste_stream_inventory(self):
        '''[tuple] All chemical inputs, fugitive gases, waste emissions, and products.'''
        return tuple(i.impact_item for i in self.lca_waste_streams)
    
    @property
    def waste_stream_impacts(self):
        '''[dict] Total impacts associated with WasteStreams (e.g., chemicals, emissions).'''
        return self._get_ws_impacts()
        
    @property
    def other_items (self):
        '''[dict] Other ImpactItems (e.g., electricity) and their quantities.'''
        return self._other_items
    @other_items.setter
    def other_items(self, item, quantity, unit=''):
        self.add_other_item(item, quantity, unit)
        
    @property
    def other_impacts(self):
        '''[dict] Total impacts associated with other ImpactItems (e.g., electricity).'''
        return self._get_add_impacts(**self.other_items)
    
    @property
    def total_impacts(self):
        '''[dict] Total impacts of the entire system (construction, transportation, and wastestream).'''
        return self._get_total_impacts()






