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
'''


# %%

import numpy as np
import pandas as pd
import qsdsan as qs
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
    **item_quantities : kwargs, ImpactItem or str = float/callable or (float/callable, unit)
        Other ImpactItems (e.g., electricity) and their quantities. Note that
        callable functions are used so that quantity of items can be updated.
    
    '''
    
    __slots__ = ('_system',  '_lifetime', '_uptime_ratio',
                 '_construction_units', '_transportation_units',
                 '_lca_streams', '_impact_indicators',
                 '_other_items', '_other_items_f')
    
    
    def __init__(self, system, lifetime, lifetime_unit='yr', uptime_ratio=1,
                 **item_quantities):
        system.simulate()
        self._construction_units = set()
        self._transportation_units = set()
        self._lca_streams = set()
        self._update_system(system)
        self._update_lifetime(lifetime, lifetime_unit)
        self.uptime_ratio = uptime_ratio
        self._other_items = {}
        self._other_items_f = {}
        for item, val in item_quantities.items():
            try:
                f_quantity, unit = val # unit provided for the quantity
            except:
                f_quantity = val
                unit = ''
            if not callable(f_quantity):
                f_quantity = lambda: f_quantity
            self.add_other_item(item, f_quantity, unit)
                
    
    def _update_system(self, system):
        for u in system.units:
            if u.construction:
                self._construction_units.add(u)
            if u.transportation:
                self._transportation_units.add(u)
        self._construction_units = sorted(self._construction_units,
                                          key=lambda u: u.ID)
        self._transportation_units = sorted(self._transportation_units,
                                            key=lambda u: u.ID)
        for s in system.streams:
            if s.impact_item:
                self._lca_streams.add(s)
        self._lca_streams = sorted(self._lca_streams, key=lambda s: s.ID)
        self._system = system
    
    def _refresh_stream_items(self, system=None):
        if not system:
            system = self.system
        for item in items.values():
            if isinstance(item, StreamImpactItem):
                ws = item.linked_stream
                ws._impact_item = item
                if ws in system.streams and ws not in self._lca_streams:
                    self._lca_streams.add(ws)
        self._lca_streams = sorted(self._lca_streams,
                                         key=lambda ws: ws.ID)
        
    def _update_lifetime(self, lifetime=0., unit='yr'):
        if not unit or unit == 'yr':
            self._lifetime = float(lifetime)
        else:
            converted = auom(unit).convert(float(lifetime), 'yr')
            self._lifetime = converted
    
    def add_other_item(self, item, f_quantity, unit='', return_dict=None):
        '''Add other ``ImpactItem`` in LCA.'''
        if isinstance(item, str):
            item = items[item]
        fu = item.functional_unit
        quantity = f_quantity()
        if unit and unit != fu:
            try:
                quantity = auom(unit).convert(quantity, fu)
            except:
                raise ValueError(f'Conversion of the given unit {unit} to '
                                 f'item functional unit {fu} is not supported.')
        self._other_items_f[item.ID] = {'item':item, 'f_quantity':f_quantity, 'unit':unit}
        if not return_dict:
            self.other_items[item.ID] = {'item':item, 'quantity':quantity}
        else:
            return_dict[item.ID] = {'item':item, 'quantity':quantity}
            return return_dict
    
    def refresh_other_items(self):
        '''Refresh quantities of other items using the given functions.'''
        for item, record in self._other_items_f.keys():
            item, quantity, unit = record.values()
            self.add_other_item(item, quantity, unit)
        
        
    def __repr__(self):
        return f'<LCA: {self.system}>'

    def show(self, lifetime_unit='yr'):
        '''Show basic information of this ``LCA`` object.'''
        lifetime = auom('yr').convert(self.lifetime, lifetime_unit)
        info = f'LCA: {self.system} (lifetime {f_num(lifetime)} {lifetime_unit})'
        info += '\nImpacts:'
        print(info)
        if len(self.indicators) == 0:
            print(' None')
        else:
            index = pd.Index((i.ID+' ('+i.unit+')' for i in self.indicators))
            df = pd.DataFrame({
                'Construction': tuple(self.total_construction_impacts.values()),
                'Transportation': tuple(self.total_transportation_impacts.values()),
                'WasteStream': tuple(self.total_stream_impacts.values()),
                'Others': tuple(self.total_other_impacts.values()),
                'Total': tuple(self.total_impacts.values())
                },
                index=index)
            # print(' '*9+df.to_string().replace('\n', '\n'+' '*9))
            print(df.to_string())
    
    _ipython_display_ = show
    
    
    def get_construction_impacts(self, units, time=None, time_unit='hr'):
        '''
        Return all construction-related impacts for the given unit,
        normalized to a certain time frame.
        '''
        try: iter(units)
        except: units = (units,)
        if not time:
            ratio = 1
        else:
            converted = auom(time_unit).convert(float(time), 'hr')
            ratio = converted/self.lifetime_hr
        impacts = dict.fromkeys((i.ID for i in self.indicators), 0.)
        for i in units:
            if i.lifetime:
                factor = self.lifetime/i.lifetime
            else:
                factor = 1
            for j in i.construction:
                impact = j.impacts
                for m, n in impact.items():
                    impacts[m] += n*ratio*factor
        return impacts
    
    def get_transportation_impacts(self, units, time=None, time_unit='hr'):
        '''
        Return all transportation-related impacts for the given unit,
        normalized to a certain time frame.
        '''
        try: iter(units)
        except: units = (units,)
        if not time:
            time = self.lifetime_hr
        else:
            time = auom(time_unit).convert(float(time), 'hr')
        impacts = dict.fromkeys((i.ID for i in self.indicators), 0.)
        for i in units:
            for j in i.transportation:
                impact = j.impacts
                for m, n in impact.items():
                    impacts[m] += n*time/j.interval
        return impacts
    
    
    def get_stream_impacts(self, stream_items=None, exclude=None,
                           kind='all', time=None, time_unit='hr'):
        '''
        Return all stream-related impacts for the given unit,
        normalized to a certain time frame.
        '''
        try: iter(stream_items)
        except: stream_items = (stream_items,)
        impacts = dict.fromkeys((i.ID for i in self.indicators), 0.)
        if not time:
            time = self.lifetime_hr
        else:
            time = auom(time_unit).convert(float(time), 'hr')
        if None in stream_items:
            self._refresh_stream_items()
            stream_items = self.stream_inventory
        for j in stream_items:
            # In case that ws instead of the item is given
            if isinstance(j, qs.WasteStream):
                ws = j
                j = ws.impact_item
            else:
                ws = j.linked_stream
            if ws is exclude: continue
            for m, n in j.CFs.items():
                if kind == 'all':
                    pass
                elif kind == 'direct_emission':
                    n = max(n, 0)
                elif kind == 'offset':
                    n = min(n, 0)
                else:
                    raise ValueError('kind can only be "all", "direct_emission", or "offset", '
                                     f'not {kind}.')
                impacts[m] += n*time*ws.F_mass
        return impacts
    
    def get_other_impacts(self, **item_quantities):
        '''
        Return all additional impacts from "other" ``ImpactItems`` objects,
        based on defined quantity.
        '''
        impacts = dict.fromkeys((i.ID for i in self.indicators), 0.)
        for j, k in item_quantities.items():
            if isinstance(j, str):
                item = items[j]
            for m, n in item.CFs.items():
                impacts[m] += n*k
        return impacts
    
    def get_total_impacts(self, exclude=None, time=None, time_unit='hr'):
        '''Return total impacts, normalized to a certain time frame.'''
        impacts = dict.fromkeys((i.ID for i in self.indicators), 0.)
        ws_impacts = self.get_stream_impacts(stream_items=self.stream_inventory,
                                                   exclude=exclude, time=time, time_unit=time_unit)
        for i in (self.total_construction_impacts,
                  self.total_transportation_impacts,
                  ws_impacts,
                  self.total_other_impacts):
            for m, n in i.items():
                impacts[m] += n
        return impacts
        
    def get_normalized_impacts(self, stream):
        '''Normalize all impacts based on the mass flow of a stream.'''
        assert stream in self.system.streams, \
               f'WasteStream {stream} not in the System {self.system}.'
        impacts = self.get_total_impacts(exclude=stream)
        for i in impacts.values():
            i /= stream.F_mass
        return impacts
    
    def get_units_impacts(self, units, time=None, time_unit='hr',
                          exclude=None, **item_quantities):
        '''Return total impacts with certain units, normalized to a certain time frame. '''
        try: iter(units)
        except: units = (units,)
        constr = self.get_construction_impacts(units, time, time_unit)
        trans = self.get_transportation_impacts(units, time, time_unit)
        ws_items = set()
        for i in sum((tuple(unit.ins+unit.outs) for unit in units), ()):
            if i in self._lca_streams:
                ws_items.add(i.impact_item)
        ws = self.get_stream_impacts(stream_items=ws_items,
                                           exclude=exclude,
                                           time=time, time_unit=time_unit)
        other = self.get_other_impacts(**item_quantities)
        tot = constr.copy()
        for m in tot.keys():
            tot[m] += trans[m] + ws[m] + other[m]
        return tot
    
    def _append_cat_sum(self, cat_table, cat, tot):
        num = len(cat_table)
        cat_table.loc[num] = ''
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
    
    def get_impact_table(self, category=None, time=None, time_unit='hr'):
        '''
        Return a ``pandas.DataFrame`` table for the given impact category,
        normalized to a certain time frame.
        '''
        if not time:
            time = self.lifetime_hr
        else:
            time = auom(time_unit).convert(float(time), 'hr')
        
        if category in ('Construction', 'Other'):
            time = time/self.lifetime_hr
        
        cat = category.lower()
        tot = getattr(self, f'total_{cat}_impacts')
        if category in ('Construction', 'Transportation'):
            cat = category.lower()
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
                for i in getattr(su, cat):
                    item_dct[i.item.ID]['SanUnit'].append(su.ID)
                    item_dct[i.item.ID]['Quantity'].append(i.quantity*time)
                    if cat == 'transportation':
                        item_dct[i.item.ID]['Quantity'][-1] /= i.interval
            dfs = []
            for item in items:
                dct = item_dct[item.ID]
                dct['SanUnit'].append('Total')
                dct['Quantity'] = np.array(dct['Quantity'])
                dct['Quantity'] = np.append(dct['Quantity'], dct['Quantity'].sum())
                dct['Item Ratio'] = dct['Quantity']/dct['Quantity'].sum()*2
                for i in self.indicators:
                    dct[f'{i.ID} [{i.unit}]'] = impact = dct['Quantity']*item.CFs[i.ID]
                    dct[f'Category {i.ID} Ratio'] = impact/tot[i.ID]
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
        
        if category == 'Stream':
            headings = ['Stream', 'Mass [kg]', *ind_head]
            item_dct = dict.fromkeys(headings)
            for key in item_dct.keys():
                item_dct[key] = []
            for ws_item in self.stream_inventory:
                ws = ws_item.linked_stream
                item_dct['Stream'].append(ws.ID)
                mass = ws.F_mass * time
                item_dct['Mass [kg]'].append(mass)
                for ind in self.indicators:
                    impact = ws_item.CFs[ind.ID]*mass
                    item_dct[f'{ind.ID} [{ind.unit}]'].append(impact)
                    item_dct[f'Category {ind.ID} Ratio'].append(impact/tot[ind.ID])
            table = pd.DataFrame.from_dict(item_dct)
            table.set_index(['Stream'], inplace=True)
            return self._append_cat_sum(table, cat, tot)

        elif category == 'Other':
            headings = ['Other', 'Quantity', *ind_head]
            item_dct = dict.fromkeys(headings)
            for key in item_dct.keys():
                item_dct[key] = []
            for other_ID in self.other_items.keys():
                other = self.other_items[other_ID]['item']
                item_dct['Other'].append(f'{other_ID} [{other.functional_unit}]')
                quantity = self.other_items[other_ID]['quantity'] * time
                item_dct['Quantity'].append(quantity)
                for ind in self.indicators:
                    impact = other.CFs[ind.ID]*quantity
                    item_dct[f'{ind.ID} [{ind.unit}]'].append(impact)
                    item_dct[f'Category {ind.ID} Ratio'].append(impact/tot[ind.ID])
            table = pd.DataFrame.from_dict(item_dct)
            table.set_index(['Other'], inplace=True)
            return self._append_cat_sum(table, cat, tot)
        
        else:
            raise ValueError(
                "category can only be 'Construction', 'Transportation', 'Stream', or 'Other', " \
                f'not {category}.')

    def save_report(self, file=None, sheet_name='LCA',
                    time=None, time_unit='hr',
                    n_row=0, row_space=2):
        '''Save all LCA tables as an Excel file.'''
        if not file:
            file = f'{self.system.ID}_lca.xlsx'
        tables = [self.get_impact_table(cat, time, time_unit)
                  for cat in ('Construction', 'Transportation',
                              'Stream', 'Other')]
        with pd.ExcelWriter(file) as writer:
            for table in tables:
                table.to_excel(writer, sheet_name=sheet_name, startrow=n_row)
                n_row += table.shape[0] + row_space + 1 # one extra line for the heading

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
    def lifetime_hr(self):
        '''[float] Lifetime of the system in hours, [hr].'''
        return self._lifetime*365*24*self.uptime_ratio
    
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
        if not self.stream_inventory:
            ws = set()
        else:
            ws = set(sum((i.indicators for i in self.stream_inventory
                          if i is not None), ()))
        if not self.other_items:
            other = set()
        else:
            other = set(sum((items[i].indicators for i in self.other_items.keys()), ()))
        tot = constr.union(trans, ws, other)
        if len(tot) == 0:
            raise ValueError('No ImpactIndicators have been added.')
        return tot
    
    @property
    def construction_units(self):
        '''[set] All units in the linked system with constrution activity.'''
        return self._construction_units
    
    @property
    def construction_inventory(self):
        '''[tuple] All construction activities.'''
        return sum((i.construction for i in self.construction_units), ())
    
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
        return sum((i.transportation for i in self.transportation_units), ())
    
    @property
    def total_transportation_impacts(self):
        '''[dict] Total impacts associated with transportation activities.'''
        return self.get_transportation_impacts(self.transportation_units)
    
    @property
    def lca_streams(self):
        '''[set] All streams in the linked system with StreamImpactItems.'''
        return self._lca_streams
    
    @property
    def stream_inventory(self):
        '''[tuple] All chemical inputs, fugitive gases, waste emissions, and products.'''
        return tuple(i.impact_item for i in self.lca_streams)
    
    @property
    def total_stream_impacts(self):
        '''[dict] Total impacts associated with WasteStreams (e.g., chemicals, emissions).'''
        return self.get_stream_impacts(stream_items=self.stream_inventory)
        
    @property
    def other_items (self):
        '''[dict] Other ImpactItems (e.g., electricity) and their quantities.'''
        return self._other_items
    @other_items.setter
    def other_items(self, item, f_quantity, unit=''):
        self.add_other_item(item, f_quantity, unit)
        
    @property
    def total_other_impacts(self):
        '''[dict] Total impacts associated with other ImpactItems (e.g., electricity).'''
        items = self.other_items.keys()
        quantities = (self.other_items[i]['quantity'] for i in self.other_items.keys())
        return self.get_other_impacts(**dict(zip(items, quantities)))
    
    @property
    def total_impacts(self):
        '''[dict] Total impacts of the entire system (construction, transportation, and wastestream).'''
        return self.get_total_impacts()






