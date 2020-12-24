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
                 '_construction_sanunits', '_transportation_sanunits',
                 '_lca_waste_streams',
                 '_impact_indicators', '_other_items')
    
    
    def __init__(self, system, lifetime, lifetime_unit='yr', uptime_ratio=1,
                 **item_quantities):

        self._construction_sanunits = set()
        self._transportation_sanunits = set()
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
        for su in system.units:
            if su.construction:
                self._construction_sanunits.add(su)
            if su.transportation:
                self._transportation_sanunits.add(su)
        self._construction_sanunits = sorted(self._construction_sanunits,
                                             key=lambda su: su.ID)
        self._transportation_sanunits = sorted(self._transportation_sanunits,
                                               key=lambda su: su.ID)
        for ws in system.streams:
            if ws.impact_item:
                self._lca_waste_streams.add(ws)
        self._lca_waste_streams = sorted(self._lca_waste_streams,
                                         key=lambda ws: ws.ID)
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
                                         key=lambda ws: ws.ID)
        
    
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
            self.other_items[item.ID] = {'item':item, 'quantity':quantity}
        else:
            return_dict[item.ID] = {'item':item, 'quantity':quantity}
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
                'Construction': tuple(self.total_construction_impacts.values()),
                'Transportation': tuple(self.total_transportation_impacts.values()),
                'WasteStream': tuple(self.total_waste_stream_impacts.values()),
                'Others': tuple(self.total_other_impacts.values()),
                'Total': tuple(self.total_impacts.values())
                },
                index=index)
            # print(' '*9+df.to_string().replace('\n', '\n'+' '*9))
            print(df.to_string())
    
    _ipython_display_ = show
    
    
    def get_construction_impacts(self, iterable, time=None, time_unit='hr'):
        if not time:
            ratio = 1
        else:
            converted = auom(time_unit).convert(float(time), 'hr')
            ratio = converted/self.lifetime_hr
        impacts = dict.fromkeys((i.ID for i in self.indicators), 0.)
        for i in iterable:
            for j in i.construction:
                impact = j.impacts
                for m, n in impact.items():
                    impacts[m] += n*ratio
        return impacts
    
    def get_transportation_impacts(self, iterable, time=None, time_unit='hr'):
        if not time:
            time = self.lifetime_hr
        else:
            time = auom(time_unit).convert(float(time), 'hr')
        impacts = dict.fromkeys((i.ID for i in self.indicators), 0.)
        for i in iterable:
            for j in i.transportation:
                impact = j.impacts
                for m, n in impact.items():
                    impacts[m] += n*time/j.interval
        return impacts
    
    
    def get_waste_stream_impacts(self, iterable=None, exclude=None, time=None, time_unit='hr'):
        impacts = dict.fromkeys((i.ID for i in self.indicators), 0.)
        if not time:
            time = self.lifetime_hr
        else:
            time = auom(time_unit).convert(float(time), 'hr')
        if not iterable:
            iterable = self.waste_stream_inventory
        if None in iterable:
            self._refresh_ws_items()
            iterable = self.waste_stream_inventory
        for j in iterable:
            ws = j.linked_ws
            if ws is exclude: continue
            for m, n in j.CFs.items():
                impacts[m] += n*time*ws.F_mass
        return impacts
    
    def get_other_impacts(self, **item_quantities):
        if not item_quantities:
            return self.total_other_impacts
        impacts = dict.fromkeys((i.ID for i in self.indicators), 0.)
        for j, k in item_quantities.items():
            if isinstance(j, str):
                item = items[j]
            for m, n in item.CFs.items():
                impacts[m] += n*k
        return impacts
    
    def get_total_impacts(self, exclude=None, time=None, time_unit='hr'):
        impacts = dict.fromkeys((i.ID for i in self.indicators), 0.)
        ws_impacts = self.get_waste_stream_impacts(exclude=exclude, time=time, time_unit=time_unit)
        for i in (self.total_construction_impacts,
                  self.total_transportation_impacts,
                  ws_impacts,
                  self.total_other_impacts):
            for m, n in i.items():
                impacts[m] += n
        return impacts
        
    def get_normalized_impacts(self, waste_stream):
        '''[dict] Normalize all impacts based on the mass flow of a WasteStream.'''
        assert waste_stream in self.system.streams, \
               f'WasteStream {waste_stream} not in the System {self.system}.'
        impacts = self.get_total_impacts(exclude=waste_stream)
        for i in impacts.values():
            i /= waste_stream.F_mass
        return impacts
    
    def get_sanunits_impacts(self, units, exclude=None, **item_quantities):
        try: iter(units)
        except: units = (units,)
        constr = self.get_construction_impacts(units)
        trans = self.get_transportation_impacts(units)
        ws_list = set()
        for i in (unit.ins+unit.outs for unit in units):
            if i in self._lca_waste_streams:
                ws_list.add(i)
        ws = self.get_waste_stream_impacts(iterable=ws_list, exclude=exclude)
        other = self.get_other_impacts(**item_quantities)
        tot = constr.copy()
        for m, n in tot.items():
            n += trans[m] + ws[m] + other[m]
        return tot
    
    def _append_cat_sum(self, cat_table, cat, tot):
        num = len(cat_table)
        cat_table.loc[num] = ''
        for i in self.indicators:
            cat_table[f'{i.ID} [{i.unit}]'][num] = tot[i.ID]
            cat_table[f'Total {i.ID} Ratio'][num] = 1
        if cat in ('construction', 'transportation'):        
            cat_table.rename(index={num: ('Sum', 'All')}, inplace=True)
            cat_table.index = \
                pd.MultiIndex.from_tuples(cat_table.index,
                                          names=[cat.capitalize(), 'SanUnit'])
        else:
            cat_table.rename(index={num: 'Sum'}, inplace=True)
        return cat_table
    
    def get_impact_table(self, category=None, time=None, time_unit='hr'):
        if not time:
            time = self.lifetime_hr
        else:
            time = auom(time_unit).convert(float(time), 'hr')
        
        if category in ('Construction', 'Other'):
            time = time/self.lifetime_hr
        
        cat = category.lower()
        if cat == 'wastestream': cat = 'waste_stream'
        tot = getattr(self, f'total_{cat}_impacts')
        if category in ('Construction', 'Transportation'):
            cat = category.lower()
            sanunits = sorted(getattr(self, f'_{cat}_sanunits'),
                              key=(lambda su: su.ID))
            items = sorted(set(i.item for i in getattr(self,  f'{cat}_inventory')),
                           key=(lambda item: item.ID))
            if len(items) == 0:
                return f'No {cat}-related impacts.'

            # Note that item_dct = dict.fromkeys([item.ID for item in items], []) won't work
            item_dct = dict.fromkeys([item.ID for item in items])
            for item_ID in item_dct.keys():
                item_dct[item_ID] = dict(SanUnit=[], Quantity=[])
            for su in sanunits:
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
                    dct[f'Total {i.ID} Ratio'] = impact/tot[i.ID]
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
                         f'Total {i.ID} Ratio'] for i in self.indicators), [])
        
        if category == 'WasteStream':
            headings = ['WasteStream', 'Mass [kg]', *ind_head]
            item_dct = dict.fromkeys(headings)
            for key in item_dct.keys():
                item_dct[key] = []
            for ws_item in self.waste_stream_inventory:
                ws = ws_item.linked_ws
                item_dct['WasteStream'].append(ws.ID)
                mass = ws.F_mass * time
                item_dct['Mass [kg]'].append(mass)
                for ind in self.indicators:
                    impact = ws_item.CFs[ind.ID]*mass
                    item_dct[f'{ind.ID} [{ind.unit}]'].append(impact)
                    item_dct[f'Total {ind.ID} Ratio'].append(impact/tot[ind.ID])
            table = pd.DataFrame.from_dict(item_dct)
            table.set_index(['WasteStream'], inplace=True)
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
                    item_dct[f'Total {ind.ID} Ratio'].append(impact/tot[ind.ID])
            table = pd.DataFrame.from_dict(item_dct)
            table.set_index(['Other'], inplace=True)
            return self._append_cat_sum(table, cat, tot)
        
        else:
            raise ValueError(
                "category can only be 'Construction', 'Transportation', 'WasteStream', or 'Other', " \
                f'not {category}.')

    def save_lca_report(self, file='lca_report.xlsx', sheet_name='LCA',
                        time=None, time_unit='hr',
                        n_row=0, row_space=2):
        
        tables = [self.get_impact_table(cat, time, time_unit)
                  for cat in ('Construction', 'Transportation',
                              'WasteStream', 'Other')]
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
        if not self.waste_stream_inventory:
            ws = set()
        else:
            ws = set(sum((i.indicators for i in self.waste_stream_inventory
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
    def construction_sanunits(self):
        '''[set] All SanUnits in the linked System with constrution activity.'''
        return self._construction_sanunits
    
    @property
    def construction_inventory(self):
        '''[tuple] All construction activities.'''
        return sum((i.construction for i in self.construction_sanunits), ())
    
    @property
    def total_construction_impacts(self):
        '''[dict] Total impacts associated with construction activities.'''
        return self.get_construction_impacts(self.construction_sanunits)
    
    @property
    def transportation_sanunits(self):
        '''[set] All SanUnits in the linked System with transportation activity.'''
        return self._transportation_sanunits
    
    @property
    def transportation_inventory(self):
        '''[tuple] All transportation activities.'''
        return sum((i.transportation for i in self.transportation_sanunits), ())
    
    @property
    def total_transportation_impacts(self):
        '''[dict] Total impacts associated with transportation activities.'''
        return self.get_transportation_impacts(self.transportation_sanunits)
    
    @property
    def lca_waste_streams(self):
        '''[set] All WasteStreams in the linked System with StreamImpactItems.'''
        return self._lca_waste_streams
    
    @property
    def waste_stream_inventory(self):
        '''[tuple] All chemical inputs, fugitive gases, waste emissions, and products.'''
        return tuple(i.impact_item for i in self.lca_waste_streams)
    
    @property
    def total_waste_stream_impacts(self):
        '''[dict] Total impacts associated with WasteStreams (e.g., chemicals, emissions).'''
        return self.get_waste_stream_impacts()
        
    @property
    def other_items (self):
        '''[dict] Other ImpactItems (e.g., electricity) and their quantities.'''
        return self._other_items
    @other_items.setter
    def other_items(self, item, quantity, unit=''):
        self.add_other_item(item, quantity, unit)
        
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






