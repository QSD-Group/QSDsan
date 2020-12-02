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

import pandas as pd
from ._units_of_measure import auom
from .utils.formatting import format_number as f_num

__all__ = ('LCA',)


class LCA:
    '''For life cycle assessment (LCA) of a System.'''
    
    __slots__ = ('_system', '_construction_units', '_transportation_units',
                 '_lca_waste_streams', '_impact_indicators',
                 '_life_time')
    
    
    def __init__(self, system, life_time, life_time_unit='hr'):
        self._construction_units = set()
        self._transportation_units = set()
        self._lca_waste_streams = set()
        self._update_system(system)
        self._update_life_time(life_time, life_time_unit)
    
    def _update_system(self, system):
        for unit in system.units:
            if unit.construction:
                self._construction_units.add(unit)
            if unit.transportation:
                self._transportation_units.add(unit)
        for ws in system.streams:
            if ws.impact_item:
                self._lca_waste_streams.add(ws)
        self._system = system
        
    
    def _update_life_time(self, life_time=0., unit='hr'):
        if not unit or unit == 'hr':
            self._life_time = float(life_time)
        else:
            converted = auom(unit).convert(float(life_time), 'hr')
            self._life_time = converted
    
    def __repr__(self):
        return f'<LCA: {self.system}>'

    def show(self, life_time_unit='yr'):
        life_time = auom('hr').convert(self.life_time, life_time_unit)
        print(f'LCA: {self.system} (life time {f_num(life_time)} {life_time_unit})')
        index = pd.Index((i.ID+' ('+i.unit+')' for i in self.indicators))
        df = pd.DataFrame({
            'Construction': tuple(self.construction_impacts.values()),
            'Transportation': tuple(self.transportation_impacts.values()),
            'WasteStream': tuple(self.waste_stream_impacts.values()),
            'Total': tuple(self.total_impacts.values())
            },
            index=index)
        print(df)
    
    _ipython_display_ = show
    
    
    def _get_ws_impacts(self, exclude=None):
        impacts = dict.fromkeys((i.ID for i in self.indicators), 0.)
        hr = self.life_time
        for j in self.waste_stream_inventory:
            ws = j.linked_ws
            if ws is exclude: continue
            for m, n in j.CFs.items():
                impacts[m] += n*hr*ws.F_mass
        return impacts
    
    def _get_total_impacts(self, exclude=None):
        impacts = dict.fromkeys((i.ID for i in self.indicators), 0.)
        ws_impacts = self._get_ws_impacts(exclude)
        for i in (self.construction_impacts,
                  self.transportation_impacts,
                  ws_impacts):
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
        
    
    @property
    def system(self):
        '''[biosteam.System] The System linked to this LCA.'''
        return self._system
    @system.setter
    def system(self, i):
        self._update_system(i)
    
    @property
    def life_time (self):
        '''[float] Life time of the system, [hr].'''
        return self._life_time
    @life_time.setter
    def life_time(self, life_time, unit='hr'):
        self._update_life_time(life_time, unit)
    
    @property
    def indicators(self):
        '''[set] All ImpactIndicators associated with this LCA.'''
        constr = set(sum((i.indicators for i in self.construction_inventory), ()))
        trans = set(sum((i.indicators for i in self.transportation_inventory), ()))
        ws = set(sum((i.indicators for i in self.waste_stream_inventory), ()))
        return constr.union(trans, ws)
    
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
        impacts = dict.fromkeys((i.ID for i in self.indicators), 0.)
        for j in self.construction_inventory:
            impact = j.impacts
            for m, n in impact.items():
                impacts[m] += n
        return impacts
    
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
        impacts = dict.fromkeys((i.ID for i in self.indicators), 0.)
        hr = self.life_time
        for j in self.transportation_inventory:
            impact = j.impacts
            for m, n in impact.items():
                impacts[m] += n*hr/j.interval
        return impacts
    
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
    def total_impacts(self):
        '''[dict] Total impacts of the entire system (construction, transportation, and wastestream).'''
        return self._get_total_impacts()






