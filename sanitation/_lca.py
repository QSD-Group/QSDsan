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

import numpy as np
import pandas as pd
from copy import copy


__all__ = ('LCA',)


class LCA:
    '''
    For life cycle assessment (LCA) of a System.

    '''
    
    def __init__(self, system, CFs, life_time):
        self._system = system
        self.units = set(i for i in system.units if i._if_LCA)
        self._CFs = CFs
        self._life_time = life_time
    
    @property
    def CFs (self):
        '''
        [pandas.DataFrame] Characterization factors for different impact categories,
        the function unit is 1 kg for SanUnit construction materials,
        1 kJ for HeatUtility, and 1 kWh for PowerUtility.

        Notes
        -----
        The value should be negative for credits.

        '''
        return self._CFs
    @CFs.setter
    def CFs(self, i):
        self._CFs = i
    
    @property
    def life_time (self):
        '''[float] Life time of the system in hours.'''
        return self._life_time
    @life_time.setter
    def life_time(self, i):
        self._life_time = i
    
    def get_ws_impacts(self, exclude=(), summarize=True):
        CF_dict = dict.fromkeys(self.CFs.columns)
        unit_IDs = []
        ws_IDs = []
        ws_impacts = []
        for unit in self.units:
            for ws in unit.get_ws_inv():
                if ws in exclude: continue
                unit_IDs.append(unit.ID)
                ws_IDs.append(ws.ID)
                CF_dict.update(ws.CFs)
                ws_impact = np.asarray(list(CF_dict.values())) * ws.F_mass * self.life_time
                ws_impacts.append(ws_impact)
                CF_dict.clear()
        index = pd.MultiIndex(['WasteStream']*len(unit_IDs), unit_IDs, ws_IDs,
                              names=('Category', 'Unit', 'Item'))
        df = pd.DataFrame(ws_impacts, index=index)
        df.fillna(0)
        if summarize:
            df.loc[('WasteStream', 'All', 'All')] = df.sum()
        return df

    def get_heating_impacts(self, summarize=True):
        CF_array = np.asarray(list(self.CFs['heating']))
        unit_IDs = []
        hu_IDs = []
        hu_impacts = []
        for unit in self.units:
            for hu in unit.heat_utilities:
                if hu.duty*hu.flow>0:
                    unit_IDs.append(unit.ID)
                    hu_IDs.append(hu.ID)
                    hu_impact = CF_array * hu.duty * self.life_time
                    hu_impacts.append(hu_impact)
        index = pd.MultiIndex(['Heating']*len(unit_IDs), unit_IDs, hu_IDs,
                              names=('Category', 'Unit', 'Item'))
        df = pd.DataFrame(hu_impacts, index=index)
        df.fillna(0)
        if summarize:
            df.loc[('Heating', 'All', 'All')] = df.sum()
        return df
    
    def get_cooling_impacts(self, summarize=True):
        CF_array = np.asarray(list(self.CFs['cooling']))
        unit_IDs = []
        hu_IDs = []
        hu_impacts = []
        for unit in self.units:
            for hu in unit.heat_utilities:
                if hu.duty*hu.flow<0:
                    unit_IDs.append(unit.ID)
                    hu_IDs.append(hu.ID)
                    hu_impact = CF_array * -hu.duty * self.life_time
                    hu_impacts.append(hu_impact)
        index = pd.MultiIndex(['Cooling']*len(unit_IDs), unit_IDs, hu_IDs,
                              names=('Category', 'Unit', 'Item'))
        df = pd.DataFrame(hu_impacts, index=index)
        df.fillna(0)
        if summarize:
            df.loc[('Cooling', 'All', 'All')] = df.sum()
        return df

    def get_power_impacts(self, summarize=True):
        CF_array = np.asarray(list(self.CFs['power']))
        unit_IDs = []
        pu_impacts = []
        for unit in self.units:
            if not bool(unit.power_utility):
                unit_IDs.append(unit.ID)
                pu_impact = CF_array * unit.power_utility.rate * self.life_time
                pu_impacts.append(pu_impact)
        index = pd.MultiIndex(['Power']*len(unit_IDs), unit_IDs,
                              ['power_utility']*len(unit_IDs),
                              names=('Category', 'Unit', 'Item'))
        df = pd.DataFrame(pu_impacts, index=index)
        df.fillna(0)
        if summarize:
            df.loc[('Power', 'All', 'All')] = df.sum()
        return df

    def get_construction_material_impacts(self, summarize=True):
        CF_dict = dict.fromkeys(self.CFs.columns)
        unit_IDs = []
        cm_IDs = []
        cm_impacts = []
        for unit in self.untis:
            if not self.construction_materials: continue
            for material, mass in unit.construction_materials.items():
                unit_IDs.append(unit.ID)
                cm_IDs.append(material)
                CF_dict.updates(self.CFs[material])
                cm_impact = np.asarray(list(CF_dict.values())) * mass
                cm_impacts.append(cm_impact)
                CF_dict.clear()
        index = pd.MultiIndex(['Construction']*len(unit_IDs), unit_IDs, cm_IDs,
                              names=('Category', 'Unit', 'Item'))
        df = pd.DataFrame(cm_impacts, index=index)
        df.fillna(0)
        if summarize:
            df.loc[('Construction', 'All', 'All')] = df.sum()
        return df

    def get_total_impacts(self):
        df = pd.concat((
            self.get_ws_impacts(summarize=False),
            self.get_heating_impacts(summarize=False),
            self.get_cooling_impacts(summarize=False),
            self.get_power_impacts(summarize=False),
            self.get_construction_material_impacts(summarize=False)
            ))
        df.loc[('All', 'All', 'All')] = df.sum()
        return df
        
    
    @staticmethod
    def like(system, model_lca):
        '''
        Create an LCA object for the system based on assumptions in model_lca

        Parameters
        ----------
        system : System
            New system for which the LCA will be peformed.
        model_lca : LCA
            Model LCA object.

        Returns
        -------
        new_lca : LCA.
            A new LCA object for the given system.

        '''
        
        new_lca = copy(model_lca)
        new_lca.units = (unit for unit in system.units if unit._if_lca)
        # new_lca.units = sorted(system._lcaunits, key=lambda x: x.line)
        new_lca.system = system
        return new_lca














