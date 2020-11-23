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

Ref:
    [1] Trimmer et al., Navigating Multidimensional Social–Ecological System
        Trade-Offs across Sanitation Alternatives in an Urban Informal Settlement.
        Environ. Sci. Technol. 2020, 54 (19), 12641–12653.
        https://doi.org/10.1021/acs.est.0c03296.

'''



# %%

import pandas as pd
from .utils.loading import data_path
data_path += 'lca_data/materials.xlsx'

class Material:
    '''
    A class containing different types of materials.
    
    '''
    
    __slots__ = ('_ID', '_functional_unit', '_price', '_CFs')
    
    def __new__(cls, ID, functional_unit='kg', price=0., CFs={}):
        self = super().__new__(cls)
        self._ID = ID
        self._functional_unit = functional_unit
        self._price = price
        self._CFs = CFs
        return self
    
    # This makes sure it won't be shown as memory location of the object
    def __repr__(self):
        return f"Material('{self}')"
    
    def show(self):
        info = f'Material: {self.ID} [per {self.functional_unit}]'
        info += f'\n price: {self.price}'
        info += '\n CFs:'
        CFs = self.CFs
        if len(CFs) == 0:
            info += ' None'
        else:
            for cat in self.CFs.keys():
                #!!! Need to add units!
                info += f'\n     {cat}: {CFs[cat]} [unit]'
        print(info)    
    
    _ipython_display_ = show

    
    #!!! Are the values GWP100 from ref [1]?
    @classmethod
    def load_default_materials(cls, path=data_path):
        # TODO: add prices and other LCA categories
        data_file = pd.ExcelFile(data_path)
        materials = {}
        for cat in data_file.sheet_names:
            data = data_file.parse(cat, index_col=0)
            for material in data.index:
                if material not in materials.keys():
                    new = cls.__new__(cls, ID=material)
                    materials[material] = new
                    materials[material].CFs[cat] = float(data.loc[material]['expected'])
                    # materials[material] = new
                else:
                    old = materials[material]
                    old.CFs[cat] = float(data.loc[material]['expected'])
        return materials
    
    
    #!!! Temporary, units will be taken care of separately
    supported_units = {
        'mass': ('kg',),
        'area': ('m2',),
        'volume': ('m3',)
        }
    
    @property
    def ID(self):
        '''Material ID.'''    
        return self._ID
    
    @property
    def functional_unit(self):
        return self._functional_unit
    @functional_unit.setter
    def functional_unit(self, unit):
        assert unit in sum((i for i in self.supported_units.values()), ()), \
            f'{unit} not supported, see Material.supported_units'
    
    @property
    def price(self):
        '''[float] Price of the material per functional unit.'''
        return self._price
    @price.setter
    def price(self, i):
        self._price = float(i)
    
    @property
    def CFs(self):
        '''[dict] Characterization factors of the material.'''
        return self._CFs
    @CFs.setter
    def CFs(self, i):
        self._CFs = i














    
    
    
    
class ConstructionMaterial:
    '''
    A class to calculate the cost and environmental impacts associated with construction materials.
    
    '''

    __slots__ = ('quantity',)
    
    def __init__(self):
        self.quantity = 0.












    
    
    
    