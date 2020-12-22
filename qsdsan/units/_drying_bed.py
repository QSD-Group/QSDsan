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

Ref:
    [1] Trimmer et al., Navigating Multidimensional Social–Ecological System
        Trade-Offs across Sanitation Alternatives in an Urban Informal Settlement.
        Environ. Sci. Technol. 2020, 54 (19), 12641–12653.
        https://doi.org/10.1021/acs.est.0c03296.

'''


# %%

import numpy as np
from warnings import warn
from .. import SanUnit, Construction
from ._decay import Decay
from ..utils.loading import load_data, data_path

__all__ = ('DryingBed',)

data_path += 'unit_data/_drying_bed.csv'


class DryingBed(SanUnit, Decay):
    '''Unplanted and planted drying bed for solids.'''
    
    def __init__(self, ID='', ins=None, outs=(), design_type='unplanted',
                 **kwargs):
        
        '''

        Parameters
        ----------
        ins : WasteStream
            Solid for drying.
        outs : WasteStream
            Dried solids, evaporated water, fugitive CH4, and fugitive N2O.
        design_type : [str]
            Can be 'unplanted' or 'planted'. The default 'unplanted' process has
            a number of 'covered', 'uncovered', and 'storage' beds. The 'storage'
            bed is similar to the 'covered' bed, but with higher wall height.

        '''        
        
        SanUnit.__init__(self, ID, ins, outs)
        N_unplanted = {'covered': 19,
                       'uncovered': 30,
                       'storage': 19,
                       'planted': 0}
        if design_type == 'unplanted':
            self._N_bed = N_unplanted
            self.design_type = 'unplanted'
        else:
            self._N_bed = dict.fromkeys(N_unplanted.keys(), 0)
            self._N_bed['planted'] = 2
            self.design_type = 'planted'
            
        data = load_data(path=data_path)
        for para in data.index:
            if para == 'N_bed': continue
            if para in ('sol_frac', 'bed_L', 'bed_W', 'bed_H'):
                value = eval(data.loc[para]['expected'])
            else:
                value = float(data.loc[para]['expected'])
            setattr(self, '_'+para, value)
        del data
        
        for attr, value in kwargs.items():
            setattr(self, attr, value)

    __doc__ += __init__.__doc__
    __init__.__doc__ = __doc__
        
    _N_ins = 1
    _N_outs = 4
    
    def _run(self):
        waste = self.ins[0]
        sol, evaporated, CH4, N2O = self.outs
        sol.copy_like(waste)
        evaporated.phase = CH4.phase = N2O.phase = 'g'
        
        # COD degradation in settled solids
        COD_loss = self.first_order_decay(k=self.decay_k_COD,
                                          t=self.tau/365,
                                          max_decay=self.COD_max_decay)

        sol._COD *= 1 - COD_loss
        sol.imass['OtherSS'] *= 1 - COD_loss
        CH4.imass['CH4'] = waste.COD/1e3*waste.F_vol*COD_loss * \
            self.max_CH4_emission*self.MCF_decay # COD in mg/L (g/m3)

        # N degradation
        N_loss = self.first_order_decay(k=self.decay_k_N,
                                        t=self.tau/365,
                                        max_decay=self.N_max_decay)
        N_loss_tot = N_loss*waste.TN/1e3*waste.F_vol
        NH3_rmd, NonNH3_rmd = \
            self.allocate_N_removal(N_loss_tot, sol.imass['NH3'])
        sol.imass ['NH3'] -=  NH3_rmd
        sol.imass['NonNH3'] -= NonNH3_rmd
        N2O.imass['N2O'] = N_loss_tot*self.N2O_EF_decay*44/28

        # Adjust water content in the dried solids
        set_sol_frac = self.sol_frac[self.design_type]
        solid_content = 1 - sol.imass['H2O']/sol.F_mass
        if solid_content > set_sol_frac:
            msg = f'Solid content of the solid after COD removal is {solid_content},'\
                'larger than the set sol_frac of {set_sol_frac} for the {self.design_type} ' \
                'process type, the set value is ignored.'
            warn(msg, source=self)
            evaporated.empty()
        else:
            sol.imass['H2O'] = (sol.F_mass-sol.imass['H2O'])/set_sol_frac
            evaporated.imass['H2O'] = waste.imass['H2O'] - sol.imass['H2O']

    _units = {
        'Single covered bed volume': 'm3',
        'Single uncovered bed volume': 'm3',
        'Single storage bed volume': 'm3',
        'Single planted bed volume': 'm3',
        'Total cover area': 'm2',
        'Total column length': 'm'
        }

    def _design(self):
        design = self.design_results
        # Bed
        L = np.fromiter(self.bed_L.values(), dtype=float)
        W = np.fromiter(self.bed_W.values(), dtype=float)
        H = np.fromiter(self.bed_H.values(), dtype=float)
        V = L * W * H
        N_bed = self.N_bed
        N = np.fromiter(N_bed.values(), dtype=int)
        for n, i in enumerate(N_bed.keys()):
            design[f'Number of {i} bed'] = N_bed[i]
            design[f'Single {i} bed volume'] = V[n]
        cover_array = np.array((1, 0, 1, 0)) # covered, uncovered, storage, planted
        design['Total cover area'] = tot_cover_area = \
            (cover_array*N*L*W/(np.cos(self.cover_slope/180*np.pi))).sum()
        design['Total column length'] = tot_column_length = \
            (cover_array*N*2*self.column_per_side*self.column_H).sum()
        
        concrete = (N*self.concrete_thickness*(L*W+2*L*H+2*W*H)).sum()
        steel = tot_cover_area*self.cover_unit_mass + \
            tot_column_length*self.column_unit_mass
        self.construction = (
            Construction(item='Concrete', quantity=concrete, unit='m3'),
            Construction(item='Steel', quantity=steel, unit='kg'),
            )
        for i in self.construction:
            self._BM[i.item.ID] = 1
        
    def _cost(self):
        pass


    @property
    def tau(self):
        '''[float] Retention time, [d].'''
        return self._tau
    @tau.setter
    def tau(self, i):
        self._tau = float(i)

    @property
    def sol_frac(self):
        '''[float] Final solid content of the dried solids.'''
        return self._sol_frac
    @sol_frac.setter
    def sol_frac(self, i):
        self._sol_frac = float(i)

    @property
    def design_type(self):
        '''[str] Drying bed type, can be either 'unplanted' or 'planted'.'''
        return self._design_type
    @design_type.setter
    def design_type(self, i):
        if i in ('unplanted', 'planted'):
            self._design_type = i
            self.line =f'{i.capitalize()} drying bed'
        else:
            if i is None:
                self.line = 'Drying bed'
            else:
                raise ValueError(f"design_type can only be 'unplanted', 'planted', "
                                 f"or 'None', not {i}.")

    @staticmethod
    def _set_bed_prop(prop, value):
        if isinstance(value, dict):
            for key in value.keys():
                prop[key] = value[key]
        else:
            raise TypeError(f'Only dict object is allowed, not {type(value).__name__}.')

    @property
    def N_bed(self):
        '''
        [dict] Number of the different types of drying beds,
        float will be converted to the smallest integer.
        '''
        #!!! Think of a better way to do this
        for i, j in self._N_bed.items():
            self._N_bed[i] = int(np.ceil(j))
        return self._N_bed
    @N_bed.setter
    def N_bed(self, i):
        self._set_bed_prop(self._N_bed, i)

    @property
    def bed_L(self):
        '''[dict] Length of the different types of drying beds, [m].'''
        return self._bed_L
    @bed_L.setter
    def bed_L(self, i):
        self._set_bed_prop(self._bed_L, i)

    @property
    def bed_W(self):
        '''[dict] Width of the different types of drying beds, [m].'''
        return self._bed_W
    @bed_W.setter
    def bed_W(self, i):
        self._set_bed_prop(self._bed_W, i)

    @property
    def bed_H(self):
        '''[dict] Wall height of the different types of drying beds, [m].'''
        return self._bed_H
    @bed_H.setter
    def bed_H(self, i):
        self._set_bed_prop(self._bed_H, i)

    @property
    def column_H(self):
        '''[float] Column height for covered bed, [m].'''
        return self._column_H
    @column_H.setter
    def column_H(self, i):
        self._column_H = float(i)

    @property
    def column_per_side(self):
        '''[int] Number of columns per side of covered bed, float will be converted to the smallest integer.'''
        return self._column_per_side
    @column_per_side.setter
    def column_per_side(self, i):
        self._column_per_side = int(np.ceil(i))

    @property
    def column_unit_mass(self):
        '''[float] Unit mass of the column, [kg/m].'''
        return self._column_unit_mass
    @column_unit_mass.setter
    def column_unit_mass(self, i):
        self._column_unit_mass = float(i)

    @property
    def concrete_thickness(self):
        '''[float] Thickness of the concrete wall.'''
        return self._concrete_thickness
    @concrete_thickness.setter
    def concrete_thickness(self, i):
        self._concrete_thickness = float(i)

    @property
    def cover_slope(self):
        '''[float] Slope of the bed cover, [°].'''
        return self._cover_slope
    @cover_slope.setter
    def cover_slope(self, i):
        self._cover_slope = float(i)

    @property
    def cover_unit_mass(self):
        '''[float] Unit mass of the bed cover, [kg/m2].'''
        return self._cover_unit_mass
    @cover_unit_mass.setter
    def cover_unit_mass(self, i):
        self._cover_unit_mass = float(i)








