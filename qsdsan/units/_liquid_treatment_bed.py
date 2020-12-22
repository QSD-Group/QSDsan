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
from .. import SanUnit, Construction
from ._decay import Decay
from ..utils.loading import load_data, data_path

__all__ = ('LiquidTreatmentBed',)

data_path += 'unit_data/_liquid_treatment_bed.csv'


class LiquidTreatmentBed(SanUnit, Decay):
    '''For secondary treatment of liquid.'''
    
    def __init__(self, ID='', ins=None, outs=(), if_N2O_emission=False, **kwargs):

        '''

        Parameters
        ----------
        ins : WasteStream
            Waste for treatment.
        outs : WasteStream
            Treated waste, fugitive CH4, and fugitive N2O.
        if_N2O_emission : [bool]
            If consider N2O emission from N degradation the process.

        '''

        SanUnit.__init__(self, ID, ins, outs)
        self.if_N2O_emission = if_N2O_emission
    
        data = load_data(path=data_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, '_'+para, value)
        del data
        
        for attr, value in kwargs.items():
            setattr(self, attr, value)
    
    __doc__ += __init__.__doc__
    __init__.__doc__ = __doc__
    
    _N_ins = 1
    _N_outs = 3

    def _run(self):
        waste = self.ins[0]
        treated, CH4, N2O = self.outs
        treated.copy_like(self.ins[0])
        CH4.phase = N2O.phase = 'g'
        
        # COD removal
        COD_loss = self.first_order_decay(k=self.decay_k_COD,
                                          t=self.tau/365,
                                          max_decay=self.COD_max_decay)
        COD_loss_tot = COD_loss*waste.COD/1e3*waste.F_vol
        
        treated._COD *= (1-COD_loss)
        #!!! Which assumption is better?
        treated.imass['OtherSS'] *= (1-COD_loss)
        # treated.mass *= (1-COD_loss)
        CH4_prcd = COD_loss_tot*self.MCF_decay*self.max_CH4_emission
        CH4.imass['CH4'] = CH4_prcd

        if self.if_N2O_emission:
            N_loss = self.first_order_decay(k=self.decay_k_N,
                                            t=self.tau/365,
                                            max_decay=self.N_max_decay)
            N_loss_tot = N_loss*waste.TN/1e3*waste.F_vol
            NH3_rmd, NonNH3_rmd = \
                self.allocate_N_removal(N_loss_tot, waste.imass['NH3'])
            treated.imass ['NH3'] = waste.imass['NH3'] - NH3_rmd
            treated.imass['NonNH3'] = waste.imass['NonNH3'] - NonNH3_rmd
            N2O.imass['N2O'] = N_loss_tot*self.N2O_EF_decay*44/28
        else:
            N2O.empty()
    
    _units = {
        'Residence time': 'd',
        'Bed length': 'm',
        'Bed width': 'm',
        'Bed height': 'm',
        'Single bed volume': 'm3'
        }

    def _design(self):
        design = self.design_results
        design['Residence time'] = self.tau
        design['Bed number'] = N = self.N_bed
        design['Bed length'] = L = self.bed_L
        design['Bed width'] = W = self.bed_W
        design['Bed height'] = H = self.bed_H
        design['Single bed volume'] = L*W*H
        
        concrete = N*self.concrete_thickness*(L*W+2*L*H+2*W*H)
        self.construction = (
            Construction(item='Concrete', quantity=concrete, unit='m3'),
            )
        self.add_construction()


    def _cost(self):
        pass


    @property
    def tau(self):
        '''[float] Residence time, [d].'''
        return self._tau
    @tau.setter
    def tau(self, i):
        self._tau = float(i)

    @property
    def N_bed(self):
        '''[int] Number of treatment beds, float will be converted to the smallest integer.'''
        return self._N_bed
    @N_bed.setter
    def N_bed(self, i):
        self._N_bed = int(np.ceil(i))

    @property
    def bed_L(self):
        '''[float] Bed length, [m].'''
        return self._bed_L
    @bed_L.setter
    def bed_L(self, i):
        self._bed_L = float(i)

    @property
    def bed_W(self):
        '''[float] Bed width, [m].'''
        return self._bed_W
    @bed_W.setter
    def bed_W(self, i):
        self._bed_W = float(i)

    @property
    def bed_H(self):
        '''[float] Bed height, [m].'''
        return self._bed_H
    @bed_H.setter
    def bed_H(self, i):
        self._bed_H = float(i)

    @property
    def concrete_thickness(self):
        '''[float] Thickness of the concrete wall.'''
        return self._concrete_thickness
    @concrete_thickness.setter
    def concrete_thickness(self, i):
        self._concrete_thickness = float(i)







