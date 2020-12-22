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

from warnings import warn
from .. import currency, SanUnit, Transportation
from .._units_of_measure import auom

__all__ = ('Trucking',)


class Trucking(SanUnit):
    '''For transportation of materials with considerations on material loss.'''
    
    def __init__(self, ID='', ins=None, outs=(),
                 load_type='mass', load=1., load_unit='tonne',
                 distance=1., distance_unit='km',
                 interval=1., interval_unit='day',
                 fee=0., fee_unit=currency,
                 if_material_loss=True, loss_ratio=0.02):
        '''

        Parameters
        ----------
        load_type : [str]
            Either 'mass' or 'volume'.
        load : [float]
            Transportation load per trip.
        load_unit : [str]
            Unit of the load.
        distance : [float]
            Transportation distance per trip.
        distance_unit : [float]
            Unit of the distance.
        interval : [float]
            Time interval between trips.
        fee : [float]
            Transportation fee per trip.
        fee : [float]
            Transportation fee per trip.
        if_material_loss : [bool]
            If material loss occurs during transportation.
        loss_ratio : [float] or [dict]
            Fractions of material losses due to transportation.

        '''
        SanUnit.__init__(self, ID, ins, outs)
        self.single_truck = \
            Transportation(item='Trucking',
                           load_type=load_type, load=load, load_unit=load_unit,
                           distance=distance, distance_unit=distance_unit,
                           interval=interval, interval_unit=interval_unit)
        self.total_truck = self.single_truck.copy()
        self._update_fee(fee, fee_unit)
        self.if_material_loss = if_material_loss
        self.loss_ratio = loss_ratio

    __doc__ += __init__.__doc__
    __init__.__doc__ = __doc__
    
    _N_ins = 1
    _N_outs = 2

    def _update_fee(self, fee=0., unit=''):
        if not unit or unit == currency:
            self._fee = float(fee)
        else:
            converted = auom(unit).convert(float(fee), currency)
            self._fee = converted

    def _run(self):
        transported, loss = self.outs
        transported.copy_like(self.ins[0])
        loss.empty()
        if self.if_material_loss:
            if self._loss_ratio_type == 'float':
                loss.copy_like(transported)
                transported.mass *= 1 - self.loss_ratio
                loss.mass = self.ins[0].mass - transported.mass
            else:
                for cmp, ratio in self.loss_ratio.items():                        
                    transported.imass[cmp] *= 1 - ratio
                    loss.imass[cmp] = self.ins[0].imass[cmp] - transported.imass[cmp]


    def _design(self):        
        self.transportation = (trans,) = (self.total_truck,)
        if trans.load_type == 'volume':        
            trans._update_value('load', self.F_vol_in*trans.interval) # interval in hr
        else:
            trans.load_type == 'mass'
            trans._update_value('load', self.F_mass_in*trans.interval)
        self.design_results['Paralle trucks'] = N = \
            self.total_truck.load / self.single_truck.load
        self._add_OPEX = N*self.fee/trans.interval/24 # add_OPEX in hr

    @property
    def fee(self):
        '''[float] Transportation fee per trip.'''
        return self._fee
    @fee.setter
    def fee(self, fee, unit=''):
        self._update_fee(fee, unit)

    @property
    def loss_ratio(self):
        '''
        [float] or [dict] Fractions of material losses due to transportation.
        If a single number is provided, then it is assumed that losses of
        all Components in the WasteStream are the same.
        
        Note
        ----
            Set state variable values (e.g., COD) will be retained if the loss
            ratio is a single number (treated like the loss stream is split
            from the original stream), but not when the ratio is a dict.

        '''
        return self._loss_ratio
    @loss_ratio.setter
    def loss_ratio(self, i):
        if not self.if_material_loss:
            msg = f'if_material_loss is False, the set value {i} is ignored.'
            warn(msg, source=self)
        try:
            self._loss_ratio = float(i)
            self._loss_ratio_type = 'float'
        except:
            if isinstance(i, dict):
                self._loss_ratio = i
                self._loss_ratio_type = 'dict'
            else:
                raise TypeError(f'Only float or dict allowed, not {type(i).__name__}.')












