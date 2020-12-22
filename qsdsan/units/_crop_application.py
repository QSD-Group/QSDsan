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
from .. import SanUnit

__all__ = ('CropApplication',)

class CropApplication(SanUnit):
    '''Recovery nutrients in the recycled excreta (energy not recovered).'''
    
    def __init__(self, ID='', ins=None, outs=(),
                 if_material_loss=True, loss_ratio=0.02):
        '''

        Parameters
        ----------
        if_material_loss : [bool] or [dict]
            If material loss occurs during application.
        loss_ratio : [float] or [dict]
            Fractions of material losses during application (if if_materiloass is True).

        '''
        SanUnit.__init__(self, ID, ins, outs)
        self.if_material_loss = if_material_loss
        self.loss_ratio = loss_ratio

    __doc__ += __init__.__doc__
    __init__.__doc__ = __doc__
    
    _N_ins = 1
    _N_outs = 2

    def _run(self):
        applied, loss = self.outs
        applied.copy_like(self.ins[0])
        loss.empty()
        if self.if_material_loss:
            if self._loss_ratio_type == 'float':
                loss.copy_like(applied)
                applied.mass *= 1 - self.loss_ratio
                loss.mass = self.ins[0].mass - applied.mass
            else:
                for cmp, ratio in self.loss_ratio.items():
                    applied.imass[cmp] *= 1 - ratio
                    loss.imass[cmp] = self.ins[0].imass[cmp] - applied.imass[cmp]


    @property
    def loss_ratio(self):
        '''
        [float] or [dict] Fractions of material losses during application.
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
        else:
            try:
                self._loss_ratio = float(i)
                self._loss_ratio_type = 'float'
            except:
                if isinstance(i, dict):
                    self._loss_ratio = i
                    self._loss_ratio_type = 'dict'
                else:
                    raise TypeError(f'Only float or dict allowed, not {type(i).__name__}.')



















