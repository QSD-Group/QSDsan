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

TODO:
    [1] Incorporate ADM, or change this to simple_AD or something

'''


# %%


from .. import SanUnit

__all__ = ('AnaerobicDigestion',)


class AnaerobicDigestion(SanUnit):
    '''Anaerobic digestion of wastes with production of biogas.'''
    
    def __init__(self, ID='', ins=None, outs=(), if_combustion=False,
                 biogas_loss=0.1, biogas_eff=0.55):
        '''

        Parameters
        ----------
        if_combustion : [bool]
            If include combusion reaction during simulation.
        biogas_loss : [float]
            Fraction of biogas loss.
        biogas_eff : [float]
            Combustion efficiency of biogas as a fraction of CH4.

        '''
        SanUnit.__init__(self, ID, ins, outs)
        self.if_combustion = if_combustion
        self._biogas_loss = biogas_loss
        self._biogas_eff = biogas_eff
    
    _N_ins = 2
    _N_outs = 3