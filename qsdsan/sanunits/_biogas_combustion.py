#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/master/LICENSE.txt
for license details.
'''


# %%

from .. import SanUnit

__all__ = ('BiogasCombustion',)


class BiogasCombustion(SanUnit):
    '''
    Use of biogas in combustion.

    Parameters
    ----------
    if_combustion : bool
        If include combusion reaction during simulation.
    biogas_loss : float
        Fraction of biogas loss.
    biogas_eff : float
        Combustion efficiency of biogas as a fraction of CH4.
    
    '''
    
    def __init__(self, ID='', ins=None, outs=(), if_combustion=False,
                 biogas_loss=0.1, biogas_eff=0.55):

        SanUnit.__init__(self, ID, ins, outs)
        self.if_combustion = if_combustion
        self._biogas_loss = biogas_loss
        self._biogas_eff = biogas_eff
    
    _N_ins = 2
    _N_outs = 3

    def _run(self):
        biogas, air = self.ins
        biogas.phase = air.phase = 'g'
        for i in self.outs:
            i.copy_like(biogas)
        used, lost, wasted = self.outs
        lost.mass *= self.biogas_loss
        used.mass -= lost.mass
        wasted.mass = used.mass * (1-self.biogas_eff)
        used.mass -= wasted.mass
        if self.if_combustion:
            rxn = self.chemicals.get_combustion_reactions()
            rxn.force_reaction(used.mol)
            air.imol['O2'] = -used.imol['O2']
            used.imol['O2'] = 0.
            air.imol['N2'] = 0.79/0.21 * air.imol['O2']
        else:
            air.empty()
        
        
    @property
    def biogas_loss(self):
        '''[float] Fraction of biogas loss.'''
        return self._biogas_loss
    @biogas_loss.setter
    def biogas_loss(self, i):
        self._biogas_loss = float(i)

    @property
    def biogas_eff(self):
        '''[float] Combustion efficiency of biogas as a fraction of CH4.'''
        return self._biogas_eff
    @biogas_eff.setter
    def biogas_eff(self, i):
        self._biogas_eff = float(i)










