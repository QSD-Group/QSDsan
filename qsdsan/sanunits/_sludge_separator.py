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

from warnings import warn
from .. import SanUnit
from ._decay import Decay
from ..utils.loading import load_data, data_path

__all__ = ('SludgeSeparator',)

data_path += 'sanunit_data/_sludge_separator.csv'

allocate_N_removal = Decay.allocate_N_removal

class SludgeSeparator(SanUnit):
    '''
    For sludge separation based on Trimmer et al. [1]_, note that no default
    cost or environmental impacts are included.
    
    Parameters
    ----------
    ins : WasteStream
        Waste for treatment.
    outs : WasteStream
        Liquid, settled solids.
    split : float or dict
        Fractions of material retention in the settled solids.
        Default values will be used if not given.
    settled_frac : float
        Fraction of influent that settles as solids.
        The default value will be used if not given.
        
    References
    ----------
    .. [1] Trimmer et al., Navigating Multidimensional Social–Ecological System
        Trade-Offs across Sanitation Alternatives in an Urban Informal Settlement.
        Environ. Sci. Technol. 2020, 54 (19), 12641–12653.
        https://doi.org/10.1021/acs.est.0c03296.
    
    '''
    
    def __init__(self, ID='', ins=None, outs=(), split=None, settled_frac=None):    
        
        SanUnit.__init__(self, ID, ins, outs)
        data = load_data(path=data_path)
        if not split:
            value = eval(data.loc['split']['expected'])
            setattr(self, 'split', value)
        if not settled_frac:
            value = float(data.loc['settled_frac']['expected'])
            setattr(self, 'settled_frac', value)
        del data
    
    _N_ins = 1
    _outs_size_is_fixed = False
    
    def _adjust_solid_water(self, influent, liq, sol):
        sol.imass['H2O'] = 0
        sol.imass['H2O'] = influent.F_mass * self.settled_frac - sol.F_mass
        if sol.imass['H2O'] < 0:
            sol.imass['H2O'] = 0
            msg = 'Negative water content calcualted for settled solids, ' \
                'try smaller split or larger settled_frac.'
            warn(msg, stacklevel=2)
        liq.imass['H2O'] = influent.imass['H2O'] - sol.imass['H2O']
        return liq, sol
        
    def _run(self):
        waste = self.ins[0]
        liq, sol = self.outs[0], self.outs[1]
        
        # Retention in the settled solids
        split = self.split
        if self._split_type == 'float':
            liq.copy_like(waste)
            sol.copy_like(waste)
            sol.mass *= self.split
            liq.mass -= sol.mass
        else:
            for var in self.split.keys():
                #!!! In the future this should be best by changing the state variable
                if var == 'TS':
                    sol.imass['OtherSS'] = split[var] * waste.imass['OtherSS']
                elif var == 'COD':
                    sol_COD = split[var] * waste._COD * waste.F_vol
                    liq_COD = waste._COD * waste.F_vol - sol_COD
                elif var == 'N':
                    N_sol = split[var]*(waste.imass['NH3']+waste.imass['NonNH3'])
                    NonNH3_rmd, NH3_rmd = \
                        allocate_N_removal(N_sol, waste.imass['NonNH3'])
                    sol.imass['NonNH3'] = NonNH3_rmd
                    sol.imass['NH3'] = NH3_rmd
                else:
                    sol.imass[var] = split[var] * waste.imass[var]
            liq.mass = waste.mass - sol.mass
        
        # Adjust total mass of of the settled solids by changing water content.
        liq, sol = self._adjust_solid_water(waste, liq, sol)
        sol._COD = sol_COD / sol.F_vol
        liq._COD = liq_COD / liq.F_vol


    @property
    def split(self):
        '''
        [float] or [dict] Fractions of material retention in the settled solids
        before degradation. If a single number is provided, then it is assumed
        that retentions of all Components in the WasteStream are the same.
        
        Note
        ----
        Set state variable values (e.g., COD) will be retained if the retention
        ratio is a single number (treated like the loss stream is split
        from the original stream), but not when the ratio is a dict.

        '''
        return self._split
    @split.setter
    def split(self, i):
        try:
            self._split = float(i)
            self._split_type = 'float'
        except:
            if isinstance(i, dict):
                self._split = i
                self._split_type = 'dict'
            else:
                raise TypeError(f'Only float or dict allowed, not {type(i).__name__}.')

    @property
    def settled_frac(self):
        '''[float] Fraction of influent that settles as solids.'''
        return self._settled_frac
    @settled_frac.setter
    def settled_frac(self, i):
        self._settled_frac = float(i)























