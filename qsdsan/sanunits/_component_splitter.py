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
'''

# %%

from biosteam._graphics import splitter_graphics
from .. import SanUnit

__all__ = ('ComponentSplitter',)


class ComponentSplitter(SanUnit):
    '''
    Split the influent into individual ``Component`` objects,
    the last effluent contains all remaining ``Component`` objects.
    '''
    
    def __init__(self, ID='', ins=None, outs=(), split_keys=()):
        
        SanUnit.__init__(self, ID, ins, outs)
        self.split_keys = split_keys
    
    _ins_size_is_fixed = False
    _outs_size_is_fixed = False
    _graphics = splitter_graphics

    def _run(self):
        last = self.outs[-1]
        last.mix_from(self.ins)
        splitted = []
        num = 0
        for cmps in self.split_keys:
            if isinstance(cmps, str):
                self.outs[num].imass[cmps] = last.imass[cmps]
                last.imass[cmps] = 0
                if cmps in splitted:
                    raise ValueError(f'The Component {cmps} appears more than once in split_dict')
                splitted.append(cmps)
            else:
                try: iter(cmps)
                except:
                    raise ValueError('Elements of the split must be str or iterable, '
                                     f'not {type(cmps).__name__}.')
                for cmp in cmps:
                    self.outs[num].imass[cmp] = last.imass[cmp]
                    last.imass[cmp] = 0
                    if cmp in splitted:
                        raise ValueError(f'The Component {cmps} appears more than once in split_dict')
                    splitted.append(cmp)

            num += 1

    @property
    def split_keys(self):
        '''
        [iterable] An iterable containing IDs of Components to be splitted to
        different effluents. Element of the item in the iterable can be str of
        another iterable containing Component IDs. If the item is also iterable,
        all Components whose ID are in the iterable will be splitted to the same
        effluent. Note that the split is 1 (i.e., all of the Component will be
        diverted to the effluent).
        '''
        return self._split_keys
    @split_keys.setter
    def split_keys(self, i):        
        self._split_keys = i











