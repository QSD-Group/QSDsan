#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

# %%

from collections.abc import Iterable
from biosteam._graphics import splitter_graphics
from .. import SanUnit

__all__ = ('ComponentSplitter',)


class ComponentSplitter(SanUnit):
    '''
    Split the influent into individual components,
    the last effluent contains all remaining components.

    Parameters
    ----------
    split_keys : iterable
        IDs of components to be splitted to different effluents.
        Element of the item in the iterable can be str or another iterable
        containing component IDs.
        If the item is also iterable, all components whose ID are in the iterable
        will be splitted to the same effluent.
        The split is always 1 for a certain component to an effluent (i.e., complete split).

        .. note::

            Length of the `split_keys()` (which determines size of the outs) \
            cannot be changed after initiation.

    Examples
    --------
    `bwaise systems <https://github.com/QSD-Group/EXPOsan/blob/main/exposan/bwaise/systems.py>`_
    '''

    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', split_keys=()):
        if not split_keys:
            raise ValueError('`split_keys` cannot be empty.')

        if isinstance(split_keys, str):
            self._N_outs = 2
        else:
            self._N_outs = len(split_keys) + 1
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)

        self._split_keys = split_keys


    _ins_size_is_fixed = False
    _outs_size_is_fixed = False
    _graphics = splitter_graphics


    def _run(self):
        last = self.outs[-1]
        last.mix_from(self.ins)

        splitted = []
        for num, cmps in enumerate(self.split_keys):
            if isinstance(cmps, str):
                cmps = (cmps,)

            elif not isinstance(cmps, Iterable):
                raise ValueError('`split_keys` must be an iterable, '
                                 f'not {type(cmps).__name__}.')

            for cmp in cmps:
                self.outs[num].imass[cmp] = last.imass[cmp]
                last.imass[cmp] = 0
                if cmp in splitted:
                    raise ValueError(f'The component {cmps} appears more than once in `split_keys`.')
                splitted.append(cmp)


    @property
    def split_keys(self):
        '''
        [iterable] IDs of components to be splitted to different effluents.
        Element of the item in the iterable can be str or another iterable
        containing component IDs.
        If the item is also iterable, all components whose ID are in the iterable
        will be splitted to the same effluent.
        The split is always 1 for a certain component to an effluent (i.e., complete split).

        .. note::

            Length of the `split_keys()` (which determines size of the outs) \
                cannot be changed after initiation.
        '''
        return self._split_keys
    @split_keys.setter
    def split_keys(self, i):
        if isinstance(i, str):
            i = (i,)

        if len(i) != len(self.outs):
            raise ValueError('Size of `split_keys` cannot be changed after initiation.')

        self._split_keys = i