#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:

    Yalin Li <mailto.yalin.li@gmail.com>

    Joy Zhang <joycheung1994@gmail.com>

    Jianan Feng <jiananf2@illinois.edu>

Part of this module is based on the biosteam package:
https://github.com/BioSTEAMDevelopmentGroup/biosteam

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

import thermosteam as tmo
from warnings import warn
from collections.abc import Iterable
from ... import SanUnit
from ..bst import Splitter

__all__ = (
    'PhaseChanger',
    'ComponentSplitter',
    )


# %%

class PhaseChanger(SanUnit):
    '''
    Change the effluent phase to the desired one; can also bridge between
    stream types (e.g., :class:`~.WasteStream` to plain
    :class:`thermosteam.Stream`) via the ``init_with`` argument.

    The outlet stream class is selected by ``init_with``. The inlet is
    copied into the outlet (creating the type conversion if those classes
    differ), then the outlet's ``phase`` is set to the requested value.

    Parameters
    ----------
    ins : Iterable(stream)
        Influent.
    outs : Iterable(stream)
        Effluent.
    init_with : str
        Stream class for the outlet: ``'WasteStream'`` (default),
        ``'SanStream'``, or ``'Stream'``. Use ``'Stream'`` to bridge from a
        ``WasteStream`` feed to a plain ``thermosteam.Stream`` outlet -- for
        example, as the inlet of a :class:`~.HeatExchangerNetwork` branch,
        which cannot consume ``WasteStream``.
    phase : str
        Desired outlet phase; one of ``'g'``, ``'l'``, or ``'s'``.

    Examples
    --------
    Setting the outlet phase only (outlet remains a ``WasteStream``):

    >>> import qsdsan as qs
    >>> from qsdsan.utils import create_example_components
    >>> qs.set_thermo(create_example_components())
    >>> ws = qs.WasteStream('pc_feed', Water=1000, Methanol=10, units='kg/hr')
    >>> pc = qs.unit_operations.PhaseChanger('PC1', ins=ws, phase='g')
    >>> pc.simulate()
    >>> pc.outs[0].phase
    'g'
    >>> type(pc.outs[0]).__name__
    'WasteStream'

    Bridging from ``WasteStream`` to plain ``Stream`` (e.g., to feed a
    :class:`~.HeatExchangerNetwork` branch that cannot accept
    ``WasteStream``):

    >>> ws2 = qs.WasteStream('pc_feed2', Water=1000, units='kg/hr')
    >>> pc2 = qs.unit_operations.PhaseChanger(
    ...     'PC2', ins=ws2, init_with='Stream', phase='l',
    ... )
    >>> pc2.simulate()
    >>> type(pc2.outs[0]).__name__
    'Stream'
    '''
    _N_ins = 1
    _N_outs = 1
    _ins_size_is_fixed = False
    _outs_size_is_fixed = False

    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', phase='l'):

        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.phase = phase


    def _run(self):
        influent = self.ins[0]
        effluent = self.outs[0]
        if isinstance(influent, tmo.MultiStream): # issue a warning
            warn(f'MultiStream {influent.ID} converted to Stream in {self.ID}')
            influent.as_stream()
        effluent.copy_like(influent)
        effluent.phase = self.phase

    @property
    def phase(self):
        return self._phase
    @phase.setter
    def phase(self, i):
        if not i in ('g', 'l', 's'):
            raise ValueError('`phase` must be one of ("g", "l", or "s"), '
                             f'not "{i}".')
        self._phase = i

# %%

class ComponentSplitter(SanUnit):
    '''
    Split the influent into individual components,
    the last effluent contains all remaining components.

    Parameters
    ----------
    split_keys : iterable
        IDs of components to be split to different effluents.
        Element of the item in the iterable can be str or another iterable
        containing component IDs.
        If the item is also iterable, all components whose ID are in the iterable
        will be split to the same effluent.
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
    _graphics = Splitter._graphics


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
        [iterable] IDs of components to be split to different effluents.
        Element of the item in the iterable can be str or another iterable
        containing component IDs.
        If the item is also iterable, all components whose ID are in the iterable
        will be split to the same effluent.
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
