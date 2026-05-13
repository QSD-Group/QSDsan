#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:

    Yalin Li <mailto.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

import biosteam as bst
from biosteam._flowsheet import Flowsheet as _BstFlowsheet, MainFlowsheet as _BstMainFlowsheet
from thermosteam.utils import Registry

__all__ = ('SanFlowsheet', 'SanMainFlowsheet')


class SanFlowsheet(_BstFlowsheet):
    '''
    A Flowsheet subclass that also holds per-flowsheet registries for
    ImpactItem, Construction, and Transportation objects.

    ImpactIndicator objects remain global and are shared across all flowsheets.
    '''

    def __new__(cls, ID=None):
        self = super().__new__(cls, ID)
        object.__setattr__(self, 'item',           Registry())
        object.__setattr__(self, 'construction',   Registry())
        object.__setattr__(self, 'transportation', Registry())
        return self

    @classmethod
    def from_registries(cls, ID, stream, unit, system,
                        item=None, construction=None, transportation=None):
        self = _BstFlowsheet.from_registries(ID, stream, unit, system)
        object.__setattr__(self, '__class__', cls)
        object.__setattr__(self, 'item',           item          if item          is not None else Registry())
        object.__setattr__(self, 'construction',   construction  if construction  is not None else Registry())
        object.__setattr__(self, 'transportation', transportation if transportation is not None else Registry())
        return self

    @property
    def registries(self):
        return (self.stream, self.unit, self.system,
                self.item, self.construction, self.transportation)

    def clear(self, reset_ticket_numbers=True):
        # Detach StreamImpactItems from their streams
        for lca_item in list(self.item.data.values()):
            linked = getattr(lca_item, '_linked_stream', None)
            if linked is not None:
                linked._stream_impact_item = None
                object.__setattr__(lca_item, '_linked_stream', None)

        # Detach Construction objects from their units
        for c in list(self.construction.data.values()):
            if getattr(c, '_linked_unit', None) is not None:
                c._linked_unit._construction = []

        # Detach Transportation objects from their units
        for t in list(self.transportation.data.values()):
            if getattr(t, '_linked_unit', None) is not None:
                t._linked_unit._transportation = []

        super().clear(reset_ticket_numbers)
        self.item.clear()
        self.construction.clear()
        self.transportation.clear()


class SanMainFlowsheet(SanFlowsheet, _BstMainFlowsheet):
    '''
    MainFlowsheet subclass that also swaps ImpactItem, Construction, and
    Transportation registries when set_flowsheet() is called.
    '''

    __slots__ = ()

    def set_flowsheet(self, flowsheet, new=False):
        if isinstance(flowsheet, str):
            if not new and flowsheet in self.flowsheet:
                flowsheet = self.flowsheet[flowsheet]
            else:
                flowsheet = SanFlowsheet(flowsheet)
        elif not isinstance(flowsheet, SanFlowsheet):
            # Upgrade a plain _BstFlowsheet to SanFlowsheet so we always
            # have item/construction/transportation registries available.
            object.__setattr__(flowsheet, '__class__', SanFlowsheet)
            if not hasattr(flowsheet, 'item'):
                object.__setattr__(flowsheet, 'item',           Registry())
            if not hasattr(flowsheet, 'construction'):
                object.__setattr__(flowsheet, 'construction',   Registry())
            if not hasattr(flowsheet, 'transportation'):
                object.__setattr__(flowsheet, 'transportation', Registry())

        _BstMainFlowsheet.set_flowsheet(self, flowsheet)

        # Swap LCA registries — lazy import to avoid circular imports
        dct = self.__dict__
        if 'item' not in dct:
            return
        from . import ImpactItem, Construction, Transportation
        ImpactItem.registry     = dct['item']
        Construction.registry   = dct['construction']
        Transportation.registry = dct['transportation']
