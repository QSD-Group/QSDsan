#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:

    Jianan Feng <jiananf2@illinois.edu>

    Yalin Li <mailto.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

import biosteam as bst, qsdsan as qs

__all__ = (
    'BinaryDistillation',
    'ShortcutColumn',
    'MESHDistillation',
    'AdiabaticMultiStageVLEColumn',
    )

_lb_to_kg = qs.utils.auom('lb').conversion_factor('kg')

# SanUnit add-on kwargs siphoned by these wrappers before forwarding to the
# BioSTEAM ``__init__``; the BST ``__init__`` does not accept them and the
# wrapper's ``_init_sanunit_addons`` consumes them after.
_SANUNIT_ADDON_KEYS = (
    'include_construction', 'construction', 'transportation', 'equipment',
    'add_OPEX', 'lifetime', 'F_BM_default', 'isdynamic', 'exogenous_vars',
)


def _split_sanunit_addons(kwargs):
    return {k: kwargs.pop(k) for k in tuple(_SANUNIT_ADDON_KEYS) if k in kwargs}


class BinaryDistillation(bst.units.BinaryDistillation, qs.SanUnit):
    '''
    Similar to biosteam.units.BinaryDistillation, but can include construction impact calculation.

    See Also
    --------
    `biosteam.units.BinaryDistillation <https://biosteam.readthedocs.io/en/latest/API/units/distillation.html>`_
    '''

    include_construction = True

    def __init__(self, *args, **kwargs):
        addons = _split_sanunit_addons(kwargs)
        super().__init__(*args, **kwargs)
        self._init_sanunit_addons(**addons)

    def _design(self):
        super()._design()
        D = self.design_results
        if self.include_construction:
            construction = getattr(self, 'construction', [])
            if construction: construction[0].quantity = (D['Rectifier weight'] + D['Stripper weight'])*_lb_to_kg
            else:
                self.construction = [
                    qs.Construction('carbon_steel', linked_unit=self, item='Carbon_steel',
                                    quantity=(D['Rectifier weight'] + D['Stripper weight'])*_lb_to_kg, quantity_unit='kg'),
                    ]




class ShortcutColumn(bst.units.ShortcutColumn, qs.SanUnit):
    '''
    biosteam.units.ShortcutColumn with QSDsan properties.

    See Also
    --------
    `biosteam.units.ShortcutColumn <https://biosteam.readthedocs.io/en/latest/API/units/distillation.html>`_
    '''

    def __init__(self, *args, **kwargs):
        addons = _split_sanunit_addons(kwargs)
        super().__init__(*args, **kwargs)
        self._init_sanunit_addons(**addons)


class MESHDistillation(bst.units.MESHDistillation, qs.SanUnit):
    '''
    biosteam.units.MESHDistillation with QSDsan properties.

    See Also
    --------
    `biosteam.units.MESHDistillation <https://biosteam.readthedocs.io/en/latest/API/units/distillation.html>`_
    '''

    def __init__(self, *args, **kwargs):
        addons = _split_sanunit_addons(kwargs)
        super().__init__(*args, **kwargs)
        self._init_sanunit_addons(**addons)


class AdiabaticMultiStageVLEColumn(bst.units.AdiabaticMultiStageVLEColumn, qs.SanUnit):
    '''
    biosteam.units.AdiabaticMultiStageVLEColumn with QSDsan properties.

    See Also
    --------
    `biosteam.units.AdiabaticMultiStageVLEColumn <https://biosteam.readthedocs.io/en/latest/API/units/distillation.html>`_
    '''

    def __init__(self, *args, **kwargs):
        addons = _split_sanunit_addons(kwargs)
        super().__init__(*args, **kwargs)
        self._init_sanunit_addons(**addons)