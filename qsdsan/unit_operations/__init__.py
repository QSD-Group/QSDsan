#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:

    Yalin Li <mailto.yalin.li@gmail.com>

    Joy Zhang <joycheung1994@gmail.com>

    Lewis Rowles <stetsonsc@gmail.com>

    Hannah Lohman <hlohman94@gmail.com>

    Tori Morgan <vlmorgan@illinois.edu>

    Shion Watabe <swatabe2@illinois.edu>

    Lane To <lane20@illinois.edu>

    Smiti Mittal <smitimittal@gmail.com>

    Anna Kogler <akogler@stanford.edu>

    Jianan Feng <jiananf2@illinois.edu>

    Saumitra Rai <raisaumitra9@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''
# %%
from numba import njit
@njit(cache=True)
def dydt_cstr(QC_ins, QC, V, _dstate):
    Q_ins = QC_ins[:, -1]
    C_ins = QC_ins[:, :-1]
    _dstate[-1] = 0
    _dstate[:-1] = (Q_ins @ C_ins - sum(Q_ins)*QC[:-1])/V
    
#%%
# **NOTE** PLEASE ORDER THE MODULES ALPHABETICALLY #

# Units that do not rely on other units
from ._bst_abstract import *
from ._dynamic_abstract import *
from ._static_abstract import *
from ._static_combustion import *
from ._bst_compressor import *
from ._static_crop_application import *
from ._dynamic_influent import *
from ._static_electrochemical_cell import *
from ._static_excretion import *
from ._bst_facilities import *
from ._bst_heat_exchanging import *
from ._dynamic_junction import *
from ._static_membrane_gas_extraction import *
from ._static_metal_dosage import *
from ._static_non_reactive import *
from ._bst_pumping import *
from ._dynamic_pumping import *
from ._static_pumping import *
from ._static_reactor import *
from ._static_screening import *
from ._static_sludge_pasteurization import *
from ._static_sludge_thickening import *
from ._dynamic_bioreactor import *
from ._bst_tank import *
from ._static_trucking import *

# Units that rely on other units
from ._static_activated_sludge_process import *
from ._dynamic_anaerobic_reactor import *
from ._static_anaerobic_reactor import *
from ._static_clarifier import *
from ._bst_distillation import *
from ._bst_flash import *
from ._static_hydroprocessing import *
from ._static_hydrothermal import *
from ._static_internal_circulation_rx import *
from ._static_lagoon import *
from ._dynamic_membrane_bioreactor import *
from ._static_membrane_distillation import *
from ._static_polishing_filter import *
from ._static_sedimentation import *
from ._static_sludge_treatment import *
from ._static_septic_tank import *
from ._static_toilet import *
from ._static_treatment_bed import *

# System-specific unit (public)
from ._static_biogenic_refinery import *
from ._static_eco_san import *
from ._static_reclaimer import *

from . import bst, static, dynamic

# From then on the order doesn't matter, listed alphabetically
from . import (
        _bst_abstract,
        _dynamic_abstract,
        _static_abstract,
        _static_activated_sludge_process,
        _dynamic_anaerobic_reactor,
        _static_anaerobic_reactor,
        _static_clarifier,
        _static_combustion,
        _bst_compressor,
        _static_crop_application,
        _bst_distillation,
        _dynamic_influent,
        _static_electrochemical_cell,
        _static_excretion,
        _bst_facilities,
        _bst_flash,
        _bst_heat_exchanging,
        _static_hydroprocessing,
        _static_hydrothermal,
        _static_internal_circulation_rx,
        _dynamic_junction,
        _static_lagoon,
        _dynamic_membrane_bioreactor,
        _static_membrane_distillation,
        _static_membrane_gas_extraction,
        _static_metal_dosage,
        _static_non_reactive,
        _static_polishing_filter,
        _bst_pumping,
        _dynamic_pumping,
        _static_pumping,
        _static_reactor,
        _static_screening,
        _static_sedimentation,
        _static_septic_tank,
        _static_sludge_pasteurization,
        _static_sludge_thickening,
        _dynamic_bioreactor,
        _bst_tank,
        _static_toilet,
        _static_treatment_bed,
        _static_trucking,

        # System-specific units (public)
        _static_biogenic_refinery,
        _static_eco_san,
        _static_reclaimer,
        _static_sludge_treatment,
        )


__all__ = (
        *_bst_abstract.__all__,
        *_dynamic_abstract.__all__,
        *_static_abstract.__all__,
        *_static_activated_sludge_process.__all__,
        *_dynamic_anaerobic_reactor.__all__,
        *_static_anaerobic_reactor.__all__,
        *_static_clarifier.__all__,
        *_static_combustion.__all__,
        *_bst_compressor.__all__,
        *_static_crop_application.__all__,
        *_bst_distillation.__all__,
        *_dynamic_influent.__all__,
        *_static_electrochemical_cell.__all__,
        *_static_excretion.__all__,
        *_bst_facilities.__all__,
        *_bst_flash.__all__,
        *_bst_heat_exchanging.__all__,
        *_static_hydroprocessing.__all__,
        *_static_hydrothermal.__all__,
        *_static_internal_circulation_rx.__all__,
        *_dynamic_junction.__all__,
        *_static_lagoon.__all__,
        *_dynamic_membrane_bioreactor.__all__,
        *_static_membrane_distillation.__all__,
        *_static_metal_dosage.__all__,
        *_static_non_reactive.__all__,
        *_static_polishing_filter.__all__,
        *_bst_pumping.__all__,
        *_dynamic_pumping.__all__,
        *_static_pumping.__all__,
        *_static_reactor.__all__,
        *_static_screening.__all__,
        *_static_sedimentation.__all__,
        *_static_septic_tank.__all__,
        *_static_sludge_pasteurization.__all__,
        *_static_sludge_thickening.__all__,
        *_dynamic_bioreactor.__all__,
        *_bst_tank.__all__,
        *_static_toilet.__all__,
        *_static_treatment_bed.__all__,
        *_static_trucking.__all__,

        # System-specific units (public)
        *_static_biogenic_refinery.__all__,
        *_static_reclaimer.__all__,
        *_static_eco_san.__all__,
        *_static_sludge_treatment.__all__,
        *_static_membrane_gas_extraction.__all__,
        'bst',
        'static',
        'dynamic',
        )
