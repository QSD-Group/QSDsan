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

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

# **NOTE** PLEASE ORDER THE MODULES ALPHABETICALLY #

# Units that do not rely on other units
from ._abstract import *
from ._clarifier import *
from ._combustion import *
from ._compressor import *
from ._crop_application import *
from ._dynamic_influent import *
from ._electrochemical_cell import *
from ._excretion import *
from ._facilities import *
from ._heat_exchanging import *
from ._junction import *
from ._non_reactive import *
from ._pumping import *
from ._reactor import *
from ._screening import *
from ._sludge_pasteurization import *
from ._sludge_thickening import *
from ._suspended_growth_bioreactor import *
from ._tank import *
from ._trucking import *

# Units that rely on other units
from ._activated_sludge_process import *
from ._anaerobic_reactor import *
from ._distillation import *
from ._flash import *
from ._hydroprocessing import *
from ._hydrothermal import *
from ._internal_circulation_rx import *
from ._lagoon import *
from ._membrane_bioreactor import *
from ._membrane_distillation import *
from ._polishing_filter import *
from ._sedimentation import *
from ._septic_tank import *
from ._toilet import *
from ._treatment_bed import *

# System-specific unit (public)
from ._biogenic_refinery import *
from ._eco_san import *
from ._reclaimer import *

# From then on the order doesn't matter, listed alphabetically
from . import (
        _abstract,
        _activated_sludge_process,
        _anaerobic_reactor,
        _clarifier,
        _combustion,
        _compressor,
        _crop_application,
        _distillation,
        _dynamic_influent,
        _electrochemical_cell,
        _excretion,
        _facilities,
        _flash,
        _heat_exchanging,
        _hydroprocessing,
        _hydrothermal,
        _internal_circulation_rx,
        _junction,
        _lagoon,
        _membrane_bioreactor,
        _membrane_distillation,
        _non_reactive,
        _polishing_filter,
        _pumping,
        _reactor,
        _screening,
        _sedimentation,
        _septic_tank,
        _sludge_pasteurization,
        _sludge_thickening,
        _suspended_growth_bioreactor,
        _tank,
        _toilet,
        _treatment_bed,
        _trucking,

        # System-specific units (public)
        _biogenic_refinery,
        _eco_san,
        _reclaimer,
        )


__all__ = (
        *_abstract.__all__,
        *_activated_sludge_process.__all__,
        *_anaerobic_reactor.__all__,
        *_clarifier.__all__,
        *_combustion.__all__,
        *_compressor.__all__,
        *_crop_application.__all__,
        *_distillation.__all__,
        *_dynamic_influent.__all__,
        *_electrochemical_cell.__all__,
        *_excretion.__all__,
        *_facilities.__all__,
        *_flash.__all__,
        *_heat_exchanging.__all__,
        *_hydroprocessing.__all__,
        *_hydrothermal.__all__,
        *_internal_circulation_rx.__all__,
        *_junction.__all__,
        *_lagoon.__all__,
        *_membrane_bioreactor.__all__,
        *_membrane_distillation.__all__,
        *_non_reactive.__all__,
        *_polishing_filter.__all__,
        *_pumping.__all__,
        *_reactor.__all__,
        *_screening.__all__,
        *_sedimentation.__all__,
        *_septic_tank.__all__,
        *_sludge_pasteurization.__all__,
        *_sludge_thickening.__all__,
        *_suspended_growth_bioreactor.__all__,
        *_tank.__all__,
        *_toilet.__all__,
        *_treatment_bed.__all__,
        *_trucking.__all__,

        # System-specific units (public)
        *_biogenic_refinery.__all__,
        *_reclaimer.__all__,
        *_eco_san.__all__,
        )