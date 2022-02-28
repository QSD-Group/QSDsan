#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

# Units that do not rely on other units
from ._abstract import *
from ._clarifier import *
from ._combustion import *
from ._crop_application import *
from ._dynamic_influent import *
from ._decay import *
from ._electrochemical_cell import *
from ._excretion import *
from ._hx import *
from ._lumped_cost import *
from ._pumping import *
from ._screening import *
from ._sludge_thickening import *
from ._suspended_growth_bioreactor import *
from ._tanks import *
from ._trucking import *
from ._encapsulation_bioreactor import *

# Units that rely on other units
from ._activated_sludge_process import *
from ._anaerobic_reactors import *
from ._internal_circulation_rx import *
from ._lagoon import *
from ._membrane_bioreactors import *
from ._polishing_filter import *
from ._sedimentation import *
from ._toilets import *
from ._treatment_beds import *


# From then on the order doesn't matter, listed alphabetically
from . import (
    _abstract,
    _activated_sludge_process,
    _anaerobic_reactors,
    _clarifier,
    _combustion,
    _crop_application,
    _decay,
    _dynamic_influent,
    _electrochemical_cell,
    _encapsulation_bioreactor,
    _excretion,
    _hx,
    _internal_circulation_rx,
    _lagoon,
    _lumped_cost,
    _membrane_bioreactors,
    _polishing_filter,
    _pumping,
    _screening,
    _sedimentation,
    _sludge_thickening,
    _suspended_growth_bioreactor,
    _tanks,
    _toilets,
    _treatment_beds,
    _trucking,
    )

__all__ = (
    *_abstract.__all__,
    *_activated_sludge_process.__all__,
    *_anaerobic_reactors.__all__,
    *_clarifier.__all__,
    *_combustion.__all__,
    *_crop_application.__all__,
    *_decay.__all__,
    *_dynamic_influent.__all__,
    *_electrochemical_cell.__all__,
    *_encapsulation_bioreactor.__all__,
    *_excretion.__all__,
    *_hx.__all__,
    *_internal_circulation_rx.__all__,
    *_lagoon.__all__,
    *_lumped_cost.__all__,
    *_membrane_bioreactors.__all__,
    *_polishing_filter.__all__,
    *_pumping.__all__,
    *_screening.__all__,
    *_sedimentation.__all__,
    *_sludge_thickening.__all__,
    *_suspended_growth_bioreactor.__all__,
    *_tanks.__all__,
    *_toilets.__all__,
    *_treatment_beds.__all__,
    *_trucking.__all__,
    )