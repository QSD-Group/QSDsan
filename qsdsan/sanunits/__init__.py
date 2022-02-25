#!/usr/bin/env python3
	# -*- coding: utf-8 -*-


'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>
    Joy Zhang <joycheung1994@gmail.com>
    Lewis Rowles <stetsonsc@gmail.com>
    Lane To <lane20@illinois.edu>
    Smiti Mittal <smitimittal@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

# **NOTE** PLEASE ORDER THE MODULES ALPHABETICALLY #

# Units that do not rely on other units
from ._abstract import *
from ._carbonizer_base import *
from ._clarifier import *
from ._combustion import *
from ._control_boxes import *
from ._crop_application import *
from ._decay import *
from ._dynamic_influent import *
from ._electrochemical_cell import *
from ._excretion import *
from ._grinder import *
from ._housing import *
from ._hxs import *
from ._ion_exchange_NH3 import *
from ._lumped_cost import *
from ._pollution_control_device import *
from ._pumping import *
from ._screening import *
from ._sludge_thickening import *
from ._struvite_precipitation import *
from ._suspended_growth_bioreactor import *
from ._tanks import *
from ._trucking import *


# Units that rely on other units
from ._activated_sludge import *
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
	    _activated_sludge,
	    _anaerobic_reactors,
        _carbonizer_base,
	    _clarifier,
	    _combustion,
        _control_boxes,
	    _crop_application,
	    _decay,
	    _dynamic_influent,
	    _electrochemical_cell,
	    _excretion,
        _housing,
	    _hxs,
	    _internal_circulation_rx,
	    _lagoon,
	    _lumped_cost,
	    _membrane_bioreactors,
	    _polishing_filter,
        _pollution_control_device,
	    _pumping,
	    _screening,
	    _sedimentation,
	    _sludge_thickening,
        _struvite_precipitation,
	    _suspended_growth_bioreactor,
	    _tanks,
	    _toilets,
	    _treatment_beds,
	    _trucking,



        _ion_exchange_NH3,
        _grinder,
	    )


__all__ = (
	    *_abstract.__all__,
	    *_activated_sludge.__all__,
	    *_anaerobic_reactors.__all__,
        *_carbonizer_base.__all__,
	    *_clarifier.__all__,
	    *_combustion.__all__,
        *_control_boxes.__all__,
	    *_crop_application.__all__,
	    *_decay.__all__,
	    *_dynamic_influent.__all__,
	    *_electrochemical_cell.__all__,
	    *_excretion.__all__,
        *_housing.__all__,
	    *_hxs.__all__,
	    *_internal_circulation_rx.__all__,
	    *_lagoon.__all__,
	    *_lumped_cost.__all__,
	    *_membrane_bioreactors.__all__,
	    *_polishing_filter.__all__,
        *_pollution_control_device.__all__,
	    *_pumping.__all__,
	    *_screening.__all__,
	    *_sedimentation.__all__,
	    *_sludge_thickening.__all__,
        *_struvite_precipitation.__all__,
	    *_suspended_growth_bioreactor.__all__,
	    *_tanks.__all__,
	    *_toilets.__all__,
	    *_treatment_beds.__all__,
	    *_trucking.__all__,




        *_ion_exchange_NH3.__all__,
        *_grinder.__all__,
	    )