#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>
    Joy Cheung

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/master/LICENSE.txt
for license details.
'''

import pkg_resources
try:
    __version__ = pkg_resources.get_distribution('qsdsan').version
except pkg_resources.DistributionNotFound:  # pragma: no cover
    __version__ = None

import biosteam as bst
CEPCI = bst.CE # Chemical Engineering Plant Cost Index
CEPCI_by_year = bst.units.design_tools.CEPCI_by_year
del bst
currency = 'USD'

from .utils import descriptors
from ._cod import *
from ._component import *
from ._components import *
from ._waste_stream import *
from ._impact_indicator import *
from ._impact_item import *
from ._construction import *
from ._transportation import *
from ._sanunit import *
from ._simple_tea import *
from ._lca import *

from . import (
    _units_of_measure, # if not included here, then need to add to setup.py
    _cod,
    _component,
    _components,
    _waste_stream,
    _impact_indicator,
    _impact_item,
    _construction,
    _transportation,
    _sanunit,
    _simple_tea,
    _lca,
    utils,
    sanunits,
    systems,
    stats,
    )

utils.secondary_importing()

__all__ = (
    *_cod.__all__,
    *_component.__all__,
    *_components.__all__,
    *_waste_stream.__all__,
    *_impact_indicator.__all__,
    *_impact_item.__all__,
    *_construction.__all__,
    *_transportation.__all__,
    *_sanunit.__all__,
    *_simple_tea.__all__,
    *_lca.__all__,
           )






