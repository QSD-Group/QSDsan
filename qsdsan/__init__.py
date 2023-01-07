#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:

    Yalin Li <mailto.yalin.li@gmail.com>

    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

import pkg_resources
try:
    __version__ = pkg_resources.get_distribution('qsdsan').version
except pkg_resources.DistributionNotFound:  # pragma: no cover
    __version__ = None
del pkg_resources


# BioSTEAM/Thermosteam APIs
import thermosteam as tmo, biosteam as bst
Chemical = tmo.Chemical
Chemicals = tmo.Chemicals
CompiledChemicals = tmo.CompiledChemicals
Stream = tmo.Stream
MultiStream = tmo.MultiStream
set_thermo = tmo.settings.set_thermo
get_components = tmo.settings.get_chemicals
get_thermo = tmo.settings.get_thermo

HeatUtility = bst.HeatUtility
PowerUtility = bst.PowerUtility
Unit = bst.Unit
System = bst.System
Scope = bst.utils.Scope
Model = bst.Model
Metric = bst.Metric
Parameter = bst.Parameter
Flowsheet = bst.Flowsheet
main_flowsheet = bst.main_flowsheet
CEPCI_by_year = bst.units.design_tools.CEPCI_by_year
del tmo, bst

# Global variables
currency = 'USD'
CHECK_IMPACT_INDICATOR_CONSISTENCY = True
CHECK_IMPACT_ITEM_CONSISTENCY = True


from . import utils
from ._component import *
from ._components import *
from ._sanstream import *
from ._waste_stream import *
from ._process import *
from ._impact_indicator import *
from ._impact_item import *
from ._construction import *
from ._equipment import *
from ._transportation import *
from ._sanunit import *
from ._tea import *
from ._lca import *


from . import (
    _component,
    _components,
    _construction,
    _equipment,
    _impact_indicator,
    _impact_item,
    _lca,
    _process,
    _sanstream,
    _sanunit,
    _tea,
    _transportation,
    _waste_stream,
    equipments,
    processes,
    sanunits,
    stats,
    )

utils._secondary_importing()
for _slot in utils.doc_examples.__all__:
    setattr(utils, _slot, getattr(utils.doc_examples, _slot))

# Add the `pump` decorator to the util module
from .sanunits import wwtpump
utils.__all__ = (*utils.__all__, 'wwtpump')
setattr(utils, 'wwtpump', wwtpump)


__all__ = (
    *_component.__all__,
    *_components.__all__,
    *_construction.__all__,
    *_equipment.__all__,
    *_impact_indicator.__all__,
    *_impact_item.__all__,
    *_lca.__all__,
    *_process.__all__,
    *_sanstream.__all__,
    *_sanunit.__all__,
    *_tea.__all__,
    *_transportation.__all__,
    *_waste_stream.__all__,
    )