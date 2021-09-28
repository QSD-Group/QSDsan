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

# Check system environment, Python 3.7 and below will have issues unpickling saved results
import sys
py_version = sys.version.split('.')
_PY_MAJOR, _PY_MINOR = int(py_version[0]), int(py_version[1])

if (_PY_MAJOR, _PY_MINOR) <= (3, 7):
    from warnings import warn
    if (_PY_MAJOR, _PY_MINOR) >= (3, 5):
        try: import pickle5 as _pk
        except ModuleNotFoundError:
            warn(f'Python version {_PY_MAJOR}.{_PY_MINOR} does not support Pickle Protocol 5, '
                 'installing `pickle5` by running `pip install pickle5` in your '
                 'command/Anaconda prompt or terminal can reduce the loading time.\n'
                 'For further information, check https://pypi.org/project/pickle5/.')
            _pk = None
    else:
        warn(f'Python version {_PY_MAJOR}.{_PY_MINOR} does not support Pickle Protocol 5, '
             'and will have slower speed in when loading the default processes.')
        _pk = None
    del warn
else:
    import pickle as _pk

del sys, py_version


import pkg_resources
try:
    __version__ = pkg_resources.get_distribution('qsdsan').version
except pkg_resources.DistributionNotFound:  # pragma: no cover
    __version__ = None


import thermosteam as tmo
import biosteam as bst
Chemical = tmo.Chemical
Chemicals = tmo.Chemicals
CompiledChemicals = tmo.CompiledChemicals
Stream = tmo.Stream
set_thermo = tmo.settings.set_thermo
get_components = tmo.settings.get_chemicals
get_thermo = tmo.settings.get_thermo

MultiStream = tmo.MultiStream
Unit = bst.Unit
System = bst.System
Flowsheet = bst.Flowsheet
CEPCI = bst.CE # Chemical Engineering Plant Cost Index
CEPCI_by_year = bst.units.design_tools.CEPCI_by_year
del tmo, bst

currency = 'USD'

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
from ._simple_tea import *
from ._lca import *


from . import (
    _component,
    _components,
    _sanstream,
    _waste_stream,
    _process,
    _impact_indicator,
    _impact_item,
    _construction,
    _equipment,
    _transportation,
    _sanunit,
    _simple_tea,
    _lca,
    processes,
    equipments,
    sanunits,
    stats,
    )

utils._secondary_importing()
for _slot in utils.doc_examples.__all__:
    setattr(utils, _slot, getattr(utils.doc_examples, _slot))

__all__ = (
    *_component.__all__,
    *_components.__all__,
    *_sanstream.__all__,
    *_waste_stream.__all__,
    *_process.__all__,
    *_impact_indicator.__all__,
    *_impact_item.__all__,
    *_construction.__all__,
    *_transportation.__all__,
    *_equipment.__all__,
    *_sanunit.__all__,
    *_simple_tea.__all__,
    *_lca.__all__,
           )