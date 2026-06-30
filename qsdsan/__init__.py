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

import importlib as _importlib
import importlib.metadata as impmeta
try:
    __version__ = impmeta.version('qsdsan')
except impmeta.PackageNotFoundError:
    __version__ = None
del impmeta

# BioSTEAM/Thermosteam APIs
import biosteam as _bst
Chemical = _bst.Chemical
Chemicals = _bst.Chemicals
CompiledChemicals = _bst.CompiledChemicals
Stream = _bst.Stream
MultiStream = _bst.MultiStream
set_thermo = _bst.settings.set_thermo
get_components = _bst.settings.get_chemicals
get_thermo = _bst.settings.get_thermo

HeatUtility = _bst.HeatUtility
PowerUtility = _bst.PowerUtility
Unit = _bst.Unit
System = _bst.System
AgileSystem = _bst.AgileSystem
Facility = _bst.Facility
Scope = _bst.utils.Scope
Model = _bst.Model
Metric = _bst.Metric
Parameter = _bst.Parameter
default_utilities = _bst.default_utilities
get_OSBL = _bst.get_OSBL

# Thermo/utility configuration objects, re-exported so users can configure QSDsan
# without ``import biosteam``. ``settings`` and ``stream_utility_prices`` are
# mutated in place (e.g. ``qs.settings.thermo.ideal()``), never rebound, so a
# plain re-export of the shared object is correct; ``test_public_api`` asserts the
# identity so a future BioSTEAM rebind would fail loudly. The settable scalar
# global ``bst.CE`` is exposed separately as the ``CEPCI`` property below.
settings = _bst.settings
# ``preferences`` controls display/diagram defaults (dark mode, stream labels,
# graphviz format, etc.); re-exported so users can tweak them via ``qs.preferences``
# and as the home for any future QSDsan-specific display preferences.
preferences = _bst.preferences
stream_utility_prices = _bst.stream_utility_prices
Thermo = _bst.Thermo
UtilityAgent = _bst.UtilityAgent
MissingStream = _bst.utils.MissingStream

# ``biosteam.report`` is imported eagerly by ``import biosteam`` above, so a plain
# re-export (no lazy machinery) is enough for ``qs.report.voc_table`` etc.
report = _bst.report

# Reaction APIs (defined in Thermosteam, re-exported through BioSTEAM) so users can
# do ``from qsdsan import Reaction`` instead of reaching into BioSTEAM/Thermosteam.
Reaction = _bst.Reaction
ReactionItem = _bst.ReactionItem
ReactionSet = _bst.ReactionSet
ParallelReaction = _bst.ParallelReaction
SeriesReaction = _bst.SeriesReaction
ReactionSystem = _bst.ReactionSystem
Rxn = _bst.Rxn
RxnI = _bst.RxnI
RxnS = _bst.RxnS
PRxn = _bst.PRxn
SRxn = _bst.SRxn
RxnSys = _bst.RxnSys

# Temporary placeholder — will be upgraded to SanMainFlowsheet at module bottom
main_flowsheet = _bst.main_flowsheet

from ._sanflowsheet import SanFlowsheet, SanMainFlowsheet

# Global variables
currency = 'USD'
CHECK_IMPACT_INDICATOR_CONSISTENCY = True
CHECK_IMPACT_ITEM_CONSISTENCY = True


from . import units_of_measure, utils
CEPCI_by_year = utils.indices.tea_indices['CEPCI_by_year']
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
    )

utils._secondary_importing()
for _slot in utils.doc_examples.__all__:
    setattr(utils, _slot, getattr(utils.doc_examples, _slot))

# Add the `pump` decorator to the util module
def wwtpump(*args, **kwargs):
    from .unit_operations import wwtpump as _wwtpump
    return _wwtpump(*args, **kwargs)

utils.__all__ = (*utils.__all__, 'wwtpump')
setattr(utils, 'wwtpump', wwtpump)

_lazy_modules = frozenset(('equipments', 'process_models', 'unit_operations', 'stats'))
_legacy_aliases = {'sanunits': 'unit_operations', 'processes': 'process_models'}

def __getattr__(name):
    target = _legacy_aliases.get(name, name)
    if target in _lazy_modules:
        module = _importlib.import_module(f'{__name__}.{target}')
        globals()[target] = module
        globals()[name] = module
        return module
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

def __dir__():
    return sorted((*globals(), *_lazy_modules, *_legacy_aliases, 'CEPCI'))

def default(utilities=True, CEPCI=True, flowsheet=True):
    '''
    Reset utilities, the chemical plant cost index (CEPCI), and/or the active
    flowsheet back to ``qsdsan``'s defaults (whichever are requested).

    Parameters
    ----------
    utilities : bool, optional
        Whether to reset the heating/cooling utility agents (and their prices)
        to BioSTEAM's defaults. The default is True.
    CEPCI : bool, optional
        Whether to reset the Chemical Engineering Plant Cost Index (:func:`qsdsan.CEPCI`)
        to BioSTEAM's default. The default is True.
    flowsheet : bool, optional
        Whether to reset the active flowsheet to a fresh one named ``'default'``,
        clearing all registered streams, units, systems, and LCA objects
        (impact indicators/items, construction, and transportation). The default
        is True.

    Notes
    -----
    Regardless of the arguments, the stream/unit/system auto-naming counters and
    the cached chemical lookups are always reset (this mirrors ``biosteam.default``).

    Unlike ``biosteam.default`` (whose ``flowsheet`` argument defaults to False),
    ``qsdsan.default`` resets the flowsheet by default, pass ``flowsheet=False`` to keep the current flowsheet.

    See Also
    --------
    `biosteam.default <https://biosteam.readthedocs.io/en/latest/API/process_tools/default.html>`_
    '''
    _bst.default(utilities=utilities, CEPCI=CEPCI, flowsheet=False)
    # `biosteam.default` resets the native stream/unit/system ticket counters; do the
    # same for the qsdsan LCA classes, whose ticket_numbers are also class-level
    # globals that `set_flowsheet` does not swap.
    for i in (ImpactIndicator, ImpactItem, Construction, Transportation):
        i.ticket_numbers.clear()
    if flowsheet:
        main_flowsheet.set_flowsheet('default', new=True)

# ── Upgrade main_flowsheet to SanMainFlowsheet in-place ──────────────────────
# _construction.py / _transportation.py imported main_flowsheet by reference
# above; upgrading the object in-place keeps all those references valid.
from biosteam._unit import AbstractUnit as _AbstractUnit
from thermosteam import AbstractStream as _AbstractStream
from biosteam._flowsheet import MainFlowsheet as _BstMainFlowsheet

_qs_default_fs = SanFlowsheet.from_registries(
    'default',
    _AbstractStream.registry,
    _AbstractUnit.registry,
    _bst.System.registry,
)
# Upgrade the existing _bst.main_flowsheet instance to SanMainFlowsheet
main_flowsheet.__class__ = SanMainFlowsheet
# set_flowsheet swaps __dict__ and LCA registries
main_flowsheet.set_flowsheet(_qs_default_fs)

Flowsheet = SanFlowsheet
F = main_flowsheet

del _qs_default_fs, _AbstractUnit, _AbstractStream, _BstMainFlowsheet


# ── Expose BioSTEAM's global cost index as a settable ``qsdsan.CEPCI`` ────────
# Lets users read/set the CEPCI without importing biosteam, and pairs consistently
# with ``qsdsan.CEPCI_by_year``. Plain module-level attributes cannot intercept
# assignment, so ``CEPCI`` is a property on a ModuleType subclass; the module-level
# ``__getattr__``/``__dir__`` above still work.
import sys as _sys

class _SanModule(_sys.modules[__name__].__class__):
    @property
    def CEPCI(self):
        '''
        [float] Chemical Engineering Plant Cost Index (CEPCI) used to scale equipment costs.

        Look up the index for a given year with ``qsdsan.CEPCI_by_year`` (e.g.,
        ``qsdsan.CEPCI = qsdsan.CEPCI_by_year[2023]`` to report costs in 2023 dollars).

        Setting ``qsdsan.CEPCI`` is the same as setting ``biosteam.CE`` (BioSTEAM
        abbreviates the index as ``CE``); this is a live view of that value.
        '''
        return _bst.CE
    @CEPCI.setter
    def CEPCI(self, value):
        _bst.CE = value

_sys.modules[__name__].__class__ = _SanModule


# BioSTEAM/Thermosteam objects re-exported for a self-contained ``qsdsan`` surface
_bst_reexports = (
    'Chemical', 'Chemicals', 'CompiledChemicals',
    'Stream', 'MultiStream', 'MissingStream',
    'set_thermo', 'get_components', 'get_thermo', 'settings', 'preferences',
    'HeatUtility', 'PowerUtility', 'stream_utility_prices',
    'Unit', 'System', 'AgileSystem', 'Facility',
    'Scope', 'Model', 'Metric', 'Parameter',
    'default_utilities', 'get_OSBL',
    'Thermo', 'UtilityAgent', 'report',
    'Reaction', 'ReactionItem', 'ReactionSet',
    'ParallelReaction', 'SeriesReaction', 'ReactionSystem',
    'Rxn', 'RxnI', 'RxnS', 'PRxn', 'SRxn', 'RxnSys',
    'CEPCI', 'CEPCI_by_year',
    'Flowsheet', 'main_flowsheet', 'F',
    )

__all__ = (
    'units_of_measure',
    'SanFlowsheet',
    'SanMainFlowsheet',
    *_bst_reexports,
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
