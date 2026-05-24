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

# Units of measure
from ..units_of_measure import ureg, auom, ruom

from . import (
    cod,
    colors,
    components,
    construction,
    formatting,
    loading,
    dynamics,
    indices,
    misc,
    model_eval,
    parsing,
    utilities,
    wwt_design,
    scope,
    )

from .cod import *
from .colors import *
from .components import *
from .construction import *
from .formatting import *
from .loading import *
from .dynamics import *
from .indices import *
from .misc import *
from .model_eval import *
from .parsing import *
from .utilities import *
from .wwt_design import *
from .scope import *

# BioSTEAM/Thermosteam helpers re-exported so users can build and evaluate units
# through ``qsdsan.utils`` without ``import biosteam``/``import thermosteam``:
#   - rho_to_V/V_to_rho: density <-> molar-volume conversion (component models)
#   - cost: the ``@cost`` decorator for adding cost/design to custom units
#   - var_columns/var_indices: extract Model parameter/metric labels from results
# These are plain re-exports; ``test_public_api`` asserts identity so a
# BioSTEAM/Thermosteam rename surfaces as a failing test rather than a user error.
from thermosteam.functional import rho_to_V, V_to_rho
from biosteam.units.decorators import cost
from biosteam.evaluation._utils import var_columns, var_indices

__all__ = (
    'ureg', 'auom', 'ruom',
    'NotImplementedMethod',
    'rho_to_V', 'V_to_rho', 'cost', 'var_columns', 'var_indices',
    *cod.__all__,
    *colors.__all__,
    *components.__all__,
    *construction.__all__,
    *formatting.__all__,
    *loading.__all__,
    *dynamics.__all__,
    *indices.__all__,
    *model_eval.__all__,
    *misc.__all__,
    *parsing.__all__,
    *utilities.__all__,
    *wwt_design.__all__,
    *scope.__all__,
    )


# This tiered importing is because some modules in utils need to be imported
# before the the main modules (e.g., _component) since the main modules depend
# on them, while other modules in utils depend on the main modules.
def _secondary_importing():
    from . import (
        doc_examples,
        )
