#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

# Units of measure
import os
from thermosteam.units_of_measure import (
    ureg,
    AbsoluteUnitsOfMeasure as auom,
    RelativeUnitsOfMeasure as ruom
    )

path = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                    'units_definition.txt')
ureg.load_definitions(path)


from . import (
    cod,
    colors,
    components,
    construction,
    formatting,
    loading,
    misc,
    model_eval,
    parsing,
    utilities,
    wwt_design,
    scope,
    )

from biosteam.utils import NotImplementedMethod

from .cod import *
from .colors import *
from .components import *
from .construction import *
from .formatting import *
from .loading import *
from .misc import *
from .model_eval import *
from .parsing import *
from .utilities import *
from .wwt_design import *
from .scope import *

__all__ = (
    'ureg', 'auom', 'ruom',
    'NotImplementedMethod',
    *cod.__all__,
    *colors.__all__,
    *components.__all__,
    *construction.__all__,
    *formatting.__all__,
    *loading.__all__,
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