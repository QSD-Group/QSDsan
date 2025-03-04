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
from thermosteam.units_of_measure import (
    ureg,
    AbsoluteUnitsOfMeasure as auom,
    RelativeUnitsOfMeasure as ruom,
    )

# Additional unit definition
ureg.define('sq_m = m2')
ureg.define('cu_m = m3')
ureg.define('sq_cm = cm2')
ureg.define('cu_cm = cm3')
ureg.define('sq_ft = ft2')
ureg.define('cu_ft = ft3')
ureg.define('cu_in = in3')
ureg.define('yd3 = yard**3 = yd3 = cu_yd')
ureg.define('cfm = cf/minute = CFM')
ureg.define('cfs = cf/second = CFS')
ureg.define('yr = year = yr = y')
ureg.define('hr = hour = hr = h')
ureg.define('d = day')
ureg.define('each = count = ea')
ureg.define('unit = count')
ureg.define('point = points')
ureg.define('MGD = 1e6 * gallon / day')
ureg.define('mgd = MGD')

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

__all__ = (
    'ureg', 'auom', 'ruom',
    'NotImplementedMethod',
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