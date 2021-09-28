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

# This tiered importing is because some modules in utils need to be imported
# before the the main modules (e.g., _component) since the main modules depend
# on them, while other modules in utils depend on the main modules.

from . import (
    cod,
    colors,
    decorators,
    # descriptors, # currently not in use
    evaluation,
    getters,
    formatting,
    loading,
    parse,
    setters,
    units_of_measure,
    )

from biosteam.utils import NotImplementedMethod

from .cod import *
from .colors import *
from .decorators import *
# from .descriptors import *
from .evaluation import *
from .getters import *
from .formatting import *
from .loading import *
from .parse import *
from .setters import *
from .units_of_measure import *


__all__ = (
    *cod.__all__,
    *colors.__all__,
    *decorators.__all__,
    # *descriptors.__all__,
    *evaluation.__all__,
    *getters.__all__,
    *formatting.__all__,
    *loading.__all__,
    *parse.__all__,
    *setters.__all__,
    *units_of_measure.__all__,
    'NotImplementedMethod',
    )


def _secondary_importing():
    from . import (
        doc_examples,
        )