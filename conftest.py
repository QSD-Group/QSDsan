# -*- coding: utf-8 -*-
'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''
import pytest
import numpy as np

@pytest.fixture(autouse=True)
def set_np_legacy_mode():
    try: np.set_printoptions(legacy='1.25')
    except: pass
