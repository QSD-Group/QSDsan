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

import biosteam as bst, qsdsan as qs

__all__ = (
    'ProcessWaterCenter',
    )

class ProcessWaterCenter(bst.facilities.ProcessWaterCenter, qs.SanUnit):
    '''
    biosteam.facilities.ProcessWaterCenter with QSDsan properties.
    
    See Also
    --------
    `biosteam.facilities.ProcessWaterCenter <https://biosteam.readthedocs.io/en/latest/API/facilities/ProcessWaterCenter.html>`_
    '''