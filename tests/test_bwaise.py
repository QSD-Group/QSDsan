#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/master/LICENSE.txt
for license details.
'''

def test_bwaise():
    from exposan import bwaise as bw
    bw.print_summaries((bw.sysA, bw.sysB, bw.sysC))
    
# If pytest runs this module, it calls the test_bwaise function and test exposan
if __name__ == '__main__':
    test_bwaise()