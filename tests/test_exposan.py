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

# Just making sure the systems can run is sufficient
# TODO: change to the new loading routine
__all__ = ('test_exposan',)

def test_exposan():
    ##### Systems without costs/impacts #####
    from exposan import adm
    # adm.load()
    adm.sys.simulate(state_reset_hook='reset_cache', t_span=(0, 200), method='BDF')

    from exposan import asm
    asm1_an = asm.create_asm_validation_system(process_model='ASM1', aerated=False)
    asm1_aer = asm.create_asm_validation_system(process_model='ASM1', aerated=True)
    asm2d_an = asm.create_asm_validation_system(process_model='ASM2d', aerated=False)
    asm2d_aer = asm.create_asm_validation_system(process_model='ASM2d', aerated=True)
    syses = (asm1_an, asm1_aer, asm2d_an, asm2d_aer)
    for sys in syses: sys.simulate(state_reset_hook='reset_cache', t_span=(0, 10), method='BDF')

    from qsdsan.utils import get_SRT
    from exposan.bsm1 import bsm1
    bsm1.reset_cache()
    bsm1.simulate(t_span=(0,10), method='BDF')
    print(get_SRT(bsm1, biomass_IDs=('X_BH', 'X_BA'))) # to test the `get_SRT` function

    from exposan import cas
    cas.sys.simulate()

    ##### Systems with costs/impacts #####
    from exposan import biogenic_refinery as br
    br.load()
    br.print_summaries((br.sysA, br.sysB, br.sysC, br.sysD))

    from exposan import bwaise as bw
    bw.load()
    bw.print_summaries((bw.sysA, bw.sysB, bw.sysC))

    from exposan import eco_san as es
    es.load()
    es.print_summaries((es.sysA, es.sysB, es.sysC))

    from exposan import reclaimer as re
    re.load()
    re.print_summaries((re.sysA, re.sysB, re.sysC, re.sysD))


if __name__ == '__main__':
    test_exposan()