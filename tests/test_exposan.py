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
__all__ = ('test_exposan',)

def test_exposan():
    ##### Systems without costs/impacts #####
    from exposan.adm import create_system as create_adm_system
    adm_sys = create_adm_system()
    adm_sys.simulate(state_reset_hook='reset_cache', t_span=(0, 200), method='BDF')

    from exposan.asm import create_system as create_asm_system
    for process_model in ('ASM1', 'ASM2d'):
        for aerated in (False, True):
            asm_sys = create_asm_system(process_model=process_model, aerated=aerated)
            asm_sys.simulate(t_span=(0, 10), method='BDF')

    from qsdsan.utils import get_SRT
    from exposan.bsm1 import create_system as create_bsm1_system, biomass_IDs
    bsm1_sys = create_bsm1_system()
    bsm1_sys.simulate(t_span=(0,10), method='BDF')
    print(get_SRT(bsm1_sys, biomass_IDs=biomass_IDs)) # to test the `get_SRT` function

    from exposan.cas import create_system as create_cas_system
    cas_sys = create_cas_system()
    cas_sys.simulate()

    try: # skip this test while GH is still using the outdated version of EXPOsan
        from exposan.interface import create_system as create_inter_system
        sys_inter = create_inter_system()
        sys_inter.simulate(t_span=(0, 3))
    except: pass


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