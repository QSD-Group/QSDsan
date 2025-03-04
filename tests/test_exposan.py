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
    from qsdsan import default
    default()
    
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
    print(get_SRT(bsm1_sys, biomass_IDs=biomass_IDs['asm1'])) # to test the `get_SRT` function

    #!!! Will use bsm2 to test the junction models
    # from exposan.interface import create_system as create_inter_system
    # sys_inter = create_inter_system()
    # sys_inter.simulate(method='BDF', t_span=(0, 3)) # the default 'RK45' method can't solve it

    from exposan.cas import create_system as create_cas_system
    cas_sys = create_cas_system()
    cas_sys.simulate()

    ##### Systems with costs/impacts #####
    from qsdsan.utils import clear_lca_registries
    
    clear_lca_registries()
    from exposan import biogenic_refinery as br
    br.load()
    br.print_summaries((br.sysA, br.sysB, br.sysC, br.sysD))

    clear_lca_registries()
    from exposan import bwaise as bw
    bw.load()
    bw.print_summaries((bw.sysA, bw.sysB, bw.sysC))

    clear_lca_registries()
    from exposan import eco_san as es
    es.load()
    es.print_summaries((es.sysA, es.sysB, es.sysC))
    
    clear_lca_registries()
    from exposan import htl
    htl.load()
    
    clear_lca_registries()
    from exposan import metab
    UASB_M = metab.create_system(n_stages=2, reactor_type='UASB', gas_extraction='M', tot_HRT=4)
    UASB_M.simulate(state_reset_hook='reset_cache', method='BDF', t_span=(0, 400))
    FB_H = metab.create_system(n_stages=2, reactor_type='FB', gas_extraction='H', tot_HRT=4)
    # # Just simulate one system to save testing time
    # # (all configurations are included in EXPOsan test)
    # FB_H.simulate(state_reset_hook='reset_cache', method='BDF', t_span=(0, 400))
    PB_P = metab.create_system(n_stages=2, reactor_type='PB', gas_extraction='P', tot_HRT=4)
    # Might fail the first time it runs, re-running will usually fix the problem
    # try: PB_P.simulate(state_reset_hook='reset_cache', method='BDF', t_span=(0, 400))
    # except: PB_P.simulate(state_reset_hook='reset_cache', method='BDF', t_span=(0, 400))

    clear_lca_registries()
    from exposan import reclaimer as re
    re.load()
    re.print_summaries((re.sysA, re.sysB, re.sysC, re.sysD))
    
    clear_lca_registries()
    from exposan import pou_disinfection as pou
    pou.load()
    pou.print_summaries((pou.sysA, pou.sysB, pou.sysC, pou.sysD))


if __name__ == '__main__':
    test_exposan()