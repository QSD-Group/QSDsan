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
__all__ = ('test_exposan', 'test_all_exposan_modules_accounted_for')

import pkgutil
import pytest

pytestmark = [pytest.mark.integration, pytest.mark.slow]


# --- per-system construction smoke checks ---------------------------------
# This is the non-blocking integration *canary*: it confirms each EXPOsan system
# still constructs/wires up against the current QSDsan. It deliberately does NOT
# simulate the dynamic systems or static systems with a long simulation time (e.g. `hap`)# -- QSDsan's own dynamic machinery is covered by
# test_dyn_sys.py and the example treatment systems, and full simulation is owned
# by EXPOsan's CI (which runs the EXPOsan suite against QSDsan main). Hence
# load(simulate=False) / create_system() without simulate(), keeping this lane fast.

def _load(module_name, **load_kwargs):
    """Import an EXPOsan module and call its load(): dynamic systems construct only
    (load defaults to simulate=False); steady-state systems solve on load."""
    import importlib
    importlib.import_module(f'exposan.{module_name}').load(**load_kwargs)


def _summarize(module_name, sys_attrs):
    """load() (steady-state systems already solve on load) and print summaries."""
    import importlib
    mod = importlib.import_module(f'exposan.{module_name}')
    mod.load()
    mod.print_summaries(tuple(getattr(mod, a) for a in sys_attrs))


def _run_asm():
    from exposan import asm
    for process_model in ('ASM1', 'ASM2d'):
        for aerated in (False, True):
            asm.load(reload=True, process_model=process_model, aerated=aerated)


def _run_bsm2():
    from exposan import bsm2
    for kind in ('bsm2', 'bsm2p'):
        bsm2.load(kind=kind)  # changing `kind` triggers reconstruction


def _run_metab():
    # metab has no single default config; construct (do not simulate) a few
    # representative reactor / gas-extraction combinations via its central load().
    from exposan import metab
    for reactor_type, gas_extraction in (('UASB', 'M'), ('FB', 'H'), ('PB', 'P')):
        metab.load(reload=True, n_stages=2, reactor_type=reactor_type,
                   gas_extraction=gas_extraction, tot_HRT=4)


def _run_werf():
    # werf has 18 configs behind one central load(); construct one representative
    # full plant (N1: mASM2d/ADM1p + AD + junctions) as the canary.
    from exposan import werf
    werf.load('N1')


# Load-only for dynamic systems/static systems with long solve times; load+print for static systems with fast solve times. They all get exercised by the EXPOsan CI, but this ensures they at least construct against the current QSDsan and don't have any immediate import-time errors.
SYSTEMS = {
    'adm': lambda: _load('adm'),
    'asm': _run_asm,
    'bsm1': lambda: _load('bsm1'),
    'bsm2': _run_bsm2,
    'cas': lambda: _load('cas'),
    'htl': lambda: _load('htl'),
    'metab': _run_metab,
    'biogenic_refinery': lambda: _summarize('biogenic_refinery', ('sysA', 'sysB', 'sysC', 'sysD')),
    'bwaise': lambda: _summarize('bwaise', ('sysA', 'sysB', 'sysC')),
    'eco_san': lambda: _summarize('eco_san', ('sysA', 'sysB', 'sysC')),
    'reclaimer': lambda: _summarize('reclaimer', ('sysA', 'sysB', 'sysC', 'sysD')),
    'pou_disinfection': lambda: _summarize('pou_disinfection', ('sysA', 'sysB', 'sysC', 'sysD')),
    'hap': lambda: _load('hap'),
    'werf': _run_werf,
    'pm2_batch': lambda: _load('pm2_batch'),
    'pm2_ecorecover': lambda: _load('pm2_ecorecover'),
    'biobinder': lambda: _load('biobinder'),            # construct only; distillation sim is CI-fragile
    'saf': lambda: _load('saf'),                        # construct only; distillation sim is CI-fragile
    'g2rt': lambda: _load('g2rt'),
}

# Modules deliberately NOT exercised here, each with a reason. A new EXPOsan
# module that is neither in SYSTEMS nor here will fail the completeness guard,
# forcing a conscious decision instead of silent omission.
# NOTE TO MAINTAINER: refine these reasons; this seeds the registry with the
# modules currently uncovered by test_exposan.
KNOWN_SKIP = {
    'new_generator': 'NDA-protected system, no public entry point',
}


def _discover_exposan_modules():
    import exposan
    return {m.name for m in pkgutil.iter_modules(exposan.__path__)
            if not m.name.startswith('_') and m.name != 'utils'}


def test_all_exposan_modules_accounted_for():
    """Every EXPOsan system module must be either tested or explicitly skipped."""
    discovered = _discover_exposan_modules()
    accounted = set(SYSTEMS) | set(KNOWN_SKIP)
    unaccounted = discovered - accounted
    assert not unaccounted, (
        f"EXPOsan modules neither tested nor in KNOWN_SKIP: {sorted(unaccounted)}. "
        f"Add run-logic to SYSTEMS or an entry (with reason) to KNOWN_SKIP."
    )
    stale = accounted - discovered
    assert not stale, (
        f"SYSTEMS/KNOWN_SKIP reference modules that no longer exist: {sorted(stale)}."
    )


def test_exposan():
    for name, run in SYSTEMS.items():
        run()


if __name__ == '__main__':
    test_all_exposan_modules_accounted_for()
    test_exposan()
