# -*- coding: utf-8 -*-
"""
Simulate all 18 WERF WRRF systems and extract stream flow data.
Outputs: wrrf_flow_data.json (same directory as this script)
"""

import warnings
warnings.filterwarnings('ignore')

import json, os, time

# ── path setup ──────────────────────────────────────────────────────────────
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT  = os.path.join(HERE, 'wrrf_flow_data.json')

import sys
sys.path.insert(0, ROOT)

from exposan.werf import create_system

SYS_LIST = [
    'B1', 'B2', 'B3',
    'C1', 'C2', 'C3',
    'E2', 'E2P',
    'F1',
    'G1', 'G2', 'G3',
    'H1',
    'I1', 'I2', 'I3',
    'N1', 'N2',
]

results = {}

for sid in SYS_LIST:
    t0 = time.time()
    print(f'\n{"="*40}')
    print(f'  {sid}')
    print(f'{"="*40}')

    try:
        sys_obj = create_system(sid)
    except Exception as e:
        print(f'  create_system failed: {e}')
        results[sid] = {'_error': str(e)}
        continue

    # simulate with fallback
    try:
        sys_obj.simulate(t_span=(0, 300), method='BDF')
    except Exception:
        try:
            sys_obj.simulate(state_reset_hook='reset_cache', t_span=(0, 300), method='BDF')
        except Exception as e2:
            print(f'  simulate failed: {e2}')
            results[sid] = {'_error': str(e2)}
            continue

    elapsed = time.time() - t0
    print(f'  Simulated in {elapsed:.1f}s')

    # ── extract stream data ──────────────────────────────────────────────────
    stream_data = {}
    for s in sys_obj.streams:
        if s is None or not s.ID:
            continue

        try:
            q_m3d   = round(s.F_vol * 24, 2)        # m³/d  (F_vol in m³/hr)
            is_gas  = (s.phase == 'g')
            if is_gas:
                tss_kgd = None
            else:
                tss_conc = s.get_TSS()               # mg/L
                tss_kgd  = round(tss_conc * q_m3d / 1000, 2)  # kg/d
            stream_data[s.ID] = {
                'id'     : s.ID,
                'q_m3d'  : q_m3d,
                'tss_kgd': tss_kgd,
                'is_gas' : is_gas,
                'phase'  : s.phase,
            }
        except Exception as e:
            stream_data[s.ID] = {'id': s.ID, 'error': str(e)}

    results[sid] = stream_data

    # print a summary
    print(f'  {len(stream_data)} streams:')
    for sid2, sd in stream_data.items():
        if 'error' in sd:
            print(f'    {sid2:30s}  ERROR: {sd["error"]}')
        elif sd['is_gas']:
            print(f'    {sid2:30s}  Q={sd["q_m3d"]:>8.1f} m³/d  [gas]')
        else:
            print(f'    {sid2:30s}  Q={sd["q_m3d"]:>8.1f} m³/d  TSS={sd["tss_kgd"]:>8.1f} kg/d')

# ── write JSON ──────────────────────────────────────────────────────────────
os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, 'w') as f:
    json.dump(results, f, indent=2)

print(f'\n\nSaved → {OUT}')
