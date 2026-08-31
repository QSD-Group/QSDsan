#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:

    Joy Zhang <joycheung1994@gmail.com>

    Yalin Li <mailto.yalin.li@gmail.com>

    Yoel Cortes-Pena <yoelcortes@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

__all__ = ('test_dyn_sys', 'test_dynamic_influent_noncyclic')

def test_dyn_sys():
    from qsdsan import process_models as pc, unit_operations as su, set_thermo, System
    import numpy as np
    from numpy.testing import assert_allclose

    cmps = pc.create_asm1_cmps()
    DI = su.DynamicInfluent('Dyn_Inf')
    S1 = su.Splitter('Split', ins=DI-0, split=0.3, init_with='WasteStream')
    M1 = su.Mixer('Mix', ins=(S1-0, S1-1), outs=('Dyn_Eff'))
    sys = System('test_sys', path=(DI, S1, M1))
    sys.set_dynamic_tracker(DI.outs[0], M1.outs[0])

    t = 1
    t_step = 0.05
    sys.simulate(t_span=(0,t),
                 t_eval=np.arange(0, t+t_step, t_step))
    dinf = sys.units[0].outs[0]
    deff = sys.units[-1].outs[0]
    assert_allclose(deff.scope.record, dinf.scope.record, rtol=1e-12)


def test_dynamic_influent_noncyclic(tmp_path):
    """Non-cyclic input (first and last rows differ) used to enter a branch that
    called ``df.append`` -- a silent no-op on pre-pandas-2.0 and an AttributeError
    on pandas >=2.0. The fix uses ``pd.concat`` and re-assigns ``self._data``;
    the phantom final point should now actually be appended.
    """
    import numpy as np
    import pandas as pd
    from qsdsan import process_models as pc, unit_operations as su, set_thermo, System

    set_thermo(pc.create_asm1_cmps())
    p = tmp_path / 'noncyclic.csv'
    t = np.linspace(0, 10, 21)
    pd.DataFrame({'t': t,
                  'S_S': np.where(t < 4.0, 50.0, 300.0),
                  'Q': 1000.0}).to_csv(p, index=False)
    DI = su.DynamicInfluent('DI_nc', data_file=str(p))
    # Phantom point appended: _data has +1 row, t_end extrapolated past 10.
    assert len(DI._data) == 22
    assert DI._t_end > 10.0
    # Construction is still functional end-to-end.
    sys = System('s_nc', path=(DI,))
    sys.set_dynamic_tracker(DI.outs[0])
    sys.simulate(t_span=(0, 10), t_eval=np.linspace(0, 10, 11))


if __name__ == '__main__':
    test_dyn_sys()