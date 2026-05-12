# -*- coding: utf-8 -*-
'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:

    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''



__all__ = ('test_process', 'test_kinetic_reaction')

def test_process():
    import pytest, os, qsdsan.process_models as pc
    from sympy import symbols, Eq
    from sympy.parsing.sympy_parser import parse_expr
    from math import isclose
    from qsdsan import set_thermo, Components, Process, Processes, CompiledProcesses

    cmps = Components.load_default()

    S_A = cmps.S_Ac.copy('S_A')
    S_ALK = cmps.S_CO3.copy('S_ALK')      # measured as g C
    S_F = cmps.S_F.copy('S_F')
    S_I = cmps.S_U_E.copy('S_I')
    S_N2 = cmps.S_N2.copy('S_N2')
    S_NH4 = cmps.S_NH4.copy('S_NH4')
    S_NO3 = cmps.S_NO3.copy('S_NO3')
    S_O2 = cmps.S_O2.copy('S_O2')
    S_PO4 = cmps.S_PO4.copy('S_PO4')
    X_AUT = cmps.X_AOO.copy('X_AUT')
    X_H = cmps.X_OHO.copy('X_H')
    X_I = cmps.X_U_OHO_E.copy('X_I')
    X_MeOH = cmps.X_FeOH.copy('X_MeOH')
    X_MeP = cmps.X_FePO4.copy('X_MeP')
    X_PAO = cmps.X_PAO.copy('X_PAO')
    X_PHA = cmps.X_PAO_PHA.copy('X_PHA')
    X_PP = cmps.X_PAO_PP_Lo.copy('X_PP')
    X_S = cmps.X_B_Subst.copy('X_S')

    S_I.i_N = 0.01
    S_F.i_N = 0.03
    X_I.i_N = 0.02
    X_S.i_N = 0.04
    X_H.i_N = X_PAO.i_N = X_AUT.i_N = 0.07

    S_I.i_P = 0.00
    S_F.i_P = 0.01
    X_I.i_P = 0.01
    X_S.i_P = 0.01
    X_H.i_P = X_PAO.i_P = X_AUT.i_P = 0.02

    cmps_asm2d = Components([S_O2, S_F, S_A, S_NH4, S_NO3, S_PO4, S_I, S_ALK, S_N2,
                             X_I, X_S, X_H, X_PAO, X_PP, X_PHA, X_AUT, X_MeOH, X_MeP])

    cmps_asm2d.compile(ignore_inaccurate_molar_weight=True)
    set_thermo(cmps_asm2d)

    p1 = Process('aero_hydrolysis',
                 'X_S -> [1-f_SI]S_F + [f_SI]S_I + [?]S_NH4 + [?]S_PO4 + [?]S_ALK',
                 ref_component='X_S',
                 rate_equation='K_h * S_O2/(K_O2+S_O2) * X_S/(K_X*X_H+X_S) * X_H',
                 parameters=('f_SI', 'K_h', 'K_O2', 'K_X'))

    f_SI = symbols('f_SI')
    assert abs(sum(p1._stoichiometry * p1._components.i_N).subs({'f_SI':1})) < 1e-8
    assert abs(sum(p1._stoichiometry * p1._components.i_N).subs({'f_SI':0})) < 1e-8
    assert abs(sum(p1._stoichiometry * p1._components.i_P).subs({'f_SI':1})) < 1e-8
    assert abs(sum(p1._stoichiometry * p1._components.i_charge).subs({'f_SI':1})) < 1e-8

    p1.set_parameters(f_SI = 0.0)
    assert p1.parameters['f_SI'] == 0.0
    assert Eq(p1._stoichiometry[p1._components._index['S_I']], parse_expr('1*f_SI'))

    p12 = Process('anox_storage_PP',
                  'S_PO4 + [Y_PHA]X_PHA + [?]S_NO3 -> X_PP + [?]S_N2 + [?]S_NH4 + [?]S_ALK',
                  ref_component='X_PP',
                  rate_equation='q_PP * S_O2/(K_O2+S_O2) * S_PO4/(K_PS+S_PO4) * S_ALK/(K_ALK+S_ALK) * (X_PHA/X_PAO)/(K_PHA+X_PHA/X_PAO) * (K_MAX-X_PP/X_PAO)/(K_PP+K_MAX-X_PP/X_PAO) * X_PAO * eta_NO3 * K_O2/S_O2 * S_NO3/(K_NO3+S_NO3)',
                  parameters=('Y_PHA', 'q_PP', 'K_O2', 'K_PS', 'K_ALK', 'K_PHA', 'eta_NO3', 'K_PP', 'K_NO3'),
                  conserved_for=('COD', 'N', 'P', 'NOD', 'charge'))

    p14 = Process('PAO_anox_growth',
                 '[1/Y_H]X_PHA + [?]S_NO3 + [?]S_PO4 -> X_PAO + [?]S_N2 + [?]S_NH4  + [?]S_ALK',
                 ref_component='X_PAO',
                 rate_equation='mu_PAO * S_O2/(K_O2 + S_O2) * S_NH4/(K_NH4 + S_NH4) * S_PO4/(K_P + S_PO4) * S_CO3/(K_ALK + S_ALK) * (X_PHA/X_PAO)/(K_PHA + X_PHA/X_PAO) * X_PAO * eta_NO3 * K_O2/S_O2 * S_NO3/(K_NO3 + S_NO3)',
                 parameters=('Y_H', 'mu_PAO', 'K_O2', 'K_NH4', 'K_P', 'K_ALK', 'K_PHA', 'eta_NO3', 'K_NO3'),
                 conserved_for=('COD', 'N', 'P', 'NOD', 'charge'))

    PAO_anox_processes = Processes([p12, p14])
    assert PAO_anox_processes.PAO_anox_growth.ref_component == X_PAO.ID
    with pytest.raises(AttributeError):
        print(PAO_anox_processes.production_rates)

    params = ('f_SI', 'Y_H', 'f_XI', 'Y_PO4', 'Y_PHA', 'Y_A',
              'K_h', 'eta_NO3', 'eta_fe', 'K_O2', 'K_NO3', 'K_X',
              'mu_H', 'q_fe', 'eta_NO3_deni', 'b_H', 'K_F', 'K_fe', 'K_A',
              'K_NH4', 'K_P', 'K_ALK', 'q_PHA', 'q_PP', 'mu_PAO', 'b_PAO',
              'b_PP', 'b_PHA', 'K_PS', 'K_PP', 'K_MAX', 'K_IPP', 'K_PHA',
              'mu_AUT', 'b_AUT', 'K_O2_AUT', 'K_NH4_AUT', 'K_ALK_2',
              'k_PRE', 'k_RED')

    path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'ASM2d_original.tsv')
    asm2d = Processes.load_from_file(path,
                                     conserved_for=('COD', 'N', 'P', 'charge'),
                                     parameters=params,
                                     compile=False)
    asm2d.extend(PAO_anox_processes)
    asm2d.compile()
    assert isinstance(asm2d, CompiledProcesses)
    assert p12 in asm2d
    assert set(asm2d.parameters.keys()) == set(params)

    pc.create_adm1_cmps()
    pc.create_asm1_cmps()
    pc.create_asm2d_cmps()
        

def test_kinetic_reaction():
    import pytest
    import qsdsan as qs
    from math import log, isclose
    from qsdsan.process_models import KineticReaction as KRxn

    gas_kwargs = dict(phase='g', particle_size='Dissolved gas',
                      degradability='Undegradable', organic=False)
    SO2Cl2 = qs.Component('SO2Cl2', **gas_kwargs)
    SO2    = qs.Component('SO2',    **gas_kwargs)
    Cl2    = qs.Component('Cl2',    **gas_kwargs)
    qs.set_thermo(qs.Components((SO2Cl2, SO2, Cl2)))

    # --- first-order reaction ---
    rxn1 = KRxn('SO2Cl2', n=1, k=2.2e-5, t=1e5, reaction='SO2Cl2 -> SO2 + Cl2')
    assert rxn1.n == 1
    assert rxn1.k == 2.2e-5
    assert rxn1.t == 1e5
    assert rxn1.C0 is None  # not yet set

    s1 = qs.Stream('s1_kr', SO2Cl2=100, SO2=10, Cl2=5)
    rxn1(s1)
    assert isclose(rxn1.X, 1 - pow(2.71828, -2.2e-5 * 1e5), rel_tol=1e-3)

    hl = rxn1.half_life
    assert isclose(hl, log(2) / 2.2e-5, rel_tol=1e-6)

    assert rxn1.rate_equation is not None
    assert rxn1.integrated_rate_equation is not None

    # --- property setters re-calculate X when C0 is already set ---
    rxn1.k = 3.0e-5   # triggers k.setter recalc branch
    rxn1.t = 2e5      # triggers t.setter recalc branch
    rxn1.n = 1        # same value; exercises n.setter without changing

    # --- second-order reaction: C0 must be positive at construction for n!=1 ---
    rxn2 = KRxn('SO2Cl2', n=2, k=1e-3, t=100,
                reaction='SO2Cl2 -> SO2 + Cl2', C0=1.0)
    s2 = qs.Stream('s2_kr', SO2Cl2=50, SO2=5, Cl2=5)
    rxn2(s2)   # __call__ updates C0 from stream concentration
    assert 0 < rxn2.X < 1
    assert rxn2.half_life > 0
    assert rxn2.integrated_rate_equation is not None

    # C0 <= 0 raises ValueError (exercises C0.setter error branch)
    with pytest.raises(ValueError, match='C0'):
        rxn2.C0 = -1

    # --- zeroth-order reaction: test properties without stream call ---
    # (stream call requires k*t < C0(stream) to stay under 100% conversion;
    #  we just need to cover the n=0 branch of integrated_rate_equation)
    rxn0 = KRxn('SO2Cl2', n=0, k=5e-4, t=500,
                reaction='SO2Cl2 -> SO2 + Cl2', C0=1.0)
    assert rxn0.integrated_rate_equation is not None  # covers n=0 branch: C0 - k*t
    assert rxn0.half_life > 0

    # half_life guard: n!=1 with falsy C0
    rxn_hl = KRxn('SO2Cl2', n=2, k=1e-3, t=100,
                  reaction='SO2Cl2 -> SO2 + Cl2', C0=2.0)
    rxn_hl._C0 = None  # bypass setter to put C0 in falsy state
    with pytest.raises(ValueError, match='C0'):
        _ = rxn_hl.half_life

    # --- constructor error cases ---
    with pytest.raises(ValueError, match='basis'):
        KRxn('SO2Cl2', n=1, k=1e-3, t=100,
             reaction='SO2Cl2 -> SO2 + Cl2', basis='mass')

    with pytest.raises(ValueError, match='single-phase'):
        KRxn('SO2Cl2', n=1, k=1e-3, t=100,
             reaction='SO2Cl2 -> SO2 + Cl2', phases=('g', 'l'))

    with pytest.raises(ValueError, match='non-negative integer'):
        KRxn('SO2Cl2', n=1.5, k=1e-3, t=100, reaction='SO2Cl2 -> SO2 + Cl2')


if __name__ == '__main__':
    test_process()
    test_kinetic_reaction()