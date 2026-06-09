# -*- coding: utf-8 -*-
"""
Self-contained mass-balance tests for QSDsan Junction process-model converters.

Each junction is exercised standalone via `Junction._run()`. There is not
reactor/system and no dependency on EXPOsan. We assert that the
conversion conserves the quantities the interface is designed to conserve:
COD and TKN for every junction, plus TP for the phosphorus-aware variants.

These are characterization tests of existing converters: they pass when the
balance holds and fail if a future change breaks it.
"""
import pytest
import qsdsan as qs
from qsdsan import process_models as pc, unit_operations as su, WasteStream

T = 273.15 + 35  # 35 °C, the temperature these converters are written for


def _rel_err(a, b):
    denom = max(abs(a), abs(b))
    return abs(a - b) / denom if denom else 0.0


def assert_conserved(inf, eff, props, rtol):
    """Assert each composite property in `props` is conserved across a junction."""
    errors = {}
    for prop in props:
        a, b = getattr(inf, prop), getattr(eff, prop)
        errors[prop] = (a, b, _rel_err(a, b))
    bad = {p: v for p, v in errors.items() if v[2] > rtol}
    assert not bad, (
        f"Non-conserved quantities (rtol={rtol}): "
        + "; ".join(f"{p}: in={a:.4g} out={b:.4g} rel_err={e:.2e}"
                    for p, (a, b, e) in bad.items())
    )


# Representative influents (literature values, from EXPOsan bsm2 test, inlined here)
_ASM1_INF = dict(
    S_I=28.0665, S_S=48.9526, X_I=10361.7101, X_S=20375.0176, X_BH=10210.0698,
    X_BA=553.2808, X_P=3204.6601, S_O=0.25225, S_NO=1.6871, S_NH=28.9098,
    S_ND=4.6834, X_ND=906.0933, S_ALK=7.1549 * 12,
)
# ASM2d influent (Flores-Alsina et al. 2016, Table 1.1), shared by the ASM2d/mASM2d pairs
_ASM2D_INF = dict(
    S_O2=0, S_F=26.44, S_A=17.66, S_I=27.23, S_NH4=18.58, S_N2=5.07, S_NO3=0.02,
    S_PO4=4.69, S_IC=78.99, X_I=10964.41, X_S=19084.76, X_H=9479.39,
    X_PAO=3862.20, X_PP=450.87, X_PHA=24.64, X_AUT=333.79, S_K=19.79,
    S_Mg=189.87, S_Na=70, S_Cl=1035, S_Ca=300,
)


def _make_stream(name, create_cmps, concentrations):
    # NOTE: cmps_factory() sets the global default thermo as a side effect, so a
    # caller that later builds a second component set must capture qs.get_thermo()
    # for this stream's side BEFORE making that call (see test_asm1_adm1_pair).
    create_cmps()
    ws = WasteStream(name, T=T)
    ws.set_flow_by_concentration(
        flow_tot=178.4674, concentrations=concentrations, units=('m3/d', 'mg/L'))
    return ws


def test_asm1_adm1_pair():
    # Forward: ASM1 -> ADM1
    inf_asm = _make_stream('inf_asm1', pc.create_asm1_cmps, _ASM1_INF)
    thermo_asm = qs.get_thermo()        # capture ASM1 thermo before switching sides
    pc.create_adm1_cmps()
    adm = pc.ADM1()
    thermo_adm = qs.get_thermo()
    fwd = su.ASMtoADM('fwd_asm1', upstream=inf_asm, downstream='adm1_mid',
                      thermo=thermo_adm, isdynamic=False, adm1_model=adm, T=T, pH=7.2631)
    fwd._run()
    assert_conserved(fwd.ins[0], fwd.outs[0], ('COD', 'TKN'), rtol=1e-2)

    # Reverse: ADM1 -> ASM1, fed by the forward junction's (valid ADM1) output.
    # Reuse the captured ASM1 thermo for the downstream side.
    rev = su.ADMtoASM('rev_asm1', upstream=fwd.outs[0], downstream='asm1_out',
                      thermo=thermo_asm, isdynamic=False, adm1_model=adm)
    rev._run()
    assert_conserved(rev.ins[0], rev.outs[0], ('COD', 'TKN'), rtol=1e-2)

    qs.main_flowsheet.clear()


def test_asm2d_adm1_pair():
    # Plain ASM2d uses S_ALK (gC/m3) rather than the S_IC/S_K/S_Mg/S_Na/S_Cl/S_Ca
    # keys that mASM2d carries.  Build the influent from _ASM2D_INF's shared keys
    # and substitute S_ALK for S_IC (same units: gC/m3 expressed as alkalinity).
    _asm2d_inf = dict(
        S_O2=_ASM2D_INF['S_O2'], S_F=_ASM2D_INF['S_F'], S_A=_ASM2D_INF['S_A'],
        S_I=_ASM2D_INF['S_I'], S_NH4=_ASM2D_INF['S_NH4'], S_N2=_ASM2D_INF['S_N2'],
        S_NO3=_ASM2D_INF['S_NO3'], S_PO4=_ASM2D_INF['S_PO4'],
        S_ALK=_ASM2D_INF['S_IC'],  # S_IC (gC/m3) maps to S_ALK in plain ASM2d
        X_I=_ASM2D_INF['X_I'], X_S=_ASM2D_INF['X_S'], X_H=_ASM2D_INF['X_H'],
        X_PAO=_ASM2D_INF['X_PAO'], X_PP=_ASM2D_INF['X_PP'],
        X_PHA=_ASM2D_INF['X_PHA'], X_AUT=_ASM2D_INF['X_AUT'],
    )

    # Forward: ASM2d -> ADM1
    inf_asm = _make_stream('inf_asm2d', pc.create_asm2d_cmps, _asm2d_inf)
    cmps_asm2d = qs.get_components()        # plain ASM2d component set
    thermo_asm = qs.get_thermo()            # capture ASM2d thermo before switching sides
    cmps_adm1 = pc.create_adm1_cmps()
    adm = pc.ADM1()
    # ASM2d X_I.i_N (0.02) differs from ADM1 default (0.06); align AFTER
    # ADM1() construction (which resets the value) to avoid the RuntimeError
    # in Step 5 of the ASM2d->ADM1 interface.  Same fix as EXPOsan's
    # test_adm1_junctions for the ASM1 pair.
    cmps_adm1.X_I.i_N = cmps_asm2d.X_I.i_N
    cmps_adm1.refresh_constants()
    thermo_adm = qs.get_thermo()
    fwd = su.ASM2dtoADM1('fwd_asm2d', upstream=inf_asm, downstream='adm1_mid2',
                         thermo=thermo_adm, isdynamic=False, adm1_model=adm,
                         T=T, pH=7.2631)
    fwd._run()
    # ASM2d carries P, but plain ADM1 does not track it -> only COD/TKN conserve.
    assert_conserved(fwd.ins[0], fwd.outs[0], ('COD', 'TKN'), rtol=1e-2)

    # Reverse: ADM1 -> ASM2d (T/pH read from the upstream stream, like ADMtoASM).
    rev = su.ADM1toASM2d('rev_asm2d', upstream=fwd.outs[0], downstream='asm2d_out',
                         thermo=thermo_asm, isdynamic=False, adm1_model=adm)
    rev._run()
    assert_conserved(rev.ins[0], rev.outs[0], ('COD', 'TKN'), rtol=1e-2)

    qs.main_flowsheet.clear()


def test_asm2d_madm1_pair():
    # Plain ASM2d <-> "mADM1" (ADM1_p_extension): the mADMjunction family
    # (ASM2dtomADM1 / mADM1toASM2d). These hard-code plain ASM2d's 19-component
    # layout, so they take create_asm2d_cmps()/ASM2d() on the ASM side and
    # create_adm1_p_extension_cmps()/ADM1_p_extension() on the ADM side. Plain
    # ASM2d still carries phosphorus (S_PO4, X_PP, X_PHA, X_PAO), so COD/TKN/TP
    # all conserve. ASM2dtomADM1 accepts T/pH kwargs; mADM1toASM2d needs neither.
    _asm2d_inf = dict(
        S_O2=_ASM2D_INF['S_O2'], S_F=_ASM2D_INF['S_F'], S_A=_ASM2D_INF['S_A'],
        S_I=_ASM2D_INF['S_I'], S_NH4=_ASM2D_INF['S_NH4'], S_N2=_ASM2D_INF['S_N2'],
        S_NO3=_ASM2D_INF['S_NO3'], S_PO4=_ASM2D_INF['S_PO4'],
        S_ALK=_ASM2D_INF['S_IC'],  # S_IC (gC/m3) maps to S_ALK in plain ASM2d
        X_I=_ASM2D_INF['X_I'], X_S=_ASM2D_INF['X_S'], X_H=_ASM2D_INF['X_H'],
        X_PAO=_ASM2D_INF['X_PAO'], X_PP=_ASM2D_INF['X_PP'],
        X_PHA=_ASM2D_INF['X_PHA'], X_AUT=_ASM2D_INF['X_AUT'],
    )

    # Forward: ASM2d -> mADM1
    inf_asm = _make_stream('inf_asm2d_m', pc.create_asm2d_cmps, _asm2d_inf)
    cmps_asm = qs.get_components()
    thermo_asm = qs.get_thermo()        # capture ASM2d thermo before switching sides
    asm = pc.ASM2d()
    cmps_adm = pc.create_adm1_p_extension_cmps()
    adm = pc.ADM1_p_extension()
    cmps_adm.X_I.i_N = cmps_asm.X_I.i_N   # align i_N to avoid the converter RuntimeError
    cmps_adm.refresh_constants()
    thermo_adm = qs.get_thermo()
    fwd = su.ASM2dtomADM1('fwd_madm1', upstream=inf_asm, downstream='madm1_mid',
                          thermo=thermo_adm, isdynamic=False,
                          adm1_model=adm, asm2d_model=asm, T=T, pH=7.0)
    fwd._run()
    assert_conserved(fwd.ins[0], fwd.outs[0], ('COD', 'TKN', 'TP'), rtol=1e-2)

    # Reverse: mADM1 -> ASM2d
    rev = su.mADM1toASM2d('rev_madm1', upstream=fwd.outs[0], downstream='asm2d_out_m',
                          thermo=thermo_asm, isdynamic=False,
                          adm1_model=adm, asm2d_model=asm)
    rev._run()
    assert_conserved(rev.ins[0], rev.outs[0], ('COD', 'TKN', 'TP'), rtol=1e-2)

    qs.main_flowsheet.clear()


def test_masm2d_adm1p_pair():
    # The correct junction classes for mASM2d <-> ADM1p (the precipitation-aware
    # extension of ADM1) are mASM2dtoADM1p (forward) and ADM1ptomASM2d (reverse).
    # They implement the A1 algorithm (Flores-Alsina et al. 2016) and handle
    # mASM2d's full 31-component set including mineral precipitation species
    # (X_CaCO3, X_struv, etc.).  The ADM side must use create_adm1p_cmps() /
    # ADM1p() -- not create_adm1_p_extension_cmps() / ADM1_p_extension() -- because
    # the A1 interface explicitly accesses precipitation stoichiometry rows that only
    # ADM1p carries.  The A1junction subclasses require no T or pH arguments (they
    # perform no charge-balance step and map components via i_C/i_N/i_P directly).

    # Forward: mASM2d -> ADM1p
    inf_asm = _make_stream('inf_masm2d_p', pc.create_masm2d_cmps, _ASM2D_INF)
    cmps_asm = qs.get_components()
    thermo_asm = qs.get_thermo()        # capture mASM2d thermo before switching sides
    # Build mASM2d model while the mASM2d thermo is still active
    asm = pc.mASM2d()
    cmps_adm = pc.create_adm1p_cmps()
    adm = pc.ADM1p()
    # Align i_N AFTER ADM model construction (construction can reset the value).
    # A1junction.check_component_properties also performs this alignment, but
    # doing it here pre-emptively makes the intent explicit.
    cmps_adm.X_I.i_N = cmps_asm.X_I.i_N
    cmps_adm.refresh_constants()
    thermo_adm = qs.get_thermo()
    fwd = su.mASM2dtoADM1p('fwd_adm1p', upstream=inf_asm, downstream='adm1p_mid',
                            thermo=thermo_adm, isdynamic=False,
                            adm1_model=adm, asm2d_model=asm)
    fwd._run()
    assert_conserved(fwd.ins[0], fwd.outs[0], ('COD', 'TKN', 'TP'), rtol=1e-2)

    # Reverse: ADM1p -> mASM2d
    rev = su.ADM1ptomASM2d('rev_adm1p', upstream=fwd.outs[0], downstream='masm2d_out_p',
                           thermo=thermo_asm, isdynamic=False,
                           adm1_model=adm, asm2d_model=asm)
    rev._run()
    assert_conserved(rev.ins[0], rev.outs[0], ('COD', 'TKN', 'TP'), rtol=1e-2)

    qs.main_flowsheet.clear()
