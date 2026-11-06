# -*- coding: utf-8 -*-
"""
Self-contained mass-balance tests for QSDsan Junction process-model converters.

Representative feeds exercise each junction standalone via `Junction._run()`;
phosphorus-deficient edge cases also call `junction.reactions()` directly.
There is no reactor, numerical integration, or dependency on EXPOsan. These
checks isolate interface allocation errors from kinetics and solver behavior.
The conserved quantities are COD and TKN for every junction, plus TP for the
phosphorus-aware variants.

The phosphorus regressions are necessary because converting between component
sets with different i_P values can require more P than the source contains.
The forward interface must limit its product allocation; the reverse interface
may dissolve only available minerals and must reject an unrecoverable deficit.
Clipping a substantial phosphate deficit to zero would create P, while reusing
struvite already dissolved for N recovery would overdraw the mineral inventory.
The tests therefore check nonnegative outputs as well as COD/N/P conservation.
Their low-P compositions are deliberate stress cases, not calibration defaults.
"""
import pytest
import numpy as np
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


# %% Representative forward/reverse interface balances
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


# %% Phosphorus-limited forward allocation
@pytest.mark.parametrize('substrate_p,inert_p', [(0.00311, 0.00361),
                                               (0.006149, 0.007139)])
def test_masm2d_adm1p_low_p_allocation(substrate_p, inert_p):
    """Restrict lipid allocation when low-P material cannot supply its P demand.

    Retain these deliberately low structural-biomass i_P values as regression
    inputs, not current plant settings. Reallocating COD to carbohydrate must
    prevent negative S_IP without losing COD, nitrogen, or phosphorus.
    """
    concentrations = dict(_ASM2D_INF)
    concentrations.update(
        S_PO4=3.3505, S_F=79.8232, X_S=36251.9857,
        X_H=22629.0355, X_AUT=1697.3254, X_I=44697.5586,
        )
    inf_asm = _make_stream(
        'inf_masm2d_low_p', pc.create_masm2d_cmps, concentrations,
        )
    cmps_asm = qs.get_components()
    for component_ID, i_N, i_P in (
            ('S_F', 0.023464, substrate_p),
            ('X_S', 0.023464, substrate_p),
            ('S_I', 0.042021, inert_p),
            ('X_I', 0.042021, inert_p),
            ('X_H', 0.060305, 0.002369),
            ('X_PAO', 0.060305, 0.002369),
            ('X_AUT', 0.060305, 0.002369),
            ):
        component = getattr(cmps_asm, component_ID)
        component.i_N = i_N
        component.i_P = i_P
    cmps_asm.refresh_constants()
    asm = pc.mASM2d()

    pc.create_adm1p_cmps()
    adm = pc.ADM1p()
    thermo_adm = qs.get_thermo()
    fwd = su.mASM2dtoADM1p(
        'fwd_adm1p_low_p', upstream=inf_asm, downstream='adm1p_low_p',
        thermo=thermo_adm, isdynamic=False,
        adm1_model=adm, asm2d_model=asm,
        )
    fwd._run()

    eff = fwd.outs[0]
    assert eff.iconc['S_IP'] >= 0.
    assert eff.iconc['X_ch'] > 0.
    assert np.all(np.asarray(eff.conc) >= -1e-8)
    assert_conserved(inf_asm, eff, ('COD', 'TKN', 'TP'), rtol=1e-6)

    qs.main_flowsheet.clear()


# %% Phosphorus-deficient reverse mapping and mineral recovery
@pytest.fixture
def reverse_phosphorus_junction():
    """Build an isolated reverse interface whose organic products require P.

    The mASM2d X_S/S_F i_P values intentionally make carbohydrate/sugar mapping
    demand more P than the dissolved feed supplies. This exercises post-mapping
    recovery, which a balanced forward/reverse round trip need not trigger.
    """
    ca = pc.create_masm2d_cmps()
    ca.X_S.i_P = ca.S_F.i_P = 0.006149
    ca.refresh_constants()
    ta = qs.get_thermo()
    asm = pc.mASM2d()
    cd = pc.create_adm1p_cmps()
    adm = pc.ADM1p()
    stream = qs.WasteStream('reverse_p_feed')
    junction = su.ADM1ptomASM2d(
        'reverse_p', upstream=stream, thermo=ta,
        adm1_model=adm, asm2d_model=asm)
    yield junction, cd, ca
    qs.main_flowsheet.clear()


def _reverse_phosphorus_feed(components):
    """Return a low-dissolved-P feed in the interface's mg/L convention."""
    values = np.zeros(len(components))
    for name, value in dict(X_ch=10000., S_su=100., S_IP=0.1,
                            S_IN=1000., S_IC=10000.).items():
        values[components.index(name)] = value
    return values


@pytest.mark.parametrize('mineral', ['X_struv', 'X_newb', 'X_ACP',
                                    'X_AlPO4', 'X_FePO4'])
def test_post_mapping_p_recovered(reverse_phosphorus_junction, mineral):
    """Recover a mapping-induced P deficit from each supported P mineral.

    The input must remain unchanged and the mapped output must remain
    nonnegative and conserve COD/N/P. This catches checking only the original
    S_IP value instead of the phosphate remaining after organic remapping.
    """
    junction, cd, ca = reverse_phosphorus_junction
    values = _reverse_phosphorus_feed(cd)
    values[cd.index(mineral)] = 100. / getattr(cd, mineral).i_P
    before = values.copy()
    out = junction.reactions(values)
    np.testing.assert_array_equal(values, before)
    assert np.min(out) >= -1e-8
    assert out[ca.index('S_PO4')] == pytest.approx(0., abs=1e-8)
    for prop in ('i_COD', 'i_N', 'i_P'):
        assert out @ getattr(ca, prop) == pytest.approx(
            values @ getattr(cd, prop), rel=1e-9, abs=1e-8)


@pytest.mark.parametrize('mineral_p', [0., 10.])
def test_unrecoverable_p_rejected(reverse_phosphorus_junction, mineral_p):
    """Reject insufficient P inventory instead of silently creating phosphate.

    Both absent and insufficient mineral supplies must raise: clipping the
    negative S_PO4 would hide an infeasible allocation and violate total P.
    """
    junction, cd, ca = reverse_phosphorus_junction
    values = _reverse_phosphorus_feed(cd)
    values[cd.index('X_struv')] = mineral_p / cd.X_struv.i_P
    with pytest.raises(RuntimeError, match='negative S_PO4'):
        junction.reactions(values)


def test_n_then_p_uses_remaining_struvite(reverse_phosphorus_junction):
    """Do not spend the original struvite inventory twice on N and P deficits.

    With no dissolved inorganic N, nitrogen correction dissolves struvite
    first. Subsequent P recovery must use the remaining inventory, keeping
    mineral concentrations nonnegative and conserving both nutrients and COD.
    """
    junction, cd, ca = reverse_phosphorus_junction
    values = _reverse_phosphorus_feed(cd)
    values[cd.index('S_IN')] = 0.
    values[cd.index('X_struv')] = 1000. / cd.X_struv.i_P
    out = junction.reactions(values)
    assert np.min(out) >= -1e-8
    for prop in ('i_COD', 'i_N', 'i_P'):
        assert out @ getattr(ca, prop) == pytest.approx(
            values @ getattr(cd, prop), rel=1e-9, abs=1e-8)
