#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:

    Joy Zhang <joycheung1994@gmail.com>

    Yalin Li <mailto.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/master/LICENSE.txt
for license details.
'''

__all__ = ('test_component',)

def test_component():
    import pytest
    from qsdsan import Chemical, Component, Components, set_thermo, \
        _waste_stream as ws_module
    from chemicals.elements import molecular_weight
    from math import isclose

    S_NH4 = Component('S_NH4', formula='NH4+', measured_as='N',
                 f_BOD5_COD=0, f_uBOD_COD=0, f_Vmass_Totmass=0,
                 description="Ammonium", particle_size="Soluble",
                 degradability="Undegradable", organic=False)
    assert S_NH4.i_N == 1
    assert S_NH4.i_NOD == molecular_weight({'O':4})/molecular_weight({'N':1})
    S_NH4.measured_as = None
    assert S_NH4.i_mass == 1

    S_Ac = Component('S_Ac', formula='CH3COO-', measured_as='COD', f_BOD5_COD=0.717,
                    f_uBOD_COD=0.863, f_Vmass_Totmass=1,
                    description="Acetate", particle_size="Soluble",
                    degradability="Readily", organic=True)
    assert S_Ac.i_COD == 1
    S_Ac.measured_as = None
    assert S_Ac.i_mass == 1
    assert S_Ac.i_COD == molecular_weight({'O':4})/molecular_weight({'C':2, 'H':3, 'O':2})

    S_HS = Component.from_chemical('S_HS', Chemical('Hydrosulfide'),
                                  particle_size="Soluble",
                                  degradability="Undegradable", organic=False)
    assert S_HS.i_charge < 0
    S_HS.measured_as = 'S'
    assert S_HS.i_mass > 1

    # Check default components
    cmps1 = Components.load_default(default_compile=False)
    H2O_chemical = Chemical('H2O')
    H2O = Component.from_chemical('H2O', H2O_chemical)
    with pytest.raises(ValueError): # H2O already in default components
        cmps1.append(H2O)
    with pytest.raises(RuntimeError): # key chemical-related properties missing
        cmps1.compile(ignore_inaccurate_molar_weight=True)
    # Can compile with default-filling those missing properties
    cmps1.default_compile(lock_state_at='', particulate_ref='NaCl', ignore_inaccurate_molar_weight=True)

    cmps2 = Components((cmp for cmp in cmps1 if cmp.ID != 'H2O'))
    H2O = Component.from_chemical('H2O', Chemical('H2O'),
                                  particle_size='Soluble',
                                  degradability='Undegradable', organic=False)
    cmps2.append(H2O)
    cmps2.default_compile(lock_state_at='', particulate_ref='NaCl', ignore_inaccurate_molar_weight=True)

    cmps3 = Components.load_default()
    assert cmps3.S_H2.measured_as == 'COD'
    assert cmps3.S_H2.i_COD == 1
    assert isclose(cmps3.S_NO2.i_COD, -3.4268, rel_tol=1e-3)
    assert isclose(cmps3.S_NO3.i_COD, -4.569, rel_tol=1e-3)
    set_thermo(cmps3)

    # Check if the default groups are up-to-date
    cached_cmp_IDs = ws_module._default_cmp_IDs
    cached_cmp_groups = ws_module._specific_groups
    assert set(cmps3.IDs) == cached_cmp_IDs
    get_IDs = lambda attr: {cmp.ID for cmp in getattr(cmps3, attr)}
    for attr, IDs in cached_cmp_groups.items():
        assert IDs == get_IDs(attr)


import pytest
from qsdsan import Component
from thermosteam import Chemical

_kw = dict(particle_size='Soluble', degradability='Undegradable', organic=False)

def test_search_ID_rejects_chemical_object():
    obj = Chemical('Water')
    with pytest.raises(TypeError):
        Component('W', search_ID=obj, **_kw)

def test_chemical_object_with_search_ID_raises():
    obj = Chemical('Water')
    with pytest.raises(ValueError):
        Component('W', chemical=obj, search_ID='Water', **_kw)

def test_chemical_string_conflicting_search_ID_raises():
    with pytest.raises(ValueError):
        Component('S_su', chemical='glucose', search_ID='fructose', **_kw)

def test_chemical_string_aliases_to_search_ID():
    a = Component('A', chemical='glucose', **_kw)
    b = Component('B', search_ID='glucose', **_kw)
    assert a.formula == b.formula
    assert abs(a.MW - b.MW) < 1e-9

def test_chemical_string_same_search_ID_ok():
    c = Component('C', chemical='glucose', search_ID='glucose', **_kw)
    assert c.formula == Chemical('glucose').formula


from qsdsan import Components, set_thermo, WasteStream

def _mix_H(comp_id, **kw):
    """Build a forced-phase component, put it in a stream, return stream.H."""
    H2O = Component.from_chemical('H2O', phase='l', **_kw)
    X = Component(comp_id, **kw)
    for c in (H2O, X):
        c.copy_models_from(H2O, names=['V']); c.default()
    set_thermo(Components((H2O, X)))
    s = WasteStream('s'); s.set_flow_by_concentration(0.5, {comp_id: 1620}, units=('L/hr', 'mg/L'))
    return s.H

def test_potassium_forced_solid_does_not_crash():
    assert _mix_H('K+', phase='s', **_kw) == 0

def test_search_id_equivalent_to_from_chemical_phase_override():
    a = Component('Na_s', search_ID='Na+', phase='s', **_kw)
    b = Component.from_chemical('Na_f', chemical='Na+', phase='s', **_kw)
    assert a.formula == b.formula
    assert abs(a.MW - b.MW) < 1e-9
    assert type(a.H).__name__ == type(b.H).__name__ == 'Enthalpy'

def test_bare_id_miss_creates_blank_component():
    c = Component('S_F', formula='C5H7O2N', particle_size='Particulate',
                  degradability='Slowly', organic=True)
    assert c.formula == 'C5H7O2N'
    # Organic component built from a formula: thermosteam estimates Hf from the
    # composition, and the inorganic guard keeps it (organic estimates are
    # in-domain for the Dulong/Boie correlations).
    assert c.Hf is not None

def test_blank_custom_component_with_phase_is_locked():
    c = Component('X_solid', formula='C5H7O2N', phase='s', particle_size='Particulate',
                  degradability='Slowly', organic=True)
    assert c.locked_state == 's'

def test_explicit_search_id_miss_raises():
    with pytest.raises(LookupError):
        Component('Bad', search_ID='definitely_not_a_chemical_xyz', **_kw)

from math import isclose

def test_component_copy_matches_thermosteam_copy():
    """Path A hardcodes a mimic of thermosteam Chemical.copy(); this guards against
    upstream drift. Compare thermodynamic state, not bookkeeping slots."""
    src = Chemical('Ethanol')
    tmo_copy = src.copy('Eth_tmo')                       # thermosteam's own copy
    qs_comp = Component('Eth_qs', chemical=src, **_kw)   # QSDsan Path A
    assert qs_comp.formula == tmo_copy.formula
    assert isclose(qs_comp.MW, tmo_copy.MW, rel_tol=1e-9)
    P = 101325.
    for phase in ('l', 'g'):
        for T in (300., 350.):
            assert isclose(qs_comp.H(phase, T, P), tmo_copy.H(phase, T, P), rel_tol=1e-6), \
                f'enthalpy drift at phase={phase}, T={T}'
            assert isclose(qs_comp.S(phase, T, P), tmo_copy.S(phase, T, P), rel_tol=1e-6), \
                f'entropy drift at phase={phase}, T={T}'


def test_component_copy_roundtrips():
    # Component.copy() calls __new__ with search_db=False (via **chemical_properties);
    # this must not collide with the constructor's internal search_db handling.
    base = Component('Glu', search_ID='glucose', **_kw)
    cp = base.copy('Glu2')
    assert cp.ID == 'Glu2'
    assert cp.formula == base.formula
    assert abs(cp.MW - base.MW) < 1e-9


# --- formula override gate (formula_override) ---
# These use the MgNH4PO4*6H2O hexahydrate to exercise formula_override. The
# mass-balanced recovery component set (qsdsan.utils.doc_examples) deliberately
# uses anhydrous NH4MgPO4 instead; the two representations are intentional.
_pkw = dict(particle_size='Particulate', degradability='Undegradable', organic=False)

def test_formula_mismatch_raises_without_override():
    # hexahydrate formula vs the anhydrous database chemical -> atoms differ -> raise
    with pytest.raises(ValueError):
        Component('Struvite', search_ID='MagnesiumAmmoniumPhosphate',
                  formula='NH4MgPO4·H12O6', phase='s', **_pkw)

def test_formula_override_allows_mismatch_and_recomputes_MW():
    c = Component('Struvite', search_ID='MagnesiumAmmoniumPhosphate',
                  formula='NH4MgPO4·H12O6', phase='s',
                  formula_override=True, **_pkw)
    assert c.formula == 'NH4MgPO4·H12O6'
    # MW must reflect the override formula (hexahydrate ~245), not the DB value (~137)
    assert c.MW > 240 and c.chem_MW > 240

def test_formula_respelling_needs_no_override():
    # same atoms, different string -> not a mismatch -> no exception
    c = Component('Propane', search_ID='Propane', formula='CH3CH2CH3', phase='g', **_pkw)
    assert get_atoms_equal(c.atoms, 'C3H8')

def test_from_chemical_mismatch_raises_without_override():
    with pytest.raises(ValueError):
        Component.from_chemical('Struvite', chemical='MagnesiumAmmoniumPhosphate',
                                formula='NH4MgPO4·H12O6', phase='s', **_pkw)

def test_from_chemical_override_works():
    c = Component.from_chemical('Struvite', chemical='MagnesiumAmmoniumPhosphate',
                                formula='NH4MgPO4·H12O6', phase='s',
                                formula_override=True, **_pkw)
    assert c.MW > 240

def test_blank_component_formula_not_gated():
    # no database chemical -> no source formula to conflict with -> no override needed
    c = Component('S_F', formula='C5H7O2N', particle_size='Particulate',
                  degradability='Slowly', organic=True)
    assert c.formula == 'C5H7O2N'

def get_atoms_equal(atoms_dict, formula):
    from chemicals.elements import get_atoms
    return dict(atoms_dict) == get_atoms(formula)


if __name__ == '__main__':
    test_component()