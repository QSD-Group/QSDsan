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

__all__ = (
    'test_aerobic_polishing_filter_no_negative_O2_in_outlets',
    'test_anaerobic_polishing_filter_no_air_injection',
    )


def _make_polishing_filter_cmps():
    '''
    Minimal QSDsan Components for PolishingFilter tests.

    Uses standard chemical IDs (CO2, O2, N2) so that get_combustion_reaction()
    and get_bmp_stoichiometry() can build balanced reaction stoichiometries.
    ``Biomass`` is a custom component not in the thermosteam database, which
    avoids the Cython-level type-slot conflict triggered by known IDs like
    ``WWTsludge``.
    '''
    import qsdsan as qs
    from qsdsan.utils import create_example_components

    ref_cmps = create_example_components()

    def _cmp(ID, particle_size, degradability, organic, **kw):
        c = qs.Component(ID, particle_size=particle_size,
                         degradability=degradability, organic=organic, **kw)
        c.copy_models_from(ref_cmps.H2O, names=('V',))
        c.default()
        return c

    AceticAcid = _cmp('AceticAcid', 'Soluble', 'Readily', True)
    CH4 = _cmp('CH4', 'Dissolved gas', 'Undegradable', False)
    CO2 = _cmp('CO2', 'Dissolved gas', 'Undegradable', False)
    H2O = _cmp('H2O', 'Soluble', 'Undegradable', False)
    O2 = _cmp('O2', 'Dissolved gas', 'Undegradable', False)
    N2 = _cmp('N2', 'Dissolved gas', 'Undegradable', False)
    # 'Biomass' is not in the thermosteam DB so Component creation works cleanly.
    Biomass = qs.Component('Biomass', phase='s', formula='CH1.8O0.5N0.2',
                            particle_size='Particulate', degradability='Slowly',
                            organic=True)
    Biomass.copy_models_from(ref_cmps.H2O, names=('V',))
    Biomass.default()

    cmps = qs.Components([AceticAcid, CH4, CO2, H2O, O2, N2, Biomass])
    cmps.compile()
    qs.set_thermo(cmps)
    return cmps


# Explicit split: solubles (except O2/N2) partially to effluent; biomass to waste.
_SPLIT = {'AceticAcid': 0.125, 'CH4': 0.15, 'CO2': 0.15, 'H2O': 0.125,
          'O2': 1.0, 'N2': 1.0, 'Biomass': 0.0}


def test_aerobic_polishing_filter_no_negative_O2_in_outlets():
    import qsdsan as qs
    qs.main_flowsheet.clear()
    _make_polishing_filter_cmps()

    # High organic load with no O2 in feed forces an O2 deficit after reactions;
    # the fix must zero out negative O2 in `inf` before split_to distributes it.
    influent = qs.WasteStream('pf_influent', AceticAcid=1, H2O=100, units='kg/hr')
    recycle = qs.WasteStream('pf_recycle')
    air = qs.WasteStream('pf_air', phase='g')
    pf = qs.sanunits.PolishingFilter(
        'PF_aerobic',
        ins=(influent, recycle, air),
        outs=('pf_biogas', 'pf_eff', 'pf_waste', 'pf_air_out'),
        filter_type='aerobic',
        biomass_ID='Biomass',
        split=_SPLIT,
        solids=('Biomass',),
    )
    pf.run()

    assert pf.outs[1].imol['O2'] >= -1e-12, (
        f"Effluent O2 is negative ({pf.outs[1].imol['O2']:.4e}): "
        "O2 deficit was not corrected before split_to"
    )
    assert pf.outs[2].imol['O2'] >= -1e-12, (
        f"Waste sludge O2 is negative ({pf.outs[2].imol['O2']:.4e}): "
        "O2 deficit was not corrected before split_to"
    )


def test_anaerobic_polishing_filter_no_air_injection():
    import qsdsan as qs
    qs.main_flowsheet.clear()
    _make_polishing_filter_cmps()

    influent = qs.WasteStream('apf_influent', AceticAcid=1, H2O=100, units='kg/hr')
    recycle = qs.WasteStream('apf_recycle')
    air = qs.WasteStream('apf_air', phase='g')
    pf = qs.sanunits.PolishingFilter(
        'PF_anaerobic',
        ins=(influent, recycle, air),
        outs=('apf_biogas', 'apf_eff', 'apf_waste', 'apf_air_out'),
        filter_type='anaerobic',
        biomass_ID='Biomass',
        split=_SPLIT,
        solids=('Biomass',),
    )
    pf.run()

    # Anaerobic filter must never inject air regardless of O2 balance.
    assert pf.ins[2].F_mol == 0, (
        "Air inlet should remain empty for anaerobic filter"
    )


if __name__ == '__main__':
    test_aerobic_polishing_filter_no_negative_O2_in_outlets()
    test_anaerobic_polishing_filter_no_air_injection()
