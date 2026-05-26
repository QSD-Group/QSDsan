#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

__all__ = (
    'test_get_impact_table_empty_categories_return_message',
    'test_write_lca_sheet_writes_impact_tables',
    'test_lca_save_report_delegates_to_system',
    'test_system_save_report_patched_once',
    'test_get_allocated_impact_table',
    'test_save_report_is_unified',
    'test_save_report_allocation_sheet_opt_in',
    )

import pytest


def _sheet_names(path):
    '''Return the set of worksheet names in an xlsx file.'''
    import openpyxl
    wb = openpyxl.load_workbook(path, read_only=True)
    try:
        return set(wb.sheetnames)
    finally:
        wb.close()


def _build(flowsheet_ID):
    '''A small self-contained SanUnit system with a cost, a construction item,
    and a stream impact item (but no "other" items, so an empty LCA category is
    exercised), plus a TEA and an LCA.'''
    import qsdsan as qs
    from qsdsan import SanUnit, WasteStream, Component, Components
    cmps = Components([Component('Water', search_ID='Water', particle_size='Soluble',
                                 degradability='Undegradable', organic=False)])
    cmps.compile()
    qs.set_thermo(cmps)
    qs.main_flowsheet.set_flowsheet(qs.Flowsheet(flowsheet_ID))

    qs.ImpactIndicator('GlobalWarming', alias='GWP', unit='kg CO2-eq')
    Steel = qs.ImpactItem('Steel', functional_unit='kg', GWP=2.)

    class Treatment(SanUnit):
        _N_ins = 1
        _N_outs = 1
        def _run(self):
            self.outs[0].copy_like(self.ins[0])
        def _cost(self):
            self.baseline_purchase_costs['equipment'] = 1e6  # F_BM defaults to 1

    ws = WasteStream('ws', Water=1000)
    u = Treatment('u', ins=ws, outs='e')
    u.construction = [qs.Construction(linked_unit=u, item=Steel,
                                      quantity=1000, quantity_unit='kg')]
    qs.StreamImpactItem(linked_stream=ws, GWP=1.)

    sys = qs.System('sys', path=(u,))
    sys.simulate()
    tea = qs.TEA(sys, discount_rate=0.05, lifetime=10, simulate_system=False)
    lca = qs.LCA(system=sys, lifetime=10, simulate_system=False)
    return sys, tea, lca


# sheets scoped to those a minimal system can produce: include 'Flowsheet'
# (BioSTEAM leaves an internal flag unbound if it is omitted) and exclude
# 'Water mass balance' (needs a ProcessWaterCenter facility).
_SHEETS = {'Flowsheet', 'Itemized costs', 'Stream table', 'Design requirements'}


def _skip_if_biosteam_report_broken(sys, tmp_path):
    '''The installed BioSTEAM's ``System.save_report`` has pandas-3.0 issues
    (positional ``to_excel``); skip the cross-package integration test when it
    cannot run here. The QSDsan-side behavior is covered by the other tests.'''
    original = type(sys).save_report.__wrapped__   # unwrapped BioSTEAM method
    try:
        original(sys, str(tmp_path / 'probe.xlsx'), sheets=set(_SHEETS))
    except Exception as e:  # pragma: no cover - depends on installed BioSTEAM
        pytest.skip(f'installed BioSTEAM save_report cannot run under this '
                    f'pandas ({type(e).__name__}: {e}); see the BioSTEAM issues '
                    f'note. QSDsan-side unification is covered by the other tests.')


# %%

def test_get_impact_table_empty_categories_return_message():
    '''Empty Stream/Other categories return a message string (consistent with
    Construction/Transportation) instead of crashing on pandas >= 3.0.'''
    sys, tea, lca = _build('test_lca_empty_cat')
    # this system has construction and a stream item, but no transportation/other
    assert hasattr(lca.get_impact_table('Construction'), 'to_dict')  # non-empty -> DataFrame
    assert hasattr(lca.get_impact_table('Stream'), 'to_dict')        # non-empty -> DataFrame
    assert lca.get_impact_table('Transportation') == 'No transportation-related impacts.'
    assert lca.get_impact_table('Other') == 'No other-related impacts.'


def test_write_lca_sheet_writes_impact_tables(tmp_path):
    '''The shared helper appends the LCA tables (with impact columns) to an
    already-open writer, skipping empty categories.'''
    import pandas as pd
    import qsdsan._lca as _lca
    sys, tea, lca = _build('test_lca_helper')
    f = tmp_path / 'r.xlsx'
    with pd.ExcelWriter(f, engine='openpyxl') as writer:
        pd.DataFrame({'a': [1]}).to_excel(writer, sheet_name='Stream table')
    with pd.ExcelWriter(f, mode='a', engine='openpyxl', if_sheet_exists='replace') as writer:
        _lca._write_lca_sheet(lca, writer, sheet_name='LCA')
    assert 'LCA' in _sheet_names(f)
    raw = pd.read_excel(f, sheet_name='LCA', header=None).astype(str)
    assert raw.apply(lambda col: col.str.contains('GlobalWarming').any()).any()


def test_lca_save_report_delegates_to_system(monkeypatch):
    '''LCA.save_report forwards to System.save_report with lca_sheet_name/annual.'''
    sys, tea, lca = _build('test_lca_delegate')
    calls = {}
    def spy(self, *args, **kwargs):
        calls['args'] = args
        calls['kwargs'] = kwargs
    monkeypatch.setattr(type(sys), 'save_report', spy)
    lca.save_report('out.xlsx', annual=True)
    assert calls['args'][0] == 'out.xlsx'
    assert calls['kwargs'].get('lca_sheet_name') == 'LCA'
    assert calls['kwargs'].get('annual') is True


def test_system_save_report_patched_once():
    '''The wrapper marks itself and is not applied twice on reload.'''
    import importlib
    import qsdsan as qs
    assert getattr(qs.System.save_report, '_qsdsan_unified', False) is True
    patched = qs.System.save_report
    importlib.reload(qs._lca)
    assert qs.System.save_report is patched   # sentinel guard prevents re-wrapping


def test_get_allocated_impact_table():
    '''The allocation table reuses get_allocated_impacts; factors sum to 1, and
    a single stream returns a message.'''
    sys, tea, lca = _build('test_lca_alloc')
    streams = list(sys.streams)
    df = lca.get_allocated_impact_table(streams=streams, allocate_by='mass')
    assert 'Allocation factor' in df.columns
    assert df['Allocation factor'].sum() == pytest.approx(1.0)
    assert df.shape[0] == len(streams)
    assert isinstance(lca.get_allocated_impact_table(streams=streams[:1]), str)


@pytest.mark.parametrize('entry', ['system', 'tea', 'lca'])
def test_save_report_is_unified(tmp_path, entry):
    '''System, TEA, and LCA all produce the same unified workbook
    (process/TEA sheets plus the QSDsan LCA sheet). Skipped if the installed
    BioSTEAM's report cannot run under the current pandas.'''
    sys, tea, lca = _build(f'test_lca_unified_{entry}')
    _skip_if_biosteam_report_broken(sys, tmp_path)
    obj = {'system': sys, 'tea': tea, 'lca': lca}[entry]
    f = tmp_path / 'report.xlsx'
    obj.save_report(str(f), sheets=set(_SHEETS))
    names = _sheet_names(f)
    assert 'LCA' in names                # QSDsan LCA tables appended
    assert 'Stream table' in names       # BioSTEAM process sheet
    assert 'Itemized costs' in names     # TEA sheet (system has a TEA)


def test_save_report_allocation_sheet_opt_in(tmp_path):
    '''Passing lca_allocate_streams adds an 'LCA allocation' sheet. Skipped if the
    installed BioSTEAM's report cannot run under the current pandas.'''
    sys, tea, lca = _build('test_lca_alloc_sheet')
    _skip_if_biosteam_report_broken(sys, tmp_path)
    f = tmp_path / 'report.xlsx'
    sys.save_report(str(f), sheets=set(_SHEETS),
                    lca_allocate_streams=list(sys.streams), lca_allocate_by='mass')
    names = _sheet_names(f)
    assert 'LCA' in names
    assert 'LCA allocation' in names


if __name__ == '__main__':
    test_get_impact_table_empty_categories_return_message()
    test_system_save_report_patched_once()
    print('non-fixture tests ok (run via pytest for the rest)')
