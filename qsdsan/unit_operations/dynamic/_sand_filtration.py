
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Zixuan Wang <wyatt4428@gmail.com>

Examples
--------
These examples use the standalone ESAP validation system.

Constant splitting
~~~~~~~~~~~~~~~~~~
The continuous model sends fixed fractions of liquid and particulate solids to
the backwash stream::

    from exposan.esap.systems import create_sand_filter_validation_system

    sys = create_sand_filter_validation_system(
        wastewater_strength='mid',
        model='constant_splitting',
    )
    sys.simulate(
        method='BDF', t_span=(0, 300),
        rtol=1e-3, atol=1e-2,
    )
    sys.diagram()

    SF = sys.flowsheet.unit.SF
    influent = SF.ins[0]
    effluent, backwash = SF.outs
    particulate_IDs = tuple(
        cmp.ID for cmp, x in zip(SF.components, SF.components.x) if x > 0
    )

    def solids(stream):
        return stream.imass[particulate_IDs].sum()

    solids_in = solids(influent)
    actual_solids_fraction = solids(backwash) / solids_in
    actual_liquid_fraction = (
        backwash.imass['H2O'] / influent.imass['H2O']
    )

    print(f'Expected solids fraction: {SF.backwash_solids_mass_fraction:.4f}')
    print(f'Actual solids fraction:   {actual_solids_fraction:.4f}')
    print(f'Expected liquid fraction: {SF.backwash_flow_fraction:.4f}')
    print(f'Actual liquid fraction:   {actual_liquid_fraction:.4f}')

Empirical fouling
~~~~~~~~~~~~~~~~~
The mass-balance model varies particulate removal over the filtration cycle,
retains captured solids, and releases that inventory with the backwash::

    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    from exposan.esap.systems import create_sand_filter_validation_system

    sys = create_sand_filter_validation_system(
        wastewater_strength='mid',
        model='empirical_fouling',
    )
    sys.simulate(
        method='BDF', t_span=(0, 300),
        rtol=1e-3, atol=1e-2,
    )
    sys.diagram()

    SF = sys.flowsheet.unit.SF
    influent = SF.ins[0]
    backwash_in = SF.ins[1]
    effluent = SF.outs[0]
    backwash_out = SF.outs[1]
    recycle = SF.outs[2]
    particulate_IDs = tuple(
        cmp.ID for cmp, x in zip(SF.components, SF.components.x) if x > 0
    )

    def solids(stream):
        return stream.imass[particulate_IDs].sum()

    def flow(stream):
        # Dynamic-state flow is in m3/d; F_vol is sufficient for static results.
        return stream.state[-1] if stream.state is not None else stream.F_vol

    solids_in = solids(influent)
    t_final = sys.scope.sol.t[-1]
    expected_capture = SF.get_removal_efficiency(t_final)

    # Effluent and recycle together comprise the filtered liquid.
    filtered_solids = solids(effluent) + solids(recycle)
    actual_capture = 1 - filtered_solids / solids_in

    # Backwash solids above its inlet loading are released filter inventory.
    released_solids = solids(backwash_out) - solids(backwash_in)
    released_fraction = released_solids / solids_in
    expected_liquid_fraction = SF.average_backwash_fraction
    actual_recycle_fraction = flow(recycle) / flow(influent)
    actual_backwash_fraction = flow(backwash_out) / flow(influent)

    print(f'Cycle phase: {SF.get_cycle_phase(t_final):.3f} h')
    print(f'Expected removal: {expected_capture:.4f}')
    print(f'Actual removal:   {actual_capture:.4f}')
    print(f'Released solids/input solids: {released_fraction:.4f}')
    print()
    print(f'Expected recycle fraction:     {expected_liquid_fraction:.4f}')
    print(f'Actual recycle fraction:       {actual_recycle_fraction:.4f}')
    print(f'Actual backwash flow fraction: {actual_backwash_fraction:.4f}')
    print(f'Retained-solids inventory: '
          f'{SF.captured_solids_inventory.sum():.3f} g')

    # The validation system connects the recycle outlet directly to the
    # backwash inlet. Recycle and effluent have identical concentrations but
    # different flow rates; backwash out is enriched with released solids.
    print('Recycle outlet is backwash inlet:', recycle is backwash_in)
    print('Recycle and effluent concentrations identical:',
          np.allclose(recycle.conc, effluent.conc))
    composition = pd.DataFrame(
        {
            'Effluent [mg/L]': effluent.conc,
            'Recycle/backwash in [mg/L]': recycle.conc,
            'Backwash out [mg/L]': backwash_out.conc,
        },
        index=SF.components.IDs,
    )
    print(composition.loc[(composition != 0).any(axis=1)])

    # Plot the prescribed removal efficiency over 25 complete cycles.
    n_cycles = 25
    cycle_days = SF.filter_cycle_time / 24
    t = np.linspace(0, n_cycles * cycle_days, n_cycles * 200 + 1)
    removal = np.array([SF.get_removal_efficiency(ti) for ti in t])

    fig, ax = plt.subplots(figsize=(12, 4))
    ax.plot(t, 100 * removal, color='black', linewidth=1.5)
    for cycle in range(n_cycles + 1):
        ax.axvline(
            cycle * cycle_days,
            color='0.8', linewidth=0.5, zorder=0,
        )
    ax.axhline(
        100 * SF.minimum_removal_efficiency,
        color='red', linestyle='--', linewidth=1,
        label='Diagnostic minimum',
    )
    ax.set(
        xlabel='Time [d]',
        ylabel='Particulate removal efficiency [%]',
        xlim=(0, n_cycles * cycle_days),
        ylim=(0, 100),
    )
    ax.legend(frameon=False)
    plt.tight_layout()
    plt.show()

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

import numpy as np

from ... import SanUnit

__all__ = ('SandFiltration',)

class SandFiltration(SanUnit):
    '''
    Zero-volume sand filtration with constant-splitting and empirical-fouling
    models.

    Parameters
    ----------
    ID : str, optional
        Unit ID.
    ins : sequence[WasteStream], optional
        For ``model='constant_splitting'``, one influent is required. For
        ``model='empirical_fouling'``, the streams are filter influent and
        continuous-average backwash water, respectively.
    outs : sequence[WasteStream], optional
        For ``model='constant_splitting'``, the streams are filtered effluent
        and backwash. For ``model='empirical_fouling'``, the streams are
        discharged effluent, backwash out, and filtered-water recycle.
    model : {'constant_splitting', 'empirical_fouling'}, optional
        Sand-filter model. ``'constant_splitting'`` splits influent into effluent
        and backwash out; ``'empirical_fouling'`` accounts for backwash in and 
        includes an exponential decline of removal effciency in the sand filter 
        until backwash. 
    backwash_flow_fraction : float, optional
        Fraction of influent flow discharged as continuous backwash by the
        constant-splitting model. Must be within [0, 1]. The GPS-X Figure 9-1
        default is 0.02.
    backwash_solids_mass_fraction : float, optional
        Fraction of influent particulate mass discharged in continuous
        backwash by the constant-splitting model. Must be within [0, 1]. The
        GPS-X Figure 9-1 default is 0.8.
    filter_cycle_time : float, optional
        Complete empirical-fouling cycle duration [h]. The GPS-X Figure 9-2
        default is 24 h.
    filter_backwash_time : float, optional
        Backwash/regeneration duration at the end of each empirical-fouling
        cycle [h]. It must be shorter than ``filter_cycle_time``. To avoid
        intermittent hydraulics, its ratio to cycle time sets the continuous-
        average filtered-water recycle fraction. It is also the retained-solids
        washout time constant. The GPS-X Figure 9-2 default is 1 h.
    timeshift_for_cycles : float, optional
        Shift applied to the empirical cycle clock [h]. The GPS-X Figure 9-2
        default is 0 h.
    fouling_coefficient : float, optional
        Exponential fouling coefficient in GPS-X Equation 9-5, applied to
        elapsed cycle time in hours. The GPS-X Figure 9-2 default is 0.03.
    best_removal_efficiency : float, optional
        Particulate removal efficiency immediately after backwash, expressed as
        a fraction. The GPS-X Figure 9-2 default is 0.75 (75%).
    minimum_removal_efficiency : float, optional
        Diagnostic minimum particulate removal efficiency, expressed as a
        fraction. As in GPS-X, this value does not clamp the calculated
        efficiency. Use ``below_minimum_removal_efficiency`` to inspect the
        current condition. The Figure 9-2 default is 0.30 (30%).
    energy_consumption : float, optional
        Specific electricity consumption [kWh/m3 of filter influent]. It is
        stored for later energy accounting and does not alter stream balances.
        The default is 0 because Figure 9-1/9-2 specifies no energy value.
    isdynamic : bool, optional
        Whether the unit participates in dynamic simulation. The default is
        True.

    Notes
    -----
    The empirical model applies Equation 9-5 to particulate component
    concentrations. Soluble components pass unchanged. During the backwash
    window, a cubic smoothstep restores the removal efficiency to its best
    value, avoiding a concentration jump at the cycle boundary. Captured
    particulate mass is retained as a dynamic inventory and released through
    the continuous-average backwash stream. This preserves component mass
    balances without reproducing intermittent pump flows.

    The empirical recycle flow is ``Q_in * filter_backwash_time /
    filter_cycle_time``. This closes the hydraulics after omitting the two GPS-X
    pump-rate parameters requested by the user; it assumes the instantaneous
    backwash hydraulic rate equals the filter loading rate. The backwash inlet
    may instead be supplied externally, but connecting the recycle outlet to it
    gives the intended closed loop.

    References
    ----------
    Hydromantis Environmental Software Solutions, Inc. GPS-X Technical
    Reference, Chapter 9: Sand Filtration Models, version 8.5.
    https://www.hydromantis.com/help/GPS-X/docs/8.5/Technical/9_CHAPTER9.htm
    '''

    _N_ins = 1
    _N_outs = 2
    _ins_size_is_fixed = True
    _outs_size_is_fixed = True

    _model_aliases = {
        'constant_splitting': 'constant_splitting',
        'constant_spliting': 'constant_splitting',
        'continuous': 'constant_splitting',
        'empirical_fouling': 'empirical_fouling',
        'mass_balance': 'empirical_fouling',
        'massbalance': 'empirical_fouling',
        }

    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 model='constant_splitting',
                 backwash_flow_fraction=0.02,
                 backwash_solids_mass_fraction=0.8,
                 filter_cycle_time=24., filter_backwash_time=1.,
                 timeshift_for_cycles=0., fouling_coefficient=0.03,
                 best_removal_efficiency=0.75,
                 minimum_removal_efficiency=0.30,
                 energy_consumption=0., isdynamic=True,
                 init_with='WasteStream', F_BM_default=None, **kwargs):
        mode = self._normalize_model(model)
        self._model = mode
        if mode == 'empirical_fouling':
            self._N_ins = 2
            self._N_outs = 3
        SanUnit.__init__(
            self, ID, ins, outs, thermo, init_with=init_with,
            F_BM_default=F_BM_default, isdynamic=isdynamic, **kwargs)

        self.backwash_flow_fraction = backwash_flow_fraction
        self.backwash_solids_mass_fraction = backwash_solids_mass_fraction
        self.filter_cycle_time = filter_cycle_time
        self.filter_backwash_time = filter_backwash_time
        self.timeshift_for_cycles = timeshift_for_cycles
        self.fouling_coefficient = fouling_coefficient
        self.best_removal_efficiency = best_removal_efficiency
        self.minimum_removal_efficiency = minimum_removal_efficiency
        self.energy_consumption = energy_consumption
        self._validate_constant_splitting()
        self._validate_cycle_times()

        n = len(self.components)
        self._initial_captured_solids = np.zeros(n)
        self._release_rate = np.zeros(n)
        self._t = 0.
        if mode == 'empirical_fouling':
            self._state_header = [
                f'retained {cmp.ID} [g]' for cmp in self.components]

    @classmethod
    def _normalize_model(cls, model):
        key = str(model).strip().lower()
        try:
            return cls._model_aliases[key]
        except KeyError:
            valid = "'constant_splitting' or 'empirical_fouling'"
            raise ValueError(f'model must be {valid}, not {model!r}') from None

    @staticmethod
    def _fraction(value, name):
        value = float(value)
        if not np.isfinite(value) or not 0. <= value <= 1.:
            raise ValueError(f'{name} must be finite and within [0, 1]')
        return value

    @property
    def model(self):
        '''[str] Canonical model name selected during construction.'''
        return self._model

    @property
    def hasode(self):
        '''Whether retained-solids inventories require ODE integration.'''
        return self._model == 'empirical_fouling'

    @property
    def backwash_flow_fraction(self):
        '''[float] Influent flow fraction sent to continuous backwash.'''
        return self._backwash_flow_fraction

    @backwash_flow_fraction.setter
    def backwash_flow_fraction(self, value):
        self._backwash_flow_fraction = self._fraction(
            value, 'backwash_flow_fraction')

    @property
    def backwash_solids_mass_fraction(self):
        '''[float] Influent particulate mass fraction sent to backwash.'''
        return self._backwash_solids_mass_fraction

    @backwash_solids_mass_fraction.setter
    def backwash_solids_mass_fraction(self, value):
        self._backwash_solids_mass_fraction = self._fraction(
            value, 'backwash_solids_mass_fraction')

    @property
    def filter_cycle_time(self):
        '''[float] Complete empirical filter cycle duration [h].'''
        return self._filter_cycle_time

    @filter_cycle_time.setter
    def filter_cycle_time(self, value):
        value = float(value)
        if not np.isfinite(value) or value <= 0.:
            raise ValueError('filter_cycle_time must be finite and positive')
        if hasattr(self, '_filter_backwash_time') \
                and value <= self._filter_backwash_time:
            raise ValueError('filter_cycle_time must exceed filter_backwash_time')
        self._filter_cycle_time = value

    @property
    def filter_backwash_time(self):
        '''[float] Backwash/regeneration duration [h].'''
        return self._filter_backwash_time

    @filter_backwash_time.setter
    def filter_backwash_time(self, value):
        value = float(value)
        if not np.isfinite(value) or value <= 0.:
            raise ValueError(
                'filter_backwash_time must be finite and positive')
        if hasattr(self, '_filter_cycle_time') \
                and value >= self._filter_cycle_time:
            raise ValueError('filter_backwash_time must be less than filter_cycle_time')
        self._filter_backwash_time = value

    @property
    def timeshift_for_cycles(self):
        '''[float] Shift of the empirical filter cycle clock [h].'''
        return self._timeshift_for_cycles

    @timeshift_for_cycles.setter
    def timeshift_for_cycles(self, value):
        value = float(value)
        if not np.isfinite(value):
            raise ValueError('timeshift_for_cycles must be finite')
        self._timeshift_for_cycles = value

    @property
    def fouling_coefficient(self):
        '''[float] Exponential fouling coefficient applied on an hourly basis.'''
        return self._fouling_coefficient

    @fouling_coefficient.setter
    def fouling_coefficient(self, value):
        value = float(value)
        if not np.isfinite(value) or value < 0.:
            raise ValueError('fouling_coefficient must be finite and nonnegative')
        self._fouling_coefficient = value

    @property
    def best_removal_efficiency(self):
        '''[float] Particulate removal efficiency immediately after backwash.'''
        return self._best_removal_efficiency

    @best_removal_efficiency.setter
    def best_removal_efficiency(self, value):
        value = self._fraction(value, 'best_removal_efficiency')
        if hasattr(self, '_minimum_removal_efficiency') \
                and value < self._minimum_removal_efficiency:
            raise ValueError('best removal efficiency cannot be below minimum')
        self._best_removal_efficiency = value

    @property
    def minimum_removal_efficiency(self):
        '''[float] Diagnostic minimum particulate removal efficiency.'''
        return self._minimum_removal_efficiency

    @minimum_removal_efficiency.setter
    def minimum_removal_efficiency(self, value):
        value = self._fraction(value, 'minimum_removal_efficiency')
        if hasattr(self, '_best_removal_efficiency') \
                and value > self._best_removal_efficiency:
            raise ValueError('minimum removal efficiency cannot exceed best')
        self._minimum_removal_efficiency = value

    @property
    def energy_consumption(self):
        '''[float] Specific electricity consumption [kWh/m3 influent].'''
        return self._energy_consumption

    @energy_consumption.setter
    def energy_consumption(self, value):
        value = float(value)
        if not np.isfinite(value) or value < 0.:
            raise ValueError('energy_consumption must be finite and nonnegative')
        self._energy_consumption = value

    @property
    def average_backwash_fraction(self):
        '''[float] Continuous-average recycle fraction for empirical fouling.'''
        return self.filter_backwash_time / self.filter_cycle_time

    @property
    def captured_solids_inventory(self):
        '''[numpy.ndarray] Retained component masses [g].'''
        if self._model != 'empirical_fouling':
            return np.zeros(len(self.components))
        state = self._state
        return (self._initial_captured_solids if state is None else state).copy()

    @property
    def below_minimum_removal_efficiency(self):
        '''[bool] Whether current empirical removal is below its diagnostic limit.'''
        if self._model != 'empirical_fouling':
            return False
        return self.get_removal_efficiency(self._t) \
            < self.minimum_removal_efficiency

    def _validate_constant_splitting(self):
        fq = self.backwash_flow_fraction
        fx = self.backwash_solids_mass_fraction
        if fq == 0. and fx != 0.:
            raise ValueError('nonzero backwash solids require backwash flow')
        if fq == 1. and fx != 1.:
            raise ValueError('zero effluent flow requires all solids in backwash')

    def _validate_cycle_times(self):
        if self.filter_backwash_time >= self.filter_cycle_time:
            raise ValueError(
                'filter_backwash_time must be less than filter_cycle_time')

    def set_init_captured_solids(self, **component_masses):
        '''Set initial retained component masses [g] for empirical fouling.'''
        if self._model != 'empirical_fouling':
            raise RuntimeError(
                'captured-solids inventory exists only for empirical_fouling')
        masses = np.zeros(len(self.components))
        for ID, value in component_masses.items():
            value = float(value)
            if not np.isfinite(value) or value < 0.:
                raise ValueError(
                    'initial captured-solids masses must be finite and nonnegative')
            masses[self.components.index(ID)] = value
        self._initial_captured_solids = masses

    def get_cycle_phase(self, t=None):
        '''Return empirical cycle phase [h] at QSDsan time ``t`` [d].'''
        if t is None:
            t = self._t
        return (24.*float(t) + self.timeshift_for_cycles) \
            % self.filter_cycle_time

    def _empirical_efficiency_and_derivative(self, t):
        phase = self.get_cycle_phase(t)
        cycle = self.filter_cycle_time
        backwash = self.filter_backwash_time
        filtration = cycle - backwash
        best = self.best_removal_efficiency
        clean_ratio = 1. - best
        k = self.fouling_coefficient

        end_ratio = clean_ratio * np.exp(k*filtration)
        end_efficiency = 1. - end_ratio
        if phase < filtration:
            ratio = clean_ratio * np.exp(k*phase)
            efficiency = 1. - ratio
            derivative = -24.*k*ratio
        else:
            u = (phase-filtration)/backwash
            smoothstep = u*u*(3.-2.*u)
            efficiency = end_efficiency \
                + (best-end_efficiency)*smoothstep
            derivative = (best-end_efficiency) * 6.*u*(1.-u) \
                * 24./backwash

        if efficiency <= 0.:
            return 0., 0.
        if efficiency >= 1.:
            return 1., 0.
        return efficiency, derivative

    def get_removal_efficiency(self, t=None):
        '''Return particulate removal efficiency at QSDsan time ``t`` [d].'''
        if self._model == 'constant_splitting':
            return self.backwash_solids_mass_fraction
        if t is None:
            t = self._t
        return self._empirical_efficiency_and_derivative(float(t))[0]

    def _constant_split_factors(self):
        self._validate_constant_splitting()
        x = self.components.x
        fq = self.backwash_flow_fraction
        fx = self.backwash_solids_mass_fraction
        component_to_backwash = (1.-x)*fq + x*fx
        if fq:
            backwash_C = component_to_backwash/fq
        else:
            backwash_C = np.zeros_like(x)
        if fq < 1.:
            effluent_C = (1.-component_to_backwash)/(1.-fq)
        else:
            effluent_C = np.zeros_like(x)
        return (
            np.append(effluent_C, 1.-fq),
            np.append(backwash_C, fq),
            )

    @staticmethod
    def _assign_state(stream, state):
        if stream.state is None:
            stream.state = np.array(state, dtype=float, copy=True)
        else:
            stream.state[:] = state

    def _run(self):
        if self._model == 'constant_splitting':
            influent, = self.ins
            effluent, backwash = self.outs
            x = self.components.x
            fq = self.backwash_flow_fraction
            fx = self.backwash_solids_mass_fraction
            split = (1.-x)*fq + x*fx
            influent.split_to(backwash, effluent, split)
            return

        influent, backwash_in = self.ins
        effluent, backwash_out, recycle = self.outs
        efficiency = self.get_removal_efficiency(0.)
        captured = influent.mol * self.components.x * efficiency
        filtered = influent.mol - captured
        recycle_fraction = self.average_backwash_fraction

        effluent.mol[:] = filtered*(1.-recycle_fraction)
        recycle.mol[:] = filtered*recycle_fraction
        backwash_out.mol[:] = backwash_in.mol + captured
        effluent.thermal_condition.copy_like(influent.thermal_condition)
        recycle.thermal_condition.copy_like(influent.thermal_condition)
        backwash_out.thermal_condition.copy_like(
            backwash_in.thermal_condition)

    def _init_state(self):
        if self._model == 'constant_splitting':
            self._state = self._ins_QC[0].copy()
        else:
            self._state = self._initial_captured_solids.copy()
        self._dstate = np.zeros_like(self._state)
        self._t = 0.

    def _update_state(self):
        if self._model == 'constant_splitting':
            effluent_factor, backwash_factor = self._constant_split_factors()
            effluent, backwash = self.outs
            self._assign_state(effluent, self._state*effluent_factor)
            self._assign_state(backwash, self._state*backwash_factor)
            return

        influent = self._ins_QC[0]
        backwash_in = self._ins_QC[1]
        C_in, Q_in = influent[:-1], influent[-1]
        C_bw, Q_bw = backwash_in[:-1], backwash_in[-1]
        efficiency = self.get_removal_efficiency(self._t)
        filtered_C = C_in*(1.-self.components.x*efficiency)
        Q_recycle = Q_in*self.average_backwash_fraction
        Q_effluent = Q_in-Q_recycle

        if Q_bw > 0.:
            release = np.maximum(self._state, 0.) \
                / (self.filter_backwash_time/24.)
            backwash_C = C_bw + release/Q_bw
        else:
            release = np.zeros_like(self._state)
            backwash_C = C_bw.copy()
        self._release_rate[:] = release

        effluent, backwash_out, recycle = self.outs
        self._assign_state(
            effluent, np.append(filtered_C, Q_effluent))
        self._assign_state(
            backwash_out, np.append(backwash_C, Q_bw))
        self._assign_state(
            recycle, np.append(filtered_C, Q_recycle))

    def _update_dstate(self):
        if self._model == 'constant_splitting':
            effluent_factor, backwash_factor = self._constant_split_factors()
            effluent, backwash = self.outs
            self._assign_state(effluent, effluent.state)
            self._assign_state(backwash, backwash.state)
            if effluent.dstate is None:
                effluent.dstate = self._dstate*effluent_factor
            else:
                effluent.dstate[:] = self._dstate*effluent_factor
            if backwash.dstate is None:
                backwash.dstate = self._dstate*backwash_factor
            else:
                backwash.dstate[:] = self._dstate*backwash_factor
            return

        influent = self._ins_QC[0]
        backwash_in = self._ins_QC[1]
        d_in = self._ins_dQC[0]
        d_bw = self._ins_dQC[1]
        C_in, Q_bw = influent[:-1], backwash_in[-1]
        dC_in, dQ_in = d_in[:-1], d_in[-1]
        dC_bw, dQ_bw = d_bw[:-1], d_bw[-1]
        efficiency, de_dt = self._empirical_efficiency_and_derivative(self._t)
        multiplier = 1.-self.components.x*efficiency
        dfiltered_C = dC_in*multiplier \
            - C_in*self.components.x*de_dt
        dQ_recycle = dQ_in*self.average_backwash_fraction
        dQ_effluent = dQ_in-dQ_recycle

        dbackwash_C = dC_bw.copy()
        if Q_bw > 0.:
            tau = self.filter_backwash_time/24.
            release = self._release_rate
            drelease = self._dstate/tau
            dbackwash_C += (drelease*Q_bw-release*dQ_bw)/Q_bw**2

        effluent, backwash_out, recycle = self.outs
        effluent.dstate = np.append(dfiltered_C, dQ_effluent)
        backwash_out.dstate = np.append(dbackwash_C, dQ_bw)
        recycle.dstate = np.append(dfiltered_C, dQ_recycle)

    @property
    def AE(self):
        if self._AE is None:
            self._compile_AE()
        return self._AE

    def _compile_AE(self):
        if self._model != 'constant_splitting':
            raise RuntimeError('empirical_fouling uses an ODE, not an AE')
        state = self._state
        dstate = self._dstate
        update_state = self._update_state
        update_dstate = self._update_dstate

        def algebraic_equations(t, QC_ins, dQC_ins):
            self._t = t
            state[:] = QC_ins[0]
            dstate[:] = dQC_ins[0]
            update_state()
            update_dstate()

        self._AE = algebraic_equations

    @property
    def ODE(self):
        if self._ODE is None:
            self._compile_ODE()
        return self._ODE

    def _compile_ODE(self):
        if self._model != 'empirical_fouling':
            raise RuntimeError('constant_splitting uses an AE, not an ODE')
        dstate = self._dstate
        x = self.components.x
        update_state = self._update_state
        update_dstate = self._update_dstate

        def ordinary_differential_equations(t, QC_ins, state, dQC_ins):
            self._t = t
            C_in = QC_ins[0, :-1]
            Q_in = QC_ins[0, -1]
            Q_bw = QC_ins[1, -1]
            efficiency = self.get_removal_efficiency(t)
            capture = Q_in*C_in*x*efficiency
            if Q_bw > 0.:
                release = np.maximum(state, 0.) \
                    / (self.filter_backwash_time/24.)
            else:
                release = np.zeros_like(state)
            dstate[:] = capture-release
            self._release_rate[:] = release
            update_state()
            update_dstate()

        self._ODE = ordinary_differential_equations
