#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>
    Joy Zhang <joycheung1994@gmail.com>

Part of this module is based on the biosteam package:
https://github.com/BioSTEAMDevelopmentGroup/biosteam

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

import numpy as np
from math import pi, ceil
from biosteam.units import Pump
from biosteam.units.design_tools.mechanical import (
    brake_efficiency as brake_eff,
    motor_efficiency as motor_eff
    )
from .. import SanUnit
from ..utils import auom, select_pipe, format_str

__all__ = ('Pump', 'HydraulicDelay', 'WWTpump', 'wwtpump')


class Pump(SanUnit, Pump):
    '''
    Similar to the :class:`biosteam.units.Pump`,
    but can be initialized with :class:`qsdsan.SanStream` and :class:`qsdsan.WasteStream`,
    and allows dynamic simulation.

    See Also
    --------
    `biosteam.units.Pump <https://biosteam.readthedocs.io/en/latest/units/Pump.html>`_
    '''
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                  P=None, pump_type='Default', material='Cast iron',
                  dP_design=405300, ignore_NPSH=True,
                  init_with='Stream', F_BM_default=None, isdynamic=False):
        SanUnit.__init__(self, ID, ins, outs, thermo,
                         init_with=init_with, F_BM_default=F_BM_default,
                         isdynamic=isdynamic)
        self.P = P
        self.pump_type = pump_type
        self.material = material
        self.dP_design = dP_design
        self.ignore_NPSH = ignore_NPSH

    @property
    def state(self):
        '''The state of the Pump, including component concentrations [mg/L] and flow rate [m^3/d].'''
        if self._state is None: return None
        else:
            return dict(zip(list(self.components.IDs) + ['Q'], self._state))

    def _init_state(self):
        self._state = self._ins_QC[0]
        self._dstate = self._state * 0.

    def _update_state(self):
        '''updates conditions of output stream based on conditions of the Mixer'''
        self._outs[0].state = self._state

    def _update_dstate(self):
        '''updates rates of change of output stream from rates of change of the Mixer'''
        self._outs[0].dstate = self._dstate

    @property
    def AE(self):
        if self._AE is None:
            self._compile_AE()
        return self._AE

    def _compile_AE(self):
        _state = self._state
        _dstate = self._dstate
        _update_state = self._update_state
        _update_dstate = self._update_dstate
        def yt(t, QC_ins, dQC_ins):
            _state[:] = QC_ins[0]
            _dstate[:] = dQC_ins[0]
            _update_state()
            _update_dstate()
        self._AE = yt


# %%

class HydraulicDelay(Pump):
    '''
    A fake unit for implementing hydraulic delay by a first-order reaction
    (i.e., a low-pass filter) with a specified time constant [d].

    See Also
    --------
    `Benchmark Simulation Model No.1 implemented in MATLAB & Simulink <https://www.cs.mcgill.ca/~hv/articles/WWTP/sim_manual.pdf>`
    '''
    def __init__(self, ID='', ins=None, outs=(), thermo=None, t_delay=1e-4, *,
                 init_with='WasteStream', F_BM_default=None, isdynamic=True):
        SanUnit.__init__(self, ID, ins, outs, thermo,
                         init_with=init_with, F_BM_default=F_BM_default,
                         isdynamic=isdynamic)
        self.t_delay = t_delay
        self._concs = None

    def set_init_conc(self, **kwargs):
        '''set the initial concentrations [mg/L].'''
        Cs = np.zeros(len(self.components))
        cmpx = self.components.index
        for k, v in kwargs.items(): Cs[cmpx(k)] = v
        self._concs = Cs

    def _init_state(self):
        '''initialize state by specifying or calculating component concentrations
        based on influents. Total flow rate is always initialized as the sum of
        influent wastestream flows.'''
        self._state = self._ins_QC[0]
        self._dstate = self._state * 0
        if self._concs is not None:
            self._state[:-1] = self._concs

    def _run(self):
        s_in, = self.ins
        s_out, = self.outs
        s_out.copy_like(s_in)

    @property
    def ODE(self):
        if self._ODE is None:
            self._compile_ODE()
        return self._ODE

    def _compile_ODE(self):
        T = self.t_delay
        _dstate = self._dstate
        _update_dstate = self._update_dstate
        def dy_dt(t, QC_ins, QC, dQC_ins):
            Q_in = QC_ins[0,-1]
            Q = QC[-1]
            C_in = QC_ins[0,:-1]
            C = QC[:-1]
            if dQC_ins[0,-1] == 0:
                _dstate[-1] = 0
                _dstate[:-1] = (Q_in*C_in - Q*C)/(Q*T)
            else:
                _dstate[-1] = (Q_in - Q)/T
                _dstate[:-1] = Q_in/Q*(C_in - C)/T
            _update_dstate()
        self._ODE = dy_dt

    def _design(self):
        pass

    def _cost(self):
        pass


# %%

_hp_to_kW = auom('hp').conversion_factor('kW')
_lb_to_kg = auom('lb').conversion_factor('kg')
_ft_to_m = auom('ft').conversion_factor('m')
_ft3_to_gal = auom('ft3').conversion_factor('gallon')
_m3_to_gal = auom('m3').conversion_factor('gallon')
F_BM_pump = 1.18*(1+0.007/100) # 0.007 is for miscellaneous costs
default_F_BM = {
        'Pump': F_BM_pump,
        'Pump building': F_BM_pump,
        }
default_equipment_lifetime = {
    'Pump': 15,
    'Pump pipe stainless steel': 15,
    'Pump stainless steel': 15,
    'Pump chemical storage HDPE': 30,
    }

class WWTpump(SanUnit):
    '''
    Generic class for pumps used in wastewater treatment, [1]_
    all pumps are assumed be made of stainless steel.

    This class is intended to be used as a part of other units
    (e.g., :class:`~.AnMBR`), but it can be used as a standalone unit.

    Note that pump building concrete usage and excavation is not included here
    as pumps are often housed together with the reactors.

    Parameters
    ----------
    prefix : str
        If provided, all keys in design and cost dicts will be prefixed with
        the provided string.
    pump_type : str
        The type of the pump that determines the design algorithms to use.
        The following types are valid:

            - "permeate_cross-flow"
            - "retentate_CSTR"
            - "retentate_AF"
            - "recirculation_CSTR"
            - "recirculation_AF"
            - "lift"
            - "sludge"
            - "chemical"
            - "" (i.e., empty)

        When left as empty, the generic algorithm will be used and the following
        values should be included in `add_inputs` (in this order):

            - N_pump: number of pumps
            - L_s: pipe length of the suction side, [ft]
            - L_d: pipe length of the discharge side, [ft]
            - H_ts: total static head, [ft]
            - H_p: pressure head, [ft]

    Q_mgd : float
        Volumetric flow rate in million gallon per day, [mgd].
        Will use total volumetric flow through the unit if not provided.
    add_inputs : Iterable
        Additional inputs that will be passed to the corresponding design algorithm.
        Check the documentation of for the corresponding pump type
        for the design algorithm of the specific input requirements.
    capacity_factor : float
        A safety factor to handle peak flows.
    include_pump_cost : bool
        Whether to include pump cost.
    include_building_cost : bool
        Whether to include the cost of the pump building.
    include_OM_cost : bool
        Whether to include the operating and maintenance cost of the pump.
    F_BM : dict
        Bare module factors of the individual equipment.
    lifetime : dict
        Lifetime of the individual equipment.
    kwargs : dict
        Other attributes to be set.

    References
    ----------
    .. [1] Shoener et al., Design of Anaerobic Membrane Bioreactors for the
        Valorization of Dilute Organic Carbon Waste Streams.
        Energy Environ. Sci. 2016, 9 (3), 1102–1112.
        https://doi.org/10.1039/C5EE03715H.

    See Also
    --------
    :class:`~.AnMBR`
    '''
    _N_ins = 1
    _N_outs = 1

    _v = 3 # fluid velocity, [ft/s]
    _C = 110 # Hazen-Williams coefficient for stainless steel (SS)

    # Pump SS (for pumps within 300-1000 gpm)
    # http://www.godwinpumps.com/images/uploads/ProductCatalog_Nov_2011_spread2.pdf
    # assume 50% of the product weight is SS
    _SS_per_pump = 725 * 0.5
    _building_unit_cost = 90 # [$/ft2]

    _units = {
        'Pump pipe stainless steel': 'kg',
        'Pump stainless steel': 'kg',
        'Pump chemical storage HDPE': 'm3',
        }

    _valid_pump_types = (
        'permeate_cross-flow',
        'retentate_CSTR',
        'retentate_AF',
        'recirculation_CSTR',
        'recirculation_AF',
        'lift',
        'sludge',
        'chemical',
        '',
        )

    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream',
                 prefix='', pump_type='', Q_mgd=None, add_inputs=(),
                 capacity_factor=1.,
                 include_pump_cost=True, include_building_cost=False,
                 include_OM_cost=False,
                 F_BM=default_F_BM,
                 lifetime=default_equipment_lifetime,
                 **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with=init_with)
        self.pump_type = pump_type
        self.Q_mgd = Q_mgd
        try: iter(add_inputs)
        except: add_inputs = (add_inputs,)
        self.add_inputs = add_inputs
        self.capacity_factor = capacity_factor
        self.include_pump_cost = include_pump_cost
        self.include_building_cost = include_building_cost
        self.include_OM_cost = include_OM_cost
        self.F_BM.update(F_BM)
        self._default_equipment_lifetime.update(lifetime)

        self.prefix = prefix
        if prefix:
            self._units = {prefix+' '+[k][0].lower()+k[1:]:v for k, v in self._units.items()}
            self.F_BM = {prefix+' '+[k][0].lower()+k[1:]:v for k, v in self.F_BM.items()}
            self._default_equipment_lifetime = \
                {prefix+' '+[k][0].lower()+k[1:]:v for k, v in self._default_equipment_lifetime.items()}

        for attr, val in kwargs.items():
            setattr(self, attr, val)

    def _run(self):
        self.outs[0].copy_like(self.ins[0])

    def _format_key_start_with_prefix(self, start):
        return start.upper() if not self.prefix else self.prefix+' '+start.lower()

    def _design(self):
        pump_type = format_str(self.pump_type)
        if not pump_type:
            pipe, pumps = self._design_generic(self.Q_mgd, *self.add_inputs)
            hdpe = 0.
        else:
            design_func = getattr(self, f'design_{pump_type}')
            pipe, pumps, hdpe = design_func()

        D = self.design_results
        start = self._format_key_start_with_prefix('P')
        D[f'{start}ump pipe stainless steel'] = pipe
        D[f'{start}ump stainless steel'] = pumps
        if hdpe:
            D[f'{start}ump chemical storage HDPE'] = hdpe
            self._units[f'{start}ump chemical storage HDPE'] = 'm3'
        else:
            try:
                self._units.pop(f'{start}ump chemical storage HDPE')
            except KeyError:
                pass
            try:
                self.design_results.pop(f'{start}ump chemical storage HDPE')
            except KeyError:
                pass


    def _cost(self):
        C = self.baseline_purchase_costs
        C.clear()
        add_OPEX = self.add_OPEX
        add_OPEX.clear()
        Q_mgd, capacity_factor = self.Q_mgd, self.capacity_factor

        start = self._format_key_start_with_prefix('P')
        C[f'{start}ump'] = C[f'{start}ump building'] = 0.
        add_OPEX[f'{start}ump operating'] = add_OPEX[f'{start}ump maintenance'] = 0.
        # Pump
        if self.include_pump_cost:
            C[f'{start}ump'] = 2.065e5 + 7.721*1e4*Q_mgd # fitted curve

        # Operations and maintenance
        if self.include_OM_cost:
            FPC = capacity_factor * Q_mgd # firm pumping capacity
            O = M = 0. # USD/yr
            if 0 < FPC <= 7:
                O = 440*25*FPC**0.1285
                M = 360*25*FPC**0.1478
            elif 7 < FPC <= 41:
                O = 294.4*25*FPC**0.3335
                M = 255.2*25*FPC**0.3247
            elif 41 < FPC <= 80:
                O = 40.5*25*FPC**0.8661
                M = 85.7*25*FPC**0.6456
            else:
                O = 21.3*25*FPC**1.012
                M = 30.6*25*FPC**0.8806
            self.add_OPEX = {
                f'{start}ump operating': O/365/24,
                f'{start}ump maintenance': M/365/24,
                }

        # Pump building
        if self.include_building_cost:
            # Design capacity of intermediate pumps, gpm,
            GPM = capacity_factor * Q_mgd * 1e6 / 24 / 60
            if GPM == 0:
                N = 0
            else:
                N = 1 # number of buildings
                GPMi = GPM
                while GPMi > 80000:
                    N += 1
                    GPMi = GPM / N
            PBA = N * (0.0284*GPM+640) # pump building area, [ft2]
            C[f'{start}ump building'] = PBA * self.building_unit_cost

        self.power_utility.consumption = self.BHP/self.motor_efficiency * _hp_to_kW


    # Generic algorithms that will be called by all design functions
    def _design_generic(self, Q_mgd, N_pump, L_s=0., L_d=0., H_ts=0., H_p=0.):
        self.Q_mgd, self._H_ts, self._H_p = Q_mgd, H_ts, H_p
        v, C, Q_cfs = self.v, self.C, self.Q_cfs # [ft/s], -, [ft3/s]

        ### Suction side ###
        # Suction pipe (permeate header) dimensions
        OD_s, t_s, ID_s = select_pipe(Q_cfs/N_pump, v) # [in]

        # Suction friction head, [ft]
        self._H_sf = 3.02 * L_s * (v**1.85) * (C**(-1.85)) * ((ID_s/12)**(-1.17))

        ### Discharge side ###
        # Discharge pipe (permeate collector) dimensions
        OD_d, t_d, ID_d = select_pipe(Q_cfs, v)

        # Discharge friction head, [ft]
        self._H_df = 3.02 * L_d * (v**1.85) * (C**(-1.85)) * ((ID_d/12)**(-1.17))

        ### Material usage ###
        # Pipe SS, assume stainless steel, density = 0.29 lbs/in3
        # SS volume for suction, [in3]
        self._N_pump = N_pump
        V_s = N_pump * pi/4*((OD_s)**2-(ID_s)**2) * (L_s*12)
        # SS volume for discharge, [in3]
        V_d = pi/4*((OD_d)**2-(ID_d)**2) * (L_d*12)

        # Total SS mass, [kg]
        M_SS_pipe = 0.29 * (V_s+V_d) * _lb_to_kg
        M_SS_pump = N_pump * self.SS_per_pump
        return M_SS_pipe, M_SS_pump


    def design_permeate_cross_flow(self, Q_mgd=None, N_pump=None, D=None,
                                   TMP=None, include_aerobic_filter=False,
                                   **kwargs):
        '''
        Design pump for the permeate stream of cross-flow membrane configuration.

        Parameters defined through the `add_inputs` argument upon initialization of
        this unit (Q_mgd listed separately) will be used if not provided
        when calling this function.

        Parameters
        ----------
        Q_mgd : float
            Volumetric flow rate in million gallon per day, [mgd].
        N_pump : int
            Number of the pumps.
        D : float
            Depth of the reactor, [ft].
        TMP : float
            Transmembrane pressure, [psi].
        include_aerobic_filter : bool
            Whether aerobic filter is included in the reactor design,
            additional head will be added if the filter is included.
        kwargs : dict
            Additional attribute values to set (e.g., `L_s`, `H_ts`),
            this will overwrite the default values.
        '''
        add_inputs = self.add_inputs
        Q_mgd = Q_mgd or self.Q_mgd
        N_pump = N_pump or add_inputs[0]
        D = D or add_inputs[1]
        TMP = TMP or add_inputs[2]
        include_aerobic_filter = include_aerobic_filter or add_inputs[3]

        H_ts_PERM = D if include_aerobic_filter else 0

        val_dct = dict(
            L_s=20, # based on a 30-module unit with a total length of 6 m, [ft]
            L_d=10*N_pump, # based on a 30-module unit with a total width of 1.6 m and extra space, [ft]
            H_ts=H_ts_PERM, #  H_ds_PERM (D_tank) - H_ss_PERM (0 or D_tank)
            H_p=TMP*2.31 # TMP in water head, [ft], comment below on 2.31
            )
        val_dct.update(kwargs)
        M_SS_IR_pipe, M_SS_IR_pump = self._design_generic(
            Q_mgd=Q_mgd, N_pump=N_pump, **val_dct)

        # # factor = 2.31 calculated by
        # factor = auom('psi').conversion_factor('Pa') # Pa is kg/m/s2, now in [Pa]
        # factor /= 9.81 # divided by the standard gravity in m/s2, now in [kg/m2]
        # factor /= 1e3 # divided by water's density in kg/m3, now in [m]
        # factor *= auom('m').conversion_factor('ft') # m to ft

        return M_SS_IR_pipe, M_SS_IR_pump, 0


    def design_retentate_CSTR(self, Q_mgd=None, N_pump=None, **kwargs):
        '''
        Design pump for the retent stream of CSTR reactors.

        Parameters defined through the `add_inputs` argument upon initialization of
        this unit (Q_mgd listed separately) will be used if not provided
        when calling this function.

        Parameters
        ----------
        Q_mgd : float
            Volumetric flow rate in million gallon per day, [mgd].
        N_pump : int
            Number of the pumps.
        kwargs : dict
            Additional attribute values to set (e.g., `L_s`, `H_ts`),
            this will overwrite the default values.
        '''
        Q_mgd = Q_mgd or self.Q_mgd
        N_pump = N_pump or self.add_inputs[0]

        val_dct = dict(
            L_s=100, # pipe length per module
            L_d=30, # pipe length per module (same as the discharge side of lift pump)
            H_ts=0., # H_ds_IR (D_tank) - H_ss_IR (D_tank)
            H_p=0. # no pressure
            )
        val_dct.update(kwargs)

        M_SS_IR_pipe, M_SS_IR_pump = self._design_generic(
            Q_mgd=Q_mgd, N_pump=N_pump, **val_dct)

        return M_SS_IR_pipe, M_SS_IR_pump, 0


    def design_retentate_AF(self, Q_mgd=None, N_pump=None, D=None, **kwargs):
        '''
        Design pump for the retentate stream of AF reactors.

        Parameters defined through the `add_inputs` argument upon initialization of
        this unit (Q_mgd listed separately) will be used if not provided
        when calling this function.

        Parameters
        ----------
        Q_mgd : float
            Volumetric flow rate in million gallon per day, [mgd].
        N_pump : int
            Number of the pumps.
        D : float
            Depth of the reactor, [ft].
        kwargs : dict
            Additional attribute values to set (e.g., `L_s`, `H_ts`),
            this will overwrite the default values.
        '''
        add_inputs = self.add_inputs
        Q_mgd = Q_mgd or self.Q_mgd
        N_pump = N_pump or add_inputs[0]
        D = D or add_inputs[1]

        val_dct = dict(
            L_s=100, # assumed pipe length per filter, [ft]
            L_d=30, # same as discharge side of lift pumping, [ft]
            H_ts=0., # H_ds_IR (D) - H_ss_IR (D)
            H_p=0. # no pressure
            )
        val_dct.update(kwargs)

        M_SS_IR_pipe, M_SS_IR_pump = self._design_generic(
            Q_mgd=Q_mgd, N_pump=N_pump, **val_dct)

        return M_SS_IR_pipe, M_SS_IR_pump, 0


    def design_recirculation_CSTR(self, Q_mgd=None, N_pump=None, L=None, **kwargs):
        '''
        Design pump for the recirculation stream of reactors.

        Parameters defined through the `add_inputs` argument upon initialization of
        this unit (Q_mgd listed separately) will be used if not provided
        when calling this function.

        Parameters
        ----------
        Q_mgd : float
            Volumetric flow rate in million gallon per day, [mgd].
        N_pump : int
            Number of the pumps.
        L : float
            Length of the reactor, [ft].
        kwargs : dict
            Additional attribute values to set (e.g., `L_s`, `H_ts`),
            this will overwrite the default values.
        '''
        Q_mgd = Q_mgd or self.Q_mgd
        L = L or self.add_inputs[0]

        val_dct = dict(
            L_s=0., # ignore suction side
            L_d=L, # pipe length per train
            H_ts=5., # H_ds_IR (5) - H_ss_IR (0)
            H_p=0. # no pressure
            )
        val_dct.update(kwargs)

        M_SS_IR_pipe, M_SS_IR_pump = self._design_generic(
            Q_mgd=Q_mgd, N_pump=N_pump, **val_dct)

        return M_SS_IR_pipe, M_SS_IR_pump, 0


    def design_recirculation_AF(self, Q_mgd=None, N_pump=None, d=None,
                                D=None, **kwargs):
        '''
        Design pump for the recirculation stream of AF reactors.

        Parameters defined through the `add_inputs` argument upon initialization of
        this unit (Q_mgd listed separately) will be used if not provided
        when calling this function.

        Parameters
        ----------
        Q_mgd : float
            Volumetric flow rate in million gallon per day, [mgd].
        N_pump : int
            Number of the pumps.
        d : float
            Diameter (or width) of the reactor, [ft].
        D : float
            Depth of the reactor, [ft].
        kwargs : dict
            Additional attribute values to set (e.g., `L_s`, `H_ts`),
            this will overwrite the default values.
        '''
        add_inputs = self.add_inputs
        Q_mgd = Q_mgd or self.Q_mgd
        N_pump = N_pump or add_inputs[0]
        d = d or add_inputs[1]
        D = D or add_inputs[2]

        val_dct = dict(
            L_s=d+D, # pipe length per filter, [ft]
            L_d=30, # same as discharge side of lift pumping, [ft]
            H_ts=0., # H_ds_IR (D) - H_ss_IR (D)
            H_p=0. # no pressure
            )
        val_dct.update(kwargs)

        M_SS_IR_pipe, M_SS_IR_pump = self._design_generic(
            Q_mgd=Q_mgd, N_pump=N_pump, **kwargs)

        return M_SS_IR_pipe, M_SS_IR_pump, 0


    def design_lift(self, Q_mgd=None, N_pump=None, D=None, **kwargs):
        '''
        Design pump for the filter tank to lift streams.

        Parameters defined through the `add_inputs` argument upon initialization of
        this unit (Q_mgd listed separately) will be used if not provided
        when calling this function.

        Parameters
        ----------
        Q_mgd : float
            Volumetric flow rate in million gallon per day, [mgd].
        N_pump : int
            Number of the pumps.
        D : float
            Depth of the filter tank, [ft].
        kwargs : dict
            Additional attribute values to set (e.g., `L_s`, `H_ts`),
            this will overwrite the default values.
        '''
        add_inputs = self.add_inputs
        Q_mgd = Q_mgd or self.Q_mgd
        N_pump = N_pump or add_inputs[0]
        D = D or add_inputs[1]

        val_dct = dict(
            L_s=150, # length of suction pipe per filter, [ft]
            L_d=30, # pipe length per filter, [ft]
            H_ts=D, # H_ds_LIFT (D) - H_ss_LIFT (0)
            H_p=0. # no pressure
            )
        val_dct.update(kwargs)

        M_SS_IR_pipe, M_SS_IR_pump = self._design_generic(
            Q_mgd=Q_mgd, N_pump=N_pump, **kwargs)

        return M_SS_IR_pipe, M_SS_IR_pump, 0


    def design_sludge(self, Q_mgd=None, N_pump=None, **kwargs):
        '''
        Design pump for handling waste sludge.

        Parameters
        ----------
        Q_mgd : float
            Volumetric flow rate in million gallon per day, [mgd].
        N_pump : int
            Number of the pumps.
        kwargs : dict
            Additional attribute values to set (e.g., `L_s`, `H_ts`),
            this will overwrite the default values.
        '''
        Q_mgd = Q_mgd or self.Q_mgd
        N_pump = N_pump or 1

        val_dct = dict(
            L_s=50, # length of suction pipe, [ft]
            L_d=50, # length of discharge pipe, [ft]
            H_ts=0., # H_ds_LIFT (D) - H_ss_LIFT (0)
            H_p=0. # no pressure
            )
        val_dct.update(kwargs)

        M_SS_IR_pipe, M_SS_IR_pump = self._design_generic(
            Q_mgd=Q_mgd, N_pump=N_pump, **kwargs)

        return M_SS_IR_pipe, M_SS_IR_pump, 0


    def design_chemical(self, Q_mgd=None, N_pump=None, **kwargs):
        '''
        Design pump for membrane cleaning chemicals (NaOCl and citric acid),
        storage containers are included, and are assumed to be cubic in shape
        and made of HDPE.

        Parameters defined through the `add_inputs` argument upon initialization of
        this unit (Q_mgd listed separately) will be used if not provided
        when calling this function.

        Parameters
        ----------
        Q_mgd : float
            Volumetric flow rate in million gallon per day, [mgd].
        N_pump : int
            Number of the pumps.
        kwargs : dict
            Additional attribute values to set (e.g., `L_s`, `H_ts`),
            this will overwrite the default values.
        '''
        if not Q_mgd:
            V_CHEM = self.ins[0].F_vol * 24 * 7 * 2 # for two weeks of storage, [m3]
            Q_CHEM_mgd = self.Q_mgd
        else:
            V_CHEM = (Q_mgd*1e6/_m3_to_gal) * 7 * 2
            Q_CHEM_mgd = Q_mgd
        N_pump = N_pump or 1

        # HDPE volume, [m3], 0.003 [m] is the thickness of the container
        V_HDPE = 0.003 * (V_CHEM**(1/3))**2*6
        # # Mass of HDPE, [m3], 950 is the density of the HDPE in [kg/m3]
        # M_HDPE = 950 * V_HDPE

        H_ss_CHEM = V_CHEM**(1/3) / _ft_to_m
        # 9'-7" is the water level in membrane trains
        # 18" is the distance from C/L of the pump to the ground
        H_ds_CHEM = 9 + 7/12 - 18/12
        H_ts_CHEM = H_ds_CHEM - H_ss_CHEM

        val_dct = dict(
            L_s=0., # no suction pipe
            L_d=30.,
            H_ts=H_ts_CHEM,
            H_p=0. # no pressure
            )
        val_dct.update(kwargs)

        M_SS_CHEM_pipe, M_SS_CHEM_pump = self._design_generic(
            Q_mgd=Q_CHEM_mgd, N_pump=N_pump, **kwargs)

        return M_SS_CHEM_pipe, M_SS_CHEM_pump, V_HDPE


    @property
    def pump_type(self):
        '''
        [str] The type of the pump that determines the design algorithms to use.
        Use `valid_pump_type` to see acceptable pump types.
        '''
        return self._pump_type
    @pump_type.setter
    def pump_type(self, i):
        i_lower = i.lower()
        i_lower = i_lower.replace('cstr', 'CSTR')
        i_lower = i_lower.replace('af', 'AF')
        if i_lower not in self.valid_pump_types:
            raise ValueError(f'The given `pump_type` "{i}" is not valid, '
                             'check `valid_pump_types` for acceptable pump types.')
        self._pump_type = i_lower

    @property
    def valid_pump_types(self):
        '''[tuple] Acceptable pump types.'''
        return self._valid_pump_types

    @property
    def Q_mgd(self):
        '''
        [float] Volumetric flow rate in million gallon per day, [mgd].
        Will use total volumetric flow through the unit if not provided.
        '''
        if self._Q_mgd:
            return self._Q_mgd
        return self.F_vol_in*_m3_to_gal*24/1e6
    @Q_mgd.setter
    def Q_mgd(self, i):
        self._Q_mgd = i

    @property
    def Q_gpm(self):
        '''[float] Volumetric flow rate in gallon per minute, [gpm].'''
        return self.Q_mgd*1e6/24/60

    @property
    def Q_cmd(self):
        '''
        [float] Volumetric flow rate in cubic meter per day, [cmd].
        '''
        return self.Q_mgd *1e6/_m3_to_gal # [m3/day]

    @property
    def Q_cfs(self):
        '''[float] Volumetric flow rate in cubic feet per second, [cfs].'''
        return self.Q_mgd*1e6/24/60/60/_ft3_to_gal

    @property
    def capacity_factor(self):
        '''[float] A safety factor to handle peak flow.'''
        return self._capacity_factor
    @capacity_factor.setter
    def capacity_factor(self, i):
        self._capacity_factor = i

    @property
    def N_pump(self):
        '''[int] Number of pumps.'''
        return self._N_pump or 1
    @N_pump.setter
    def N_pump(self, i):
        self._N_pump = ceil(i)

    @property
    def v(self):
        '''[float] Fluid velocity, [ft/s].'''
        return self._v
    @v.setter
    def v(self, i):
        self._v = i

    @property
    def C(self):
        '''[float] Hazen-Williams coefficient to calculate fluid friction.'''
        return self._C
    @C.setter
    def C(self, i):
        self._C = i

    @property
    def SS_per_pump(self):
        '''[float] Quantity of stainless steel per pump, [kg/ea].'''
        return self._SS_per_pump
    @SS_per_pump.setter
    def SS_per_pump(self, i):
        self._SS_per_pump = i

    @property
    def H_sf(self):
        '''[float] Suction friction head, [ft].'''
        return self._H_sf

    @property
    def H_df(self):
        '''[float] Discharge friction head, [ft].'''
        return self._H_df

    @property
    def H_ts(self):
        '''[float] Total static head, [ft].'''
        return self._H_ts

    @property
    def H_p(self):
        '''[float] Pressure head, [ft].'''
        return self._H_p

    @property
    def TDH(self):
        '''[float] Total dynamic head, [ft].'''
        return self.H_ts+self.H_sf+self.H_df+self.H_p

    @property
    def BHP(self):
        '''[float] Brake horsepower, [hp].'''
        return (self.TDH*self.Q_gpm)/3960/self.brake_efficiency

    @property
    def brake_efficiency(self):
        '''[float] Brake efficiency.'''
        return brake_eff(self.Q_gpm)

    @property
    def motor_efficiency(self):
        '''[float] Motor efficiency.'''
        return motor_eff(self.BHP)

    @property
    def building_unit_cost(self):
        '''[float] Unit cost of the pump building, [USD/ft2].'''
        return self._building_unit_cost if self.include_cost else 0.
    @building_unit_cost.setter
    def building_unit_cost(self, i):
        self._building_unit_cost = i


# %%

# =============================================================================
# Decorator
# =============================================================================

def wwtpump(ID, ins=(), prefix='', pump_type='', Q_mgd=None, add_inputs=(),
            capacity_factor=1., include_pump_cost=True, include_building_cost=False,
            include_OM_cost=False, F_BM=F_BM_pump, lifetime=default_equipment_lifetime,
            **kwargs):
    '''
    Handy decorator to add a :class:`~.WWTpump` as an attribute of
    a :class:`qsdsan.SanUnit`.

    Refer to class:`WWTpump` for the parameters needed for using this decorator.

    See Also
    --------
    :class:`~.WWTpump`

    :class:`~.AnMBR`

    References
    ----------
    [1] Shoener et al., Design of Anaerobic Membrane Bioreactors for the
    Valorization of Dilute Organic Carbon Waste Streams.
    Energy Environ. Sci. 2016, 9 (3), 1102–1112.
    https://doi.org/10.1039/C5EE03715H.
    '''
    return lambda cls: add_pump(cls, ID, ins, prefix, pump_type, Q_mgd, add_inputs,
                                capacity_factor, include_pump_cost, include_building_cost,
                                include_OM_cost, F_BM, lifetime, **kwargs)


def add_pump(cls, ID, ins, prefix, pump_type, Q_mgd, add_inputs,
             capacity_factor, include_pump_cost, include_building_cost,
             include_OM_cost, F_BM, lifetime, **kwargs):
    pump = WWTpump(
        ID, ins=ins,
        prefix=prefix, pump_type=pump_type, Q_mgd=Q_mgd, add_inputs=add_inputs,
        capacity_factor=capacity_factor,
        include_pump_cost=include_pump_cost,
        include_building_cost=include_building_cost,
        include_OM_cost=include_OM_cost,
        F_BM=F_BM, lifetime=lifetime, **kwargs)
    setattr(cls, f'{ID}_pump', pump)
    return cls