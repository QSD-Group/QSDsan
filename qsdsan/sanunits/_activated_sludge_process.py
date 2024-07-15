#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:

    Yalin Li <mailto.yalin.li@gmail.com>
    
    Saumitra Rai <raisaumitra9@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

from warnings import warn
from math import ceil
from biosteam import Stream
from .. import SanUnit, WasteStream
from ..sanunits import HXutility, WWTpump, Mixer

from ..sanunits._pumping import (
    default_F_BM as default_WWTpump_F_BM,
    default_equipment_lifetime as default_WWTpump_equipment_lifetime,
    )

from ..equipments import Blower, GasPiping
from ..utils import auom, calculate_excavation_volume

from qsdsan.utils import ospath, time_printer, load_data, get_SRT, get_oxygen_heterotrophs, get_oxygen_autotrophs, get_airflow, get_P_blower


__all__ = ('ActivatedSludgeProcess','TreatmentTrains',)

_ft2_to_m2 = auom('ft2').conversion_factor('m2')
F_BM_pump = 1.18*(1+0.007/100) # 0.007 is for miscellaneous costs
default_F_BM = {
        'Pumps': F_BM_pump,
        'Pump building': F_BM_pump,
        }
default_equipment_lifetime = {
    'Pumps': 15,
    'Pump pipe stainless steel': 15,
    'Pump stainless steel': 15,
    }

class ActivatedSludgeProcess(SanUnit):
    '''
    A steady-state activated sludge process containing rectangular
    continuously stirred tank reactors and clarifiers for settling
    based on [1]_ (largely followed Example 6.2) and code scripts for [2]_ .

    Single tank and clarifier are assumed to have the same width;
    pump/blower buildings are assumed to have the same width as the
    total tanks/clarifiers.

    Parameters
    ----------
    ins : Iterable
        Influent, air for aeration.
    outs : Iterable
        Effluent, waste activated sludge sludge, gas emission.
    N_train : int
        Number of treatment train, should be at least two in case one failing.
    T : float
        Temperature in the aeration tank/anaerobic digester, [K].
    X_i0 : float
        Inert biomass concentration in the influent, [mg VSS/L],
        will calculated based on `inert_biomass` of the current `CompiledComponents` if not provided.
    X_v : float
        Reactor volatile suspended solids (biomass concentration in the aeration tank),
        referred to as mixed volatile liquor suspended solids
        (MLVSS, assumed to be the same as MLSS in [2]_), [mg/L].
        X_v = X_i + X_a (inert + active)
    X_e : float
        Biomass concentration in the effluent, [mg VSS/L].
    X_w : float
        Biomass concentration in the waste activated sludge, [mg VSS/L].
    SLR : float
        Solids loading rate for the clarifier, [lb/ft2/d].
    SF : float
        Safety factor to scale up the minimum solids retention time (SRT), should be larger than 1.
    aeration_power : float
        Unit power usage for aeration, [kWh/kg O2].
    q_hat: float
        Maximum specific rate of substrate utilization, [mg BOD/mg VSS/d],
        other important definitions are:

            - q_hat = mu_hat * Y (mu_hat is maximum specific growth rate, [1/time])
            - r_ut = -q_hat*(S/(K+S))*X_a (r_ut is rate of substrate utilization, [substrate mass/volume/time])

    K : float
        Substrate concentration giving one-half the maximum rate, [mg COD/L].
    Y : float
        Biomass yield, [mg VSS/mg BOD].
    b : float
        Endogenous decay coefficient, [1/d].
    f_d : float
        Fraction of the active biomass that is biodegradable.
    COD_factor : float
        Biomass-to-COD conversion factor, [mg COD/mg VSS].
    q_UAP : float
        Maximum specific rate of utilization-associated products degradation,
        [mg COD/mg VSS/d].

        SMP = UAP + BAP

        - SMP: soluble microbial products
        - UAP: utilization-associated products, produced directly from substrate metabolism
        - BAP: biomass-associated products, produced by basic metabolism

    q_BAP: float
        Maximum specific rate of biomass-associated products degradation,
        [mg COD/mg VSS/d].
    k1 : float
        UAP-formation coefficient, [mg COD/mg BOD].
    k2 : float
        BAP-formation coefficient, [mg COD/mg VSS/d].
    K_UAP : float
        Half-maximum rate concentrations for UAP, [mg COD/L].
    K_BAP : float
        Half-maximum rate concentrations for BAP, [mg COD/L].
    kwargs : dict
        Other attributes to be set.

    References
    ----------
    .. [1] Rittmann, B.; McCarty, P.; McCarty, P. L.; Bruce, R.
        Environmental Biotechnology: Principles and Applications;
        McGraw-Hill Companies,Incorporated, 2001.
    .. [2] Shoener, B. D.; Zhong, C.; Greiner, A. D.; Khunjar, W. O.; Hong, P.-Y.; Guest, J. S.
        Design of Anaerobic Membrane Bioreactors for the Valorization
        of Dilute Organic Carbon Waste Streams.
        Energy Environ. Sci. 2016, 9 (3), 1102–1112.
        https://doi.org/10.1039/C5EE03715H.

    See Also
    --------
    `MATLAB codes <https://github.com/QSD-Group/AnMBR>`_ used in ref 1,
    especially the system layout `diagrams <https://github.com/QSD-Group/AnMBR/blob/main/System_Design.pdf>`_.

    '''

    _N_ins = 2
    _N_outs = 3

    # Aeration tank, all in ft
    _W_tank = 21
    _D_tank = 12
    _freeboard = 2
    _W_dist = 4.5
    _W_eff = 5

    # Pump/blower buildings, all in ft
    _W_PB = 38 + 3/12
    _W_BB = 38 + 3/12
    _L_BB = 76 + 8/12
    _t_wall = None
    _t_slab = None

    # Wet well for mixed liquor storage, all in ft
    _L_WW = 8
    _W_WW = 8
    _D_WW = 12

    # Excavation
    _excav_slope = 1.5 # horizontal/vertical
    _constr_access = 3 # ft

    # Costs
    excav_unit_cost = (8+0.3) / 27 # $/ft3, 27 is to convert from $/yd3
    wall_concrete_unit_cost = 650 / 27 # $/ft3
    slab_concrete_unit_cost = 350 / 27 # $/ft3

    pumps = ('inf', 'recir')
    auxiliary_unit_names = ('heat_exchanger',)

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 F_BM_default=1, include_construction=False,
                 N_train=2, T=35+273.15, X_i0=None, X_v=1210, X_e=15, X_w=10000,
                 SLR=20, SF=4, aeration_power=1, # p337 of ref 1
                 q_hat=12, K=20, Y=0.5, b=0.396, f_d=0.8, # change from 0.2 of ref 2 based on p171 of ref 1
                 COD_factor=1.42, # p8 of Complex Systems in ref 1
                 q_UAP=1.8, q_BAP=0.1, k1=0.12, k2=0.09, K_UAP=100, K_BAP=85,
                 F_BM=default_F_BM, lifetime=default_equipment_lifetime,
                 **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with, F_BM_default=1,
                         include_construction=include_construction)
        ID, ins0 = self.ID, self.ins[0]
        self._inf = ins0.copy(f'{ID}_inf')
        self._ras = ins0.copy(f'{ID}_ras')
        self.N_train = N_train
        self.T = T
        self.X_i0 = X_i0
        self.X_v = X_v
        self.X_e = X_e
        self.X_w = X_w
        self.SLR = SLR
        self.SF = SF
        self.aeration_power = aeration_power
        self.q_hat = q_hat
        self.K = K
        self.Y = Y
        self.b = b
        self.f_d = f_d
        self.COD_factor = COD_factor
        self.q_UAP = q_UAP
        self.q_BAP = q_BAP
        self.k1 = k1
        self.k2 = k2
        self.K_UAP = K_UAP
        self.K_BAP = K_BAP
        ID = self.ID
        hx_in = Stream(f'{ID}_hx_in')
        hx_out = Stream(f'{ID}_hx_out')
        # Add '.' in ID for auxiliary units
        self.heat_exchanger = HXutility(ID=f'.{ID}_hx', include_construction=include_construction,
                                        ins=hx_in, outs=hx_out, T=T)
        self.blower = blower = Blower('blower', linked_unit=self, N_reactor=N_train)
        self.air_piping = air_piping = GasPiping('air_piping', linked_unit=self, N_reactor=N_train)
        self.equipments = (blower, air_piping)
        self.F_BM.update(F_BM)
        self._default_equipment_lifetime.update(lifetime)

        for attr, value in kwargs.items(): setattr(self, attr, value)


    # Equation/page numbers are noted for the 2001 Rittmann and McCarty book
    def _run(self):
        inf, air = self.ins
        self._inf.copy_like(inf)
        eff, was, emission = self.outs
        air.phase = emission.phase = 'g'
        Q = auom('m3/hr').convert(inf.F_vol, 'L/d')
        self._Q = Q/1e3 # convert to m3/d
        S_0 = inf.COD # influent substrate concentration
        cmps = self.components
        X_i0 = self.X_i0
        if not X_i0:
            try:
                inf_iconc = inf.iconc['inert_biomass']
            except AttributeError: # group not defined
                raise AttributeError('The `inert_biomass` group of current `CompiledComponents` '
                                     'has not been defined. '
                                     'Please either define it using the `define_group` function of '
                                     '`CompiledComponents` or provide `X_i0` for this '
                                     f'`{self.__class__.__name__}` unit.')
            X_i0 = inf_iconc if isinstance(inf_iconc, (int, float)) else inf_iconc.sum()

        X_v, X_e, X_w, SF = self.X_v, self.X_e, self.X_w, self.SF
        K, Y, b, f_d, COD_factor, q_hat = \
            self.K, self.Y, self.b, self.f_d, self.COD_factor, self.q_hat
        q_UAP, q_BAP, k1, k2, K_UAP, K_BAP = \
            self.q_UAP, self.q_BAP, self.k1, self.k2, self.K_UAP, self.K_BAP

        # Maximum washout solids retention time, [d]
        # (i.e., SRT below which washout occurs, Figure 3.3/eq. 3.26 when S0>>K)
        SRT_min = 1/(Y*q_hat-b)
        SRT = round(SRT_min*SF) # [d]

        # Steady state substrate concentration in the reactor, [mg/L]
        # eq. 5.39 (S=Se, effluent substrate conc.)
        S = K*((1+b*SRT)/(SRT*(Y*q_hat-b)-1))

        # Hydraulic retention time, [d], eq. 5.47
        HRT = self._HRT = SRT/X_v*(X_i0+Y*(S_0-S)*(1+(1-f_d)*b*SRT)/(1+b*SRT))
        V = Q * HRT # reactor volume, L

        # Rate of substrate utilization [mg/L/d], eq. 5.42
        r_ut = -(S_0-S)/HRT

        # Active biomass in the reactor, [mg VSS/L], eq. 5.43
        X_a = (SRT/HRT) * Y*(S_0-S)/(1+b*SRT)

        # Use mass balance on X within the entire control volume
        # (aeration tank and the clarifier)
        # to solve waste activated sludge flow, [L/d]
        dX_vdt = X_v * V / SRT # [mg VSS/d]
        Q_was = (dX_vdt-(Q*X_e))/(X_w-X_e)
        self._Q_was = Q_was/1e3 # convert to m3/d

        # Use mass balance on X for the clarifier to solve the return activated sludge flow, [L/d]
        Q_ras = ((X_w*Q_was)+X_e*(Q-Q_was)-(X_v*Q))/(X_v-X_w)
        self._Q_ras = Q_ras/1e3 # convert to m3/d

        # Soluble microbial products, eqs. 3.38 and 3.39
        # equations for chemostat holds as solubles cannot be settled out
        UAP = -(q_UAP*X_a*HRT+K_UAP+k1*r_ut*HRT)/2 + \
            ((q_UAP*X_a*HRT+K_UAP+k1*r_ut*HRT)**2-4*K_UAP*k1*r_ut*HRT)**0.5/2
        BAP = -(K_BAP+(q_BAP-k2)*X_a*HRT)/2 + \
            ((K_BAP+(q_BAP-k2)*X_a*HRT)**2+4*K_BAP*k2*X_a*HRT)**0.5/2
        SMP = UAP + BAP

        # Aeration requirements
        inf_O2_demand = Q*S_0 + Q*X_i0*COD_factor # [mg O2/d]
        eff_O2_demand = Q*(S+SMP) + dX_vdt*COD_factor # [mg O2/d]
        O2_uptake = inf_O2_demand - eff_O2_demand # [mg O2/d]

        # Update stream flows
        eff.copy_flow(inf)
        was.copy_flow(inf)
        Q_eff = Q - Q_was

        # Non-particulate components
        substrates = cmps.substrates
        S_conc = eff.iconc['substrates'] * S/S_0
        eff_dct = dict.fromkeys(cmps.IDs, 0)
        eff_dct.pop('H2O')
        if len(substrates) > 1:
            eff_dct.update({cmp.ID: conc for cmp, conc in zip(substrates, S_conc)})
        else: # Only one component in `substrates`
            eff_dct[substrates[0].ID] = S_conc
        # Non-substrate and non-particulate (soluble, colloidal, dissolved gas) components
        # will be split in the same way as flow rate mass-wise
        was_dct = eff_dct.copy()

        # Particulate components
        try:
            active_biomass = cmps.active_biomass
            X_a_inf = inf.iconc['active_biomass'].sum()
        except AttributeError:
            warn('The `active_biomass` group of current `CompiledComponents` '
                 'has not been defined and is assumed to be 0.')
            X_a_inf = 0
        if X_a_inf != 0:
            X_a_e_conc = eff.iconc['active_biomass'] * X_e/X_a_inf
            X_a_w_conc = eff.iconc['active_biomass'] * X_w/X_a_inf
        else:
            X_a_e_conc = X_a_w_conc = eff.iconc['active_biomass']
        if len(active_biomass) > 1:
            eff_dct.update({cmp.ID:conc for cmp, conc in zip(active_biomass, X_a_e_conc)})
            was_dct.update({cmp.ID:conc for cmp, conc in zip(active_biomass, X_a_w_conc)})
        else:
            eff_dct[active_biomass[0].ID] = X_a_e_conc
            was_dct[active_biomass[0].ID] = X_a_w_conc
        eff.set_flow_by_concentration(
            flow_tot=Q_eff, concentrations=eff_dct, units=('L/d', 'mg/L'))
        was.set_flow_by_concentration(
            flow_tot=Q_was, concentrations=was_dct, units=('L/d', 'mg/L'))

        # Assume all non-active solids (including inorganic solids) go to sludge
        X_i = set(cmps.solids).difference(set(active_biomass))
        X_i_IDs = tuple((cmp.ID for cmp in X_i))
        eff.imass[X_i_IDs] = 0
        was.imass[X_i_IDs] = inf.imass[X_i_IDs]
        ras = self._ras
        ras.copy_like(was)
        ras.mass *= Q_ras/Q_was

        air.empty()
        emission.empty()
        air.imass['O2'] = auom('mg/d').convert(O2_uptake, 'kg/hr')
        emission.imass['N2'] = air.imass['N2'] = air.imol['O2']/0.21*0.79

        eff.T = was.T = emission.T = self.T


    _units = {
        'N_train': '',
        'HRT': 'd',
        'Aeration tank volume': 'ft3',
        'Aeration tank length': 'ft',
        'Aeration tank width': 'ft',
        'Aeration tank depth': 'ft',
        'Q_was': 'm3/d',
        'Q_ras': 'm3/d',
        'Clarifier area': 'ft2',
        'Wall concrete': 'ft3',
        'Slab concrete': 'ft3',
        'Excavation': 'ft3',
        'Pump pipe stainless steel': 'kg',
        'Pump stainless steel': 'kg',
        }
    
    def _design(self):
        D = self.design_results
        D['HRT'] = self.HRT
        D['Q_was'] = self.Q_was
        D['Q_ras'] = self.Q_ras

        # Calculate the number of trains so that
        # the aeration tank length is within 23-30 m
        V_tot = auom('m3').convert((self.Q*self.HRT), 'ft3') # m3 to ft3
        N_train = 2
        V_tank = V_tot / N_train
        W_tank = D['Aeration tank width'] = self.W_tank
        D_tank = D['Aeration tank depth'] = self.D_tank
        A_tank = W_tank * D_tank
        L_tank = V_tank / A_tank
        while L_tank > 30:
            N_train += 1
            V_tank = V_tot / N_train
            L_tank = V_tank / A_tank

        D['N_train'] = self.N_train = N_train
        D['Aeration tank volume'] = self._V_tank = V_tank
        D['Aeration tank length'] = self.L_tank

        # Clarifier
        Q, Q_ras, X_v, SLR = self.Q, self.Q_ras, self.X_v, self.SLR
        X_mass = (Q+Q_ras)*X_v
        X_mass = auom('mg').convert(X_mass, 'lb')
        A_clarifier = D['Clarifier area'] = self._A_clarifier = X_mass/SLR/N_train

        # Blower and gas piping,
        # note that original code in ref 2 only considers O2,
        # here includes N2
        air_cfm = auom('m3/hr').convert(self.ins[1].F_vol, 'cfm')
        blower, piping = self.equipments
        blower.N_reactor = piping.N_reactor = N_train
        blower.gas_demand_per_reactor = piping.gas_demand_per_reactor = air_cfm
        self.add_equipment_design()

        # Concrete, modified from the MATLAB CAS codes to be consistent with
        # the membrane bioreactor units
        t_wall, t_slab = self.t_wall, self.t_slab
        W_N_trains = (W_tank+2*t_wall)*N_train - t_wall*(N_train-1)

        D_tot = D_tank + self.freeboard
        t = t_wall + t_slab

        get_VWC = lambda L1, N: N * t_wall * L1 * D_tot # get volume of wall concrete
        get_VSC = lambda L2: t * L2 * W_N_trains # get volume of slab concrete

        # Distribution channel, [ft3],
        # double for both the tank and the clarifier
        W_dist, W_eff = self.W_dist, self.W_eff
        VWC = 2 * get_VWC(L1=(W_N_trains+W_dist), N=2)
        VSC = 2 * get_VSC(L2=(W_dist+2*t_wall))

        # Aeration tanks, [ft3]
        VWC += get_VWC(L1=L_tank, N=(N_train+1))
        VSC += get_VSC(L2=L_tank)

        # Clarifiers, [ft3]
        L_clarifier = A_clarifier / W_tank
        VWC += get_VWC(L1=L_clarifier, N=(N_train+1))
        VSC += get_VSC(L2=L_clarifier)

        # Effluent channel, [ft3]
        VWC += 2 * get_VWC(L1=(W_N_trains+W_eff), N=2)
        VSC += 2 * get_VSC(L2=(W_eff+2*t_wall))

        # Pump/blower building, [ft3]
        W_PB, W_BB = self.W_PB, self.W_BB
        VWC += get_VWC(L1=(W_N_trains+W_PB+W_BB), N=2)
        VSC += get_VSC(L2=(W_PB+t_wall+W_BB))

        # Wet wells for mixed liquor storage, [ft3]
        L_WW, W_WW, D_WW = self.L_WW, self.W_WW, self.D_WW
        VWC += D_WW * (2*L_WW*t_wall+2*W_WW*t_wall+4*(t_wall**2))
        VSC += t * (L_WW+2*t_wall) * (W_WW+2*t_wall)

        D['Wall concrete'] = VWC
        D['Slab concrete'] = VSC

        # Excavation
        excav_slope, constr_access = self.excav_slope, self.constr_access
        # Aeration tank and clarifier
        VEX = calculate_excavation_volume(
            L=(W_dist+L_tank+L_clarifier), W=W_N_trains, D=D_tank,
            excav_slope=excav_slope, constr_access=constr_access)
        # Pump/blower building
        VEX += calculate_excavation_volume(
            L=(W_PB+W_BB), W=W_N_trains, D=D_tank,
            excav_slope=excav_slope, constr_access=constr_access)
        D['Excavation'] = VEX

        # Pumps
        pipe, pumps = self._design_pump()
        D['Pump pipe stainless steel'] = pipe
        D['Pump stainless steel'] = pumps


    def _design_pump(self):
        ID, pumps = self.ID, self.pumps
        inf, ras = self._inf, self._ras
        ins_dct = {
            'inf': inf,
            'recir': ras,
            }
        type_dct = dict.fromkeys(pumps, '')
        inputs_dct = dict.fromkeys(pumps, (self.N_train,))
        for i in pumps:
            if hasattr(self, f'{i}_pump'):
                p = getattr(self, f'{i}_pump')
                setattr(p, 'add_inputs', inputs_dct[i])
            else:
                ID = f'{ID}_{i}'
                capacity_factor=2. if i=='inf' else ras.F_vol/inf.F_vol
                pump = WWTpump(
                    ID=ID, ins=ins_dct[i], pump_type=type_dct[i],
                    Q_mgd=None, add_inputs=inputs_dct[i],
                    capacity_factor=capacity_factor,
                    include_pump_cost=True,
                    include_building_cost=False,
                    include_OM_cost=False,
                    )
                setattr(self, f'{i}_pump', pump)

        pipe_ss, pump_ss = 0., 0.
        for i in pumps:
            p = getattr(self, f'{i}_pump')
            p.simulate()
            p_design = p.design_results
            pipe_ss += p_design['Pump pipe stainless steel']
            pump_ss += p_design['Pump stainless steel']
        return pipe_ss, pump_ss


    def _cost(self):
        D, C = self.design_results, self.baseline_purchase_costs
        ### Capital ###
        # Concrete and excavation
        VEX, VWC, VSC = \
            D['Excavation'], D['Wall concrete'], D['Slab concrete']
        C['Reactor excavation'] = VEX * self.excav_unit_cost
        C['Wall concrete'] = VWC * self.wall_concrete_unit_cost
        C['Slab concrete'] = VSC * self.slab_concrete_unit_cost

        # Pump
        pumps, add_OPEX = self.pumps, self.add_OPEX
        pump_cost, building_cost, opex_o, opex_m = 0., 0., 0., 0.
        for i in pumps:
            p = getattr(self, f'{i}_pump')
            p_cost, p_add_opex = p.baseline_purchase_costs, p.add_OPEX
            pump_cost += p_cost['Pump']
            building_cost += p_cost['Pump building']
            opex_o += p_add_opex['Pump operating']
            opex_m += p_add_opex['Pump maintenance']

        C['Pumps'] = pump_cost
        C['Pump building'] = building_cost
        add_OPEX['Pump operating'] = opex_o
        add_OPEX['Pump maintenance'] = opex_m

        # Blower
        self.add_equipment_cost()

        ### Heat and power ###
        # Fluid heating
        hx = self.heat_exchanger
        ins0 = self.ins[0]
        # duty = ins0.thermo.mixture.H(mol=ins0.mol, phase='l', T=self.T, P=ins0.P) # save for later reference
        hx.ins[0].copy_flow(ins0)
        hx.outs[0].copy_flow(ins0)
        hx.ins[0].T = ins0.T
        hx.outs[0].T = self.T
        hx.ins[0].P = hx.outs[0].P = ins0.P
        hx.simulate_as_auxiliary_exchanger(ins=hx.ins, outs=hx.outs)
        # Power
        pumping = 0.
        for ID in self.pumps:
            p = getattr(self, f'{ID}_pump')
            if p is None:
                continue
            pumping += p.power_utility.rate
        self.power_utility.rate = self.blower.design_results['Total blower power'] + pumping


    @property
    def N_train(self):
        '''
        [int] Number of treatment train, should be at least two in case one failing.
        '''
        return self._N_train
    @N_train.setter
    def N_train(self, i):
        i = ceil(i)
        if i < 2:
            raise ValueError('`N_train` should be at least 2.')
        self._N_train = i

    @property
    def R(self):
        '''[float] Sludge recycle ratio.'''
        return self.X_v*(1-self.HRT/self.SRT)/(self.Xw-self.X_v) # eq. 6.14a

    @property
    def HRT(self):
        '''[float] Hydraulic retention time, [d].'''
        return self._HRT

    @property
    def A_clarifier(self):
        '''[float] Area of the clarifier, [ft2].'''
        return self._A_clarifier

    @property
    def Q(self):
        '''[float] Influent flow rate, [m3/d].'''
        return self._Q

    @property
    def Q_was(self):
        '''[float] Waste activated sludge flow rate, [m3/d].'''
        return self._Q_was

    @property
    def Q_ras(self):
        '''[float] Return activated sludge flow rate, [m3/d].'''
        return self._Q_ras

    @property
    def V_tank(self):
        '''[float] Volume of the aeration tank, [ft3].'''
        return self._V_tank

    @property
    def L_tank(self):
        '''[float] Length of the aeration tank, [ft].'''
        return self.V_tank/self.W_tank/self.D_tank

    @property
    def W_tank(self):
        '''[float] Width of the aeration tank, [ft].'''
        return self._W_tank
    @W_tank.setter
    def W_tank(self, i):
        self._W_tank = i

    @property
    def D_tank(self):
        '''[float] Depth of the aeration tank, [ft].'''
        return self._D_tank
    @D_tank.setter
    def D_tank(self, i):
        self._D_tank = i

    @property
    def W_dist(self):
        '''[float] Width of the distribution channel, [ft].'''
        return self._W_dist
    @W_dist.setter
    def W_dist(self, i):
        self._W_dist = i

    @property
    def W_eff(self):
        '''[float] Width of the effluent channel, [ft].'''
        return self._W_eff
    @W_eff.setter
    def W_eff(self, i):
        self._W_eff = i

    @property
    def freeboard(self):
        '''[float] Freeboard added to the depth of the reactor tank, [ft].'''
        return self._freeboard
    @freeboard.setter
    def freeboard(self, i):
        self._freeboard = i

    @property
    def W_PB(self):
        '''[float] Width of the pump building, [ft].'''
        return self._W_PB
    @W_PB.setter
    def W_PB(self, i):
        self._W_PB = i

    @property
    def L_BB(self):
        '''[float] Length of the blower building, [ft].'''
        return self._L_BB
    @L_BB.setter
    def L_BB(self, i):
        self._L_BB = i

    @property
    def W_BB(self):
        '''[float] Width of the blower building, [ft].'''
        return self._W_BB
    @W_BB.setter
    def W_BB(self, i):
        self._W_BB = i

    @property
    def L_WW(self):
        '''[float] Length of the wet well for mixed liquor storage, [ft].'''
        return self._L_WW
    @L_WW.setter
    def L_WW(self, i):
        self._L_WW = i

    @property
    def W_WW(self):
        '''[float] Width of the wet well for mixed liquor storage, [ft].'''
        return self._W_WW
    @W_WW.setter
    def W_WW(self, i):
        self._W_WW = i

    @property
    def D_WW(self):
        '''[float] Depth of the wet well for mixed liquor storage, [ft].'''
        return self._D_WW
    @D_WW.setter
    def D_WW(self, i):
        self._D_WW = i

    @property
    def t_wall(self):
        '''
        [float] Thickness of the wall concrete, [ft].
        default to be minimum of 1 ft with 1 in added for every ft of depth over 12 ft.
        '''
        return self._t_wall or (1 + max(self.D_tank-12, 0)/12)
    @t_wall.setter
    def t_wall(self, i):
        self._t_wall = i

    @property
    def t_slab(self):
        '''
        [float] Concrete slab thickness, [ft],
        default to be 2 in thicker than the wall thickness.
        '''
        return self._t_slab or self.t_wall + 2/12
    @t_slab.setter
    def t_slab(self, i):
        self._t_slab = i

    @property
    def excav_slope(self):
        '''[float] Slope for excavation (horizontal/vertical).'''
        return self._excav_slope
    @excav_slope.setter
    def excav_slope(self, i):
        self._excav_slope = i

    @property
    def constr_access(self):
        '''[float] Extra room for construction access, [ft].'''
        return self._constr_access
    @constr_access.setter
    def constr_access(self, i):
        self._constr_access = i
        

class TreatmentTrains(Mixer):
    '''
    Dummy unit with no run function of its own. To be used to calculate the 
    design and cost of treatment trains in ASP. Code largely dereived from code scripts 
    for [1]_.

    Parameters
    ----------
    ins : Iterable
        None expected (dummy unit).
    outs : Iterable
        None expected (dummy unit).
    N_train : int
        Number of treatment train, should be at least two in case one failing.
    V_tank : float
        Volume of tank (represents one treatment chain), [m3]. Default is 1000 m3. 
    W_tank : float
        Width of tank, [m]. Default is 6.4 m. [1]
    D_tank : float
        Depth of tank, [m]. Default is 3.66 m. [1]
    freeboard : float
        Freeboard added to the depth of the reactor tank, [m]. Default is 0.61 m. [1]
    W_dist : float
        Width of distribution channel, [m]. Default is 1.37 m. [1]
    W_eff: float
        Width of effluent channel, [m]. Default is 1.5 m. [1]
    W_recirculation: float
        Width of recirculation channel, [m]. Default is 1.5 m.  
    excav_slope: float
        Slope for excavation (horizontal/vertical). Default is 1.5. [1] 
    constr_access: float
        Extra room for construction access, [m]. Default is 0.92 m. [1]
    Q_air_design: float
        Used to calculate the power of blower. 
        Air flow required for one treatment train [m3/hr].
    Q_recirculation: float
        Used to calculate pumping power. 
        Design recirculated flow for one treatment train [m3/day]. ** confirm the logic of 'one' **
    kwargs : dict
        Other attributes to be set.

    References
    ----------
    .. [1] Shoener, B. D.; Zhong, C.; Greiner, A. D.; Khunjar, W. O.; Hong, P.-Y.; Guest, J. S.
        Design of Anaerobic Membrane Bioreactors for the Valorization
        of Dilute Organic Carbon Waste Streams.
        Energy Environ. Sci. 2016, 9 (3), 1102–1112.
        https://doi.org/10.1039/C5EE03715H.

    See Also
    --------
    `MATLAB codes <https://github.com/QSD-Group/AnMBR>`_ used in ref 1,
    especially the system layout `diagrams <https://github.com/QSD-Group/AnMBR/blob/main/System_Design.pdf>`_.

    '''
    
    # Costs
    wall_concrete_unit_cost = 1081.73 # $/m3 (Hydromantis. CapdetWorks 4.0. https://www.hydromantis.com/CapdetWorks.html)
    slab_concrete_unit_cost = 582.48 # $/m3 (Hydromantis. CapdetWorks 4.0. https://www.hydromantis.com/CapdetWorks.html)
    excav_unit_cost = (8 + 0.3) / 0.765 # $/m3, 0.765 is to convert from $/yd3 **NOT UPDATED** (taken from Shoener et al.)
    
    _t_wall = None
    _t_slab = None
    
    pumps = ('recirculation_CSTR',)
        
    def __init__(self, ID='', ins= None, outs= (), thermo=None, init_with='WasteStream',
                 F_BM_default=1, isdynamic=False, N_train=2, V_tank=1000, 
                 W_tank = 6.4, D_tank = 3.66, freeboard = 0.61, W_dist = 1.37, W_eff = 1.5, # all in meter  (converted from feet to m, Shoener et al.2016)
                 W_recirculation = 1.5, # in m (assumed same as W_eff)
                 excav_slope = 1.5, # horizontal/vertical (Shoener et al.2016)
                 constr_access = 0.92, # in meter  (converted from feet to m, Shoener et al. 2016)
                 Q_air_design  = 1000, # in m3/hr  **NO SOURCE FOR DEFAULT VALUE YET**
                 Q_recirculation = 1000, # in m3/day  **NO SOURCE FOR DEFAULT VALUE YET**
                 F_BM=default_F_BM, lifetime=default_equipment_lifetime, rigorous=False,
                 **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with,
                         F_BM_default=F_BM_default, isdynamic=isdynamic)
        
        self._effluent = self.outs[0].copy(f'{ID}_effluent')
        self.rigorous = rigorous
        self._mixed = WasteStream(f'{ID}_mixed')      
        self.N_train = N_train
        self.V_tank = V_tank
        self.W_tank = W_tank
        self.D_tank = D_tank
        self.freeboard = freeboard
        self.W_dist = W_dist
        self.W_eff = W_eff
        self.W_recirculation = W_recirculation
        self.excav_slope = excav_slope
        self.constr_access = constr_access
        self.Q_air_design = Q_air_design # **NO SOURCE FOR DEFAULT VALUE YET**
        self.Q_recirculation = Q_recirculation
        
        self.blower = blower = Blower('blower', linked_unit=self, N_reactor=N_train)
        self.air_piping = air_piping = GasPiping('air_piping', linked_unit=self, N_reactor=N_train)
        self.equipments = (blower, air_piping)
        self.F_BM.update(F_BM)
        
    @property
    def N_train(self):
        '''
        [int] Number of treatment train, should be at least two in case one failing.
        '''
        return self._N_train
    
    @N_train.setter
    def N_train(self, i):
        i = ceil(i)
        if i < 2:
            raise ValueError('`N_train` should be at least 2.')
        self._N_train = i
        
    @property
    def V_tank(self):
        '''[float] Volume of tank, [m3].'''
        return self._V_tank
    
    @V_tank.setter
    def V_tank(self, i):
        self._V_tank = i

    @property
    def W_tank(self):
        '''[float] Width of one tank, [m].'''
        return self._W_tank
    
    @W_tank.setter
    def W_tank(self, i):
        self._W_tank = i

    @property
    def D_tank(self):
        '''[float] Depth of one tank, [m].'''
        return self._D_tank
    
    @D_tank.setter
    def D_tank(self, i):
        self._D_tank = i

    @property
    def W_dist(self):
        '''[float] Width of the distribution channel, [m].'''
        return self._W_dist
    
    @W_dist.setter
    def W_dist(self, i):
        self._W_dist = i

    @property
    def W_eff(self):
        '''[float] Width of the effluent channel, [m].'''
        return self._W_eff
    
    @W_eff.setter
    def W_eff(self, i):
        self._W_eff = i
        
    @property
    def W_recirculation(self):
        '''[float] Width of the recirculation channel, [m].'''
        return self._W_recirculation
    
    @W_recirculation.setter
    def W_recirculation(self, i):
        self._W_recirculation = i

    @property
    def freeboard(self):
        '''[float] Freeboard added to the depth of the reactor tank, [m].'''
        return self._freeboard
    
    @freeboard.setter
    def freeboard(self, i):
        self._freeboard = i
        
    @property
    def t_wall(self):
        '''
        [float] Thickness of the wall concrete, [m].
        default to be minimum of 1 ft with 1 in added for every ft of depth over 12 ft.
        '''
        D_tank = self.D_tank*39.37 # m to inches 
        return self._t_wall or (1 + max(D_tank - 12, 0)/12)*0.3048 # from feet to m
    
    @t_wall.setter
    def t_wall(self, i):
        self._t_wall = i

    @property
    def t_slab(self):
        '''
        [float] Concrete slab thickness, [m],
        default to be 2 in thicker than the wall thickness.
        '''
        return self._t_slab or (self.t_wall + 2/12)*0.3048 # from feet to m
    
    @t_slab.setter
    def t_slab(self, i):
        self._t_slab = i

    @property
    def excav_slope(self):
        '''[float] Slope for excavation (horizontal/vertical).'''
        return self._excav_slope
    
    @excav_slope.setter
    def excav_slope(self, i):
        self._excav_slope = i
        
    @property
    def constr_access(self):
        '''[float] Extra room for construction access, [m].'''
        return self._constr_access
    
    @constr_access.setter
    def constr_access(self, i):
        self._constr_access = i
        
    @property
    def Q_air_design(self):
        '''[float] Air flow required for one treatment train, [m3/hr].'''
        return self._Q_air_design
    
    @Q_air_design.setter
    def Q_air_design(self, i):
        self._Q_air_design = i
    
    @property
    def Q_recirculation(self):
        '''[float] Design recirculated flow in the treatment train, [m3/day].'''
        return self._Q_recirculation
    
    @Q_recirculation.setter
    def Q_recirculation(self, i):
        self._Q_recirculation = i
    
    def _design_pump(self):
        
        ID, pumps = self.ID, self.pumps
        D = self.design_results
        
        self._effluent.copy_like(self.outs[0])
        effluent = self._effluent
        
        ins_dct = {
            'recirculation_CSTR': effluent,
            }
        
        Q_recirculation  = D['Q recirculation'] 
        meter_to_feet = 3.28
        Tank_length = D['Tank length']*meter_to_feet # in ft
        
        Q_recirculation_mgd = Q_recirculation*0.000264 #m3/day to MGD
        
        Q_mgd = {
            'recirculation_CSTR': Q_recirculation_mgd,
            }
      
        type_dct = dict.fromkeys(pumps, 'recirculation_CSTR')
        inputs_dct = dict.fromkeys(pumps, (self.N_train, Tank_length,))
       
        for i in pumps:
            if hasattr(self, f'{i}_pump'):
                p = getattr(self, f'{i}_pump')
                setattr(p, 'add_inputs', inputs_dct[i])
            else:
                ID = f'{ID}_{i}'
                capacity_factor=1
                pump = WWTpump(
                    ID=ID, ins= ins_dct[i], pump_type=type_dct[i],
                    Q_mgd=Q_mgd[i], add_inputs=inputs_dct[i],
                    capacity_factor=capacity_factor,
                    include_pump_cost=True,
                    include_building_cost=False,
                    include_OM_cost=True,
                    )
                setattr(self, f'{i}_pump', pump)

        pipe_ss, pump_ss = 0., 0.
        for i in pumps:
            p = getattr(self, f'{i}_pump')
            p.simulate()
            p_design = p.design_results
            pipe_ss += p_design['Pump pipe stainless steel']
            pump_ss += p_design['Pump stainless steel']
            
        return pipe_ss, pump_ss
    
    _units = {
        'Number of trains': 'ea',
        'Tank volume': 'm3',
        'HRT': 'hr',
        'Tank width': 'm',
        'Tank depth': 'm',
        'Tank length': 'm',
        'Wall concrete': 'm3',
        'Slab concrete': 'm3',
        'Excavation': 'm3',
        'Q recirculation': 'm3/day', 
        'Pump pipe stainless steel': 'kg',
        'Pump stainless steel': 'kg',
        }
    
    def _design(self):
        self._mixed.mix_from(self.ins)
        mixed = self._mixed
        
        D = self.design_results 
        
        D['Number of trains'] = self.N_train
        D['Tank volume'] = self.V_tank # in m3
        D['HRT'] = D['Tank volume']/mixed.get_total_flow('m3/hr') # in hr
        D['Tank width'] = self.W_tank # in m
        D['Tank depth'] = self.D_tank # in m 
        D['Tank length'] = D['Tank volume']/(D['Tank width']*D['Tank depth']) # in m 
        
        
        t_wall, t_slab = self.t_wall, self.t_slab
        W_N_trains = (D['Tank width'] + 2*t_wall)*D['Number of trains'] - t_wall*(D['Number of trains'] - 1)

        D_tot = D['Tank depth'] + self.freeboard
        t = t_wall + t_slab

        get_VWC = lambda L1, N: N * t_wall * L1 * D_tot # get volume of wall concrete
        get_VSC = lambda L2: t * L2 * W_N_trains # get volume of slab concrete
        
        # Aeration tanks, [m3]
        VWC = get_VWC(L1= D['Tank length'], N=(D['Number of trains'] + 1))
        VSC = get_VSC(L2= D['Tank length'])

        # Distribution channel, [m3]
        W_dist, W_eff, W_recirculation = self.W_dist, self.W_eff, self.W_recirculation
        VWC += get_VWC(L1=(W_N_trains+W_dist), N=2) # N =2 for two walls
        VSC += get_VSC(L2=(W_dist + 2*t_wall))

        # Effluent channel, [m3]
        VWC += get_VWC(L1=(W_N_trains + W_eff), N=2) # N =2 for two walls
        VSC += get_VSC(L2=(W_eff + 2*t_wall))

        # RAS channel, [m3]
        VWC += get_VWC(L1=(W_N_trains + W_recirculation), N=2) # N =2 for two walls
        VSC += get_VSC(L2=(W_recirculation + 2*t_wall))

        D['Wall concrete'] = VWC
        D['Slab concrete'] = VSC

        # Excavation
        excav_slope, constr_access = self.excav_slope, self.constr_access
        # Aeration tank
        VEX = calculate_excavation_volume(
            L=(W_dist + D['Tank length']), W = W_N_trains, D = D['Tank depth'],
            excav_slope=excav_slope, constr_access=constr_access)
        
        D['Excavation'] = VEX
        
        D['Q recirculation'] = self.Q_recirculation
        
        # Blower and gas piping (taken from 'ActivatedSludgeProcess' SanUnit)
        
        
        # oxygen_autotroph = get_oxygen_autotrophs(inf.F_vol*24, inf.COD, eff_sCOD, 
        #                                          inf_TKN, SRT=srt, Y_H=0.625, U_AUT=1.0, 
        #                                          b_H=0.4, b_AUT=0.15, f_d=0.1, SF_DO=1.25)
        
        # oxygen_autotroph = get_oxygen_autotrophs(inf.F_vol*24, inf.COD, eff_sCOD, 
        #                                          inf_TKN, SRT=srt, Y_H=0.625, U_AUT=1.0, 
        #                                          b_H=0.4, b_AUT=0.15, f_d=0.1, SF_DO=1.25)
        
        # airflow = get_airflow(oxygen_heterotrophs = oxygen_heterotroph, oxygen_autotrophs = oxygen_autotroph, 
        #                       oxygen_transfer_efficiency = 15) # in m3/min
        
        
        
        Q_air_design = self.Q_air_design # in m3/min
        air_cfm = auom('m3/min').convert(Q_air_design, 'cfm')
        blower, piping = self.equipments
        
        blower.N_reactor = piping.N_reactor = D['Number of trains']
        blower.gas_demand_per_reactor = piping.gas_demand_per_reactor = air_cfm
        self.add_equipment_design()
        
        # Pumps
        pipe, pumps = self._design_pump()
        D['Pump pipe stainless steel'] = pipe
        D['Pump stainless steel'] = pumps
        
    def _cost(self):
        D = self.design_results
        C = self.baseline_purchase_costs
        
        ### Capital ###
        # Concrete and excavation
        C['Wall concrete'] = D['Wall concrete'] * self.wall_concrete_unit_cost
        C['Slab concrete'] = D['Slab concrete'] * self.slab_concrete_unit_cost
        C['Reactor excavation'] = D['Excavation'] * self.excav_unit_cost

        # Pump
        pumps, add_OPEX = self.pumps, self.add_OPEX
        pump_cost, building_cost, opex_o, opex_m = 0., 0., 0., 0.
        for i in pumps:
            p = getattr(self, f'{i}_pump')
            p_cost, p_add_opex = p.baseline_purchase_costs, p.add_OPEX
            pump_cost += p_cost['Pump']
            building_cost += p_cost['Pump building']
            opex_o += p_add_opex['Pump operating']
            opex_m += p_add_opex['Pump maintenance']

        C['Pumps'] = pump_cost
        C['Pump building'] = building_cost
        add_OPEX['Pump operating'] = opex_o
        add_OPEX['Pump maintenance'] = opex_m

        # Blower
        self.add_equipment_cost()
        
        # Power
        pumping = 0.
        for ID in self.pumps:
            p = getattr(self, f'{ID}_pump')
            if p is None:
                continue
            pumping += p.power_utility.rate
            
        self.power_utility.rate = self.blower.design_results['Total blower power'] + pumping