#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

from math import ceil
from .. import SanUnit
from ..equipments import Blower, GasPiping
from ..utils import auom, calculate_excavation_volume


__all__ = ('ActivatedSludgeProcess',)

class ActivatedSludgeProcess(SanUnit):
    '''
    A steady-state activated sludge process containing rectangular
    continuously stirred tank reactors and clarifiers for settling
    based on [1]_ (largely followed Example 6.2) and code scripts for [2]_ .

    Single tank and clarifier are assumed to have the same width;
    pump/blower buildings are assumed to have the same width as the
    total tanks/clarifers.

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
        Solids loading rate for the clarifer, [lb/ft2/d].
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


    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 N_train=2, T=35+273.15, X_i0=None, X_v=1210, X_e=15, X_w=10000,
                 SLR=20, SF=4, aeration_power=1, # p337 of ref 1
                 q_hat=12, K=20, Y=0.5, b=0.396, f_d=0.8, # change from 0.2 of ref 2 based on p171 of ref 1
                 COD_factor=1.42, # p8 of Complex Systems in ref 1
                 q_UAP=1.8, q_BAP=0.1, k1=0.12, k2=0.09, K_UAP=100, K_BAP=85,
                 **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
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
        self.k_BAP = K_BAP

        self.equipments = (
            Blower('blower', linked_unit=self, N_reactor=N_train),
            GasPiping('gas piping', linked_unit=self, N_reactor=N_train),
            )

        for attr, value in kwargs.items():
            setattr(self, attr, value)


    # Equation/page numbers are noted for the 2001 Rittmann and McCarty book
    def _run(self):
        inf, air = self.ins
        eff, was, emission = self.outs
        air.phase = emission.phase = 'g'
        Q = auom('m3/hr').convert(inf.F_vol, 'L/d')
        self._Q = Q/1e3 # convert to m3/d
        S_0 = inf.COD # influent substrate concentration

        cmps = self.components
        X_i0 = self.X_i0
        if not X_i0:
            try:
                X_i0 = inf.concs[cmps.inert_biomass].sum()
            except AttributeError: # group not defined
                raise AttributeError('The `inert_biomass` group of current `CompiledComponents` '
                                     'has not bee defined. '
                                     'Please either define it using the `define_group` function of '
                                     '`CompiledComponents` or provide `X_i0` for this '
                                     f'`{self.__class__.__name__}` unit.')

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
        # (aeration tank and the clarifer)
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
        S_concs = eff.concs[substrates] * S/S_0
        eff_dct = dict.fromkeys(cmps.IDs, 0)
        eff_dct.pop('H2O')
        eff_dct.update({ID:S_concs[n] for n, ID in enumerate(substrates)})
        # Non-substrate and non-particulate (soluble, colloial, dissolved gas) components
        # will be splitted in the same way as flow rate mass-wise
        was_dct = eff_dct.copy()

        # Particulate components
        try: active_biomass = cmps.active_biomass
        except AttributeError:
            raise AttributeError('The `active_biomass` group of current `CompiledComponents` '
                                 'has not bee defined. '
                                 'Please either define it using the `define_group` function of '
                                 '`CompiledComponents`.')
        X_a_inf = inf.concs[active_biomass].sum()
        X_a_e_concs = eff.concs[active_biomass] * X_e/X_a_inf
        eff_dct.update({ID:X_a_e_concs[n] for n, ID in enumerate(active_biomass)})
        X_a_w_concs = eff.concs[active_biomass] * X_w/X_a_inf
        was_dct.update({ID:X_a_w_concs[n] for n, ID in enumerate(active_biomass)})

        eff.set_flow_by_concentration(
            flow_tot=Q_eff, concentrations=eff_dct, units=('L/d', 'mg/L'))
        was.set_flow_by_concentration(
            flow_tot=Q_was, concentrations=was_dct, units=('L/d', 'mg/L'))

        # Assume all non-active solids (including inorganic solids) go to sludge
        X_i = set(cmps.solids).difference(set(active_biomass))
        eff.imass[X_i] = 0
        was.imass[X_i] = inf.imass[X_i]

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
        D['Aeration tank volume'] = self.V_tank = V_tank
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

        D = D_tank + self.freeboard
        t = t_wall + t_slab

        get_VWC = lambda L1, N: N * t_wall * L1 * D # get volume of wall concrete
        get_VSC = lambda L2: t * L2 * W_N_trains # get volume of slab concrete

        # Distribution channel, [ft3],
        # double for both the tank and the clarifer
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


    def _cost(self):
        self.add_equipment_cost()
        #!!! Consider including pumps as equipment or similar,
        # note that the `capacity_factor` should be 2 for the intermediate pumps
        # and return sludge ratio

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
        '''[float] Returnactivated sludge flow rate, [m3/d].'''
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
        '''[float] Freeboard added to the depth of the reactor/membrane tank, [ft].'''
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
    def L_well(self):
        '''[float] Length of the wet well for mixed liquor storage, [ft].'''
        return self._L_well
    @L_well.setter
    def L_well(self, i):
        self._L_well = i

    @property
    def W_well(self):
        '''[float] Width of the wet well for mixed liquor storage, [ft].'''
        return self._W_well
    @W_well.setter
    def W_well(self, i):
        self._W_well = i

    @property
    def D_well(self):
        '''[float] Depth of the wet well for mixed liquor storage, [ft].'''
        return self._D_well
    @D_well.setter
    def D_well(self, i):
        self._D_well = i

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