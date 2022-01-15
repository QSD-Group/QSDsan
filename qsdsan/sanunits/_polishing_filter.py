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

from warnings import warn
from math import pi, ceil
from . import HXutility
from ._pumping import WWTpump
from .. import SanStream, SanUnit
from ..utils import (
    auom,
    compute_stream_COD,
    get_digestion_rxns,
    default_component_dict,
    calculate_excavation_volume,
    )

__all__ = ('PolishingFilter',)

degassing = SanStream.degassing

_ft_to_m = auom('ft').conversion_factor('m')
_ft2_to_m2 = auom('ft2').conversion_factor('m2')
_ft3_to_m3 = auom('ft3').conversion_factor('m3')
_ft3_to_gal = auom('ft3').conversion_factor('gallon')
_m3_to_gal = auom('m3').conversion_factor('gal')
_cmh_to_mgd = _m3_to_gal * 24 / 1e6 # cubic meter per hour to million gallon per day
_lb_to_kg = auom('lb').conversion_factor('kg')

_d_to_A = lambda d: pi/4*(d**2)
_A_to_d = lambda A: ((4*A)/pi)**0.5

F_BM_pump = 1.18*(1+0.007/100) # 0.007 is for miscellaneous costs
default_F_BM = {
        'Membrane': 1+0.15, # assume 15% for replacement labor
        'Pumps': F_BM_pump,
        'Pump building': F_BM_pump,
        }
default_equipment_lifetime = {
    'Pumps': 15,
    'Pump pipe stainless steel': 15,
    'Pump stainless steel': 15,
    }


# %%

class PolishingFilter(SanUnit):
    '''
    A superclass for anaerobic and aerobic polishing as in
    Shoener et al. [1]_ Some assumptions adopted from Humbird et al. [2]_

    Parameters
    ----------
    filter_type : str
        Can either be "anaerobic" or "aerobic".
    OLR : float
        Organic loading rate of influent, [kg COD/m3/hr].
    HLR : float
        Hydraulic loading rate of influent, [m3/m2/hr].
    X_decomp : float
        Fraction of the influent COD converted to biogas (`filter_type` == "anaerobic")
        or CO2 (`filter_type` == "aerobic").
    X_growth : float
        Fraction of the influent COD converted to biomass growth.
    solids : Iterable(str)
        IDs of the solid components.
        If not provided, will be set to the default `solids` attribute of the components.
    split : dict
        Component-wise split to the treated water.
        E.g., {'NaCl':1, 'WWTsludge':0} indicates all of the NaCl goes to
        the treated water and all of the WWTsludge goes to the wasted sludge.
        Default splits (based on the membrane bioreactor in [2]_) will be used
        if not provided.
        Note that the split for `Water` will be ignored as it will be adjusted
        to satisfy the `solids_conc` setting.
    biodegradability : float or dict
        Biodegradability of components,
        when shown as a float, all biodegradable components are assumed to have
        the same degradability.
    solids_conc : float
        Concentration of the biomass in the waste sludge, [g/L].
    T : float
        Temperature of the filter tank.
        Will not control temperature if provided as None.
    include_degassing_membrane : bool
        If to include a degassing membrane to enhance methane
        (generated through the digestion reaction) recovery.
        No degassing membrane will be added if `filter_type` is "aerobic".
    kwargs : dict
        Other keyword arguments (e.g., t_wall, t_slab).

    References
    ----------
    .. [1] Shoener et al., Design of Anaerobic Membrane Bioreactors for the
        Valorization of Dilute Organic Carbon Waste Streams.
        Energy Environ. Sci. 2016, 9 (3), 1102–1112.
        https://doi.org/10.1039/C5EE03715H.
    .. [2] Humbird et al., Process Design and Economics for Biochemical Conversion of
        Lignocellulosic Biomass to Ethanol: Dilute-Acid Pretreatment and Enzymatic
        Hydrolysis of Corn Stover; Technical Report NREL/TP-5100-47764;
        National Renewable Energy Lab (NREL), 2011.
        https://www.nrel.gov/docs/fy11osti/47764.pdf
    '''
    _N_ins = 3 # influent, recycle, air (optional)
    _N_outs = 4 # biogas (optional), effluent, waste sludge, air (optional)

    _N_filter_min = 2
    _d_max = 12
    _L_B = 50
    _W_B = 30
    _D_B = 10
    _t_wall = 8/12
    _t_slab = 1
    _excav_slope = 1.5
    _constr_access = 3

    # Heating
    T_air = 17 + 273.15
    T_earth = 10 + 273.15
    # Heat transfer coefficients, all in W/m2/°C
    H_wall = 0.7
    H_floor = 1.7
    H_ceilling = 0.95

    # Costs
    excav_unit_cost = (8+0.3) / 27 # $/ft3, 27 is to convert from $/yd3
    wall_concrete_unit_cost = 650 / 27 # $/ft3
    slab_concrete_unit_cost = 350 / 27 # $/ft3
    LDPE_unit_cost = 195 # $/m3
    HDPE_unit_cost = 195 # $/m3

    # Other equipment
    auxiliary_unit_names = ('heat_exchanger',)
    pumps =  ('lift', 'recir', 'eff', 'sludge')


    _units = {
        'Total volume': 'ft3',
        'Wall concrete': 'ft3',
        'Slab concrete': 'ft3',
        'Excavation': 'ft3',
        'Pump pipe stainless steel': 'kg',
        'Pump stainless steel': 'kg',
        'Packing LDPE': 'm3',
        'Packing HDPE': 'm3',
        }

    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream',
                 filter_type='aerobic',
                 OLR=(0.5+4)/24, # from the 0.5-4 kg/m3/d uniform range in ref [1]
                 HLR=(0.11+0.44)/2, # from the 0.11-0.44  uniform range in ref [1]
                 X_decomp=0.74, X_growth=0.22, # X_decomp & X_growth from ref [2]
                 biodegradability=1.,
                 split={}, solids=(), solids_conc=10.5, T=30+273.15,
                 include_degassing_membrane=False,
                 F_BM=default_F_BM, lifetime=default_equipment_lifetime,
                 **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with=init_with, F_BM_default=1)
        self.filter_type = filter_type
        self.OLR = OLR
        self.HLR = HLR
        self.X_decomp = X_decomp
        self.X_growth = X_growth
        self.biodegradability = biodegradability
        cmps = self.components
        self.split = split if split else default_component_dict(
            cmps=cmps, gas=0.15, solubles=0.125, solids=0) # ref[2]
        self.solids = solids or cmps.solids
        self.solids_conc = solids_conc
        self.T = T
        self.include_degassing_membrane = include_degassing_membrane
        self.F_BM.update(F_BM)
        self._default_equipment_lifetime.update(lifetime)

        # Initiate the attributes
        self.heat_exchanger = hx = HXutility(None, None, None, T=T)
        self.heat_utilities = hx.heat_utilities
        self._refresh_rxns()

        for k, v in kwargs.items():
            setattr(self, k, v)


    def _refresh_rxns(self, X_decomp=None, X_growth=None):
        X_decomp = X_decomp if X_decomp else self.X_decomp
        X_growth = X_growth if X_growth else self.X_growth

        xcmp_ID = self._xcmp.ID
        self._decomp_rxns = get_digestion_rxns(self.ins[0], 1.,
                                               X_decomp, 0., xcmp_ID)
        self._growth_rxns = get_digestion_rxns(self.ins[0], 1.,
                                               0., X_growth, xcmp_ID)
        self._i_rm = self._decomp_rxns.X + self._growth_rxns.X


    def _run(self):
        raw, recycled, air_in = self.ins
        biogas, eff, waste, air_out = self.outs

        # Initiate the streams
        biogas.phase = 'g'
        biogas.empty()

        inf = raw.copy()
        inf.mix_from((raw, recycled))
        self._inf = inf.copy() # this stream will be preserved (i.e., no reaction)

        self.growth_rxns(inf.mol)
        self.decomp_rxns.force_reaction(inf.mol)
        inf.split_to(eff, waste, self._isplit.data)
        # tmo.separations.split(inf, eff, waste, self._isplit.data)

        solids_conc = self._solids_conc
        m_solids = waste.imass[self.solids].sum()
        if m_solids/waste.F_vol <= solids_conc:
            diff = waste.ivol['Water'] - m_solids/solids_conc
            waste.ivol['Water'] = m_solids/solids_conc
            eff.ivol['Water'] += diff

        biogas.phase = air_in.phase = air_out.phase = 'g'

        if inf.imol['O2'] < 0:
            air_in.imol['O2'] = - inf.imol['O2']
            air_in.imol['N2'] = - 0.79/0.21 * inf.imol['O2']
            inf.imol['O2'] = 0
        else:
            air_in.empty()

        if self.filter_type == 'anaerobic':
            degassing(eff, biogas)
            degassing(waste, biogas)
            air_in.empty()
            air_out.empty()
        else:
            biogas.empty()
            degassing(eff, air_out)
            degassing(waste, air_out)
            air_out.imol['N2'] += air_in.imol['N2']
            self._recir_ratio = None

        if self.T is not None:
            biogas.T = eff.T = waste.T = air_out.T = self.T


    def _design(self):
        D = self.design_results
        func = self._design_anaerobic if self.filter_type=='anaerobic' \
            else self._design_aerobic

        ### Concrete and excavation ###
        V, VWC, VSC, VEX = func()
        D['Total volume'] = V
        D['Wall concrete'] = VWC
        D['Slab concrete'] = VSC
        D['Excavation'] = VEX

        ### Pumps ###
        pipe, pumps = self._design_pump()
        D['Pump pipe stainless steel'] = pipe
        D['Pump stainless steel'] = pumps

        ### Packing ###
        # Assume 50%/50% vol/vol LDPE/HDPE
        # 0.9 is void fraction, usually 85% - 95% for plastic packing media
        # 925 is density of LDPE (910-940), [kg/m3] (not used)
        # 950 is density of LDPE (930-970), [kg/m3] (not used)
        # M_LDPE_kg = 0.5 * (1-0.9) * 925 * V_m (not used)
        # M_HDPE_kg = 0.5 * (1-0.9) * 950 * V_m (not used)
        D['Packing LDPE'] = D['Packing HDPE'] = 0.05 * V

        ### Degassing ###
        D['Degassing membrane'] = self.N_degasser


    def _design_aerobic(self):
        inf, N = self._inf, self._N_filter_min
        Q = inf.F_vol

        ### Concrete ###
        get_V = lambda N: ((Q/N)*compute_stream_COD(inf, 'kg/m3')) / self.OLR # m3
        get_A = lambda N: Q/N/self.HLR

        V, A = get_V(N), get_A(N)
        d = _A_to_d(A)
        D = V / d # D is depth

        # Check if more than one filter is needed
        while d > self.d_max:
            N += 1
            V, A = get_V(N), get_A(N)
            d = _A_to_d(A)
            D = V / d

        self._OLR = ((Q/N)*compute_stream_COD(inf, 'kg/m3')) / V
        V_ft3 = V / _ft3_to_m3 * N
        d_ft = d / _ft_to_m
        D_ft = D / _ft_to_m
        self._N_filter, self._d, self._D = N, d, D_ft

        # Volume of wall/slab concrete, [ft3]
        VWC = self.t_wall * pi * d_ft * (D_ft+self.freeboard)
        VWC *= N
        VSC = 2 * self.t_slab * _d_to_A(d_ft)
        VSC *= N

        ### Excavation ###
        VEX = calculate_excavation_volume(
            self.L_B, self.W_B, self.D_B, self.excav_slope, self.constr_access)

        return V_ft3, VWC, VSC, VEX


    # def _design_anaerobic_filter(
    #         self, Q_mgd,
    #         Ss, # readily biodegradable (soluble) substrate concentration, [kg COD/m3]
    #         Sp, # slowly biodegradable (particulate) substrate concentration, [kg COD/m3]
    #         OLR_AF, # organic loading rate, [kg-COD/m3/day]
    #         HL_AF, # hydraulic loading rate, [m3/m2/hr]
    #         R_AF # recirculation ratio
    #         ):

    #     ### Filter material ###
    #     N_AF = 2
    #     Q_cmd = self.Q_cmd
    #     # Volume of the filter packing media in each filter, [m3]
    #     V_m_AF = (Q_cmd/N_AF) * (Ss+Sp) / OLR_AF
    #     # Diameter (d) / depth (D) of each filter, [m]
    #     d_AF, D_AF = _get_d_AF(Q_cmd, R_AF, N_AF, HL_AF, V_m_AF)

    #     while D_AF > 6: # assumed maximum depth assumption, [m]
    #         R_AF = R_AF + 0.1;
    #         d_AF, D_AF = _get_d_AF(Q_cmd, R_AF, N_AF, HL_AF, V_m_AF)

    #         while d_AF > 12: # assumed maximum diameter, [m]
    #             N_AF = N_AF + 1;
    #             d_AF, D_AF = _get_d_AF(Q_cmd, R_AF, N_AF, HL_AF)

    #     # Unit conversion
    #     d_AF /= _ft_to_m # [ft]
    #     D_AF /= _ft_to_m # [ft]
    #     V_m_AF /= _ft3_to_m3 # [ft3]

    #     ### Concrete material ###
    #     # External wall concrete, [ft3]
    #     # 6/12 is wall thickness and 3 is freeboard
    #     VWC_AF = N_AF * 6/12 * math.pi * d_AF * (D_AF+self.freeboard)
    #     VWC_AF *= N_AF
    #     # Floor slab concrete, [ft3]
    #     # 8/12 is slab thickness
    #     VSC_AF = _d_to_A(d_AF)+ 8/12 * _d_to_A(d_AF)
    #     VSC_AF *= N_AF

    #     ### Excavation ###
    #     SL = 1.5 # slope = horizontal/vertical
    #     CA = 3 # construction Access, [ft]
    #     #  Excavation of pump building
    #     PBL, PBW, PBD = 50, 30, 10 # pump building length, width, depth, [ft]
    #     Area_B_P = (PBL+2*CA) * (PBW+2*CA) # bottom area of frustum, [ft2]
    #     Area_T_P = (PBL+2*CA+PBW*SL) * (PBW+2*CA+PBD*SL) # top area of frustum, [ft2]
    #     VEX_PB = 0.5 * (Area_B_P+Area_T_P) * PBD # total volume of excavaion, [ft3]

    #     return N_AF, d_AF, D_AF, V_m_AF, VWC_AF, VWC_AF, VEX_PB

    def _design_pump(self):
        ID, ins, outs = self.ID, self.ins, self.outs
        N_filter, D, d, pumps = self.N_filter, self.D, self.d, self.pumps
        ins_dct = {
            'lift': ins[0].proxy(),
            'recir': ins[1].proxy(),
            'eff': outs[1].proxy(),
            'sludge': outs[2].proxy(),
            }
        type_dct = {
            'lift': 'lift',
            'recir': 'recirculation_AF',
            'eff': 'retentate_AF',
            'sludge': 'sludge',
            }
        inputs_dct = {
            'lift': (N_filter, D),
            'recir': (N_filter, d, D),
            'eff': (N_filter, D),
            'sludge': (1,),
            }
        for i in pumps:
            if hasattr(self, f'{i}_pump'):
                p = getattr(self, f'{i}_pump')
                setattr(p, 'add_inputs', inputs_dct[i])
            else:
                ID = f'{ID}_{i}'
                capacity_factor=1. if i!='recir' else self.recir_ratio
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
            pipe_ss += p.design_results['Pump pipe stainless steel']
            pump_ss += p.design_results['Pump stainless steel']
        return pipe_ss, pump_ss


    def _cost(self):
        D, C = self.design_results, self.baseline_purchase_costs

        ### Capital ###
        # Concrete and excavaction
        VEX, VWC, VSC = \
            D['Excavation'], D['Wall concrete'], D['Slab concrete']
        # 27 is to convert the VEX from ft3 to yard3
        C['Filter tank excavation'] = VEX * self.excav_unit_cost
        C['Wall concrete'] = VWC * self.wall_concrete_unit_cost
        C['Slab concrete'] = VSC * self.slab_concrete_unit_cost

        # Packing material
        C['Packing LDPE'] = D['Packing LDPE'] * self.LDPE_unit_cost
        C['Packing HDPE'] = D['Packing HDPE'] * self.HDPE_unit_cost

        # Pump
        pumps, add_OPEX = self.pumps, self.add_OPEX
        pump_cost, building_cost, opex_o, opex_m, pumping_power = 0., 0., 0., 0., 0.
        for i in pumps:
            p = getattr(self, f'{i}_pump')
            p_cost, p_add_opex = p.baseline_purchase_costs, p.add_OPEX
            pump_cost += p_cost['Pump']
            building_cost += p_cost['Pump building']
            opex_o += p_add_opex['Pump operating']
            opex_m += p_add_opex['Pump maintenance']
            pumping_power += p.power_utility.rate

        C['Pumps'] = pump_cost
        C['Pump building'] = building_cost
        add_OPEX['Pump operating'] = opex_o
        add_OPEX['Pump maintenance'] = opex_m

        # Degassing membrane
        C['Degassing membrane'] = 10000 * D['Degassing membrane']

        ### Heat and power ###
        T = self.T
        if T is None:
            loss = 0.
        else:
            N_filter, d, D = self.N_filter, self.d, self.D
            A_W = pi * d * D
            A_F = _d_to_A(d)
            A_W *= N_filter * _ft2_to_m2
            A_F *= N_filter * _ft2_to_m2

            loss = self.H_wall * (T-self.T_air) * A_W / 1e3
            loss += self.H_floor * (T-self.T_earth) * A_F / 1e3
            loss += self.H_ceiling * (T-self.T_air) * A_F / 1e3
        self._heat_loss = loss

        # Fluid heating
        inf = self._inf
        if T:
            H_at_T = inf.thermo.mixture.H(mol=inf.mol, phase='l', T=T, P=101325)
            duty = -(inf.H - H_at_T)
        else:
            duty = 0
        self.heat_exchanger.simulate_as_auxiliary_exchanger(duty, inf)

        # Degassing
        degassing = 3 * self.N_degasser # assume each uses 3 kW

        self.power_utility.rate = loss + pumping_power + degassing


    @property
    def filter_type(self):
        '''[str] Can either be "anaerobic" or "aerobic".'''
        return self._filter_type
    @filter_type.setter
    def filter_type(self, i):
        if i.lower() in ('anaerobic', 'aerobic'):
            self._filter_type = i.lower()
        else:
            raise ValueError('`filter_type` can only be "anaerobic" or "aerobic", '
                             f'not "{i}".')

    @property
    def OLR(self):
        '''[float] Organic loading rate, [kg COD/m3/hr].'''
        return self._OLR
    @OLR.setter
    def OLR(self, i):
        if i < 0:
            raise ValueError('`OLR` should be >=0, '
                             f'the input value {i} is outside the range.')
        self._OLR = i

    @property
    def HLR(self):
        '''[float] Hydraulic loading rate, [m3/m2/hr].'''
        return self._HLR
    @HLR.setter
    def HLR(self, i):
        self._HLR = i

    @property
    def d_max(self):
        '''[float] Maximum filter diameter, [m].'''
        return self._d_max
    @d_max.setter
    def d_max(self, i):
        self._d_max = i

    @property
    def d(self):
        '''[float] Diameter of the filter tank, [ft].'''
        return self._d

    @property
    def D(self):
        '''[float] Depth of the filter tank, [ft].'''
        return self._D

    @property
    def freeboard(self):
        '''[float] Freeboard added to the depth of the reactor tank, [ft].'''
        return self._freeboard
    @freeboard.setter
    def freeboard(self, i):
        self._freeboard = i

    @property
    def excav_slope(self):
        '''[float] Slope for excavation (horizontal/vertical).'''
        return self._excav_slope
    @excav_slope.setter
    def excav_slope(self, i):
        self._excav_slope = float(i)

    @property
    def constr_access(self):
        '''[float] Extra room for construction access, [ft].'''
        return self._constr_access
    @constr_access.setter
    def constr_access(self, i):
        self._constr_access = float(i)

    @property
    def N_filter(self):
        '''[int] Number of filter tanks.'''
        return self._N_filter

    @property
    def L_B(self):
        '''[float] Length of the reactor building, [ft].'''
        return max(self._L_B, self._d*self.N_filter)
    @L_B.setter
    def L_B(self, i):
        self._L_B = i

    @property
    def W_B(self):
        '''[float] Width of the reactor building, [ft].'''
        return max(self._W_B, self._d)
    @W_B.setter
    def W_B(self, i):
        self._W_B = i

    @property
    def D_B(self):
        '''[float] Depth of the reactor building, [ft].'''
        return max(self._D_B, self._D)
    @D_B.setter
    def D_B(self, i):
        self._D_B = i

    @property
    def t_wall(self):
        '''[float] Concrete wall thickness, [ft].'''
        return self._t_wall
    @t_wall.setter
    def t_wall(self, i):
        self._t_wall = float(i)

    @property
    def t_slab(self):
        '''[float] Concrete slab thickness, [ft].'''
        return self._t_slab
    @t_slab.setter
    def t_slab(self, i):
        self._t_slab = float(i)

    @property
    def N_degasser(self):
        '''
        [int] Number of degassing membrane needed for dissolved biogas removal
        (optional).
        '''
        if self.include_degassing_membrane:
            if self.filter_type=='aerobic':
                warn('No degassing membrane needed for when `filter_type` is "aerobic".')
                return 0
            return ceil(self.Q_cmd/24/30) # assume each can hand 30 m3/d of influent
        return 0

    @property
    def Q_mgd(self):
        '''
        [float] Influent volumetric flow rate in million gallon per day, [mgd].
        '''
        return self.ins[0].F_vol*_m3_to_gal*24/1e6

    # @property
    # def Q_gpm(self):
    #     '''[float] Influent volumetric flow rate in gallon per minute, [gpm].'''
    #     return self.Q_mgd*1e6/24/60

    @property
    def Q_cmd(self):
        '''
        [float] Influent volumetric flow rate in cubic meter per day, [cmd].
        '''
        return self.Q_mgd *1e6/_m3_to_gal # [m3/day]

    @property
    def recir_ratio(self):
        '''[float] Internal recirculation ratio.'''
        return self._recir_ratio or self.ins[1].F_vol/self.ins[0].F_vol
    @recir_ratio.setter
    def recir_ratio(self, i):
        self._recir_ratio = float(i)

    @property
    def i_rm(self):
        '''[:class:`np.array`] Removal of each component in this filter tank.'''
        return self._i_rm

    @property
    def biodegradability(self):
        '''
        [float of dict] Biodegradability of components,
        when shown as a float, all biodegradable component are assumed to have
        the same degradability.
        '''
        return self._biodegradability
    @biodegradability.setter
    def biodegradability(self, i):
        if isinstance(i, float):
            if not 0<=i<=1:
                raise ValueError('`biodegradability` should be within [0, 1], '
                                 f'the input value {i} is outside the range.')
            self._biodegradability = i
            return

        for k, v in i.items():
            if not 0<=v<=1:
                raise ValueError('`biodegradability` should be within [0, 1], '
                                 f'the input value for component "{k}" is '
                                 'outside the range.')
        self._biodegradability = i

    @property
    def split(self):
        '''Component-wise split to the treated water.'''
        return self._split
    @split.setter
    def split(self, i):
        self._split = i
        self._isplit = self.components.isplit(i, order=None)

    @property
    def solids_conc(self):
        '''Concentration of solids in the waste sludge, [g/L].'''
        return self._solids_conc
    @solids_conc.setter
    def solids_conc(self, i):
        self._solids_conc = i

    @property
    def X_decomp(self):
        '''
        [float] Fraction of the influent COD converted to biogas
        (`filter_type` == "anaerobic") or CO2 (`filter_type` == "aerobic").
        '''
        return self._X_decomp
    @X_decomp.setter
    def X_decomp(self, i):
        if not 0 <= i <= 1:
            raise ValueError('`X_decomp` should be within [0, 1], '
                             f'the input value {i} is outside the range.')
        self._X_decomp = i

    @property
    def X_growth(self):
        '''
        [float] Fraction of the influent COD converted to biomass growth.
        '''
        return self._X_growth
    @X_growth.setter
    def X_growth(self, i):
        if not 0 <= i <= 1:
            raise ValueError('`X_growth` should be within [0, 1], '
                             f'the input value {i} is outside the range.')
        self._X_growth = i

    @property
    def decomp_rxns(self):
        '''
        [:class:`tmo.ParallelReaction`] Organics to biogas (`filter_type` == "anaerobic")
        or CO2 (`filter_type` == "aerobic") reactions.
        '''
        return self._decomp_rxns

    @property
    def growth_rxns(self):
        '''
        [:class:`tmo.ParallelReaction`] Biomass (WWTsludge) growth reactions.
        '''
        return self._growth_rxns

    @property
    def organic_rm(self):
        '''[float] Overall organic (COD) removal rate.'''
        Qi, Qe = self._inf.F_vol, self.outs[1].F_vol
        Si = compute_stream_COD(self._inf, 'kg/m3')
        Se = compute_stream_COD(self.outs[1], 'kg/m3')
        return 1 - Qe*Se/(Qi*Si)