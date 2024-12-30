#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:

    Yalin Li <mailto.yalin.li@gmail.com>

    Joy Zhang <joycheung1994@gmail.com>


This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

from math import ceil, floor
from biosteam import Stream
from biosteam.exceptions import DesignError
from . import HXutility, WWTpump, InternalCirculationRx
from .. import SanStream, WasteStream, SanUnit, Process, CompiledProcesses
from ..sanunits import CSTR, PFR, dydt_cstr
from warnings import warn
import numpy as np, pandas as pd
from ..equipments import Blower
from ..utils import (
    auom,
    compute_stream_COD,
    format_str,
    default_component_dict,
    calculate_excavation_volume,
    )

__all__ = ('AnMBR', 
           'CompletelyMixedMBR',
           # 'PlugFlowMBR',
           )

#%%
degassing = SanStream.degassing

_ft_to_m = auom('ft').conversion_factor('m')
_ft2_to_m2 = auom('ft2').conversion_factor('m2')
_ft3_to_m3 = auom('ft3').conversion_factor('m3')
_ft3_to_gal = auom('ft3').conversion_factor('gallon')
_m3_to_gal = auom('m3').conversion_factor('gal')
_cmh_to_mgd = _m3_to_gal * 24 / 1e6 # cubic meter per hour to million gallon per day
_lb_to_kg = auom('lb').conversion_factor('kg')

F_BM_pump = 1.18*(1+0.007/100) # 0.007 is for miscellaneous costs
default_F_BM = {
        'Membrane': 1+0.15, # assume 15% for replacement labor
        'Pumps': F_BM_pump,
        'Pump building': F_BM_pump,
        }
default_equipment_lifetime = {
    'Membrane': 10,
    'Pumps': 15,
    'Pump pipe stainless steel': 15,
    'Pump stainless steel': 15,
    'Pump chemical storage HDPE': 30,
    }


# %%

class AnMBR(SanUnit):
    '''
    Anaerobic membrane bioreactor (AnMBR) for wastewater treatment as in
    Shoener et al. [1]_ Some assumptions adopted from Humbird et al. [2]_

    In addition to the anaerobic treatment, an optional second stage can be added,
    which can be aerobic filter or granular activated carbon (GAC).

    Parameters
    ----------
    ins : Inlets(obj)
        Influent, recycle (optional), naocl, citric acid, bisulfite, air (optional).
    outs : Outlets(obj)
        Biogas, effluent, waste sludge, air (optional).
    reactor_type : str
        Can either be "CSTR" for continuous stirred tank reactor
        or "AF" for anaerobic filter.
    N_train : int
        Number of treatment train, should be at least two in case one failing.
    membrane_configuration : str
        Can either be "cross-flow" or "submerged".
    membrane_type : str
        Can be "hollow fiber" ("submerged" configuration only),
        "flat sheet" (either "cross-flow" or "submerged" configuration),
        or "multi-tube" ("cross-flow" configuration only).
    membrane_material : str
        Can be any of the plastics ("PES", "PVDF", "PET", "PTFE")
        for any of the membrane types ("hollow fiber", "flat sheet", "multi-tube"),
        or "sintered steel" for "flat sheet",
        or "ceramic" for "multi-tube".
    membrane_unit_cost : float
        Cost of membrane, [$/ft2]
    include_aerobic_filter : bool
        Whether to include an aerobic filtration process in this AnMBR,
        can only be True in "AF" (not "CSTR") reactor.
    add_GAC : bool
        If to add granular activated carbon to enhance biomass retention,
        can only be True for the "submerged" configuration.
    include_degassing_membrane : bool
        If to include a degassing membrane to enhance methane
        (generated through the digestion reaction) recovery.
    Y_biogas : float
        Biogas yield, [kg biogas/kg consumed COD].
    Y_biomass : float
        Biomass yield, [kg biomass/kg consumed COD].
    biodegradability : float or dict
        Biodegradability of components,
        when shown as a float, all biodegradable components are assumed to have
        the same degradability.
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
    biomass_ID: str
        ID of the Component that represents the biomass.
    solids_conc : float
        Concentration of the biomass in the waste sludge, [g/L].
    T : float
        Temperature of the reactor.
        Will not control temperature if provided as None.
    kwargs : dict
        Other keyword arguments (e.g., J_max, SGD).

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

    .. note::

        TODO
            - Add algorithms for other configurations (AF, submerged, sparging, GAC, flat sheet, hollow fiber)
            - Maybe add AeMBR as well (and make an MBR superclass)
            
                - AeMBR can use higher flux and allows for lower transmembrane pressure

    See Also
    --------
    `MATLAB codes <https://github.com/QSD-Group/AnMBR>`_ used in ref 1,
    especially the system layout `diagrams <https://github.com/QSD-Group/AnMBR/blob/main/System_Design.pdf>`_.
    '''
    _N_ins = 6 # influent, recycle (optional), naocl, citric acid, bisulfite, air (optional)
    _N_outs = 4 # biogas, effluent, waste sludge, air (optional)

    # Equipment-related parameters
    _cas_per_tank_spare = 2

    _mod_surface_area = {
        'hollow fiber': 370,
        'flat sheet': 1.45/_ft2_to_m2,
        'multi-tube': 32/_ft2_to_m2
        }

    _mod_per_cas = None
    _mod_per_cas_range = {
        'hollow fiber': (30, 48), # min, max
        'flat sheet': (150, 200),
        'multi-tube': (44, 48)
        }

    _cas_per_tank = None
    _cas_per_tank_range = (16, 22)

    _N_blower = 0

    _W_tank = 21 # ft
    _D_tank = 12 # ft
    _freeboard = 2 # ft
    _W_dist = 4.5 # ft
    _W_eff = 4.5 # ft

    _L_well = 8 # ft
    _W_well = 8 # ft
    _D_well = 12 # ft

    _t_wall = None
    _t_slab = None

    _excav_slope = 1.5
    _constr_access = 3 # ft

    # Operation-related parameters
    _HRT = 10 # hr
    _J_max = 12
    _TMP_dct = {
        'cross-flow': 2.5,
        'submerged': 2.5,
        }
    _TMP_aerobic = None
    _recir_ratio = 2.25 # from the 0.5-4 uniform range in ref [1]
    _v_cross_flow = 1.2 # from the 0.4-2 uniform range in ref [1]
    _v_GAC = 8
    _SGD = 0.625 # from the 0.05-1.2 uniform range in ref [1]
    _AFF = 3.33

    # Heating
    T_air = 17 + 273.15
    T_earth = 10 + 273.15
    # Heat transfer coefficients, all in W/m2/°C
    heat_transfer_coeff=dict(wall=0.7, floor=1.7, ceiling=0.95),

    # Costs
    excav_unit_cost = (8+0.3) / 27 # $/ft3, 27 is to convert from $/yd3
    wall_concrete_unit_cost = 650 / 27 # $/ft3
    slab_concrete_unit_cost = 350 / 27 # $/ft3
    GAC_price = 13.78 # $/kg

    _refresh_rxns = InternalCirculationRx._refresh_rxns

    # Other equipment
    pumps =  ('perm', 'retent', 'recir', 'sludge', 'naocl', 'citric', 'bisulfite',
              'AF', 'AeF')
    auxiliary_unit_names = ('heat_exchanger',)

    _units = {
        'Total volume': 'ft3',
        'Wall concrete': 'ft3',
        'Slab concrete': 'ft3',
        'Excavation': 'ft3',
        'Membrane': 'm3',
        'Pump pipe stainless steel': 'kg',
        'Pump stainless steel': 'kg',
        'Pump chemical storage HDPE': 'm3',
        'Total air flow': 'CFM',
        'Blower capacity': 'CFM',
        'Packing LDPE': 'm3',
        'Packing HDPE': 'm3',
        'GAC': 'kg',
        }

    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', *,
                 reactor_type='CSTR',
                 N_train=2,
                 membrane_configuration='cross-flow',
                 membrane_type='multi-tube',
                 membrane_material='ceramic',
                 membrane_unit_cost=8,
                 include_aerobic_filter=False,
                 add_GAC=False,
                 include_degassing_membrane=True,
                 Y_biogas = 0.86,
                 Y_biomass=0.05, # from the 0.02-0.08 uniform range in ref [1]
                 biodegradability=1.0,
                 solids=(), split={},
                 biomass_ID='WWTsludge', solids_conc=10.5,
                 T=35+273.15,
                 F_BM=default_F_BM, lifetime=default_equipment_lifetime,
                 F_BM_default=1, **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with=init_with, F_BM_default=1)
        self._inf = inf = self.ins[0].copy() # this stream will be preserved (i.e., no reaction)
        # For pump design
        self._retent = inf.copy(f'{ID}_rentent')
        self._recir = inf.copy(f'{ID}_recir')
        self.reactor_type = reactor_type
        self.N_train = N_train
        self.include_aerobic_filter = include_aerobic_filter
        self.membrane_configuration = membrane_configuration
        self.membrane_type = membrane_type
        self.membrane_material = membrane_material
        self.membrane_unit_cost = membrane_unit_cost
        self.add_GAC = add_GAC
        self.include_degassing_membrane = include_degassing_membrane
        self.biodegradability = biodegradability
        self.Y_biogas = Y_biogas
        self.Y_biomass = Y_biomass
        cmps = self.components
        self.split = split if split else default_component_dict(
            cmps=cmps, gas=0.15, solubles=0.125, solids=0) # ref[2]
        self.solids = solids or cmps.solids
        self._xcmp = getattr(self.components, biomass_ID)
        self.solids_conc = solids_conc
        self.T = T
        self.F_BM.update(F_BM)
        self._default_equipment_lifetime.update(lifetime)

        # Initialize the attributes
        self.AF = self.AeF = None
        hx_in = Stream(f'{ID}_hx_in')
        hx_out = Stream(f'{ID}_hx_out')
        self.heat_exchanger = HXutility(ID=f'{ID}_hx', ins=hx_in, outs=hx_out)
        self._refresh_rxns()

        for k, v in kwargs.items(): setattr(self, k, v)

        blower = self.blower = Blower(ID+'_blower', linked_unit=self)
        self.equipments = (blower,)
        self._check_design()


    def _check_design(self):
        reactor_type = self.reactor_type
        m_config = self.membrane_configuration
        m_type = self.membrane_type
        m_material = self.membrane_material

        if reactor_type == 'CSTR':
            if self.include_aerobic_filter:
                raise DesignError('Aerobic filtration cannot be used in CSTR.')

        if m_config == 'submerged':
            if not m_type in ('hollow fiber', 'flat sheet'):
                raise DesignError('Only "hollow fiber" or "flat sheet" is allowed '
                                  'for "submerged" membrane, not {m_type}.')
        else: # cross-flow
            if not m_type in ('flat sheet', 'multi-tube'):
                raise DesignError('Only "flat sheet" or "multi-tube" is allowed '
                                  'for "cross-flow" membrane, not {m_type}.')
            if self.add_GAC:
                raise DesignError('No GAC should be added '
                                  '(i.e., `add_GAC` can only be False'
                                  'for "cross-flow" membrane.')

        plastics = ('PES', 'PVDF', 'PET', 'PTFE')
        if m_type == 'hollow fiber':
            if not m_material in plastics:
                raise DesignError(f'Only plastic materials {plastics} '
                                  'allowed for "hollow fiber" membrane',
                                  f'not "{m_material}".')
        elif m_type == 'flat sheet':
            if not m_material in (*plastics, 'sintered steel'):
                raise DesignError(f'Only plastic materials {plastics} and "sintered steel"'
                                  'allowed for "flat sheet" membrane',
                                  f'not "{m_material}".')
        else: # multi-tube
            if not m_material in (*plastics, 'ceramic'):
                raise DesignError(f'Only plastic materials {plastics} and "ceramic"'
                                  'allowed for "multi-tube" membrane',
                                  f'not "{m_material}".')


    # =========================================================================
    # _run
    # =========================================================================
    def _run(self):
        raw, recycled, naocl, citric, bisulfite, air_in = self.ins
        biogas, perm, sludge, air_out = self.outs

        # Initialize the streams
        biogas.phase = 'g'
        biogas.empty()

        inf = self._inf
        inf.mix_from((raw, recycled))

        # Chemicals for cleaning, assume all chemicals will be used up
        # 2.2 L/yr/cmd of 12.5 wt% solution (15% vol)
        naocl.empty()
        naocl.imass['NaOCl', 'Water'] = [0.125, 1-0.125]
        naocl.F_vol = (2.2/1e3/365/24) * (inf.F_vol*24) # m3/hr solution

        # 0.6 L/yr/cmd of 100 wt% solution, 13.8 lb/kg
        citric.empty()
        citric.ivol['CitricAcid'] = (0.6/1e3/365/24) * (inf.F_vol*24) # m3/hr pure

        # 0.35 L/yr/cmd of 38% solution, 3.5 lb/gal
        bisulfite.empty()
        bisulfite.imass['Bisulfite', 'Water'] = [0.38, 1-0.38]
        bisulfite.F_vol = (0.35/1e3/365/24) * (inf.F_vol*24) # m3/hr solution

        # For pump design
        self._compute_mod_case_tank_N()
        Q_R_mgd, Q_IR_mgd = self._compute_liq_flows()
        retent, recir = self._retent, self._recir
        retent.F_mass *= Q_R_mgd / self.Q_mgd
        recir.F_mass *= Q_IR_mgd / self.Q_mgd
        self._retent, self._recir = retent, recir

        # Effluents
        self.growth_rxns(inf.mol)
        self.biogas_rxns(inf.mol)
        inf.split_to(perm, sludge, self._isplit.data)

        solids_conc = self._solids_conc
        m_solids = sludge.imass[self.solids].sum()
        if m_solids/sludge.F_vol <= solids_conc:
            diff = sludge.ivol['Water'] - m_solids/solids_conc
            sludge.ivol['Water'] = m_solids/solids_conc
            perm.ivol['Water'] += diff

        degassing(perm, biogas)
        degassing(sludge, biogas)

        # Gas for sparging, no sparging needed if submerged or using GAC
        air_out.link_with(air_in)
        air_in.T = 17 + 273.15

        # self._design_blower()
        self.add_equipment_design()

        if self.T is not None: perm.T = sludge.T = biogas.T = air_out.T = self.T


    # Called by _run
    def _compute_mod_case_tank_N(self):
        N_mod_min, N_mod_max = self.mod_per_cas_range[self.membrane_type]
        N_cas_min, N_cas_max = self.cas_per_tank_range

        mod_per_cas, cas_per_tank = N_mod_min, N_cas_min

        J = self.J
        while J > self.J_max:
            mod_per_cas += 1
            if mod_per_cas == N_mod_max + 1:
                if cas_per_tank == N_cas_max + 1:
                    self.N_train += 1
                    mod_per_cas, cas_per_tank = N_mod_min, N_cas_min
                else:
                    cas_per_tank += 1
                    mod_per_cas = N_mod_min

        self._mod_per_cas, self._cas_per_tank = mod_per_cas, cas_per_tank


    # Called by _run
    def _compute_liq_flows(self):
        m_type = self.membrane_type
        if m_type == 'multi-tube':
            # Cross-flow flow rate per module,
            # based on manufacture specifications for compact 33, [m3/hr]
            Q_cross_flow = 53.5 * self.v_cross_flow
            Q_R_cmh = self.N_mod_tot * Q_cross_flow # total retentate flow rate, [m3/hr]
            Q_R_mgd = Q_R_cmh * _cmh_to_mgd # [mgd]

        Q_mgd, recir_ratio = self.Q_mgd, self.recir_ratio
        if Q_mgd*recir_ratio >= Q_R_mgd:
            Q_IR_mgd = Q_mgd*recir_ratio - Q_R_mgd
        else:
            Q_IR_mgd = 0
            # # Gives really large recir_ration,
            # # probably shouldn't back-calculate this way
            # self._recir_ratio = Q_R_mgd / Q_mgd

        if self.add_GAC:
            Q_upflow_req = (self.v_GAC/_ft_to_m) * \
                self.L_membrane_tank*self.W_tank*self.N_train * 24 * _ft3_to_gal / 1e6
            Q_IR_add_mgd = max(0, (Q_mgd+Q_IR_mgd)-Q_upflow_req)
            Q_IR_mgd += Q_IR_add_mgd
            self._recir_ratio = Q_IR_mgd / Q_mgd

        return Q_R_mgd, Q_IR_mgd


    # =========================================================================
    # _design
    # =========================================================================
    def _design(self):
        D = self.design_results
        D['Treatment train'] = self.N_train
        D['Cassette per train'] = self.cas_per_tank
        D['Module per cassette'] = self.mod_per_cas
        D['Total membrane modules'] = self.N_mod_tot

        # Step A: Reactor and membrane tanks
        # Call the corresponding design function
        # (_design_CSTR or _design_AF)
        func = getattr(self, f'_design_{self.reactor_type}')

        wall, slab, excavation = func()
        D['Wall concrete'] = wall
        D['Slab concrete'] = slab
        D['Excavation'] = excavation

        # Optional addition of packing media (used in filters)
        ldpe, hdpe = 0., 0.
        for i in (self.AF, self.AeF):
            if i is None:
                continue
            ldpe += i.design_results['Packing LDPE']
            hdpe += i.design_results['Packing HDPE']

        # Optional addition of GAC
        D['GAC'] = self._design_GAC()

        # Step B: Membrane
        # Call the corresponding design function
        # (_design_hollow_fiber, _design_flat_sheet, or _design_multi_tube)
        m_type = format_str(self.membrane_type)
        func = getattr(self, f'_design_{m_type}')
        D['Membrane'] = func()

        # Step C: Pumps
        pipe, pumps, hdpe = self._design_pump()
        D['Pump pipe stainless steel'] = pipe
        D['Pump stainless steel'] = pumps
        D['Pump chemical storage HDPE'] = hdpe

        # Step D: Degassing membrane
        D['Degassing membrane'] = self.N_degasser

        # Total volume
        D['Total volume'] = self.V_tot


    ### Step A functions ###
    # Called by _design
    def _design_CSTR(self):
        N_train = self.N_train
        W_dist, L_CSTR , W_eff, L_membrane_tank = \
            self.W_dist, self.L_CSTR, self.W_eff, self.L_membrane_tank
        W_PB, W_BB, D_tank = self.W_PB, self.W_BB, self.D_tank
        t_wall, t_slab = self.t_wall, self.t_slab
        SL, CA = self.excav_slope, self.constr_access

        ### Concrete calculation ###
        W_N_trains = (self.W_tank+2*t_wall)*N_train - t_wall*(N_train-1)

        D = D_tank + self.freeboard
        t = t_wall + t_slab

        get_VWC = lambda L1, N: N * t_wall * L1 * D
        get_VSC = lambda L2: t * L2 * W_N_trains

        # Concrete for distribution channel, [ft3]
        VWC_dist = get_VWC(L1=(W_N_trains+W_dist), N=2)
        VSC_dist = get_VSC(L2=(W_dist+2*t_wall))

        # Concrete for CSTR tanks, [ft3]
        VWC_CSTR = get_VWC(L1=L_CSTR, N=(N_train+1))
        VSC_CSTR = get_VSC(L2=L_CSTR)

        # Concrete for effluent channel, [ft3]
        VWC_eff = get_VWC(L1=(W_N_trains+W_eff), N=2)
        VSC_eff = get_VSC(L2=(W_eff+2*t_wall))

        # Concrete for the pump/blower building, [ft3]
        VWC_PBB = get_VWC(L1=(W_N_trains+W_PB+W_BB), N=2)
        VSC_PBB = get_VSC(L2=(W_PB+t_wall+W_BB))

        if self.membrane_configuration == 'submerged':
            VWC_membrane_tank, VSC_membrane_tank, VWC_well, VSC_well = \
                self._design_membrane_tank(self, D, N_train, W_N_trains,
                                           self.L_membrane_tank)
        else:
            VWC_membrane_tank, VSC_membrane_tank, VWC_well, VSC_well = 0., 0., 0., 0.

        # Total volume of wall concrete, [ft3]
        VWC = VWC_dist + VWC_CSTR + VWC_eff + VWC_PBB + VWC_membrane_tank + VWC_well

        # Total volume of slab concrete [ft3]
        VSC = VSC_dist + VSC_CSTR + VSC_eff + VSC_PBB + VSC_membrane_tank + VSC_well

        ### Excavation calculation ###
        # Excavation volume for the reactor and membrane tanks, [ft3]
        VEX_tanks = calculate_excavation_volume(
            L=(W_dist+L_CSTR+W_eff+L_membrane_tank),
            W=W_N_trains, D=D_tank, excav_slope=SL, constr_acess=CA)

        # Excavation volume for pump/blower building, [ft3]
        VEX_PBB = calculate_excavation_volume(
            L=(W_PB+W_BB),
            W=W_N_trains, D=D_tank, excav_slope=SL, constr_acess=CA)

        VEX = VEX_tanks + VEX_PBB

        #!!! Need to add wet wells for submerged configurations

        return VWC, VSC, VEX


    # Called by _design
    def _design_AF(self):
        '''NOT READY YET.'''
        # Use FilterTank


    # Called by _design_CSTR/_design_AF
    def _design_membrane_tank(self, D, N_train, W_N_trains, L_membrane_tank,
                              t_wall, t_slab):
        L_well, W_well, D_well = self.L_well, self.W_well, self.D_well

        # Concrete for membrane tanks, [ft3]
        t = t_wall + t_slab
        VWC_membrane_tank = (N_train+1) * t_wall * L_membrane_tank * D
        VSC_membrane_tank = t * L_membrane_tank * W_N_trains

        # Concrete for wet well (mixed liquor storage), [ft3]
        L = L_well + 2*t_wall
        W = W_well + 2*t_wall
        VWC_well = 2 * t_wall * (L_well+W) * D_well
        VSC_well = (t_slab+t_wall) * L * W

        return VWC_membrane_tank, VSC_membrane_tank, VWC_well, VSC_well


    # Called by _design
    def _design_GAC(self):
        '''NOT READY YET.'''
        if not self.add_GAC:
            return 0

        M_GAC = 1
        return M_GAC


    ### Step B functions ###
    # Called by _design
    def _design_hollow_fiber(self):
        '''NOT READY YET.'''


    # Called by _design
    def _design_flat_sheet(self):
        '''NOT READY YET.'''


    # Called by _design
    def _design_multi_tube(self):
        # # 0.01478 is volume of material for each membrane tube, [m3]
        # #      L_tube               OD        ID
        # V_tube = 3 * math.pi/4 * ((6e-3)**2-(5.2e-3)**2)
        # V_SU = 700 * V_tube # V for each small unit [m3]
        # M_SU = 1.78*1e3 * V_SU # mass = density*volume, [kg/m3], not used
        return self.N_mod_tot*0.01478


    ### Step C function ###
    # Called by _design
    def _design_pump(self):
        pumps = self.pumps
        ID, ins, outs = self.ID, self.ins, self.outs
        rx_type, m_config, pumps = \
            self.reactor_type, self.membrane_configuration, self.pumps
        self.AF_pump = self.AF.lift_pump if self.AF else None
        self.AeF_pump = self.AeF.lift_pump if self.AeF else None
        #!!! Maybe move `ins_dct` to `__init__` so it won't be repeated
        ins_dct = {
            'perm': outs[1].proxy(f'{ID}_perm'),
            'retent': self._retent,
            'recir': self._recir,
            'sludge': outs[2].proxy(f'{ID}_sludge'),
            'naocl': ins[2].proxy(f'{ID}_NaOCl'),
            'citric': ins[3].proxy(f'{ID}_citric'),
            'bisulfite': ins[4].proxy(f'{ID}_bisulfite'),
            }
        type_dct = {
            'perm': f'permeate_{m_config}',
            'retent': f'retentate_{rx_type}',
            'recir': f'recirculation_{rx_type}',
            'sludge': 'sludge',
            'naocl': 'chemical',
            'citric': 'chemical',
            'bisulfite': 'chemical',
            }
        inputs_dct = {
            'perm': (self.cas_per_tank, self.D_tank, self.TMP_anaerobic,
                     self.include_aerobic_filter),
            'retent': (self.cas_per_tank,),
            'recir': (1, self.L_CSTR,),
            'sludge': (1,),
            'naocl': (1,),
            'citric': (1,),
            'bisulfite': (1,),
            }

        for i in pumps[:-2]:
            if hasattr(self, f'{i}_pump'):
                p = getattr(self, f'{i}_pump')
                setattr(p, 'add_inputs', inputs_dct[i])
            else:
                # Add '.' in ID for auxiliary units
                ID = f'.{ID}_{i}'
                capacity_factor=2. if i=='perm' else self.recir_ratio if i=='recir' else 1.
                pump = WWTpump(
                    ID=ID, ins=ins_dct[i], pump_type=type_dct[i],
                    Q_mgd=None, add_inputs=inputs_dct[i],
                    capacity_factor=capacity_factor,
                    include_pump_cost=True,
                    include_building_cost=False,
                    include_OM_cost=False,
                    )
                setattr(self, f'{i}_pump', pump)

        pipe_ss, pump_ss, hdpe = 0., 0., 0.
        for i in (*pumps, 'AF', 'AeF'):
            p = getattr(self, f'{i}_pump')
            if p == None: continue
            p.simulate()
            p_design = p.design_results
            pipe_ss += p_design['Pump pipe stainless steel']
            pump_ss += p_design['Pump stainless steel']
            hdpe += p_design['Chemical storage HDPE']
        return pipe_ss, pump_ss, hdpe


    # =========================================================================
    # _cost
    # =========================================================================
    def _cost(self):
        D, C = self.design_results, self.baseline_purchase_costs
        ### Capital ###
        # Concrete and excavation
        VEX, VWC, VSC = \
            D['Excavation'], D['Wall concrete'], D['Slab concrete']
        C['Reactor excavation'] = VEX * self.excav_unit_cost
        C['Wall concrete'] = VWC * self.wall_concrete_unit_cost
        C['Slab concrete'] = VSC * self.slab_concrete_unit_cost

        # Membrane
        C['Membrane'] = self.membrane_unit_cost * D['Membrane'] / _ft2_to_m2

        # GAC
        C['GAC'] = self.GAC_price * D['GAC']

        # Packing material
        ldpe, hdpe = 0., 0.
        for i in (self.AF, self.AeF):
            if i is None:
                continue
            ldpe += i.baseline_purchase_costs['Packing LDPE']
            hdpe += i.baseline_purchase_costs['Packing HDPE']

        # Pump
        pumps, add_OPEX = self.pumps, self.add_OPEX
        pump_cost, building_cost, opex_o, opex_m = 0., 0., 0., 0.
        for i in pumps:
            p = getattr(self, f'{i}_pump')
            if p == None:
                continue
            p_cost, p_add_opex = p.baseline_purchase_costs, p.add_OPEX
            pump_cost += p_cost['Pump']
            building_cost += p_cost['Pump building']
            opex_o += p_add_opex['Pump operating']
            opex_m += p_add_opex['Pump maintenance']

        C['Pumps'] = pump_cost
        C['Pump building'] = building_cost
        add_OPEX['Pump operating'] = opex_o
        add_OPEX['Pump maintenance'] = opex_m

        # Degassing membrane
        C['Degassing membrane'] = 10000 * D['Degassing membrane']

        # Blower
        self.add_equipment_cost()

        ### Heat and power ###
        # Heat loss
        T = self.T
        coeff = self.heat_transfer_coeff
        
        if T is None: loss = 0.
        else:
            N_train, L_CSTR, W_tank, D_tank = \
                self.N_train, self.L_CSTR, self.W_tank, self.D_tank
            A_W = 2 * (L_CSTR+W_tank) * D_tank
            A_F = L_CSTR * W_tank
            A_W *= N_train * _ft2_to_m2
            A_F *= N_train * _ft2_to_m2
            loss = coeff['wall'] * (T-self.T_air) * A_W # [W]
            loss += coeff['floor'] * (T-self.T_earth) # [W]
            loss += coeff['ceiling'] * (T-self.T_air) # [W]
            loss *= 60*60/1e3 # W (J/s) to kJ/hr
        
        # Stream heating
        hx = self.heat_exchanger
        inf = self._inf
        hx_ins0, hx_outs0 = hx.ins[0], hx.outs[0]
        hx_ins0.copy_flow(inf)
        hx_outs0.copy_flow(inf)
        hx_ins0.T = inf.T
        hx_outs0.T = T
        hx.H = hx_outs0.H + loss # stream heating and heat loss
        hx.simulate_as_auxiliary_exchanger(ins=hx.ins, outs=hx.outs)
        
        # Power for pumping and gas
        pumping = 0.
        for ID in self.pumps: #!!! check if cost/power of AF_pump/AeF_pump included in AF/AeF
            p = getattr(self, f'{ID}_pump')
            if p is None: continue
            pumping += p.power_utility.rate
        sparging = 0. #!!! output from submerge design
        degassing = 3 * self.N_degasser # assume each uses 3 kW
        self.power_utility.rate = self.blower.power_utility.rate + \
            sparging + degassing + pumping


    ### Reactor configuration ###
    @property
    def reactor_type(self):
        '''
        [str] Can either be "CSTR" for continuous stirred tank reactor
        or "AF" for anaerobic filter.
        '''
        return self._reactor_type
    @reactor_type.setter
    def reactor_type(self, i):
        if not i.upper() in ('CSTR', 'AF'):
            raise ValueError('`reactor_type` can only be "CSTR", or "AF", '
                             f'not "{i}".')
        self._reactor_type = i.upper()

    @property
    def membrane_configuration(self):
        '''[str] Can either be "cross-flow" or "submerged".'''
        return self._membrane_configuration
    @membrane_configuration.setter
    def membrane_configuration(self, i):
        i = 'cross-flow' if i.lower() in ('cross flow', 'crossflow') else i
        if not i.lower() in ('cross-flow', 'submerged'):
            raise ValueError('`membrane_configuration` can only be "cross-flow", '
                             f'or "submerged", not "{i}".')
        self._membrane_configuration = i.lower()

    @property
    def membrane_type(self):
        '''
        [str] Can be "hollow fiber" ("submerged" configuration only),
        "flat sheet" (either "cross-flow" or "submerged" configuration),
        or "multi-tube" ("cross-flow" configuration only).
        '''
        return self._membrane_type
    @membrane_type.setter
    def membrane_type(self, i):
        i = 'multi-tube' if i.lower() in ('multi tube', 'multitube') else i
        if not i.lower() in ('hollow fiber', 'flat sheet', 'multi-tube'):
            raise ValueError('`membrane_type` can only be "hollow fiber", '
                             f'"flat sheet", or "multi-tube", not "{i}".')
        self._membrane_type = i.lower()

    @property
    def membrane_material(self):
        '''
        [str] Can be any of the plastics ("PES", "PVDF", "PET", "PTFE")
        for any of the membrane types ("hollow fiber", "flat sheet", "multi-tube"),
        or "sintered steel" for "flat sheet",
        or "ceramic" for "multi-tube".
        '''
        return self._membrane_material
    @membrane_material.setter
    def membrane_material(self, i):
        plastics = ('PES', 'PVDF', 'PET', 'PTFE')
        if i.upper() in plastics:
            self._membrane_material = i.upper()
        elif i.lower() in ('sintered steel', 'ceramic'):
            self._membrane_material = i.lower()
        else:
            raise ValueError(f'`membrane_material` can only be plastics materials '
                             f'{plastics}, "sintered steel", or "ceramic", not {i}.')


    ### Reactor/membrane tank ###
    @property
    def AF(self):
        '''[:class:`~.FilterTank`] Anaerobic filter tank.'''
        if self.reactor_type == 'CSTR':
            return None
        return self._AF

    @property
    def AeF(self):
        '''[:class:`~.FilterTank`] Aerobic filter tank.'''
        if not self.include_aerobic_filter:
            return None
        return self._AeF

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
    def cas_per_tank_spare(self):
        '''[int] Number of spare cassettes per train.'''
        return self._cas_per_tank_spare
    @cas_per_tank_spare.setter
    def cas_per_tank_spare(self, i):
        self._cas_per_tank_spare = ceil(i)

    @property
    def mod_per_cas_range(self):
        '''
        [tuple] Range (min, max) of the number of membrane modules per cassette
        for the current membrane type.
        '''
        return self._mod_per_cas_range
    @mod_per_cas_range.setter
    def mod_per_cas_range(self, i):
        self._mod_per_cas_range[self.membrane_type] = \
            tuple(floor(i[0]), floor(i[1]))

    @property
    def mod_per_cas(self):
        '''
        [float] Number of membrane modules per cassette for the current membrane type.
        '''
        return self._mod_per_cas or self._mod_per_cas_range[self.membrane_type][0]

    @property
    def cas_per_tank_range(self):
        '''
        [tuple] Range (min, max) of the number of membrane cassette per tank
        (same for all membrane types).
        '''
        return self._cas_per_tank_range
    @cas_per_tank_range.setter
    def cas_per_tank_range(self, i):
        self._cas_per_tank_range = tuple(floor(i[0]), floor(i[1]))

    @property
    def cas_per_tank(self):
        '''
        [float] Number of membrane cassettes per tank for the current membrane type.
        '''
        return self._cas_per_tank or self._cas_per_tank_range[0]

    @property
    def N_mod_tot(self):
        '''[int] Total number of membrane modules.'''
        return self.N_train * self.cas_per_tank * self.mod_per_cas

    @property
    def mod_surface_area(self):
        '''
        [float] Surface area of the membrane for the current membrane type, [m2/module].
        Note that one module is one sheet for plat sheet and one tube for multi-tube.
        '''
        return self._mod_surface_area[self.membrane_type]
    @mod_surface_area.setter
    def mod_surface_area(self, i):
        self._mod_surface_area[self.membrane_type] = i

    @property
    def L_CSTR(self):
        '''[float] Length of the CSTR tank, [ft].'''
        if self.reactor_type == 'AF':
            return 0
        return self._inf.F_vol/_ft3_to_m3*self.HRT/(self.N_train*self.W_tank*self.D_tank)

    @property
    def L_membrane_tank(self):
        '''[float] Length of the membrane tank, [ft].'''
        if self.membrane_configuration=='cross-flow':
            return 0.
        return ceil((self.cas_per_tank+self.cas_per_tank_spare)*3.4)

    @property
    def W_tank(self):
        '''[float] Width of the reactor/membrane tank (same value), [ft].'''
        return self._W_tank
    @W_tank.setter
    def W_tank(self, i):
        self._W_tank = i

    @property
    def D_tank(self):
        '''[float] Depth of the reactor/membrane tank (same value), [ft].'''
        return self._D_tank
    @D_tank.setter
    def D_tank(self, i):
        self._D_tank = i

    @property
    def freeboard(self):
        '''[float] Freeboard added to the depth of the reactor/membrane tank, [ft].'''
        return self._freeboard
    @freeboard.setter
    def freeboard(self, i):
        self._freeboard = i

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
    def V_tot(self):
        '''[float] Total volume of the unit, [ft3].'''
        return  self.D_tank*self.W_tank*self.L_CSTR*self.N_train

    @property
    def OLR(self):
        '''[float] Organic loading rate, [kg COD/m3/hr].'''
        return compute_stream_COD(self.ins[0], 'kg/m3')*self.ins[0].F_vol/(self.V_tot*_ft3_to_m3)


    ### Pump/blower ###
    @property
    def N_blower(self):
        '''
        [int] Number of blowers needed for gas sparging
        (not needed for some designs).
        Note that this is not used in costing
        (the cost is estimated based on the total sparging gas need).
        '''
        if not self.add_GAC and self.membrane_configuration=='submerged':
            return self._N_blower
        return 0

    @property
    def N_degasser(self):
        '''
        [int] Number of degassing membrane needed for dissolved biogas removal
        (optional).
        '''
        if self.include_degassing_membrane:
            return ceil(self.Q_cmd/24/30) # assume each can hand 30 m3/d of influent
        return 0

    @property
    def W_PB(self):
        '''[float] Width of the pump building, [ft].'''
        if self.membrane_configuration == 'submerged':
            N = self.cas_per_tank
        else: # cross-flow
            N = ceil(self.L_CSTR/((1+8/12)+(3+4/12)))

        if 0 <= N <= 10:
            W_PB = 27 + 4/12
        elif 11 <= N <= 16:
            W_PB = 29 + 6/12
        elif 17 <= N <= 22:
            W_PB = 31 + 8/12
        elif 23 <= N <= 28:
            W_PB = 35
        elif N >= 29:
            W_PB = 38 + 4/12
        else:
            W_PB = 0

        return W_PB

    @property
    def L_BB(self):
        '''[float] Length of the blower building, [ft].'''
        if self.membrane_configuration == 'submerged':
            return (69+6/12) if self.cas_per_tank<=18 else (76+8/12)
        return 0

    @property
    def W_BB(self):
        '''[float] Width of the blower building, [ft].'''
        if self.membrane_configuration == 'submerged':
            return (18+8/12) if self.cas_per_tank<=18 else 22
        return 0


    ### Wet well (submerged only) ###
    @property
    def L_well(self):
        '''
        [float] Length of the wet well, [ft].
        Only needed for submerged configuration.
        '''
        return self._L_well if self.membrane_configuration == 'submerged' else 0
    @L_well.setter
    def L_well(self, i):
        self._L_well = i

    @property
    def W_well(self):
        '''
        [float] Width of the wet well, [ft].
        Only needed for submerged configuration.
        '''
        return self._W_well if self.membrane_configuration == 'submerged' else 0
    @W_well.setter
    def W_well(self, i):
        self._W_well = i

    @property
    def D_well(self):
        '''
        [float] Depth of the wet well, [ft].
        Only needed for submerged configuration.
        '''
        return self._D_well if self.membrane_configuration == 'submerged' else 0
    @D_well.setter
    def D_well(self, i):
        self._D_well = i

    @property
    def t_wall(self):
        '''
        [float] Concrete wall thickness, [ft].
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
        return self._t_slab or self.t_wall+2/12
    @t_slab.setter
    def t_slab(self, i):
        self._t_slab = i


    ### Excavation ###
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


    ### Operation-related parameters ###
    @property
    def Q_mgd(self):
        '''
        [float] Influent volumetric flow rate in million gallon per day, [mgd].
        '''
        return self.ins[0].F_vol*_m3_to_gal*24/1e6

    @property
    def Q_gpm(self):
        '''[float] Influent volumetric flow rate in gallon per minute, [gpm].'''
        return self.Q_mgd*1e6/24/60

    @property
    def Q_cmd(self):
        '''
        [float] Influent volumetric flow rate in cubic meter per day, [cmd].
        '''
        return self.Q_mgd *1e6/_m3_to_gal # [m3/day]

    @property
    def Q_cfs(self):
        '''[float] Influent volumetric flow rate in cubic feet per second, [cfs].'''
        return self.Q_mgd*1e6/24/60/60/_ft3_to_gal

    @property
    def HRT(self):
        '''
        [float] Hydraulic retention time, [hr].
        '''
        return self._HRT
    @HRT.setter
    def HRT(self, i):
        self._HRT = i

    @property
    def recir_ratio(self):
        '''
        [float] Internal recirculation ratio, will be updated in simulation
        if the originally set ratio is not adequate for the desired flow
        required by GAC (if applicable).
        '''
        return self._recir_ratio
    @recir_ratio.setter
    def recir_ratio(self, i):
        self._recir_ratio = i

    @property
    def J_max(self):
        '''[float] Maximum membrane flux, [L/m2/hr].'''
        return self._J_max
    @J_max.setter
    def J_max(self, i):
        self._J_max = i

    @property
    def J(self):
        '''[float] Membrane flux, [L/m2/hr].'''
        # Based on the flux of one train being offline
        SA = (self.N_train-1) * self.cas_per_tank * self.mod_per_cas * self.mod_surface_area
        return self._inf.F_vol*1e3/SA # 1e3 is conversion from m3 to L

    @property
    def TMP_anaerobic(self):
        '''[float] Transmembrane pressure in the anaerobic reactor, [psi].'''
        return self._TMP_dct[self.membrane_configuration]
    @TMP_anaerobic.setter
    def TMP_anaerobic(self, i):
        self._TMP_dct[self.membrane_configuration] = i

    @property
    def TMP_aerobic(self):
        '''
        [float] Transmembrane pressure in the aerobic filter, [psi].
        Defaulted to half of the reactor TMP.
        '''
        if not self._include_aerobic_filter:
            return 0.
        else:
            return self._TMP_aerobic or self._TMP_dct[self.membrane_configuration]/2
    @TMP_aerobic.setter
    def TMP_aerobic(self, i):
        self._TMP_aerobic = i

    @property
    def SGD(self):
        '''[float] Specific gas demand, [m3 gas/m2 membrane area/h].'''
        return self._SGD
    @SGD.setter
    def SGD(self, i):
        self._SGD = i

    @property
    def AFF(self):
        '''
        [float] Air flow fraction, used in air pipe costing.
        The default value is calculated as STE/6
        (STE stands for standard oxygen transfer efficiency, and default STE is 20).
        If using different STE value, AFF should be 1 if STE/6<1
        and 3.33 if STE/6>1.
        '''
        return self._AFF
    @AFF.setter
    def AFF(self, i):
        self._AFF = i

    @property
    def v_cross_flow(self):
        '''[float] Cross-flow velocity, [m/s].'''
        return self._v_cross_flow if self.membrane_configuration=='cross-flow' else 0
    @v_cross_flow.setter
    def v_cross_flow(self, i):
        self._v_cross_flow = i

    @property
    def v_GAC(self):
        '''
        [float] Upflow velocity for GAC bed expansion, [m/hr].
        '''
        return self._v_GAC if self.add_GAC==True else 0
    @v_GAC.setter
    def v_GAC(self, i):
        self._v_GAC = i

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
        if not isinstance(i, dict):
            if not 0<=i<=1:
                raise ValueError('`biodegradability` should be within [0, 1], '
                                 f'the input value {i} is outside the range.')
            self._biodegradability = dict.fromkeys(self.chemicals.IDs, i)
        else:
            for k, v in i.items():
                if not 0<=v<=1:
                    raise ValueError('`biodegradability` should be within [0, 1], '
                                     f'the input value for chemical "{k}" is '
                                     'outside the range.')
            self._biodegradability = dict.fromkeys(self.chemicals.IDs, i).update(i)
        self._refresh_rxns()

    @property
    def biomass_ID(self):
        '''[str] ID of the Component that represents the biomass.'''
        return self._xcmp.ID
    @biomass_ID.setter
    def biomass_ID(self, i):
        self._xcmp = getattr(self.components, i)

    @property
    def solids_conc(self):
        '''Concentration of solids in the waste sludge, [g/L].'''
        return self._solids_conc
    @solids_conc.setter
    def solids_conc(self, i):
        self._solids_conc = i

    @property
    def i_rm(self):
        '''[:class:`np.array`] Removal of each chemical in this reactor.'''
        return self._i_rm

    @property
    def split(self):
        '''Component-wise split to the treated water.'''
        return self._split
    @split.setter
    def split(self, i):
        self._split = i
        self._isplit = self.chemicals.isplit(i, order=None)

    @property
    def biogas_rxns(self):
        '''
        [:class:`tmo.ParallelReaction`] Biogas production reactions.
        '''
        return self._biogas_rxns

    @property
    def growth_rxns(self):
        '''
        [:class:`tmo.ParallelReaction`] Biomass growth reactions.
        '''
        return self._growth_rxns

    @property
    def organic_rm(self):
        '''[float] Overall organic (COD) removal rate.'''
        Qi, Qe = self._inf.F_vol, self.outs[1].F_vol
        Si = compute_stream_COD(self._inf, 'kg/m3')
        Se = compute_stream_COD(self.outs[1], 'kg/m3')
        return 1 - Qe*Se/(Qi*Si)
    
#%%
from numba import njit
@njit(cache=True)
def dydt_mbr(QC_ins, QC, V, Qp, f_rtn, xarr, _dstate):
    Q_ins = QC_ins[:, -1]
    C_ins = QC_ins[:, :-1]
    Q = sum(Q_ins)
    C = QC[:-1]
    _dstate[-1] = 0
    _dstate[:-1] = (Q_ins @ C_ins - (Q*(1-xarr) + (Qp+(Q-Qp)*(1-f_rtn))*xarr)*C)/V

class CompletelyMixedMBR(CSTR):
    '''
    Completely mixed membrane bioreactor, equivalent to a CSTR with ideal 
    membrane filtration at the outlet.

    See Also
    --------
    :class:`qsdsan.processes.DiffusedAeration`
    :class:`qsdsan.sanunits.CSTR`

    Examples
    --------
    >>> from qsdsan import WasteStream, processes as pc, sanunits as su
    >>> cmps = pc.create_asm1_cmps()
    >>> ww = WasteStream('ww')
    >>> ww.set_flow_by_concentration(
    ...     flow_tot=2000, 
    ...     concentrations=dict(S_I=21.6, S_S=86.4, X_I=32.4, X_S=129.6, 
    ...                         S_NH=25, S_ND=2.78, X_ND=6.28, S_ALK=84), 
    ...     units=('m3/d', 'mg/L'))
    >>> asm = pc.ASM1()
    >>> M1 = su.CompletelyMixedMBR('M1', ins=ww, outs=['filtrate', 'pumped'],                                   
    ...                            V_max=400, pumped_flow=50, solids_capture_rate=0.999, 
    ...                            aeration=2.0, DO_ID='S_O', 
    ...                            suspended_growth_model=asm)
    >>> M1.set_init_conc(X_I=1000, S_I=30, S_S=5.0, X_S=100, X_BH=500, 
    ...                  X_BA=100, X_P=100, S_O=2, S_NH=2, S_ND=1, 
    ...                  X_ND=1, S_NO=20, S_ALK=84)
    >>> M1.simulate(t_span=(0,400), method='BDF')
    >>> M1.show()
    CompletelyMixedMBR: M1
    ins...
    [0] ww
    phase: 'l', T: 298.15 K, P: 101325 Pa
    flow (g/hr): S_I    1.8e+03
                    S_S    7.2e+03
                    X_I    2.7e+03
                    X_S    1.08e+04
                    S_NH   2.08e+03
                    S_ND   232
                    X_ND   523
                    S_ALK  7e+03
                    H2O    8.31e+07
        WasteStream-specific properties:
         pH         : 7.0
         Alkalinity : 2.5 mmol/L
         COD        : 270.0 mg/L
         BOD        : 137.1 mg/L
         TC         : 170.4 mg/L
         TOC        : 86.4 mg/L
         TN         : 36.0 mg/L
         TP         : 2.7 mg/L
         TSS        : 121.5 mg/L
    outs...
    [0] filtrate
    phase: 'l', T: 298.15 K, P: 101325 Pa
    flow (g/hr): S_I    1.75e+03
                    S_S    277
                    X_I    101
                    X_S    4.83
                    X_BH   236
                    X_BA   13.7
                    X_P    44.1
                    S_O    162
                    S_NO   1.8e+03
                    S_NH   61.7
                    S_ND   60
                    X_ND   0.323
                    S_ALK  3.6e+03
                    S_N2   254
                    H2O    8.1e+07
        WasteStream-specific properties:
         pH         : 7.0
         COD        : 29.9 mg/L
         BOD        : 4.2 mg/L
         TC         : 54.0 mg/L
         TOC        : 9.7 mg/L
         TN         : 24.0 mg/L
         TP         : 0.3 mg/L
         TK         : 0.0 mg/L
         TSS        : 3.7 mg/L
    [1] pumped
    phase: 'l', T: 298.15 K, P: 101325 Pa
    flow (g/hr): S_I    45
                    S_S    7.09
                    X_I    2.6e+03
                    X_S    124
                    X_BH   6.06e+03
                    X_BA   351
                    X_P    1.13e+03
                    S_O    4.17
                    S_NO   46
                    S_NH   1.58
                    S_ND   1.54
                    X_ND   8.29
                    S_ALK  92.2
                    S_N2   6.51
                    H2O    2.07e+06
        WasteStream-specific properties:
         pH         : 7.0
         COD        : 4950.0 mg/L
         BOD        : 1779.8 mg/L
         TC         : 1795.1 mg/L
         TOC        : 1750.8 mg/L
         TN         : 381.0 mg/L
         TP         : 77.2 mg/L
         TK         : 17.3 mg/L
         TSS        : 3693.8 mg/L
    
    >>> flt, rtn = M1.outs
    >>> flt.get_TSS() / rtn.get_TSS() # doctest: +ELLIPSIS
    0.001...
    
    '''
    _N_ins = 1
    _N_outs = 2  # [0] filtrate, [1] pumped flow
    _outs_size_is_fixed = True
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', isdynamic=True, 
                 pumped_flow=50, solids_capture_rate=0.999, 
                 V_max=1000, crossflow_air=None,
                 **kwargs):
        super().__init__(ID, ins, outs, split=None, thermo=thermo,
                         init_with=init_with, V_max=V_max, isdynamic=isdynamic, 
                         **kwargs)
        self.pumped_flow = pumped_flow
        self.solids_capture_rate = solids_capture_rate
        self.crossflow_air = crossflow_air
    
    @property
    def pumped_flow(self):
        '''[float] Pumped flow rate, in m3/d'''
        return self._Q_pump
    @pumped_flow.setter
    def pumped_flow(self, Q):
        self._Q_pump = Q
    
    @property
    def solids_capture_rate(self):
        '''[float] Membrane solid capture rate, i.e., 
        filtrate-to-internal solids concentration ratio, unitless.'''
        return self._f_rtn
    @solids_capture_rate.setter
    def solids_capture_rate(self, f):
        if f < 0 or f > 1:
            raise ValueError(f'membrane solids capture rate must be within [0,1], not {f}')
        self._f_rtn = f
        cmps = self._mixed.components
        self._flt2in_conc_ratio = (1-cmps.x) + cmps.x * (1-f)
    
    @property
    def crossflow_air(self):
        '''[:class:`qsdsan.Process` or NoneType]
        Membrane cross flow air specified for process modeling, such as `qsdsan.processes.DiffusedAeration`. 
        Ignored if DO setpoint is specified by the `aeration` attribute.
        '''
        return self._cfa
    @crossflow_air.setter
    def crossflow_air(self, cfa):
        if cfa is None or isinstance(cfa, Process):
            self._cfa = cfa
        else:
            raise TypeError('crossflow_air must be a `Process` object or None, '
                            f'not {type(cfa)}')

    split = None
        
    def _run(self):
        '''Only to converge volumetric flows.'''
        mixed = self._mixed
        mixed.mix_from(self.ins)
        cmps = mixed.components
        Q = mixed.F_vol*24 # m3/d
        Qp = self._Q_pump
        f_rtn = self._f_rtn
        xsplit = Qp / ((1-f_rtn)*(Q-Qp) + Qp) # mass split of solids to pumped flow
        qsplit = Qp / Q
        flt, rtn = self.outs
        mixed.split_to(rtn, flt, xsplit*cmps.x + qsplit*(1-cmps.x))
    
    def _compile_ODE(self):
        aer = self._aeration
        cfa = self._cfa
        isa = isinstance
        cmps = self.components
        if self._model is None:
            warn(f'{self.ID} was initialized without a suspended growth model, '
                 f'and thus run as a non-reactive unit')
            r = lambda state_arr: np.zeros(cmps.size)
        else:
            r = self._model.production_rates_eval

        _dstate = self._dstate
        _update_dstate = self._update_dstate
        V = self._V_max
        Qp = self.pumped_flow
        f_rtn = self.solids_capture_rate
        xarr = cmps.x
        gstrip = self.gas_stripping
        if gstrip:
            gas_idx = self.components.indices(self.gas_IDs)
            if isa(aer, Process): kLa = aer.kLa
            else: kLa = 0.
            if cfa: kLa += cfa.kLa
            S_gas_air = np.asarray(self.K_Henry)*np.asarray(self.p_gas_atm)
            kLa_stripping = np.maximum(kLa*self.D_gas/self._D_O2, self.stripping_kLa_min)
        hasexo = bool(len(self._exovars))
        f_exovars = self.eval_exo_dynamic_vars
         
        if isa(aer, (float, int)):
            i = cmps.index(self._DO_ID)
            def dy_dt(t, QC_ins, QC, dQC_ins):
                QC[i] = aer
                dydt_mbr(QC_ins, QC, V, Qp, f_rtn, xarr, _dstate)
                if hasexo: QC = np.append(QC, f_exovars(t))
                _dstate[:-1] += r(QC)
                if gstrip: _dstate[gas_idx] -= kLa_stripping * (QC[gas_idx] - S_gas_air)
                _dstate[i] = 0
                _update_dstate()
        else:        
            if cfa:
                cfa_stoi = cfa._stoichiometry
                cfa_frho = cfa.rate_function
                dy_cfa = lambda QC: cfa_stoi * cfa_frho(QC)
            else:
                dy_cfa = lambda QC: 0.
            
            if isa(aer, Process):
                aer_stoi = aer._stoichiometry
                aer_frho = aer.rate_function
                dy_aer = lambda QC: aer_stoi * aer_frho(QC)
            else:
                dy_aer = lambda QC: 0.
                
            def dy_dt(t, QC_ins, QC, dQC_ins):
                dydt_mbr(QC_ins, QC, V, Qp, f_rtn, xarr, _dstate)
                if hasexo: QC = np.append(QC, f_exovars(t))
                _dstate[:-1] += r(QC) + dy_aer(QC) + dy_cfa(QC)
                if gstrip: _dstate[gas_idx] -= kLa_stripping * (QC[gas_idx] - S_gas_air)
                _update_dstate()
        self._ODE = dy_dt
    
    def _update_state(self):
        arr = self._state
        arr[arr < 1e-16] = 0.
        arr[-1] = sum(ws.state[-1] for ws in self.ins)
        for ws in self.outs:
            if ws.state is None: 
                ws.state = np.zeros_like(arr)
                ws.dstate = np.zeros_like(arr)
        flt, rtn = self.outs
        Qp = self.pumped_flow
        flt.state[:-1] = arr[:-1] * self._flt2in_conc_ratio
        flt.state[-1] = arr[-1] - Qp
        rtn.state[:-1] = arr[:-1]
        rtn.state[-1] = Qp
        
    def _update_dstate(self):
        arr = self._dstate
        arr[-1] = sum(ws.dstate[-1] for ws in self.ins)
        flt, rtn = self.outs
        flt.dstate[:-1] = arr[:-1] * self._flt2in_conc_ratio
        flt.dstate[-1] = arr[-1]
        rtn.dstate[:-1] = arr[:-1]
        rtn.dstate[-1] = 0


#%%

class PlugFlowMBR(PFR):
    pass