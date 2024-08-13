#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    
    Yalin Li <mailto.yalin.li@gmail.com>

    Joy Zhang <joycheung1994@gmail.com>

    Jianan Feng <jiananf2@illinois.edu>

Part of this module is based on the biosteam package:
https://github.com/BioSTEAMDevelopmentGroup/biosteam

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

from warnings import warn
from math import ceil, pi
from biosteam import HeatExchangerNetwork as HXN, HXprocess as HXP, HXutility as HXU
from biosteam.units.design_tools.specification_factors import material_densities_lb_per_ft3
from biosteam.exceptions import bounds_warning, DesignWarning
from biosteam.units.design_tools import flash_vessel_design
from .. import SanUnit, Construction
from ..utils import auom

__all__ = ('HeatExchangerNetwork', 'HXprocess', 'HXutility',)

_in_to_ft = auom('in').conversion_factor('ft')
_lb_to_kg = auom('lb').conversion_factor('kg')
_Pa_to_psi = auom('Pa').conversion_factor('psi')

        
class HeatExchangerNetwork(SanUnit, HXN):
    '''
    Similar to the :class:`biosteam.units.facilities.HeatExchangerNetwork`,
    but also a subclass of :class:`qsdsan.SanUnit`

    Examples
    --------
    >>> import qsdsan as qs
    >>> from qsdsan import Stream, sanunits as su
    >>> from qsdsan.utils import create_example_components
    >>> qs.set_thermo(create_example_components())
    >>> salt_water = Stream('salt_water', Water=2000, NaCl=50, units='kg/hr')
    >>> methanol = Stream('methanol', Methanol=20, units='kg/hr')
    >>> ethanol = Stream('ethanol', Ethanol=10, units='kg/hr')
    >>> M1 = su.MixTank('M1', ins=(salt_water, 'recycled_brine', methanol, ethanol), init_with='Stream')
    >>> P1 = su.Pump('P1', ins=M1-0, init_with='Stream')
    >>> H1 = su.HXutility('H1', ins=P1-0, T=350, init_with='Stream')
    >>> S1 = su.ComponentSplitter('S1', ins=H1-0, split_keys=('Methanol', 'Ethanol'), init_with='Stream')
    >>> M2 = su.Mixer('M2', ins=(S1-0, S1-1), outs='alcohols', init_with='Stream')
    >>> S2 = su.Splitter('S2', ins=S1-2, outs=(1-M1, 'waste_brine'), split=0.2, init_with='Stream')
    >>> H2 = su.HXutility('H2', ins=S2.outs[1], T=280, init_with='Stream')
    >>> HXN = su.HeatExchangerNetwork('HXN')
    >>> sys = qs.System('sys', path=(M1, P1, H1, S1, M2, S2, H2), facilities=(HXN,))
    >>> sys.simulate()
    >>> # The actual utility usage is just 30% of the original one (i.e., without HXN)
    >>> round(HXN.actual_heat_util_load/HXN.original_heat_util_load, 2)
    0.28
    >>> HXN.stream_life_cycles # doctest: +SKIP
    [<StreamLifeCycle: Stream_0, cold
     	life_cycle = [
     		<LifeStage: <HXprocess: HX_0_1_hs>, H_in = 0 kJ, H_out = 3.15e+05 kJ>
     		<LifeStage: <HXutility: Util_0_hs>, H_in = 3.15e+05 kJ, H_out = 4.4e+05 kJ>
     	]>,
     <StreamLifeCycle: Stream_1, hot
     	life_cycle = [
     		<LifeStage: <HXprocess: HX_0_1_hs>, H_in = 3.49e+05 kJ, H_out = 3.36e+04 kJ>
     		<LifeStage: <HXutility: Util_1_cs>, H_in = 3.36e+04 kJ, H_out = -1.22e+05 kJ>
     	]>]

    See Also
    --------
    `biosteam.units.facilities.HeatExchangerNetwork <https://biosteam.readthedocs.io/en/latest/API/units/facilities/HeatExchangerNetwork.html>`_
    '''
    ticket_name = HXN.ticket_name
    acceptable_energy_balance_error = HXN.acceptable_energy_balance_error
    raise_energy_balance_error = HXN.raise_energy_balance_error
    network_priority = HXN.network_priority
    _N_ins = HXN._N_ins
    _N_outs = HXN._N_outs
    _units= HXN._units
    __init__ = HXN.__init__
    _init_ins = HXN._init_ins
    _init_outs = HXN._init_outs


class HXprocess(SanUnit, HXP):
    '''
    Similar to :class:`biosteam.units.HXprocess`,
    but can be initialized with :class:`qsdsan.SanStream` and :class:`qsdsan.WasteStream`.

    See Also
    --------
    `biosteam.units.HXprocess <https://biosteam.readthedocs.io/en/latest/API/units/heat_exchange.html>`_
    '''

    line = HXP.line
    _graphics = HXP._graphics
    _N_ins = HXP._N_ins
    _N_outs = HXP._N_outs

    def __init__(
            self, ID='', ins=None, outs=(), thermo=None,
            init_with='Stream', F_BM_default=None, *,
            U=None, dT=5., T_lim0=None, T_lim1=None,
            material="Carbon steel/carbon steel",
            heat_exchanger_type="Floating head",
            N_shells=2, ft=None, 
            phase0=None,
            phase1=None,
            H_lim0=None,
            H_lim1=None,
            inner_fluid_pressure_drop=None,
            outer_fluid_pressure_drop=None,
            neglect_pressure_drop=True,
        ):
        SanUnit.__init__(self, ID, ins, outs, thermo,
                         init_with=init_with, F_BM_default=F_BM_default)

        #: [float] Enforced overall heat transfer coefficient (kW/m^2/K)
        self.U = U

        #: [float] Total heat transferred in kW (not including losses).
        self.total_heat_transfer = None

        #: Number of shells for LMTD correction factor method.
        self.N_shells = N_shells

        #: User imposed correction factor.
        self.ft = ft

        #: [float] Pinch temperature difference.
        self.dT = dT

        #: [float] Temperature limit of outlet stream at index 0.
        self.T_lim0 = T_lim0

        #: [float] Temperature limit of outlet stream at index 1.
        self.T_lim1 = T_lim1

        #: [float] Temperature limit of outlet stream at index 0.
        self.H_lim0 = H_lim0

        #: [float] Temperature limit of outlet stream at index 1.
        self.H_lim1 = H_lim1

        #: Enforced phase of outlet at index 0
        self.phase0 = phase0

        #: Enforced phase of outlet at index 1
        self.phase1 = phase1

        self.material = material
        self.heat_exchanger_type = heat_exchanger_type
        self.reset_streams_at_setup = False
        
        #: Optional[float] Pressure drop along the inner fluid.
        self.inner_fluid_pressure_drop = inner_fluid_pressure_drop
        
        #: Optional[float] Pressure drop along the outer fluid.
        self.outer_fluid_pressure_drop = outer_fluid_pressure_drop
        
        #: [bool] Whether to assume a negligible pressure drop.
        self.neglect_pressure_drop = neglect_pressure_drop


class HXutility(SanUnit, HXU):
    '''
    Similar to :class:`biosteam.units.HXutility`,
    but can be initialized with :class:`qsdsan.SanStream` and :class:`qsdsan.WasteStream`
    and calculate material usage based on [1]_.
    
    References
    ----------
    .. [1] Seider, W. D., Lewin, D. R., Seader, J. D., Widagdo, S., Gani, R., &
           Ng, M. K. (2017). Product and Process Design Principles. Wiley.
           Chapter 12: Heat Exchanger Design.

    See Also
    --------
    `biosteam.units.HXutility <https://biosteam.readthedocs.io/en/latest/API/units/heat_exchange.html>`_
    '''
    
    line = HXU.line
    _graphics = HXU._graphics
    
    _units = {'Area': 'ft^2',
              'Total tube length': 'ft',
              'Inner pipe weight': 'kg',
              'Outer pipe weight': 'kg',
              'Total steel weight': 'kg',
              'Shell length': 'ft',
              'Shell diameter': 'ft',
              'Shell steel weight': 'kg',
              'Tube weight': 'kg'}
    _bounds = {'Vertical vessel weight': (4200, 1e6),
               'Horizontal vessel weight': (1e3, 9.2e5),
               'Horizontal vessel diameter': (3, 21),
               'Vertical vessel length': (12, 40)}
    
    def __init__(
            self, ID='', ins=None, outs=(), thermo=None,
            init_with='Stream', F_BM_default=None,
            include_construction=True,
            T=None, V=None, rigorous=False, U=None, H=None,
            heat_exchanger_type="Floating head",
            material="Carbon steel/carbon steel",
            N_shells=2,
            ft=None,
            heat_only=None,
            cool_only=None,
            heat_transfer_efficiency=None,
            inner_fluid_pressure_drop=None,
            outer_fluid_pressure_drop=None,
            neglect_pressure_drop=True,
            ):
            SanUnit.__init__(self, ID, ins, outs, thermo,
                             init_with=init_with, F_BM_default=F_BM_default,
                             include_construction=include_construction,)
            self.T = T #: [float] Temperature of outlet stream (K).
            self.V = V #: [float] Vapor fraction of outlet stream.
            self.H = H #: [float] Enthalpy of outlet stream.

            #: [bool] If true, calculate vapor liquid equilibrium
            self.rigorous = rigorous

            #: [float] Enforced overall heat transfer coefficient (kW/m^2/K)
            self.U = U

            #: Number of shells for LMTD correction factor method.
            self.N_shells = N_shells

            #: User imposed correction factor.
            self.ft = ft

            #: [bool] If True, heat exchanger can only heat.
            self.heat_only = heat_only

            #: [bool] If True, heat exchanger can only cool.
            self.cool_only = cool_only

            self.material = material
            self.heat_exchanger_type = heat_exchanger_type

            #: Optional[float] Pressure drop along the inner fluid.
            self.inner_fluid_pressure_drop = inner_fluid_pressure_drop
            
            #: Optional[float] Pressure drop along the outer fluid.
            self.outer_fluid_pressure_drop = outer_fluid_pressure_drop
            
            #: [bool] Whether to assume a negligible pressure drop.
            self.neglect_pressure_drop = neglect_pressure_drop

            #: [bool] User enforced heat transfer efficiency. A value less than 1
            #: means that a fraction of heat transferred is lost to the environment.
            #: If value is None, it defaults to the heat transfer efficiency of the 
            #: heat utility.
            self.heat_transfer_efficiency = heat_transfer_efficiency

    def _design(self, duty=None):
        HXU._design(self)
        
        D = self.design_results
        if not D.get('Area', 0): total_steel = 0. # no HX needs
        elif D['Area'] < 150: # double pipe
            # Assume use 1 1/4 nominal size of inner tube, based on [1] page 365
            # Table 12.3, when use Schedule 40, surface area per foot is 0.435 ft2
            # and weight is 2.28 lb steel per foot
            D['Total tube length'] = D['Area']/0.435
            D['Inner pipe weight'] = D['Total tube length']*2.28*_lb_to_kg
            # Assume use 2 nominal size of outer tube, same length as inner tube
            # the weight is 3.66 lb steel per foot
            D['Outer pipe weight'] = D['Total tube length']*3.66*_lb_to_kg
            total_steel = D['Total steel weight'] = D['Inner pipe weight'] + D['Outer pipe weight']
        else: # shell and tube
            # assume all tubes are 16 ft long and 3/4 in O.D., 16 BWG
            D['Shell length'] = 16
            # D['Total tube length'] = 'N/A'
            single_shell_area = D['Area']/self.N_shells
            if single_shell_area <= 100:
                D['Shell diameter'] = 1
            elif single_shell_area <= 400:
                D['Shell diameter'] = 2
            elif single_shell_area < 1100:
                D['Shell diameter'] = 3
            else:
                D['Shell diameter'] = 3
                self.N_shells = ceil(D['Area']/1100) # ceil shell number
                # if increase the number of N_shells, ft (correction factor) will increase (max = 1),
                # then the required area will decrease, so the calculation here is conservative.
                single_shell_area = D['Area']/self.N_shells

            Shell_design = self._horizontal_vessel_design(self.ins[0].P*_Pa_to_psi, D['Shell diameter'], D['Shell length'])
            D['Shell steel weight'] = Shell_design['Weight']*self.N_shells*_lb_to_kg
            
            single_tube_area = pi*(3/4)*_in_to_ft*D['Shell length']
            
            D['Total tube length'] = D['Shell length']*single_shell_area/single_tube_area*self.N_shells
            
            # according to [1] page 367, Table 12.4, 3/4 in O.D., 16 BWG tube: 0.520 lb/ft
            D['Tube weight'] = D['Total tube length']*0.520*_lb_to_kg
            
            total_steel = D['Total steel weight'] = D['Shell steel weight'] + D['Tube weight']

        if self.include_construction:
            construction = getattr(self, 'construction', [])
            if construction: construction[0].quantity = total_steel
            else:
                self.construction = [
                    Construction('stainless_steel', linked_unit=self, item='Stainless_steel', 
                                 quantity=total_steel, quantity_unit='kg'),
                    ]
            
    def _horizontal_vessel_design(self, pressure, diameter, length) -> dict:
        pressure = pressure
        diameter = diameter
        length = length
        # Calculate vessel weight and wall thickness
        rho_M = material_densities_lb_per_ft3['Carbon steel']
        if pressure < 14.68:
            warn('vacuum pressure vessel ASME codes not implemented yet; '
                 'wall thickness may be inaccurate and stiffening rings may be '
                 'required', category=DesignWarning)
        VW, VWT = flash_vessel_design.compute_vessel_weight_and_wall_thickness(
            pressure, diameter, length, rho_M)
        bounds_warning(self, 'Horizontal vessel weight', VW, 'lb',
                       self._bounds['Horizontal vessel weight'], 'cost')
        bounds_warning(self, 'Horizontal vessel diameter', diameter, 'ft',
                       self._bounds['Horizontal vessel diameter'], 'cost')
        Design = {}
        Design['Vessel type'] = 'Horizontal'
        Design['Length'] = length # ft
        Design['Diameter'] = diameter # ft
        Design['Weight'] = VW # lb
        Design['Wall thickness'] = VWT # in
        return Design