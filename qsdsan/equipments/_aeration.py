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
from .. import Equipment
from ..utils import select_pipe, calculate_pipe_material

__all__ = ('Blower', 'GasPiping',)


# %%

class Blower(Equipment):
    '''
    Design and cost blowers based on `Shoener et al. <https://doi.org/10.1039/C5EE03715H>`_.

    Note that blower building concrete usage and excavation is not included here
    as blowers are often housed together with the reactors.

    Refer to :class:`~.sanunits.AnMBR` or :class:`~.sanunits.ActivatedSludgeProcess`
    for examples.

    Parameters
    ----------
    N_reactor : int
        Number of the reactors where the gas sparging modules will be installed.
    gas_demand_per_reactor : float
        Gas demand per reactor, [cfm] (cubic ft per minute).
    TDH : float
        Total dynamic head for the blower, [psi].
    eff_blower : float
        Efficiency of the blower in fraction (i.e., 0.7 for 70%).
    eff_motor : float
        Efficiency of the motor in fraction (i.e., 0.7 for 70%).
    AFF : float
        Air flow fraction.
        The default value is calculated as STE/6
        (STE stands for standard oxygen transfer efficiency, and default STE is 20).
        If using different STE value, AFF should be 1 if STE/6<1
        and 3.33 if STE/6>1.
    building_unit_cost : float
        Unit cost of the blower building, [USD/ft2].

    References
    ----------
    [1] Shoener, B. D.; Zhong, C.; Greiner, A. D.; Khunjar, W. O.; Hong, P.-Y.; Guest, J. S.
    Design of Anaerobic Membrane Bioreactors for the Valorization
    of Dilute Organic Carbon Waste Streams.
    Energy Environ. Sci. 2016, 9 (3), 1102–1112.
    https://doi.org/10.1039/C5EE03715H.

    See Also
    --------
    :class:`~.sanunits.AnMBR`

    :class:`~.sanunits.ActivatedSludgeProcess`

    '''

    def __init__(self, ID=None, linked_unit=None,
                 F_BM={
                     'Blowers': 2.22,
                     'Blower piping': 1,
                     'Blower building': 1.11,
                     },
                 lifetime=15, lifetime_unit='yr',
                 units={
                     'Total gas flow': 'CFM',
                     'Blower capacity': 'CFM',
                     'Number of blowers': '',
                     'Total blower power': 'kW',
                     },
                 N_reactor=2, gas_demand_per_reactor=1,
                 TDH=6, eff_blower=0.7, eff_motor=0.7, AFF=3.33,
                 building_unit_cost=90,
                 ):
        Equipment.__init__(self=self, ID=ID, linked_unit=linked_unit, F_BM=F_BM,
                           lifetime=lifetime, lifetime_unit=lifetime_unit,
                           units=units)
        self.N_reactor = N_reactor
        self.gas_demand_per_reactor = gas_demand_per_reactor
        self.TDH = TDH
        self.eff_blower = eff_blower
        self.eff_motor = eff_motor
        self.AFF = AFF
        self.building_unit_cost = building_unit_cost


    def _design(self):
        N_reactor, gas_demand_per_reactor  =  \
            self.N_reactor, self.gas_demand_per_reactor
        gas_tot = N_reactor * gas_demand_per_reactor
        TDH, eff_blower, eff_motor = self.TDH, self.eff_blower, self.eff_motor

        # Calculate brake horsepower, 14.7 is atmospheric pressure in psi
        BHP = (gas_tot*0.23)*(((14.7+TDH)/14.7)**0.283-1)/eff_blower
        # 0.746 is horsepower to kW
        D = {'Total blower power': BHP*0.746/eff_motor}

        # Calculate the number of blowers
        TCFM = ceil(gas_demand_per_reactor) # total cubic ft per min
        N = 1
        if TCFM <= 30000:
            CFMB = TCFM / N # cubic ft per min per blower
            while CFMB > 7500:
                N += 1
                CFMB = TCFM / N
        elif 30000 < TCFM <= 72000:
            CFMB = TCFM / N
            while CFMB > 18000:
                N += 1
                CFMB = TCFM / N
        else:
            CFMB = TCFM / N
            while CFMB > 100000:
                N += 1
                CFMB = TCFM / N

        D['Total gas flow'] = TCFM
        D['Blower capacity'] = CFMB
        D['Number of blowers'] = N
        return D


    def _cost(self):
        D = self.design_results
        TCFM, CFMB = D['Total gas flow'], D['Blower capacity']
        N_reactor, AFF = self.N_reactor, self.AFF
        # Air pipes
        # Note that the original codes use CFMD instead of TCFM for air pipes,
        # but based on the coding they are equivalent
        # (i.e., just used an alternative name)
        if TCFM <= 1000:
            piping = 617.2 * AFF * (TCFM**0.2553)
        elif 1000 < TCFM <= 10000:
            piping = 1.43 * AFF * (TCFM**1.1337)
        else:
            piping = 28.59 * AFF * (TCFM**0.8085)
        C = {'Blower air piping': piping}

        # Blowers
        if TCFM <= 30000:
            ratio = 0.7 * (CFMB**0.6169)
            blowers = 58000*ratio / 100
        elif 30000 < TCFM <= 72000:
            ratio = 0.377 * (CFMB**0.5928)
            blowers = 218000*ratio / 100
        else:
            ratio = 0.964 * (CFMB**0.4286)
            blowers  = 480000*ratio / 100
        C['Blowers'] = blowers * N_reactor

        # Blower building
        area = 128 * (TCFM**0.256) # building area, [ft2]
        C['Blower building'] = area * self.building_unit_cost

        return C


# %%

class GasPiping(Equipment):
    '''
    Design and cost reactor gas header pipes and manifold based on
    `Shoener et al. <https://doi.org/10.1039/C5EE03715H>`_.

    The gas pipes will be layed along the length of the reactor with
    manifold along the width of th reactor
    (i.e., gas will be pumped from the manifold to the header then into the reactor).

    Refer to :class:`~.sanunits.ActivatedSludgeProcess` for usage.

    Parameters
    ----------
    N_reactor : int
        Number of the reactor where the gas sparging modules will be installed.
    N_pipe_per_reactor : int
        Number of the pipes per reactor.
    gas_demand_per_reactor : float
        Gas demand per reactor, [cfm] (cubic ft per minute).
    v_header : float
        Velocity of gas in the header pipe
        (layed along the length of the reactor), [ft/s].
    v_manifold : float
        Velocity of gas in the manifold pipe
        (layed along the width of the reactor), [ft/s].
    L_reactor : float
        Length of the reactor, [ft].
    L_extra : float
        Extra length to be included in piping for each of the reactor, [ft].
    W_reactor : float
        Width of the reactor, [ft].
    W_extra : float
        Extra width to be included in piping for each of the reactor, [ft].
    pipe_density : float
        Density of the pipe, [kg/ft3].
    pipe_unit_cost : float
        Unit cost of the pipe, [USD/kg].
    TDH : float
        Total dynamic head for the blower, [psi].
    eff_blower : float
        Efficiency of the blower in fraction (i.e., 0.7 for 70%).
    eff_motor : float
        Efficiency of the motor in fraction (i.e., 0.7 for 70%).

    References
    ----------
    [1] Shoener, B. D.; Zhong, C.; Greiner, A. D.; Khunjar, W. O.; Hong, P.-Y.; Guest, J. S.
    Design of Anaerobic Membrane Bioreactors for the Valorization
    of Dilute Organic Carbon Waste Streams.
    Energy Environ. Sci. 2016, 9 (3), 1102–1112.
    https://doi.org/10.1039/C5EE03715H.

    See Also
    --------
    :class:`~.sanunits.ActivatedSludgeProcess`
    '''

    def __init__(self, ID=None, linked_unit=None,
                 F_BM=1., lifetime=15, lifetime_unit='yr',
                 units={'Gas pipe material': 'kg',},
                 N_reactor=2, N_pipe_per_reactor=1,
                 gas_demand_per_reactor=1, v_header=70, v_manifold=70,
                 L_reactor=12, L_extra=0, W_reactor=21, W_extra=2,
                 pipe_density=227.3, # from 0.29 lb/in3
                 pipe_unit_cost=0,):
        Equipment.__init__(self=self, ID=ID, linked_unit=linked_unit, F_BM=F_BM,
                           lifetime=lifetime, lifetime_unit=lifetime_unit,
                           units=units)
        self.N_reactor = N_reactor
        self.N_pipe_per_reactor = N_pipe_per_reactor
        self.gas_demand_per_reactor = gas_demand_per_reactor
        self.v_header = v_header
        self.v_manifold = v_manifold
        self.L_reactor = L_reactor
        self.L_extra = L_extra
        self.W_reactor = W_reactor
        self.W_extra = W_extra
        self.pipe_density = pipe_density
        self.pipe_unit_cost = pipe_unit_cost


    def _design(self):
        # Gas piping
        N_reactor, N_pipe_per_reactor, gas_demand_per_reactor, pipe_density =  \
            self.N_reactor, self.N_pipe_per_reactor, self.gas_demand_per_reactor, self.pipe_density
        L_reactor, L_extra, W_reactor, W_extra = \
            self.L_reactor, self.L_extra, self.W_reactor, self.W_extra
        # Gas header
        L_gh = L_reactor * N_pipe_per_reactor + L_extra
        gas_demand_per_pipe = gas_demand_per_reactor / N_pipe_per_reactor
        OD_gh, t_gh, ID_gh = select_pipe(gas_demand_per_pipe, self.v_header)
        M_gh = N_reactor * \
            calculate_pipe_material(OD_gh, t_gh, ID_gh, L_gh, pipe_density)
        # Gas supply manifold, used more conservative assumption than the ref
        L_gsm = W_reactor*N_reactor + W_extra
        gas_tot = N_reactor * gas_demand_per_reactor
        OD_gsm, t_gsm, ID_gsm = select_pipe(gas_tot, self.v_manifold)
        M_gsm = calculate_pipe_material(OD_gsm, t_gsm, ID_gsm, L_gsm, pipe_density) # kg
        return {'Gas pipe material': M_gh+M_gsm}


    def _cost(self):
        return self.pipe_unit_cost*self.design_results['Gas pipe material']