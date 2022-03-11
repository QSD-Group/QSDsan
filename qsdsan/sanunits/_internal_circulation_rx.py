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

'''
TODO:
    - Set a maximum on Xw (or Xe?), then increase the Qw if Xw is too high
    - C/N ratio, George's analysis shows it's very high
        - About 16, looks good
    - Consider sulfate and sulfide
    - Check with Brian's AnMBR paper and see the COD<1300 mg/L not preferable thing
'''

import sympy as sp
from math import pi
from biosteam.exceptions import DesignError
from biosteam.utils import ExponentialFunctor
from biosteam.units.design_tools.tank_design import (
    mix_tank_purchase_cost_algorithms,
    TankPurchaseCostAlgorithm
    )
from . import HXutility, MixTank, Pump
from .. import SanStream, SanUnit
from ..utils import (
    compute_stream_COD,
    get_digestion_rxns,
    )

__all__ = ('InternalCirculationRx',)

degassing = SanStream.degassing


# %%

# Tank cost algorithms
IC_purchase_cost_algorithms = mix_tank_purchase_cost_algorithms.copy()
conventional = IC_purchase_cost_algorithms['Conventional']
#!!! Need to check if the cost correlation still holds for the ranges beyond
ic = TankPurchaseCostAlgorithm(
    ExponentialFunctor(A=conventional.f_Cp.A,
                       n=conventional.f_Cp.n),
    V_min=pi/4*1.5**2*16, # 1.5 and 16 are the lower bounds of the width and height ranges in ref [1]
    V_max=pi/4*12**2*25, # 12 and 25 are the lower bounds of the width and height ranges in ref [1]
    V_units='m^3',
    CE=conventional.CE,
    material='Stainless steel')

IC_purchase_cost_algorithms['IC'] = ic


# %%

class InternalCirculationRx(MixTank):
    '''
    Internal circulation (IC) reactor for anaerobic digestion (AD),
    including a high-rate bottom reactor for rapid organic removal and
    a low-rate top reactor for polishing.
    Both reactors are similar to upflow anaerobic blanket reactor (UASB).

    Design of the reactor follows steps described in [1]_
    (assuming steady state and pseudo-zeroth-order kinetics),
    where two methods are used based on Irizar et al.[2]_ and
    Tchobanoglous et al.[3]_.

    Parameters
    ----------
    method : str
        Either "separate" to design the bottom and top reactors separately as in [2]_ ,
        or "lumped" to design the entire IC reactor as a black box following [3]_ .

        In "separate" method, design parameters include:
            - OLRall, biodegradability, Y, q_Qw, mu_max, b, Fxt, and Fxb

        In "lumped" method, design parameters include:
            - OLRall, biodegradability, Y, q_Qw, and q_Xw
    biomass_ID: str
        ID of the Component that represents the biomass.
    OLRall : float
        Overall organic loading rate, [kg COD/m3/hr].
    biodegradability : float or dict
        Biodegradability of components,
        when shown as a float, all biodegradable components are assumed to have
        the same degradability.
    Y : float
        Biomass yield, [kg biomass/kg consumed COD].
    q_Qw : float
        Ratio between the bottom reactor waste flow and the influent.
    q_Xw : float
        Ratio between the biomass concentration in the reactor and the waste flow.
    mu_max : float
        Maximum specific growth rate, [/hr].
    b : float
        Specific endogenous decay coefficient, [/hr].
    V_wf : float
        Fraction of working volume over total volume.
    vessel_type : str
        Can be "IC" to use the reactor size constraints according to [1]_,
        or "Conventional" based on :class:`biosteam.MixTank`
        (much smaller tank size, not recommended).
    vessel_material : str
        Vessel material.
    kW_per_m3 : float
        Electricity requirement per unit volume, [kW/m^3].
        Default to 0 as IC reactors realizes mixing through internal circulation
        caused by the rising force of the generated biogas.
    T : float
        Temperature of the reactor.
        Will not control temperature if provided as None.
    kwargs : dict
        Other keyword arguments (e.g., Fxb, Fxt).

    References
    ----------
    .. [1] Kontos, G. A. Advanced Anaerobic Treatment for Energy Recovery and
        Improved Process Economics in the Management of Biorefinery Wastewaters,
        University of Illinois at Urbana-Champaign, Champaign, IL, 2021.
    .. [2] Irizar et al., Model-Based Design of a Software Sensor for Real-Time
        Diagnosis of the Stability Conditions in High-Rate Anaerobic Reactors –
        Full-Scale Application to Internal Circulation Technology.
        Water Research 2018, 143, 479–491.
        `<https://doi.org/10.1016/j.watres.2018.06.055>`_.
    .. [3] Tchobanoglous et al., Wastewater Engineering: Treatment and Resource Recovery,
        5th ed.; McGraw-Hill Education: New York, 2013.
    '''
    _N_ins = 1
    _N_outs = 3 # biogas, effluent, waste sludge

    # Assumptions
    _q_Qw = 0.01
    _q_Xw = 1.5
    _mu_max = 0.01
    _b = 0.00083
    _Fxb = 0.0032
    _Fxt = 0.0281

    # Related to cost algorithm
    _default_vessel_type = 'IC'
    _default_vessel_material = 'Stainless steel'
    purchase_cost_algorithms = IC_purchase_cost_algorithms

    # Other equipment, only capital cost will be automatically accounted for
    auxiliary_unit_names = ('heat_exchanger', 'effluent_pump', 'sludge_pump')


    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                 method='lumped', biomass_ID='WWTsludge',
                 OLRall=1.25, biodegradability=1., Y=0.07,
                 vessel_type='IC', vessel_material='Stainless steel',
                 V_wf=0.8, kW_per_m3=0., T=35+273.15, init_with='WasteStream',
                 **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with=init_with)
        ID = self.ID
        self._inf = self.ins[0].copy(f'{ID}_inf')
        self._gas = self.outs[0].copy(f'{ID}_gas')
        self.method = method
        self.OLRall = OLRall
        self.biodegradability = biodegradability
        self.Y = Y
        self.V_wf = V_wf or self._default_V_wf
        self.vessel_type = 'IC'
        self.vessel_material = vessel_material
        self.kW_per_m3 = kW_per_m3
        self.T = T
        # Initialize the attributes
        self.heat_exchanger = hx = HXutility(None, None, None, T=T)
        self.heat_utilities = hx.heat_utilities
        self._refresh_rxns()
        # Conversion will be adjusted in the _run function
        self._xcmp = xcmp = getattr(self.components, biomass_ID)
        self._decay_rxn = xcmp.get_combustion_reaction(conversion=0.)
        #!!! Double-check if this would interfere
        eff = self._eff = self.outs[1].proxy(f'{ID}_eff')
        sludge = self._sludge = self.outs[2].proxy(f'{ID}_sludge')
        self.effluent_pump = Pump(f'{self.ID}_eff', ins=eff, init_with=init_with)
        self.sludge_pump = Pump(f'{self.ID}_sludge', ins=sludge, init_with=init_with)

        for k, v in kwargs.items():
            setattr(self, k, v)


    def _refresh_rxns(self, X_biogas=None, X_growth=None):
        X_biogas = X_biogas if X_biogas else 1 - self.Y
        X_growth = X_growth if X_growth else self.Y
        xcmp_ID = self._xcmp.ID
        self._biogas_rxns = get_digestion_rxns(self.ins[0], self.biodegradability,
                                               X_biogas, 0., xcmp_ID)
        self._growth_rxns = get_digestion_rxns(self.ins[0], self.biodegradability,
                                               0., X_growth, xcmp_ID)
        self._i_rm = self._biogas_rxns.X + self._growth_rxns.X


    def _run(self):
        inf = self._inf
        biogas, eff, waste  = self.outs
        # Initialize the streams
        inf.copy_like(self.ins[0])
        biogas.phase = 'g'
        biogas.empty()
        gas = self._gas
        gas.copy_like(biogas)

        inf.split_to(waste, eff, self.q_Qw)
        biogas_rxns = self.biogas_rxns
        growth_rxns = self.growth_rxns

        growth_rxns(inf.mol)
        biogas_rxns(inf.mol)

        degassing(inf, gas)
        Se = compute_stream_COD(inf, 'kg/m3')

        Qi, Si, Xi, Qe, Y = self.Qi, self.Si, self.Xi, self.Qe, self.Y
        method = self.method.lower()
        if method == 'separate':
            run_inputs = (Qi, Si, Xi, Qe, Se, self.Vliq, Y,
                          self.mu_max, self.b, self.Fxb, self.Fxt)
            Xw, Xe = self._run_separate(run_inputs)

        else:
            Xe = (Qi*Xi+Qi*(Si-Se)*Y)/(Qe+(Qi-Qe)*self.q_Xw)
            Xw = Xe * self.q_Xw
            for rxns in (growth_rxns, biogas_rxns):
                rxns(waste.mol)
                rxns(eff.mol)

        degassing(eff, biogas)
        degassing(waste, biogas)

        eff.imass['WWTsludge'] = Xe*self.Qe
        waste.imass['WWTsludge'] = Xw*self.Qw

        diff = sum(i.imol['WWTsludge'] for i in self.outs) - inf.imol['WWTsludge']
        if diff >0:
            decay_rxn = self.decay_rxn
            decay_rxn._X = diff / inf.imol['WWTsludge']

            for i in (eff, waste):
                decay_rxn.force_reaction(i.mol)
                i.imol['O2'] = max(0, i.imol['O2'])
                degassing(i, biogas)

        if self.T:
            biogas.T = eff.T = waste.T = self.T


    def _run_separate(self, run_inputs):
        Qi, Si, Xi, Qe, Se, Vliq, Y, mu_max, b, Fxb, Fxt = run_inputs
        Qw = Qi - Qe
        Xb, Xe, Sb, Vb = sp.symbols('Xb, Xe, Sb, Vb', real=True)

        # Mass balances based on biomass/substrate changes in the bottom/top rx,
        # (0 at steady state)
        biomass_b = Qi*Xi - (Qe*Xb*Fxb+Qw*Xb) + Xb*Vb*(mu_max-b)
        biomass_t = Qe*(Fxb*Xb-Fxt*Xe) + Xe*(Vliq-Vb)*(mu_max-b)
        substrate_b = Qi*(Si-Sb) - mu_max*(Xb*Vb/Y)
        substrate_t = Qe*(Sb-Se) - mu_max*((Vliq-Vb)*Xe/Y)

        parameters = (Qi, Qe, Si, Se, Vliq)

        results = sp.solve(
            (sp.Eq(biomass_b, 0),
             sp.Eq(biomass_t, 0),
             sp.Eq(substrate_b, 0),
             sp.Eq(substrate_t, 0)), (Xb, Xe, Sb, Vb))

        Xb, Xe, Sb, Vb = self._filter_results('separate', parameters, results)

        Vt = Vliq - Vb # volume of the top rx, m3
        self._Vb, self._Vt = Vb, Vt
        return Xb, Xe


    @staticmethod
    def _filter_results(method, parameters, results):
        '''Check if the solution satisfies the design constraints.'''
        Qi, Qe, Si, Se, Vliq = parameters
        solutions = []
        for result in results:
            Xb, Xe, Sb, Vb = result
            Vt = Vliq - Vb
            OLRt = Qe*Sb / Vt
            OLRb = Qi*Si / Vb

            if (
                    0 <= OLRt <= OLRb and
                    0 <= Se <= Sb <= Si and
                    0 <= Xe <= Xb and
                    0 <= Vb <= Vliq
                ):
                solutions.append(result)

        if len(solutions) == 0 :
            raise DesignError('No feasible design found for the given parameters.')

        elif len(solutions) >1: # find more than one solution
            Xbs = [i[1] for i in solutions]
            index = Xbs.index(min(Xbs)) # choose the one with lowest effluent biomass
            return solutions[index]


    _units = {
        'HRT': 'hr',
        'SRT': 'hr',
        'Single reactor liquid volume': 'm3',
        'Bottom reactor volume': 'm3',
        'Top reactor volume': 'm3',
        'Gas chamber volume': 'm3'
        }
    def _design(self):
        D = self.design_results
        D['HRT'] = D['Residence time'] = self.HRT
        D['SRT'] = self.SRT
        D['Total volume'] = self.Vtot
        D['Total liquid volume'] = self.Vliq
        if self.method == 'separate':
            D['Bottom reactor volume'] = self.Vb
            D['Top reactor volume'] = self.Vt


    def _cost(self):
        MixTank._cost(self)

        power_utility = self.power_utility
        pumps = (self.effluent_pump, self.sludge_pump)
        for p in pumps:
            p.simulate()
            power_utility.rate += p.power_utility.rate

        inf, T = self.ins[0], self.T
        if T:
            H_at_T = inf.thermo.mixture.H(mol=inf.mol, phase='l', T=T, P=101325)
            duty = -(inf.H - H_at_T) if self.T else 0.
        else:
            duty = 0.
        self.heat_exchanger.simulate_as_auxiliary_exchanger(duty, inf)


    @property
    def method(self):
        '''[str] Design method, can be "separate" or "lumped".'''
        return self._method
    @method.setter
    def method(self, i):
        if not i.lower() in ('separate', 'lumped'):
            raise ValueError('`method` can only be "separated", or "lumped", '
                             f'not "{i}".')
        self._method = i.lower()

    @property
    def biomass_ID(self):
        '''[str] ID of the Component that represents the biomass.'''
        return self._xcmp.ID
    @biomass_ID.setter
    def biomass_ID(self, i):
        self._xcmp = getattr(self.components, i)

    @property
    def OLRall(self):
        '''[float] Overall organic loading rate, [kg COD/m3/hr].'''
        return self._OLRall
    @OLRall.setter
    def OLRall(self, i):
        if i < 0:
            raise ValueError('`OLRall` should be >=0, '
                             f'the input value {i} is outside the range.')
        self._OLRall = i

    @property
    def biodegradability(self):
        '''
        [float of dict] Biodegradability of components,
        when shown as a float, all biodegradable components are assumed to have
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
    def i_rm(self):
        '''[:class:`np.array`] Removal of each component in this reactor.'''
        return self._i_rm

    @property
    def Y(self):
        '''[float] Biomass yield, [kg biomass/kg consumed COD].'''
        return self._Y
    @Y.setter
    def Y(self, i):
        if not 0 <= i <= 1:
            raise ValueError('`Y` should be within [0, 1], '
                             f'the input value {i} is outside the range.')
        self._Y = i

    @property
    def q_Qw(self):
        '''[float] Ratio between the bottom reactor waste flow and the influent.'''
        return self._q_Qw
    @q_Qw.setter
    def q_Qw(self, i):
        if not 0<=i<=1:
            raise ValueError('`q_Qw` should be within [0, 1], '
                             f'the input value {i} is outside the range.')
        self._q_Qw = i

    @property
    def q_Xw(self):
        '''
        [float] Ratio between the biomass concentration in the reactor and the waste flow,
        only relevant when the "lumped" method is used.
        '''
        return self._q_Xw if self.method=='lumped' else None
    @q_Xw.setter
    def q_Xw(self, i):
        if not i>=1:
            raise ValueError('`q_Xw` should be >=1, '
                             f'the input value {i} is outside the range.')
        self._q_Xw = i

    @property
    def mu_max(self):
        '''
        [float] Maximum specific growth rate, [/hr],
        only relevant when the "separate" method is used.
        '''
        return self._mu_max if self.method=='separate' else None
    @mu_max.setter
    def mu_max(self, i):
        if i < 0:
            raise ValueError('`mu_max` should be >= 0, '
                             f'the input value {i} is outside the range.')
        self._mu_max = i

    @property
    def b(self):
        '''
        [float] Specific endogenous decay coefficient, [/hr],
        only relevant when the "separate" method is used.
        '''
        return self._b if self.method=='separate' else None
    @b.setter
    def b(self, i):
        if i < 0:
            raise ValueError('`b` should be >= 0, '
                             f'the input value {i} is outside the range.')
        self._b = i

    @property
    def Fxb(self):
        '''
        [float] Biomass transfer ratio from the bottom reacor to the top reactor,
        should be within [0, 1] (ideal to no retention),
        only relevant when the "separate" method is used.
        '''
        return self._Fxb if self.method=='separate' else None
    @Fxb.setter
    def Fxb(self, i):
        if not 0<=i<=1:
            raise ValueError('`Fxb` should be within [0, 1], '
                             f'the input value {i} is outside the range.')
        self._Fxb = i

    @property
    def Fxt(self):
        '''
        [float] Biomass transfer ratio from the top reacor to the effluent,
        should be within [0, 1] (ideal to no retention),
        only relevant when the "separate" method is used.
        '''
        return self._Fxt if self.method=='separate' else None
    @Fxt.setter
    def Fxt(self, i):
        if not 0<=i<=1:
            raise ValueError('`Fxt` should be within [0, 1], '
                             f'the input value {i} is outside the range.')
        self._Fxt = i

    @property
    def Vb(self):
        '''
        [float] Volume of the bottom reactor, [m3],
        only relevant when the "separate" method is used.
        '''
        return self._Vb if self.method=='separate' else None

    @property
    def Vt(self):
        '''
        [float] Volume of the top reactor, [m3],
        only relevant when the "separate" method is used.
        '''
        return self._Vt if self.method=='separate' else None

    @property
    def Vliq(self):
        '''
        [float] Total volume of the liquid, not considering gas headspace
        and `V_wf`, [m3].
        '''
        return self.Qi*self.Si / self.OLRall

    @property
    def Vtot(self):
        '''
        [float] Total volume considering `V_wf`, [m3].
        '''
        return self.Vliq / self.V_wf

    @property
    def HRT(self):
        '''[float] Hydraulic retention time [hr].'''
        return self.Vliq / self.Qi

    @property
    def tau(self):
        '''
        [float] Reactor residence time, [hr]
        (same as the hydraulic retention time, HRT).
        '''
        return self.HRT

    @property
    def SRT(self):
        '''[float] Solid residence time [hr].'''
        if self.method == 'separate':
            return (self.Xb*self.Vb+self.Xe*self.Vt)/(self.q_Qw*self.Xb+self.Fxt*self.Qe*self.Xe)

        return  self.Xe*self.Vliq / (self.q_Qw*self.q_Xw+self.Qe*self.Xe)

    @property
    def biogas_rxns(self):
        '''
        [:class:`tmo.ParallelReaction`] Biogas production reactions.
        '''
        return self._biogas_rxns

    @property
    def growth_rxns(self):
        '''
        [:class:`tmo.ParallelReaction`] Biomass (WWTsludge) growth reactions.
        '''
        return self._growth_rxns

    @property
    def decay_rxn(self):
        '''
        [:class:`tmo.Reaction`] Biomass endogenous decay.

        .. note::
            Conversion is adjusted in the _run function.
        '''
        return self._decay_rxn

    @property
    def Qi(self):
        '''[float] Influent volumetric flow rate, [m3/hr].'''
        return self.ins[0].F_vol

    @property
    def Qe(self):
        '''[float] Effluent volumetric flow rate, [m3/hr].'''
        return self.outs[1].F_vol

    @property
    def Qw(self):
        '''[float] Waste flow volumetric flow rate, [m3/hr].'''
        return self.outs[2].F_vol

    @property
    def Si(self):
        '''
        [float] Influent substrate (i.e., biodegradable components)
        concentration, [kg/m3].
        '''
        return compute_stream_COD(self.ins[0], 'kg/m3')

    @property
    def Se(self):
        '''
        [float] Effluent substrate (i.e., biodegradable components)
        concentration, [kg/m3].
        '''
        return compute_stream_COD(self.outs[1], 'kg/m3')

    @property
    def Sw(self):
        '''
        [float] Waste flow substrate (i.e., biodegradable components)
        concentration, [kg/m3].
        '''
        return compute_stream_COD(self.outs[2], 'kg/m3')

    @property
    def Xi(self):
        '''[float] Influent biomass (i.e., `WWTsludge`) concentration, [kg/m3].'''
        return self.ins[0].imass['WWTsludge']/self.ins[0].F_vol

    @property
    def Xe(self):
        '''[float] Effluent  biomass (i.e., `WWTsludge`) concentration, [kg/m3].'''
        return self.outs[1].imass['WWTsludge']/self.outs[1].F_vol

    @property
    def Xw(self):
        '''[float] Waste flow biomass (i.e., `WWTsludge`) concentration, [kg/m3].'''
        return self.outs[2].imass['WWTsludge']/self.outs[2].F_vol

    @property
    def organic_rm(self):
        '''[float] Overall organic (COD) removal rate.'''
        return 1 - self.Qe*self.Se/(self.Qi*self.Si)