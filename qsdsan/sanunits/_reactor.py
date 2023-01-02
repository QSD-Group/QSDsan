#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:

    Jianan Feng <jiananf2@illinois.edu>

    Yalin Li <mailto.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

from math import pi, ceil
from biosteam.units.design_tools import PressureVessel
from biosteam.exceptions import DesignError
from qsdsan import SanUnit, Stream, Construction
from qsdsan.utils import auom

__all__ = ('Reactor',)

_lb_to_kg = auom('lb').conversion_factor('kg')
_m_to_ft = auom('m').conversion_factor('ft')
_Pa_to_psi = auom('Pa').conversion_factor('psi')

class Reactor(SanUnit, PressureVessel, isabstract=True):
    '''
    Create an abstract class for reactor unit, purchase cost of the reactor
    is based on volume calculated by residence time.
    
    Parameters
    ----------
    ins : Iterable(stream)
        Inlet.
    outs : Iterable(stream)
        Outlet.
    tau: float
        Residence time, [hr].
    V_wf: float
        Fraction of working volume over total volume.
    length_to_diameter : float
        Reactor length to diameter ratio.
    N: int
        Number of reactor.
    V: float
        Volume of reactor, [m3].
    auxiliary: bool
        Whether or not the reactor is an auxiliary unit.      
    mixing_intensity: float
        Mechanical mixing intensity, [/s].
    kW_per_m3: float
        Power usage of agitator
        (converted from 0.5 hp/1000 gal as in [1]).
        If mixing_intensity is provided, this will be calculated based on
        the mixing_intensity and viscosity of the influent mixture as in [2]_
    wall_thickness_factor=1: float
        A safety factor to scale up the calculated minimum wall thickness.
    vessel_material : str, optional
        Vessel material. Default to 'Stainless steel 316'.
    vessel_type : str, optional
        Vessel type. Can only be 'Horizontal' or 'Vertical'.
        
    References
    ----------
    .. [1] Seider, W. D.; Lewin, D. R.; Seader, J. D.; Widagdo, S.; Gani, R.;
        Ng, M. K. Cost Accounting and Capital Cost Estimation. In Product
        and Process Design Principles; Wiley, 2017; pp 470.
    .. [2] Shoener et al. Energy Positive Domestic Wastewater Treatment:
        The Roles of Anaerobic and Phototrophic Technologies.
        Environ. Sci.: Processes Impacts 2014, 16 (6), 1204–1222.
        https://doi.org/10.1039/C3EM00711A.
    '''
    _N_ins = 2
    _N_outs = 1
    _ins_size_is_fixed = False
    _outs_size_is_fixed = False

    _units = {**PressureVessel._units,
              'Residence time': 'hr',
              'Total volume': 'm3',
              'Single reactor volume': 'm3',
              'Reactor volume': 'm3'}

    # For a single reactor, based on diameter and length from PressureVessel._bounds,
    # converted from ft3 to m3
    _Vmax = pi/4*(20**2)*40/35.3147
    
    _F_BM_default = PressureVessel._F_BM_default

    _vessel_material = 'Stainless steel 316'

    def __init__(self, ID='', ins=None, outs=(), *,
                 P=101325, tau=0.5, V_wf=0.8,
                 length_to_diameter=2, N=None, V=None, auxiliary=False,
                 mixing_intensity=None, kW_per_m3=0.0985,
                 wall_thickness_factor=1,
                 vessel_material='Stainless steel 316',
                 vessel_type='Vertical'):
        SanUnit.__init__(self, ID, ins, outs)
        self.P = P
        self.tau = tau
        self.V_wf = V_wf
        self.length_to_diameter = length_to_diameter
        self.N = N
        self.V = V
        self.auxiliary = auxiliary
        self.mixing_intensity = mixing_intensity
        self._mixture = Stream(f'{self.ID}_mixture')
        self.kW_per_m3 = kW_per_m3
        self.wall_thickness_factor = wall_thickness_factor
        self.vessel_material = vessel_material
        self.vessel_type = vessel_type

        
    def _init_lca(self):
        for i in self.construction: i.registry.discard(i)
        item_name = self._vessel_material.replace(' ', '_').rstrip('_316').rstrip('_304')
        self.construction = [
            Construction(item_name.lower(), linked_unit=self, item=item_name, quantity_unit='kg'),
            ]
        
    def _design(self):
        Design = self.design_results
        if not self.auxiliary:
        # if auxiliary, in our system, only K/O drum whose N, and V are provided
        # do not need to deal with self.F_vol_in (auxiliary unit has trouble doing this)
            ins_F_vol = self.F_vol_in
            for i in range(len(self.ins)):
                ins_F_vol -= (self.ins[i].ivol['H2'] +\
                              self.ins[i].ivol['CHG_catalyst'] +\
                              self.ins[i].ivol['HT_catalyst'] +\
                              self.ins[i].ivol['HC_catalyst'])
            # not include gas (e.g. H2)
            V_total = ins_F_vol * self.tau / self.V_wf
        P = self.P * _Pa_to_psi # Pa to psi
        length_to_diameter = self.length_to_diameter
        wall_thickness_factor = self.wall_thickness_factor

        if self.N:
            if self.V:
                D = (4*self.V/pi/length_to_diameter)**(1/3)
                D *= _m_to_ft # convert from m to ft
                L = D * length_to_diameter

                Design['Total volume'] = self.V*self.N
                Design['Single reactor volume'] = self.V
                Design['Number of reactors'] = self.N
            else:
                V_reactor = V_total/self.N
                D = (4*V_reactor/pi/length_to_diameter)**(1/3)
                D *= _m_to_ft # convert from m to ft
                L = D * length_to_diameter

                Design['Residence time'] = self.tau
                Design['Total volume'] = V_total
                Design['Single reactor volume'] = V_reactor
                Design['Number of reactors'] = self.N
        else:
            N = ceil(V_total/self._Vmax)
            if N == 0:
                V_reactor = 0
                D = 0
                L = 0
            else:
                V_reactor = V_total / N
                D = (4*V_reactor/pi/length_to_diameter)**(1/3)
                D *= _m_to_ft # convert from m to ft
                L = D * length_to_diameter

            Design['Residence time'] = self.tau
            Design['Total volume'] = V_total
            Design['Single reactor volume'] = V_reactor
            Design['Number of reactors'] = N

        Design.update(self._vessel_design(P, D, L))
        if wall_thickness_factor == 1: pass
        elif wall_thickness_factor < 1:
            raise DesignError('wall_thickness_factor must be larger than 1')
        else:
             Design['Wall thickness'] *= wall_thickness_factor
             # Weight is proportional to wall thickness in PressureVessel design
             Design['Weight'] = round(Design['Weight']*wall_thickness_factor, 2)
             
        self.construction[0].quantity = Design['Weight']*Design['Number of reactors']*_lb_to_kg


    def _cost(self):
        Design = self.design_results
        purchase_costs = self.baseline_purchase_costs

        if Design['Total volume'] == 0:
            purchase_costs.clear()

        else:
            purchase_costs.update(self._vessel_purchase_cost(
                Design['Weight'], Design['Diameter'], Design['Length']))
            for i, j in purchase_costs.items():
                purchase_costs[i] *= Design['Number of reactors']

            self.power_utility(self.kW_per_m3*Design['Total volume'])
        
    @property
    def vessel_material(self):
        return self._vessel_material
    @vessel_material.setter
    def vessel_material(self, i):
        exist_material = getattr(self, '_vessel_material', None)
        PressureVessel.vessel_material.fset(self, i)
        if i and exist_material == i: return # type doesn't change, no need to reload construction items
        self._init_lca()

    @property
    def kW_per_m3(self):
        G = self.mixing_intensity
        if G is None:
            return self._kW_per_m3
        else:
            mixture = self._mixture
            mixture.mix_from(self.ins)
            kW_per_m3 = mixture.mu*(G**2)/1e3
            return kW_per_m3
        
    @kW_per_m3.setter
    def kW_per_m3(self, i):
        if self.mixing_intensity and i is not None:
            raise AttributeError('`mixing_intensity` is provided, kw_per_m3 will be calculated.')
        else:
            self._kW_per_m3 = i