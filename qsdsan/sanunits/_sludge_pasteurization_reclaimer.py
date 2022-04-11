#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>
    Shion Watabe <shionwatabe@gmail.com>
    Tori Morgan <tvlmorgan@gmail.com>
    Hannah Lohman <hlohman94@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''


from .. import SanUnit, Construction
from ..utils import ospath, load_data, data_path, price_ratio

__all__ = ('SludgePasteurizationReclaimer',)

sludge_data_path = ospath.join(data_path, 'sanunit_data/_sludge_pasteurization_reclaimer.tsv')


@price_ratio(default_price_ratio=1)
class SludgePasteurizationReclaimer(SanUnit):

    '''
    Sludge Pasteurization of sludge from AnMBR and using biogas from AnMBR,
    and/or LPG. Biogas is combusted when pasteurization occurs.

    The following impact items should be pre-constructed for life cycle assessment:
    Steel.

    Parameters
    ----------
    if_combustion : bool
        If include combusion reaction during simulation.
    temp_pasteurization : float
        Pasteurization temperature is 70C or 343.15 Kelvin
    sludge_temp : float
        Temperature of sludge is 10C or 283.15K
    target_MC : float
        Target moisture content is 10%
    heat_loss : float
        Heat loss during pasteurization process is assumed to be 10%
    lhv_lpg : float
        Lower heating value of Liquid Petroleum Gas at 25C/298.15K is
        46-51 MJ/kg from World Nuclear Organization.
        The lower heating value (also known gross calorific value or gross
        energy) of a fuel is defined as the amount of heat released
        by a specified quantity (initially at 25°C) once it is combusted and
        the products remain evaporated and into atmosphere.
    ppl: int
        Total number of users for scaling of costs.
    if_sludge_service: bool
        If share sludge pasteurization unit among multiple Reclaimer septic tanks (assume 1,000 users per sludge
        pasteurization unit or 10 septic tanks serving a population of 100 users per septic tank)

    ins:
        air: From atmosphere
        sludge: Sludge produced from the AnMBR
        lpg: Purchased LPG for pasteurization

    outs:
        treated sludge: Sludge treated from pasteurization

     References
    ----------
    [1] Shoener et al., Design of Anaerobic Membrane Bioreactors for the
    Valorization of Dilute Organic Carbon Waste Streams.
    Energy Environ. Sci. 2016, 9 (3), 1102–1112.
    https://doi.org/10.1039/C5EE03715H.

    [2] Turek et al., Proposed EU Legislation to Force Changes in Sewage
    Sludge Disposal: A Case Study.
    Front. Chem. Sci. Eng. 2018, 12 (4), 660–669.
    https://doi.org/10.1007/s11705-018-1773-0.

    Examples
    --------
    `bwaise systems <https://github.com/QSD-Group/EXPOsan/blob/main/exposan/bwaise/systems.py>`_
    '''

    # Constants
    # Specific Heat capacity of water
    Cp_w = 4.184  # kJ kg^-1 C^-1
    # Specific Heat capacity of dry matter (sludge)
    Cp_dm = 1.231  # kJ kg^-1 C^-1
    # Specific latent heat of vaporization of water
    l_w = 2260  # kJ kg^-1
    # Baseline population served by a single septic tank
    baseline_ppl = 100
    # Exponential scaling constant for scaling cost and LCA with change in users
    exponent_scale = 0.6  # exponential scaling constant

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream', 
                 heat_loss=0.1, target_MC=0.1, sludge_temp=283.15,
                 temp_pasteurization=343.15, lhv_lpg=48.5, ppl=1, if_sludge_service=True, **kwargs):

        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with, F_BM_default=1)
        self._heat_loss = heat_loss  
        self._target_MC = target_MC      
        self._sludge_temp = sludge_temp
        self._temp_pasteurization = temp_pasteurization
        self.lhv_lpg = lhv_lpg
        self.ppl = ppl
        self.if_sludge_service = if_sludge_service

        if self.ppl < self.baseline_ppl:
            self.user_scale_up = 1  # users are within the capacity of a septic tank (100 users)
        else:
            self.user_scale_up = self.ppl / self.baseline_ppl  # users exceed the capacity of a standard septic tank

        
        data = load_data(path=sludge_data_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data
        
        for attr, value in kwargs.items():
            setattr(self, attr, value)
            
    _N_ins = 3
    _N_outs = 1

    def _run(self):
        air, sludge, lpg = self.ins
        lpg.phase = 'l'
        sludge.phase = 's'
        treated_sludge = self.outs[0] 
        treated_sludge.copy_like(sludge) 

        # Mass calculations
        # total amount of water in sludge
        M_w = sludge.imass['H2O']  # kg/hr
        # total amount of dry matter in sludge
        M_dm = (sludge.F_mass - sludge.imass['H2O'])  # kg/hr
        # total amount of water to be evaporated to reach target dry matter content
        M_we = (sludge.imass['H2O'] - sludge.F_mass * self.target_MC)  # kg/hr
        
        # Overall heat required for pasteurization
        temp_diff = self.temp_pasteurization - self.sludge_temp
        Q_d = (M_w * self.Cp_w + M_dm * self.Cp_dm) * temp_diff + M_we * self.l_w  # kJ/hr
        Q_d *= (1 + self.heat_loss)/1000  # MJ/hr

        lhv_lpg = self.lhv_lpg
        lpg.imass['LPG'] = Q_d / lhv_lpg  # kg-LPG/hr = (MJ/hr)/(MJ/kg-LPG)
           
        treated_sludge.imass['H2O'] = sludge.F_mass * self.target_MC
        treated_sludge.imass['OtherSS'] = sludge.F_mass - M_we
        treated_sludge.imass['P'] = sludge.imass['P']     
        treated_sludge.imass['N'] = sludge.imass['N']
        
    def _design(self):
        design = self.design_results

        if self.if_sludge_service:
            # assume 10 septic tanks sharing 1 sludge pasteurization unit
            design['Steel'] = S_quant = (self.sludge_dryer_weight + self.sludge_barrel_weight) / 10 * self.user_scale_up
        else:
            design['Steel'] = S_quant = (self.sludge_dryer_weight + self.sludge_barrel_weight) * self.user_scale_up

        self.construction = Construction(item='Steel', quantity=S_quant, quantity_unit='kg')
        self.add_construction(add_cost=False)
        
    def _cost(self):
        C = self.baseline_purchase_costs

        if self.if_sludge_service:
            C['Dryer'] = self.sludge_dryer / 10 * (self.user_scale_up ** self.exponent_scale)
            C['Barrel'] = self.sludge_barrel / 10 * (self.user_scale_up ** self.exponent_scale)
        else:
            C['Dryer'] = self.sludge_dryer * (self.user_scale_up ** self.exponent_scale)
            C['Barrel'] = self.sludge_barrel * (self.user_scale_up ** self.exponent_scale)

        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio
        
        self.add_OPEX = self._calc_replacement_cost() + self._calc_labor_cost()  # USD/hr
        self.power_demand = 0
        self.power_utility(self.power_demand)  # kWh

    def _calc_replacement_cost(self):
        sludge_replacement_cost = 0  # Assume dryer and barrel have 25 year lifetime so will not need to be replaced
        return sludge_replacement_cost / (365 * 24) * self.price_ratio  # USD/hr

    def _calc_labor_cost(self):
        labor_cost = (self.wages * self.sludge_labor_maintenance) * (self.user_scale_up ** self.exponent_scale)
        return labor_cost/(365 * 24)  # USD/hr

    @property
    def heat_loss(self):
        '''[float] Pasteurization heat loss as a fraction.'''
        return self._heat_loss
    @heat_loss.setter
    def heat_loss(self, i):
        self._heat_loss = i

    @property
    def sludge_temp(self):
        '''[float] Temperature of sludge, [K].'''
        return self._sludge_temp
    @sludge_temp.setter
    def sludge_temp(self, i):
        self._sludge_temp = i        
        
    @property
    def temp_pasteurization(self):
        '''[float] Temperature of sludge pasteurization, [K].'''
        return self._temp_pasteurization
    @temp_pasteurization.setter
    def temp_pasteurization(self, i):
        self._temp_pasteurization = i

    @property
    def lhv_lpg(self):
        '''[float] Liquid Petroleum higher heating value, [MJ/kg].'''
        return self._lhv_lpg
    @lhv_lpg.setter
    def lhv_lpg(self, i):
        self._lhv_lpg = i
