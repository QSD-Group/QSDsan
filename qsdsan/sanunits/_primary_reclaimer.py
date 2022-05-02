#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Tori Morgan <tvlmorgan@gmail.com>
    Hannah Lohman <hlohman94@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''


from qsdsan import SanUnit, Construction
from ._decay import Decay
from ..utils import ospath, load_data, data_path, price_ratio

__all__ = ('PrimaryReclaimer',)

primary_reclaimer_path = ospath.join(data_path, 'sanunit_data/_primary_reclaimer.csv')


@price_ratio(default_price_ratio=1)
class PrimaryReclaimer(SanUnit, Decay):
    '''
    Anaerobic digestion of waste in septic tank.

    The following impact items should be pre-constructed for life cycle assessment:
    FRP, Pump.

    Parameters
    ----------
    if_include_front_end: bool
        If front end is included in analysis.
    ppl: int
        Total number of users for scaling of costs.

    ins:
        waste: liquid waste stream to be treated by septic tank unit
        MgOH2: input Mg(OH)2 for struvite precipitation

    outs:
        treated: treated liquid leaving septic tank
        CH4: fugitive CH4 emissions
        N2O: fugitive N2O emissions
        sludge: solid waste to be sent to sludge pasteurization
        struvite: precipitated struvite recovered

    References
    ----------
    [1] 2019.06 Technical report for BMGF V3 _ CC 2019.06.13.pdf

    See Also
    --------
    :ref:`qsdsan.sanunits.Decay <sanunits_Decay>`
    
    '''

    # Constants
    # Baseline population served by a single septic tank
    baseline_ppl = 100
    # Exponential scaling constant for scaling cost and LCA with change in users
    exponent_scale = 0.6  # exponential scaling constant
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 if_include_front_end=True, ppl=1, **kwargs):

        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with, F_BM_default=1)

        self.if_include_front_end = if_include_front_end
        self.ppl = ppl

        if self.ppl < 25:
            self.user_scale_up = 0.25  # don't scale smaller than 1/4 original septic tank
        else:
            self.user_scale_up = self.ppl / self.baseline_ppl  # users exceed the capacity of a standard septic tank

        data = load_data(path=primary_reclaimer_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data
        
        for attr, value in kwargs.items():
            setattr(self, attr, value)
    
    _N_ins = 2
    _N_outs = 5

    def _run(self):
        waste, MgOH2 = self.ins
        self.ins[1].imass['MagnesiumHydroxide'] = (self.orthoP_post / self.MW_P / self.Mg_dose * self.MW_MgOH2)/24  # kg Mg(OH)2 per hr;

        treated, CH4, N2O, sludge, struvite = self.outs
        treated.copy_like(self.ins[0])
        CH4.phase = N2O.phase = 'g'
        
        P_precipitated = self.orthoP_post * self.P_recovery  # kg precipitated-P/hr
        
        struvite_production_time = P_precipitated / self.MW_P * self.MW_struvite  # kg (NH4)MgPO4•6(H2O) / hr
        struvite.imass['Struvite'] = struvite_production_time
        
        # COD removal
        COD_deg = treated.COD*treated.F_vol/1e3*self.COD_removal  # kg/hr
        treated._COD = waste.COD * (1-self.COD_removal)
        
        CH4_prcd = COD_deg * self.MCF_decay * self.max_CH4_emission
        CH4.imass['CH4'] = CH4_prcd
        N_loss = self.first_order_decay(k=self.decay_k_N, t=self.tau/365, max_decay=self.N_max_decay)
        
        N_loss_tot = N_loss*waste.TN/1e3*waste.F_vol 
        NH3_rmd, NonNH3_rmd = \
            self.allocate_N_removal(N_loss_tot, waste.imass['NH3'])
        treated.imass['NH3'] = waste.imass['NH3'] - NH3_rmd
        treated.imass['NonNH3'] = waste.imass['NonNH3'] - NonNH3_rmd
        N2O.imass['N2O'] = N_loss_tot * self.N2O_EF_decay * 44/28
        
        # sludge production
        self.sludge_TN = waste.imass['N'] * self.N_max_decay
        self.sludge_TN_F_mass = self.sludge_TN * waste.F_vol * 1e-3
        sludge.imass['N'] = self.sludge_TN
        
        self.sludge_COD = waste._COD * self.COD_removal 
        self.sludge_COD_F_mass = self.sludge_COD * waste.F_vol * 1e-3
        sludge._COD = self.sludge_COD 
        
        sludge.imass['H2O'] = sludge.F_mass * 0.5

    def _design(self):
        design = self.design_results
        design['FRP'] = FRP_quant = self.FRP_per_tank * self.user_scale_up  # Fibre-reinforced plastic material
        design['Pump'] = pump_quant = self.qty_pump * self.user_scale_up

        self.construction = (Construction(item='FRP', quantity=FRP_quant, quantity_unit='kg'),
                             Construction(item='Pump', quantity=pump_quant, quantity_unit='each'))

        self.add_construction(add_cost=False)
 
    def _cost(self):
        if self.if_include_front_end:
            C = self.baseline_purchase_costs
            C['Tanks'] = self.FRP_tank_cost * (self.user_scale_up ** self.exponent_scale)
            C['Pump'] = self.pump_cost * (self.user_scale_up ** self.exponent_scale)
            ratio = self.price_ratio
            for equipment, cost in C.items():
                C[equipment] = cost * ratio
        else:
            self.baseline_purchase_costs.clear()
