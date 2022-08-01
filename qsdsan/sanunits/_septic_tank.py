#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems
This module is developed by:
    Yalin Li <mailto.yalin.li@gmail.com>
    Hannah Lohman <hlohman94@gmail.com>
    Tori Morgan <vlmorgan@illinois.edu>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

from qsdsan import SanUnit, Construction
from ._decay import Decay
from ..utils import ospath, load_data, data_path, price_ratio

__all__ = ('SepticTank',)


septic_tank_path = ospath.join(data_path, 'sanunit_data/_septic_tank.csv')

@price_ratio()
class SepticTank(SanUnit, Decay):
    '''
    Septic tank that anaerobically treat the influent waste stream,
    often used as a primary treatment unit.

    Designed based on the Reclaimer system
    as described in http://washaid.pratt.duke.edu/sites/washaid.pratt.duke.edu/files/u71/Reclaimer_July2021.pdf
    and `Trotochaud et al. <https://doi.org/10.1021/acs.est.0c02755>`_

    The following impact items should be pre-constructed for life cycle assessment:
    FRP, Pump.

    Parameters
    ----------
    ins : Iterable(stream)
        waste: liquid waste stream to be treated by septic tank unit.
        MgOH2: input Mg(OH)2 for struvite precipitation.
    outs : Iterable(stream)
        treated: treated liquid leaving septic tank.
        CH4: fugitive CH4 emissions.
        N2O: fugitive N2O emissions.
        sludge: solid waste to be sent to sludge pasteurization, could include the precipitated struvite.
        struvite: a separate struvite stream when `if_generate_struvite` is True.
    if_include_front_end : bool
        If the front end is included in the analysis.
    if_generate_struvite : bool
        If generating struvite.
    if_struvite_in_sludge : bool
        If the generated struvite is in sludge.
    ppl : int
        Total number of users for scaling of costs.
    sludge_moisture_content : float
        Moisture content of the sludge, assumed to be 0.95 based on Tchobanoglous et al.
        (sludge leaving anaerobic treatment 2-5% solids).

    References
    ----------
    [1] 2019.06 Technical report for BMGF V3 _ CC 2019.06.13.pdf
    [2] Tchobanoglous et al., Wastewater Engineering: Treatment and Resource Recovery,
    McGraw-Hill Education, New York, 5th edn., 2013.

    See Also
    --------
    :ref:`qsdsan.sanunits.Decay <sanunits_Decay>`
    '''

    # Constants
    # Baseline population served by a single septic tank, 100 instead of the 30 of other units
    baseline_ppl = 100
    # Exponential scaling constant for scaling cost and LCA with change in users
    exponent_scale = 0.6  # exponential scaling constant

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 if_include_front_end=True, if_generate_struvite=True, if_struvite_in_sludge=True,
                 ppl=1, sludge_moisture_content=0.95, **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with, F_BM_default=1)
        self.if_include_front_end = if_include_front_end
        self.if_generate_struvite = if_generate_struvite
        self.if_struvite_in_sludge = if_struvite_in_sludge
        self.ppl = ppl
        self.sludge_moisture_content = sludge_moisture_content

        data = load_data(path=septic_tank_path)
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
        treated, CH4, N2O, sludge, struvite = self.outs

        treated.copy_like(waste)
        CH4.phase = N2O.phase = 'g'

        # COD decay
        COD_deg = waste.COD*waste.F_vol/1e3*self.COD_removal  # kg/hr
        remaining_COD = COD_deg * (1-self.COD_removal) # remaining COD in the treated liquid

        CH4_prcd = COD_deg * self.MCF_decay * self.max_CH4_emission
        CH4.imass['CH4'] = CH4_prcd

        # N decay
        N_loss = self.first_order_decay(k=self.decay_k_N, t=self.tau/365, max_decay=self.N_max_decay)
        N_loss_tot = N_loss*waste.TN/1e3*waste.F_vol
        NH3_rmd, NonNH3_rmd = \
            self.allocate_N_removal(N_loss_tot, waste.imass['NH3'])
        treated.imass['NH3'] = waste.imass['NH3'] - NH3_rmd
        treated.imass['NonNH3'] = waste.imass['NonNH3'] - NonNH3_rmd
        N2O.imass['N2O'] = N_loss_tot * self.N2O_EF_decay * 44/28

        # P recovery
        P_recovery = self.P_recovery
        treated.imass['P'] *= (1-P_recovery)

        sludge.empty()
        struvite.empty()

        if self.if_generate_struvite:
            cmps = self.components
            MW_P = 30.97
            MgOH2.imol['MagnesiumHydroxide'] = (self.orthoP_post/MW_P/self.Mg_dose)/24  # mol Mg(OH)2 per hr
            struvite_production_time = waste.imass['P'] * self.P_recovery / cmps.Struvite.i_P  # kg (NH4)MgPO4•6(H2O) / hr
            struvite_stream = sludge if self.if_struvite_in_sludge else struvite
            struvite_stream.imass['Struvite'] = struvite_production_time
        else: # no MgOH2 is added if not generating struvite
            MgOH2.empty()

        # Assuming `Mg`, `Ca`, and `OtherSS` all go to sludge
        sludge.copy_flow(treated, IDs=('Ca', 'Mg', 'OtherSS'), remove=True)
        sludge_mc = self.sludge_moisture_content
        sludge_water = sludge.F_mass/(1-sludge_mc)  # sludge total mass == sludge dry mass at this stage (no water has been added)
        treated_water = max(0, waste.imass['H2O']-sludge_water) # all mass goes to sludge if the target moisture content is higher than the incoming moisture content

        # Allocate solubles based on the water content
        ratio = treated_water/waste.imass['H2O']  # ratio that remains in the treated liquid
        treated_prior_allocation = treated.mass
        treated_after_allocation = treated.mass * ratio
        sludge.mass += treated_prior_allocation - treated_after_allocation
        sludge.imass['H2O'] = sludge_water
        treated.mass = treated_after_allocation
        treated.imass['H2O'] = treated_water

        # Update the COD content of treated liquid and sludge
        treated._COD = remaining_COD * ratio * 1e3 / treated.F_vol
        sludge._COD = remaining_COD * (1-ratio) * 1e3 / sludge.F_vol


    def _design(self):
        design = self.design_results
        design['FRP'] = FRP_quant = self.FRP_per_tank * self.user_scale_up  # fiber-reinforced plastic material
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


    @property
    def user_scale_up(self):
        '''[float] Scaling factor based on the user number.'''
        if self.ppl and self.baseline_ppl:
            # Don't scale smaller than 1/4 original septic tank
            return max(0.25, self.ppl/self.baseline_ppl)
        return 1 # No scaling if no information on the user number is provided