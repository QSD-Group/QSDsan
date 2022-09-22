#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems
This module is developed by:
    Yalin Li <mailto.yalin.li@gmail.com>
    Shion Watabe <shionwatabe@gmail.com>
    Hannah Lohman <hlohman94@gmail.com>
    Tori Morgan <vlmorgan@illinois.edu>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

from math import ceil
from thermosteam.reaction import ParallelReaction
from .. import SanUnit, Construction
from ..utils import ospath, load_data, data_path, price_ratio

__all__ = ('SludgePasteurization',)

pasteurization_path = ospath.join(data_path, 'sanunit_data/_sludge_pasteurization.tsv')

@price_ratio()
class SludgePasteurization(SanUnit):

    '''
    Unit operation for the pasteurization of sludge using liquid petroleum gas (LPG)
    biogas or biogas.

    Parameters
    ----------
    ins : Iterable(stream)
        biogas: biogas input from the AnMBR and used for pasteurization (biogas combusted).
        air: used for combustion, from atmosphere.
        sludge: sludge produced from the AnMBR in the NEWgenerator.
        LPG: purchased LPG to supplement heat requirement for pasteurization.
    outs : Iterable(stream)
        used: biogas used in combustion.
        lost: biogas lost in combustion (e.g., leaked as fugitive biogas).
        wasted: biogas wasted in combustion.
        treated sludge: sludge treated from pasteurization.
    if_combustion : bool
        If include combustion reaction during simulation.
    biogas_loss : float
        Fraction of biogas loss is 0.1 (e.g., leaked).
    biogas_eff : float
        Combustion efficiency of biogas as a fraction of CH4.
    temp_pasteurization : float
        Pasteurization temperature is 70°C or 343.15 K.
    sludge_temp : float
        Temperature of sludge is 10°C or 283.15 K.
    target_MC : float
        Target moisture content is 10%
    heat_loss : float
        Heat loss during pasteurization process is assumed to be 10%
    if_biogas : bool
        If biogas is used for sludge pasteurization, otherwise LPG is used
        and biogas combusted.
    lhv_lpg : float
        Lower heating value of LPG at 298.15K is 46-51 MJ/kg
        based on World Nuclear Organization.
    lhv_methane : float
        Lower heating value of methane at 298.15K is 50-55 MJ/kg
        based on World Nuclear Organization.
    ppl: int
        Total number of users for scaling of costs.
    baseline_ppl : int
        Baseline capacity of the unit for scaling of costs.
    user_scale_up : float
        Scaling factor (calculated based on the number of users as in `ppl` compared
        to the capacity as in `baseline_ppl`)
        for consumables, electricity demand, capital parts, and replacement parts.
        If not given (or set to None), will be calculated based on the user number
        and the capacity of this unit (i.e., the `baseline_ppl` attr).
    exponent_scale : float
        Exponential factor for the scaling up of capital costs.
    if_sludge_service: bool
        If share sludge pasteurization unit among multiple septic tanks
        (assume 1,000 users per sludge pasteurization unit,
         or 10 septic tanks serving a population of 100 users per septic tank).

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
    '''

    # Specific Heat capacity of water
    Cp_w = 4.184 # kJ kg^-1 K^-1
    # Specific Heat capacity of dry matter (sludge)
    Cp_dm = 1.231 # kJ kg^-1 K^-1
    # Specific latent heat of vaporization of water
    l_w = 2260 # kJ kg^-1

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                  if_biogas=True, heat_loss=0.1, target_MC=0.1, sludge_temp=283.15,
                  temp_pasteurization=343.15, if_combustion=False, biogas_loss=0.1,
                  biogas_eff=0.55,lhv_lpg = 48.5, lhv_methane=52.5,
                  ppl=100, baseline_ppl=100, user_scale_up=1, exponent_scale=0.6,
                  if_sludge_service=True, **kwargs):

        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with,
                          F_BM_default=1)
        self.if_combustion = if_combustion
        self.biogas_loss = biogas_loss
        self.biogas_eff = biogas_eff
        self.if_biogas = if_biogas
        self.heat_loss = heat_loss
        self.target_MC = target_MC
        self.sludge_temp = sludge_temp
        self.temp_pasteurization = temp_pasteurization
        self.lhv_methane = lhv_methane
        self.lhv_lpg = lhv_lpg
        self.ppl = ppl
        self.baseline_ppl = baseline_ppl
        self.user_scale_up = user_scale_up
        self.exponent_scale = exponent_scale
        self.if_sludge_service = if_sludge_service

        data = load_data(path=pasteurization_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)

    _N_ins = 4
    _N_outs = 4


    def _run(self):
        biogas, air, sludge, lpg = self.ins
        biogas.phase = air.phase = 'g'
        lpg.phase = 'l'
        used, lost, wasted, treated_sludge = self.outs
        used.copy_like(biogas)
        lost.copy_like(biogas)
        wasted.copy_like(biogas)
        treated_sludge.copy_like(sludge)
        lost.mass *= self.biogas_loss
        used.mass -= lost.mass
        wasted.mass = used.mass * (1-self.biogas_eff)
        used.mass -= wasted.mass

        if self.if_combustion:
            rxns = []
            for component in self.components:
                try:
                    rxn = component.get_combustion_reaction()
                    if rxn is not None: rxns.append(rxn)
                except:
                    continue
            rxns = ParallelReaction(rxns)
            rxns.force_reaction(used.mol)
            air.imol['O2'] = -used.imol['O2']
            used.imol['O2'] = 0.
            air.imol['N2'] = 0.79/0.21 * air.imol['O2']
        else:
            air.empty()

        # Mass calculations
        # total amount of water in sludge
        M_w = sludge.imass['H2O']  # kg/hr
        # total amount of dry matter in sludge
        M_dm = sludge.F_mass - M_w  # kg/hr
        # total amount of water to be evaporated to reach target dry matter content
        target_mc = self.target_MC
        treated_sludge.imass['H2O'] = M_w_sludge = M_dm / (1 - target_mc) - M_dm
        M_we = M_w - M_w_sludge  # kg/hr

        # Overall heat required for pasteurization
        temp_diff = self.temp_pasteurization - self.sludge_temp
        Q_d = (M_w * self.Cp_w + M_dm * self.Cp_dm) * temp_diff + M_we * self.l_w  # kJ/hr
        Q_tot = Q_d/(1-self.heat_loss)/1e3 # MJ/hr, 10% of the total generated is lost

        lhv_methane = self.lhv_methane
        lhv_lpg = self.lhv_lpg
        if self.if_biogas:
            methane_kg_reqd = Q_tot / lhv_methane # kg-CH4/hr = (MJ/hr)/(MJ/kg-CH4)
            if biogas.imass['CH4'] >= methane_kg_reqd:
                lpg.empty()
            else: # extra biogas, converted to LPG-equivalent
                lpg_Q_d = (methane_kg_reqd-biogas.imass['CH4'])*lhv_methane #MJ/hr = (kg-CH4/hr)*(MJ/kg-CH4)
                lpg.imass['LPG'] = lpg_Q_d / lhv_lpg # kg-LPG/hr = (MJ/hr)/(MJ/kg-LPG)
        else:
            lpg.imass['LPG'] = Q_tot / lhv_lpg # kg-LPG/hr = (MJ/hr)/(MJ/kg-LPG)
            for biogas_stream in (used, lost, wasted): biogas_stream.empty()


    def _design(self):
        design = self.design_results
        if self.if_sludge_service:
            design['Steel'] = S_quant = (self.sludge_dryer_weight + self.sludge_barrel_weight)/10*self.user_scale_up
        else:
            design['Steel'] = S_quant = (self.sludge_dryer_weight + self.sludge_barrel_weight)*self.user_scale_up

        self.construction = (
            Construction(item='Steel', quantity = S_quant, quantity_unit = 'kg'),
            )
        self.add_construction(add_cost=False)

    def _cost(self):
        C = self.baseline_purchase_costs
        factor = self.user_scale_up ** self.exponent_scale
        if self.if_sludge_service:
            C['Dryer'] = self.sludge_dryer / 10 * factor
            C['Barrel'] = self.sludge_barrel / 10 * factor
        else:
            C['Dryer'] = self.sludge_dryer * factor
            C['Barrel'] = self.sludge_barrel * factor
        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio

        self.add_OPEX = self._calc_replacement_cost() + self._calc_maintenance_labor_cost()


    # Assume dryer and barrel have 25 year lifetime so will not need to be replaced
    def _calc_replacement_cost(self):
        return 0

    def _calc_maintenance_labor_cost(self):
        sludge_maintenance_labor_cost = self.sludge_labor_maintenance * self.wages * (self.user_scale_up**self.exponent_scale)
        return sludge_maintenance_labor_cost/(365 * 24)  # USD/hr


    @property
    def user_scale_up(self):
        '''
        [float] Scaling factor (calculated based on the number of users
        as in `ppl` compared to the capacity as in `baseline_ppl`)
        for consumables, electricity demand, capital parts, and replacement parts.
        If not given (or set to None), will be calculated based on the user number
        and the capacity of this unit (i.e., the `baseline_ppl` attr).
        '''
        if self._user_scale_up: return self._user_scale_up
        return ceil(self.ppl / self.baseline_ppl)
    @user_scale_up.setter
    def user_scale_up(self, i):
        self._user_scale_up = i