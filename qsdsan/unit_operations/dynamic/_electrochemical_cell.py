#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:

    Zixuan Wang <wyatt4428@gmail.com>

    Smiti Mittal <smitimittal@gmail.com>

    Yalin Li <mailto.yalin.li@gmail.com>

    Anna Kogler <akogler@stanford.edu>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.

Note: `unit_operations/static/_electrochemical_cell.py` is a separate module
holding `ElectrochemicalCell`/`ElectrochemicalStrippingAdsorptionPrecipitation`,
the original monolithic (steady-state only) design that `ESAP`/`ESAPRecovery`/
`ESAPEffluent` here later replaced with a modular, dynamic-capable pipeline.
'''

# %%

from ... import SanUnit
from ._abstract import Splitter
from ...equipments import Electrode, Machine, Membrane
from math import ceil
import numpy as np

__all__ = ('ESAPRecovery', 'ESAPEffluent', 'ESAP',)


class ESAP(SanUnit):
    '''
    Splits an electrochemical stripping influent into a recovered-product
    stream, a loss stream, and an effluent, by independent component-wise
    recovery/loss mass fractions. This unit is the front-end split step of
    the modular ESAP process (paired with :class:`ESAPRecovery` and
    :class:`ESAPEffluent` downstream, which carry the equipment design and
    costing); it is able to perform dynamic simulation.

    Parameters
    ----------
    ins :
        Wastewater stream.
    outs :
        * [0] recovery product
        * [1] loss stream
        * [2] effluent
    recovery : dict
        Keys refer to chemical component IDs. Values refer to recovery fractions (with 1 being 100%) for the respective chemicals.
    loss : dict
        Keys refer to chemical component IDs. Values refer to loss fractions (with 1 being 100%) for the respective chemicals.
    '''
    _N_outs = 3
    _ins_size_is_fixed = False
    _outs_size_is_fixed = False
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                 recovery={'NH3':0.7}, loss={'NH3':0.06}, order=None,
                 init_with ='WasteStream', F_BM_default=None, isdynamic=False,
                  OPEX_over_CAPEX=0.2, component_ID_NH3 = 'NH3',
                  component_ID_P ='Phosphate', component_ID_Mg = "Mg2+", **kwargs): #0.015 kg N/h, product 2 mol/L
        SanUnit.__init__(self, ID=ID, ins=ins, outs=outs, thermo=thermo,
                         init_with = init_with, **kwargs)
        self.recovery = recovery
        self.loss = loss

    def _run(self):
        influent, = self.ins
        recovery_product, loss, effluent = self.outs
        mass = influent.mass.copy()

        # Initialize matrices
        components = influent.components
        recovery_matrix = np.zeros_like(mass)
        loss_matrix = np.zeros_like(mass)
        remaining_matrix = np.ones_like(mass)

        # Fill matrices based on recovery and loss dictionaries
        for component, recovery_fraction in self.recovery.items():
            idx = components.index(component)
            recovery_matrix[idx] = recovery_fraction
            remaining_matrix[idx] -= recovery_fraction

        for component, loss_fraction in self.loss.items():
            idx = components.index(component)
            loss_matrix[idx] = loss_fraction
            remaining_matrix[idx] -= loss_fraction

        # Calculate mass splits
        recovery_product.mass = mass * recovery_matrix
        loss.mass = mass * loss_matrix
        effluent.mass = mass * remaining_matrix

    @property
    def state(self):
        '''Component mass flow rate.'''
        if self._state is None: return None
        else:
            return dict(zip(list(self.components.IDs), self._state))

    def _init_state(self):
        influent = self.ins[0]
        components = self.ins[0].components
        self._state = influent.mass.copy()
        self._dstate = self._state * 0.
        self._recovery_matrix = np.zeros_like(self._state)
        self._loss_matrix = np.zeros_like(self._state)
        self._remaining_matrix = np.zeros_like(self._state)

        for component, recovery_fraction in self.recovery.items():
            idx = components.index(component)
            self._recovery_matrix[idx] = recovery_fraction
            self._remaining_matrix[idx] -= recovery_fraction

        for component, loss_fraction in self.loss.items():
            idx = components.index(component)
            self._loss_matrix[idx] = loss_fraction
            self._remaining_matrix[idx] -= loss_fraction

    def _update_state(self):
        arr = self._state
        self._outs[0].state = self._recovery_matrix * arr #TODO: why underscore here?
        self._outs[1].state = self._loss_matrix * arr
        self._outs[2].state = self._remaining_matrix * arr

    def _update_dstate(self):
        arr = self._dstate
        self._outs[0].dstate = self._recovery_matrix * arr
        self._outs[1].dstate = self._loss_matrix * arr
        self._outs[2].dstate = self._remaining_matrix * arr

    @property
    def AE(self):
        if self._AE is None:
            self._compile_AE()
        return self._AE

    def _compile_AE(self):
        _state = self._state
        _dstate = self._dstate
        _update_state = self._update_state
        _update_dstate = self._update_dstate
        def yt(t, mass_ins, dmass_ins):
            _state[:] = mass_ins[0]
            _dstate[:] = dmass_ins[0]
            _update_state()
            _update_dstate()
        self._AE = yt


#%%
class ESAPRecovery(Splitter):
    '''
    Electrochemical stripping, adsorption, and precipitation for nutrient recovery.
    This unit is able to perform dynamic simulation.

    This unit has the following equipment:
        - :class:`~.equipments.Column`
        - :class:`~.equipments.Machine`
        - :class:`~.equipments.Electrode`
        - :class:`~.equipments.Membrane`

    Parameters
    ----------
    ins:
        Wastewater stream
    outs:
        * [0] recovery product
        * [1] Remainder stream
    recovery : dict
        Keys refer to chemical component IDs. Values refer to recovery fractions (with 1 being 100%) for the respective chemicals.
    equipment : list(obj)
        List of Equipment objects part of the Electrochemical Cell.
    N_treatment_capacity: float
        kg N/h designed treatment capacity for 1 tower.
    OPEX_over_CAPEX : float
        Ratio with which operating costs are calculated as a fraction of capital costs
    component_ID_NH3: string
        The ID for ammonia/ammonium in the influent wastestream
    component_ID_P: string
        The ID for dissolved phosphate in the influent wastestream
    component_ID_Mg: string
        The ID for dissolved magnesium in the influent wastestream
    N_prodcut_concentration: float
        molar concentration (mol/L) of N in the final product
    '''
    _ins_size_is_fixed = False

    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                 recovery={'NH3':0.8,'Mg2+':0.7,'Phosphate':0.95}, order=None,
                 init_with ='WasteStream', F_BM_default=None, isdynamic=False,
                 N_treatment_capacity=0.015, OPEX_over_CAPEX=0.2, component_ID_NH3 = 'NH3',
                 component_ID_P ='Phosphate', component_ID_Mg = "Mg2+", N_prodcut_concentration = 2): #0.015 kg N/h, product 2 mol/L
        Splitter.__init__(self=self, ID=ID, ins=ins, outs=outs,thermo = thermo, split = recovery,
                          order = order, init_with = init_with, F_BM_default = F_BM_default,
                          isdynamic = isdynamic
                          )
        self.recovery = recovery
        self.OPEX_over_CAPEX = OPEX_over_CAPEX
        self.component_ID_NH3 = component_ID_NH3
        self.component_ID_Mg = component_ID_Mg
        self.N_treatment_capacity = N_treatment_capacity

        self.equipment = [
            Electrode('Anode', linked_unit=self, N=1, electrode_type='anode',
                      material='Ti MMO mesh',surface_area=0.359, unit_cost=2576.19),
            Electrode('Cathode', linked_unit=self, N=1, electrode_type='cathode',
                      material='SS mesh', surface_area=0.359, unit_cost=46.67),
            Membrane('Cation_Exchange_Membrane', linked_unit=self, N=1,
                     material='Selemion CMVN', unit_cost=1016.67, surface_area= 0.479), #$1016.67/m2
            Membrane('Gas_Permeable_Membrane', linked_unit=self, N=1,
                     material='Omniphobic', unit_cost=279.45, surface_area=0.2394),
            Machine('Pumps',linked_unit=self, N=3, lifetime= 43800, unit_cost =652.36), #pump lifetime 5 years
            ]

    def _design(self):
        self.add_equipment_design()
        D = self.design_results
        D['Number of ECS towers'] = ceil(self.outs[0].imass['NH3']/17*14/self.N_treatment_capacity)
        D['H2SO4'] = self.outs[0].imol[self.component_ID_NH3] * 0.5 * 98/0.96 #assume 1 mol NH3 requires 0.5 mol H2SO4, kg/hr 96% H2SO4
        D['NaOH'] = self.outs[0].imol[self.component_ID_Mg] * 2* 40/0.9 #assume 1 mol Mg requires 2 mol NaOH, kg/hr 90% NaOH

    def _cost(self):
        self.add_equipment_cost()
        C = self.baseline_purchase_costs
        C['Flanges'] = (47.83+#https://www.mcmaster.com/4881K241
                        113.79+#https://www.mcmaster.com/95665K324
                        37.5+#https://www.mcmaster.com/4881K239
                        78.98+#https://www.mcmaster.com/95665K322
                        19.15+#https://www.mcmaster.com/4881K236
                        37.65+#https://www.mcmaster.com/95665K217
                        61.53)*2#https://www.mcmaster.com/4881K967

        C['TowerWall'] =  (548.55+ #https://www.mcmaster.com/4740K32
                           287.17+#https://www.mcmaster.com/4740K31
                           118.64)#https://www.mcmaster.com/4740K26

        self.equip_costs = self.baseline_purchase_costs.values()
        add_OPEX = sum(self.equip_costs)*self.OPEX_over_CAPEX
        recovered = self.outs[0]

        self.power_utility.rate = recovered.imass['NH3']*0.67577
        # steady state value derived from 17.57 kWh used over 26 hrs
        self._add_OPEX = {'Additional OPEX': add_OPEX}

#%%
class ESAPEffluent(Splitter):
    '''
    Splitting the ESAP effluent. This unit is able to perform dynamic simulation.

    This unit has the following equipment:
        - :class:`~.equipments.Column`
        - :class:`~.equipments.Machine`
        - :class:`~.equipments.Electrode`
        - :class:`~.equipments.Membrane`

    Parameters
    ----------
    ins:
        Inlet fluid to be split
    outs:
        * [0] loss stream during the process
        * [1] effluent stream
    loss : dict
        Keys refer to chemical component IDs. Values refer to loss fractions (with 1 being 100%)
        for the respective chemicals. The fraction is respective to the influent of ESAP_effluent unit.
    equipment : list(obj)
        List of Equipment objects part of the Electrochemical Cell.
    OPEX_over_CAPEX : float
        Ratio with which operating costs are calculated as a fraction of capital costs
    '''
    _ins_size_is_fixed = False

    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                 loss={'NH3':0.5,'Mg2+':0.3,'Phosphate':0.8}, order=None,
                 init_with ='WasteStream', F_BM_default=None, isdynamic=False):
        Splitter.__init__(self=self, ID=ID, ins=ins, outs=outs,thermo = thermo, split = loss,
                          order = order, init_with = init_with, F_BM_default = F_BM_default,
                          isdynamic = isdynamic
                          )
        self.loss = loss
