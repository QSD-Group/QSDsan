#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Smiti Mittal <smitimittal@gmail.com>
    Yalin Li <zoe.yalin.li@gmail.com>
    Anna Kogler

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

# %%

from .. import SanUnit, WasteStream
from ..equipments import Column, Electrode, Machine, Membrane

__all__ = ('ElectrochemicalCell',)


class ElectrochemicalCell(SanUnit):
    '''
    Electrochemical cell for nutrient recovery.

    This unit has the following equipment:
        - :class:`~.equipments.Column`
        - :class:`~.equipments.Machine`
        - :class:`~.equipments.Electrode`
        - :class:`~.equipments.Membrane`

    Parameters
    ----------
    recovery : dict
        Keys refer to chemical component IDs. Values refer to recovery fractions (with 1 being 100%) for the respective chemicals.
    removal : dict
        Keys refer to chemical component IDs. Values refer to removal fractions (with 1 being 100%) for the respective chemicals.
    equipments : list
        List of Equipment objects part of the Electrochemical Cell.
    OPEX_over_CAPEX : float
        Ratio with which operating costs are calculated as a fraction of capital costs

    Example
    -------
    >>> # Set components
    >>> import qsdsan as qs
    >>> kwargs = dict(particle_size='Soluble',
    ...               degradability='Undegradable',
    ...               organic=False)
    >>> H2O = qs.Component.from_chemical('H2O', phase='l', **kwargs)
    >>> NH3 = qs.Component.from_chemical('NH3', phase='g', **kwargs)
    >>> NH3.particle_size = 'Dissolved gas'
    >>> NH4OH = qs.Component.from_chemical('NH4OH', phase='l', **kwargs)
    >>> H2SO4 = qs.Component.from_chemical('H2SO4', phase='l', **kwargs)
    >>> AmmoniumSulfate = qs.Component.from_chemical('AmmoniumSulfate', phase='l',
    ...                                              **kwargs)
    >>> CleaningAgent = qs.Component('CleaningAgent', MW=1, phase='l', **kwargs)
    >>> cmps = qs.Components((H2O, NH3, NH4OH, H2SO4, AmmoniumSulfate, CleaningAgent))
    >>> # Assuming all has the same molar volume as water for demonstration purpose
    >>> for cmp in cmps:
    ...     cmp.copy_models_from(H2O, names=['V'])
    ...     defaulted_properties = cmp.default()
    >>> qs.set_thermo(cmps)
    >>> # Set waste streams
    >>> influent = qs.WasteStream('influent', H2O=1000, NH4OH=50)
    >>> cleaning_agent = qs.WasteStream('cleaning_agent', price=5)
    >>> # Set the unit
    >>> U1 = qs.sanunits.ElectrochemicalCell('U1', ins=(influent, cleaning_agent),
    ...                                 outs=('rec', 'rem', 'leftover'),
    ...                                 recovery={'NH4OH':0.6}, removal={'NH4OH':0.2},
    ...                                 OPEX_over_CAPEX = 0.2)
    >>> # Simulate and look at the results
    >>> U1.simulate()
    >>> # U1.diagram() # have a look at the diagram
    >>> U1.show() # doctest: +SKIP
    ElectrochemicalCell: U1
    ins...
    [0] influent
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (g/hr): H2O    1e+06
                     NH4OH  5e+04
        WasteStream-specific properties:
         pH         : 7.0
         TN         : 19424.7 mg/L
         TKN        : 19424.7 mg/L
    [1] cleaning_agent
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow: 0
        WasteStream-specific properties: None for empty waste streams
    outs...
    [0] rec
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (g/hr): NH4OH  3e+04
        WasteStream-specific properties:
         pH         : 7.0
         TN         : 775176.8 mg/L
         TKN        : 775176.8 mg/L
    [1] rem
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (g/hr): NH4OH  1e+04
        WasteStream-specific properties:
         pH         : 7.0
         TN         : 775176.8 mg/L
         TKN        : 775176.8 mg/L
    [2] leftover
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (g/hr): H2O    1e+06
                     NH4OH  1e+04
        WasteStream-specific properties:
         pH         : 7.0
         TN         : 3964.4 mg/L
         TKN        : 3964.4 mg/L
    >>> U1.results() # doctest: +SKIP
    Electrochemical cell                             Units            U1
    Design              Number of columns                              3
                        Material of the column                     resin
                        Surface area of columns         m2            20
                        Number of anode               None             1
                        Material of anode             None      graphite
                        Surface area of anode           m2            10
                        Number of cathode             None             1
                        Material of cathode           None        carbon
                        Surface area of cathode         m2            10
                        Number of machines                             1
                        Number of membrane modules                     2
                        Material of membrane                polyethylene
                        Surface area of membrane        m2             1
    Purchase cost       U1_column                      USD           120
                        U1_anode                       USD           0.1
                        U1_cathode                     USD             1
                        U1_fan                         USD             3
                        U1_membrane                    USD           0.4
    Total purchase cost                                USD           124
    Utility cost                                    USD/hr             0
    Additional OPEX                                 USD/hr          24.9
    '''

    _N_ins = 2
    _N_outs = 3


    def __init__(self, ID='', ins=None, outs=(),
                 recovery={'NH3':0.6}, removal={'NH3':0.2}, OPEX_over_CAPEX=0):
        SanUnit.__init__(self=self, ID=ID, ins=ins, outs=outs)
        self.recovery = recovery
        self.removal = removal
        self.OPEX_over_CAPEX = OPEX_over_CAPEX
        self.equipments = (
            Column('column', linked_unit=self, N=3,
                   material='resin', unit_cost=2, surface_area=20),
            Electrode('anode', linked_unit=self, N=1, electrode_type='anode',
                      material='graphite', surface_area=10),
            Electrode('cathode', linked_unit=self, N=1, electrode_type='cathode',
                      material='carbon', surface_area=10, unit_cost=1),
            Machine('fan', linked_unit=self, N=1, unit_cost=3),
            Membrane('membrane', linked_unit=self, N=2,
                     material='polyethylene', unit_cost=0.2, surface_area=1),
            )

    def _run(self):
        influent, cleaner = self.ins
        recovered, removed, left = self.outs[0], self.outs[1], self.outs[2]

        mixture = WasteStream()
        mixture.mix_from(self.ins)
        left.copy_like(mixture)

        for chemical, ratio in self.recovery.items():
            recovered.imass[chemical] = mixture.imass[chemical]*ratio
            left.imass[chemical] = left.imass[chemical]-mixture.imass[chemical]*ratio

        for chemical, ratio in self.removal.items():
            removed.imass[chemical] = mixture.imass[chemical]*ratio
            left.imass[chemical] = left.imass[chemical]-mixture.imass[chemical]*ratio


    def _design(self):
        self.add_equipment_design()


    def _cost(self):
        self.add_equipment_cost()
        self.equip_costs = self.baseline_purchase_costs.values()
        add_OPEX = sum(self.equip_costs)*self.OPEX_over_CAPEX
        self._add_OPEX = {'Additional OPEX': add_OPEX}