#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Smiti Mittal <smitimittal@gmail.com>
    Yalin Li <mailto.yalin.li@gmail.com>
    Anna Kogler <akogler@stanford.edu>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.

Reference for the default electrochemical cell modelled below: 
    Evaluating Membrane Performance in Electrochemical Stripping Reactor for Nitrogen Removal
    Julia Simon, Department of Chemical Engineering, Stanford University
    Principal Investigator: Professor William Tarpeh
    Second Reader: Professor Gerald Fuller
'''

# %%

from qsdsan import SanUnit, WasteStream
from qsdsan.equipments import Column, Electrode, Machine, Membrane

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
    equipment : list(obj)
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
    >>> H2SO4 = qs.Component.from_chemical('H2SO4', phase='l', **kwargs)
    >>> AmmoniumSulfate = qs.Component.from_chemical('AmmoniumSulfate', phase='l',
    ...                                              **kwargs)
    >>> Na_ion = qs.Component.from_chemical('Na+', phase='s', **kwargs)
    >>> K_ion = qs.Component.from_chemical('K+', phase='s', **kwargs)
    >>> Cl_ion = qs.Component.from_chemical('Cl-', phase='s', **kwargs)
    >>> PO4_ion = qs.Component.from_chemical('Phosphate', phase='s', **kwargs)
    >>> SO4_ion = qs.Component.from_chemical('Sulfate', phase='s', **kwargs)
    >>> C = qs.Component.from_chemical('Carbon', phase='s', **kwargs)
    >>> COD = qs.Component.from_chemical('O2', phase='s', **kwargs)
    >>> cmps = qs.Components((H2O, NH3, H2SO4, AmmoniumSulfate, Na_ion, K_ion, Cl_ion, PO4_ion, SO4_ion, C, COD))
    >>> # Assuming all has the same molar volume as water for demonstration purpose
    >>> for cmp in cmps:
    ...     cmp.copy_models_from(H2O, names=['V'])
    ...     defaulted_properties = cmp.default()
    >>> qs.set_thermo(cmps)
    >>> # Set waste streams
    >>> influent = qs.WasteStream('influent')
    >>> influent.set_flow_by_concentration(flow_tot=0.5, concentrations={'NH3':3820,
    ...                                     'Na+':1620, 'K+':1470, 'Cl-':3060, 'Phosphate':169,
    ...                                     'Sulfate':1680, 'Carbon':1860, 'O2':3460}, units=('mL/min', 'mg/L'))
    >>> influent.show()
    WasteStream: influent
     phase: 'l', T: 298.15 K, P: 101325 Pa
     flow (g/hr): H2O        29.5
                  NH3        0.115
                  Na+        0.0486
                  K+         0.0441
                  Cl-        0.0918
                  Phosphate  0.00507
                  Sulfate    0.0504
                  Carbon     0.0558
                  O2         0.104
     WasteStream-specific properties:
      pH         : 7.0
      Alkalinity : 2.5 mg/L
      TC         : 1860.0 mg/L
      TP         : 31.9 mg/L
      TK         : 1470.0 mg/L
     Component concentrations (mg/L):
      H2O              984474.3
      NH3              3820.0
      Na+              1620.0
      K+               1470.0
      Cl-              3060.0
      Phosphate        169.0
      Sulfate          1680.0
      Carbon           1860.0
      O2               3460.0
    >>> catalysts = qs.WasteStream('catalysts', H2SO4=0.00054697, units=('kg/hr'))
    >>> # kg/hr imass value derived from 0.145 moles used on average for one, 26 hour experiment in the model cell
    >>> catalysts.price = 0.12
    >>> # electricity price for ECS experiments in the default model tested by the Tarpeh Lab
        power_utility.price = 0.1741
    
        # Set the unit
    >>> U1 = qs.sanunits.ElectrochemicalCell('U1', ins=(influent, catalysts), outs=('recovered', 'removed', 'residual'))
    >>> # Simulate and look at the results
    >>> U1.simulate()
    >>> U1.diagram() # have a look at the diagram
    >>> U1.results()
    Electrochemical cell                                      Units            U1
    Design              Column U1_column - Number of col...                     3
                        Column U1_column - Material of t...                 resin
                        Column U1_column - Surface area ...      m2            20
                        Electrode U1_anode - Number of a...    None             1
                        Electrode U1_anode - Material of...    None      graphite
                        Electrode U1_anode - Surface are...      m2            10
                        Electrode U1_cathode - Number of...    None             1
                        Electrode U1_cathode - Material ...    None        carbon
                        Electrode U1_cathode - Surface a...      m2            10
                        Machine U1_fan - Number of machines                     1
                        Membrane U1_membrane - Number of...                     2
                        Membrane U1_membrane - Material ...          polyethylene
                        Membrane U1_membrane - Surface a...      m2             1
    Purchase cost       U1_column                               USD           120
                        U1_anode                                USD           0.1
                        U1_cathode                              USD             1
                        U1_fan                                  USD             3
                        U1_membrane                             USD           0.4
    Total purchase cost                                         USD           124
    Utility cost                                             USD/hr             0
    Additional OPEX                                          USD/hr             0
    >>> U1.show()
    ElectrochemicalCell: U1
    ins...
    [0] influent
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (g/hr): H2O        29.5
                     NH3        0.115
                     Na+        0.0486
                     K+         0.0441
                     Cl-        0.0918
                     Phosphate  0.00507
                     Sulfate    0.0504
                     Carbon     0.0558
                     O2         0.104
        WasteStream-specific properties:
         pH         : 7.0
         TC         : 1860.0 mg/L
         TP         : 31.9 mg/L
         TK         : 1470.0 mg/L
    [1] catalysts
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (g/hr): H2SO4  0.547
        WasteStream-specific properties:
         pH         : 7.0
    outs...
    [0] recovered
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (g/hr): NH3  0.0688
        WasteStream-specific properties:
         pH         : 7.0
    [1] removed
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (g/hr): NH3  0.0229
        WasteStream-specific properties:
         pH         : 7.0
    [2] residual
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (g/hr): H2O        29.5
                     NH3        0.0229
                     H2SO4      0.547
                     Na+        0.0486
                     K+         0.0441
                     Cl-        0.0918
                     Phosphate  0.00507
                     Sulfate    0.0504
                     Carbon     0.0558
                     O2         0.104
        WasteStream-specific properties:
         pH         : 7.0
         TC         : 1859.8 mg/L
         TP         : 31.9 mg/L
         TK         : 1469.8 mg/L
    '''

    _N_ins = 2
    _N_outs = 3


    def __init__(self, ID='', ins=(), outs=(),
                 recovery={'NH3':0.7}, removal={'NH3':0.83, 'K+':0.83, "Na+":0.8}, OPEX_over_CAPEX=0.2):
        SanUnit.__init__(self=self, ID=ID, ins=ins, outs=outs)
        self.recovery = recovery
        self.removal = removal
        self.OPEX_over_CAPEX = OPEX_over_CAPEX


        self.equipment = [
            Column('column', linked_unit=self, N=3,
                   material='resin', unit_cost=2, surface_area=20),
            Electrode('Main Anode', linked_unit=self, N=1, electrode_type='anode',
                      material='iridium tantalum oxide coated titanium grid', surface_area=1, unit_cost=288), #in m^2
            Electrode('Main Cathode', linked_unit=self, N=1, electrode_type='cathode',
                      material='stainless steel wire mesh', surface_area=30.25, unit_cost=0.182), #in in^2
            Electrode('Current Collector Cathode', linked_unit=self, N=1, electrode_type='cathode',
                      material='stainless steel sheet', surface_area=38.5, unit_cost=1.037), #in in^2
            Electrode('Reference Electrode', linked_unit=self, N=1, electrode_type='reference',
                      material='Ag/AgCl with ceramic', surface_area=1, unit_cost=94), #in unknown units (94/unit, 1 unit)
            Machine('Pump', linked_unit=self, N=3, unit_cost=2800), #$1075 per pump drive, $1095 per pump head
            Membrane('Cation Exchange Membrane', linked_unit=self, N=1,
                     material='0.45 mm Gel Polystyrene cross linked with divinylbenzene with sulphonic acid functional groups',
                     unit_cost=0.182, surface_area=30.25), # in in^2
            Membrane('Gas Permeable Membrane', linked_unit=self, N=1,
                     material='Polyethylene 0.3 micron', unit_cost=0.032, surface_area=30.25), #in in^2
            ]


    def _run(self):
        influent, catalysts = self.ins
        recovered, removed, residual = self.outs[0], self.outs[1], self.outs[2]

        mixture = WasteStream()
        mixture.mix_from(self.ins)
        residual.copy_like(mixture)

        for chemical, recovery_ratio in self.recovery.items():
            recovered.imass[chemical] = mixture.imass[chemical]*recovery_ratio
            
        for chemical, removal_ratio in self.removal.items():
            recovery_ratio = 0 if recovered.imass[chemical] is None else recovered.imass[chemical]
            removed.imass[chemical] = mixture.imass[chemical]*removal_ratio - mixture.imass[chemical]*recovery_ratio
            residual.imass[chemical] = residual.imass[chemical]-residual.imass[chemical]*removal_ratio

    def _design(self):
        self.add_equipment_design()

    def _cost(self):
        self.add_equipment_cost()
        self.baseline_purchase_costs['Exterior'] = 60.52
        '''
        TOTAL CELL_EXTERIOR_COST = 60.52 USD
        Breakdown:
        Exterior frame (7'' x 7'' x 11/16'')	$21.46
        Interior half-cells (5.5'' x 5.5'' x 11/16'')	$19.87
        Rubber Sheets	$9.196
        Threaded Rods	$2.9696
        Wingnuts	$2.5504
        Flat Washers	$0.728
        Nylon Cable Glands	$3.744
        '''
        self.equip_costs = self.baseline_purchase_costs.values()
        add_OPEX = sum(self.equip_costs)*self.OPEX_over_CAPEX
        recovered, removed = self.outs[0], self.outs[1]

        self.power_utility.rate = recovered.imass['NH3']*0.67577
    # steady state value derived from 17.57 kWh used over 26 hrs
        self._add_OPEX = {'Additional OPEX': add_OPEX}