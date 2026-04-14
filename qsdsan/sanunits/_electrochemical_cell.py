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

from .._sanunit import SanUnit
from .._waste_stream import WasteStream
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
      Alkalinity : 2.5 mmol/L
      TC         : 1860.0 mg/L
      TP         : 31.9 mg/L
      TK         : 1470.0 mg/L
     Component concentrations (mg/L):
      H2O              984514.0
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
    >>> catalysts.price = 8.4 #in $/kg
    >>> # electricity price for ECS experiments in the default model tested by the Tarpeh Lab
    >>> # if want to change the price of electricity,
    >>> # use qs.PowerUtility, e.g., qs.PowerUtility.price = 0.1741
    
    >>> # Set the unit
    >>> U1 = qs.sanunits.ElectrochemicalCell('unit_1', ins=(influent, catalysts), outs=('recovered', 'removed', 'residual'))
    >>> # Simulate and look at the results
    >>> U1.simulate()
    >>> U1.results()
    Electrochemical cell                                      Units                              unit_1
    Electricity         Power                                    kW                            5.42e-05
                        Cost                                 USD/hr                            4.24e-06
    Design              Electrode unit_1_Main_Anode - Nu...    None                                   1
                        Electrode unit_1_Main_Anode - Ma...    None  titanium grid catalyst welded t...
                        Electrode unit_1_Main_Anode - Su...      m2                                   1
                        Electrode unit_1_Main_Cathode - ...    None                                   1
                        Electrode unit_1_Main_Cathode - ...    None  timesetl 3pcs stainless steel w...
                        Electrode unit_1_Main_Cathode - ...      m2                                30.2
                        Electrode unit_1_Current_Collect...    None                                   1
                        Electrode unit_1_Current_Collect...    None    stainless steel 26 gauge 5.5 x 7
                        Electrode unit_1_Current_Collect...      m2                                38.5
                        Electrode unit_1_Reference_Elect...    None                                   1
                        Electrode unit_1_Reference_Elect...    None  re-5b ag/agcl, 7.5 cm long, wit...
                        Electrode unit_1_Reference_Elect...      m2                                   1
                        Membrane unit_1_Cation_Exchange_...                                           1
                        Membrane unit_1_Cation_Exchange_...          CMI-7000S, polystyrene 0.45mm t...
                        Membrane unit_1_Cation_Exchange_...      m2                                30.2
                        Membrane unit_1_Gas_Permeable_Me...                                           1
                        Membrane unit_1_Gas_Permeable_Me...          Aquastill 0.3-micron polyethyle...
                        Membrane unit_1_Gas_Permeable_Me...      m2                                30.2
    Purchase cost       unit_1_Main_Anode                       USD                                 288
                        unit_1_Main_Cathode                     USD                               0.847
                        unit_1_Current_Collector_Cathode        USD                                39.9
                        unit_1_Reference_Electrode              USD                                  94
                        unit_1_Cation_Exchange_Membrane         USD                                 167
                        unit_1_Gas_Permeable_Membrane           USD                                29.3
                        Exterior                                USD                                60.5
    Total purchase cost                                         USD                                 679
    Utility cost                                             USD/hr                            4.24e-06
    Additional OPEX                                          USD/hr                                 136
    >>> U1.show() # doctest: +ELLIPSIS
    ElectrochemicalCell: unit_1
    ins...
    [0] influent
    phase: 'l', T: 298.15 K, P: 101325 Pa
    flow (g/hr): H2O        29.5
    ...
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
            Electrode('Main_Anode', linked_unit=self, N=1, electrode_type='anode',
                      material='Titanium grid catalyst welded to current collector tab both coated in iridium tantalum mixed metal oxide', surface_area=1, unit_cost=288), #288/unit, 1 unit
            Electrode('Main_Cathode', linked_unit=self, N=1, electrode_type='cathode',
                      material='TIMESETL 3pcs Stainless Steel Woven Wire 20 Mesh - 12"x8"(30x21cm) Metal Mesh Sheet 1mm Hole Great for Air Ventilation - A4', surface_area=30.25, unit_cost=0.847), #in in^2
            Electrode('Current_Collector_Cathode', linked_unit=self, N=1, electrode_type='cathode',
                      material='Stainless Steel 26 gauge 5.5'' x 7''', surface_area=38.5, unit_cost=39.9245), #in unknown units (94/unit, 1 unit)
            Electrode('Reference_Electrode', linked_unit=self, N=1, electrode_type='reference',
                      material='RE-5B Ag/AgCl, 7.5 cm long, with ceramic (MF-2056)', surface_area=1, unit_cost=94), #in unknown units (94/unit, 1 unit)
            Membrane('Cation_Exchange_Membrane', linked_unit=self, N=1,
                     material='CMI-7000S, polystyrene 0.45mm thick [48'' x 20'']',
                     unit_cost=5.5055, surface_area=30.25), # in in^2
            Membrane('Gas_Permeable_Membrane', linked_unit=self, N=1,
                     material='Aquastill 0.3-micron polyethylene membrane', unit_cost=0.968, surface_area=30.25), #in in^2
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