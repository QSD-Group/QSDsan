# # <<<<<<< Updated upstream
# # # #!/usr/bin/env python3
# # # # -*- coding: utf-8 -*-

# # # '''
# # # QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

# # # This module is developed by:
# # #     Zixuan Wang <wyatt4428@gmail.com>
# # #     Smiti Mittal <smitimittal@gmail.com>
# # #     Yalin Li <mailto.yalin.li@gmail.com>
# # #     Anna Kogler <akogler@stanford.edu>

# # # This module is under the University of Illinois/NCSA Open Source License.
# # # Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
# # # for license details.

# # # Reference for the default electrochemical cell modelled below: 
# # #     Evaluating Membrane Performance in Electrochemical Stripping Reactor for Nitrogen Removal
# # #     Julia Simon, Department of Chemical Engineering, Stanford University
# # #     Principal Investigator: Professor William Tarpeh
# # #     Second Reader: Professor Gerald Fuller
# # # '''

# # # # %%

# # from .._sanunit import SanUnit
# # from ._abstract import Splitter
# =======
# # #!/usr/bin/env python3
# # # -*- coding: utf-8 -*-

# # '''
# # QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

# # This module is developed by:
# #     Smiti Mittal <smitimittal@gmail.com>
# #     Yalin Li <mailto.yalin.li@gmail.com>
# #     Anna Kogler <akogler@stanford.edu>

# # This module is under the University of Illinois/NCSA Open Source License.
# # Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
# # for license details.

# # Reference for the default electrochemical cell modelled below: 
# #     Evaluating Membrane Performance in Electrochemical Stripping Reactor for Nitrogen Removal
# #     Julia Simon, Department of Chemical Engineering, Stanford University
# #     Principal Investigator: Professor William Tarpeh
# #     Second Reader: Professor Gerald Fuller
# # '''

# # # %%

# # from .._sanunit import SanUnit
# >>>>>>> Stashed changes
# # from .._waste_stream import WasteStream
# # from ..equipments import Electrode, Machine, Membrane
# # from math import ceil
# # import numpy as np

# # __all__ = ('ElectrochemicalCell',
# #            'ESAPRecovery',
# #            'ESAPEffluent',
# #            'ESAP',
# #            'ElectrochemicalStrippingAdsorptionPrecipitation',)

# <<<<<<< Updated upstream
# # # class ElectrochemicalCell(SanUnit):
# # #     '''
# # #     Electrochemical cell for nutrient recovery.

# # #     This unit has the following equipment:
# # #         - :class:`~.equipments.Column`
# # #         - :class:`~.equipments.Machine`
# # #         - :class:`~.equipments.Electrode`
# # #         - :class:`~.equipments.Membrane`

# # #     Parameters
# # #     ----------
# # #     recovery : dict
# # #         Keys refer to chemical component IDs. Values refer to recovery fractions (with 1 being 100%) for the respective chemicals.
# # #     removal : dict
# # #         Keys refer to chemical component IDs. Values refer to removal fractions (with 1 being 100%) for the respective chemicals.
# # #     equipment : list(obj)
# # #         List of Equipment objects part of the Electrochemical Cell.
# # #     OPEX_over_CAPEX : float
# # #         Ratio with which operating costs are calculated as a fraction of capital costs

# # #     Example
# # #     -------
# # #     >>> # Set components
# # #     >>> import qsdsan as qs
# # #     >>> kwargs = dict(particle_size='Soluble',
# # #     ...               degradability='Undegradable',
# # #     ...               organic=False)
# # #     >>> H2O = qs.Component.from_chemical('H2O', phase='l', **kwargs)
# # #     >>> NH3 = qs.Component.from_chemical('NH3', phase='g', **kwargs)
# # #     >>> NH3.particle_size = 'Dissolved gas'
# # #     >>> H2SO4 = qs.Component.from_chemical('H2SO4', phase='l', **kwargs)
# # #     >>> AmmoniumSulfate = qs.Component.from_chemical('AmmoniumSulfate', phase='l',
# # #     ...                                              **kwargs)
# # #     >>> Na_ion = qs.Component.from_chemical('Na+', phase='s', **kwargs)
# # #     >>> K_ion = qs.Component.from_chemical('K+', phase='s', **kwargs)
# # #     >>> Cl_ion = qs.Component.from_chemical('Cl-', phase='s', **kwargs)
# # #     >>> PO4_ion = qs.Component.from_chemical('Phosphate', phase='s', **kwargs)
# # #     >>> SO4_ion = qs.Component.from_chemical('Sulfate', phase='s', **kwargs)
# # #     >>> C = qs.Component.from_chemical('Carbon', phase='s', **kwargs)
# # #     >>> COD = qs.Component.from_chemical('O2', phase='s', **kwargs)
# # #     >>> cmps = qs.Components((H2O, NH3, H2SO4, AmmoniumSulfate, Na_ion, K_ion, Cl_ion, PO4_ion, SO4_ion, C, COD))
# # #     >>> # Assuming all has the same molar volume as water for demonstration purpose
# # #     >>> for cmp in cmps:
# # #     ...     cmp.copy_models_from(H2O, names=['V'])
# # #     ...     defaulted_properties = cmp.default()
# # #     >>> qs.set_thermo(cmps)
# # #     >>> # Set waste streams
# # #     >>> influent = qs.WasteStream('influent')
# # #     >>> influent.set_flow_by_concentration(flow_tot=0.5, concentrations={'NH3':3820,
# # #     ...                                     'Na+':1620, 'K+':1470, 'Cl-':3060, 'Phosphate':169,
# # #     ...                                     'Sulfate':1680, 'Carbon':1860, 'O2':3460}, units=('mL/min', 'mg/L'))
# # #     >>> influent.show()
# # #     WasteStream: influent
# # #      phase: 'l', T: 298.15 K, P: 101325 Pa
# # #      flow (g/hr): H2O        29.5
# # #                   NH3        0.115
# # #                   Na+        0.0486
# # #                   K+         0.0441
# # #                   Cl-        0.0918
# # #                   Phosphate  0.00507
# # #                   Sulfate    0.0504
# # #                   Carbon     0.0558
# # #                   O2         0.104
# # #      WasteStream-specific properties:
# # #       pH         : 7.0
# # #       Alkalinity : 2.5 mmol/L
# # #       TC         : 1860.0 mg/L
# # #       TP         : 31.9 mg/L
# # #       TK         : 1470.0 mg/L
# # #      Component concentrations (mg/L):
# # #       H2O              984514.0
# # #       NH3              3820.0
# # #       Na+              1620.0
# # #       K+               1470.0
# # #       Cl-              3060.0
# # #       Phosphate        169.0
# # #       Sulfate          1680.0
# # #       Carbon           1860.0
# # #       O2               3460.0
# # #     >>> catalysts = qs.WasteStream('catalysts', H2SO4=0.00054697, units=('kg/hr'))
# # #     >>> # kg/hr imass value derived from 0.145 moles used on average for one, 26 hour experiment in the model cell
# # #     >>> catalysts.price = 8.4 #in $/kg
# # #     >>> # electricity price for ECS experiments in the default model tested by the Tarpeh Lab
# # #     >>> # if want to change the price of electricity,
# # #     >>> # use qs.PowerUtility, e.g., qs.PowerUtility.price = 0.1741
    
# # #     >>> # Set the unit
# # #     >>> U1 = qs.sanunits.ElectrochemicalCell('unit_1', ins=(influent, catalysts), outs=('recovered', 'removed', 'residual'))
# # #     >>> # Simulate and look at the results
# # #     >>> U1.simulate()
# # #     >>> U1.results()
# # #     Electrochemical cell                                      Units                              unit_1
# # #     Electricity         Power                                    kW                            5.42e-05
# # #                         Cost                                 USD/hr                            4.24e-06
# # #     Design              Electrode unit_1_Main_Anode - Nu...    None                                   1
# # #                         Electrode unit_1_Main_Anode - Ma...    None  titanium grid catalyst welded t...
# # #                         Electrode unit_1_Main_Anode - Su...      m2                                   1
# # #                         Electrode unit_1_Main_Cathode - ...    None                                   1
# # #                         Electrode unit_1_Main_Cathode - ...    None  timesetl 3pcs stainless steel w...
# # #                         Electrode unit_1_Main_Cathode - ...      m2                                30.2
# # #                         Electrode unit_1_Current_Collect...    None                                   1
# # #                         Electrode unit_1_Current_Collect...    None    stainless steel 26 gauge 5.5 x 7
# # #                         Electrode unit_1_Current_Collect...      m2                                38.5
# # #                         Electrode unit_1_Reference_Elect...    None                                   1
# # #                         Electrode unit_1_Reference_Elect...    None  re-5b ag/agcl, 7.5 cm long, wit...
# # #                         Electrode unit_1_Reference_Elect...      m2                                   1
# # #                         Membrane unit_1_Cation_Exchange_...                                           1
# # #                         Membrane unit_1_Cation_Exchange_...          CMI-7000S, polystyrene 0.45mm t...
# # #                         Membrane unit_1_Cation_Exchange_...      m2                                30.2
# # #                         Membrane unit_1_Gas_Permeable_Me...                                           1
# # #                         Membrane unit_1_Gas_Permeable_Me...          Aquastill 0.3-micron polyethyle...
# # #                         Membrane unit_1_Gas_Permeable_Me...      m2                                30.2
# # #     Purchase cost       unit_1_Main_Anode                       USD                                 288
# # #                         unit_1_Main_Cathode                     USD                               0.847
# # #                         unit_1_Current_Collector_Cathode        USD                                39.9
# # #                         unit_1_Reference_Electrode              USD                                  94
# # #                         unit_1_Cation_Exchange_Membrane         USD                                 167
# # #                         unit_1_Gas_Permeable_Membrane           USD                                29.3
# # #                         Exterior                                USD                                60.5
# # #     Total purchase cost                                         USD                                 679
# # #     Utility cost                                             USD/hr                            4.24e-06
# # #     Additional OPEX                                          USD/hr                                 136
# # #     >>> U1.show() # doctest: +ELLIPSIS
# # #     ElectrochemicalCell: unit_1
# # #     ins...
# # #     [0] influent
# # #     phase: 'l', T: 298.15 K, P: 101325 Pa
# # #     flow (g/hr): H2O        29.5
# # #     ...
# # #     '''

# # #     _N_ins = 2
# # #     _N_outs = 3


# # #     def __init__(self, ID='', ins=(), outs=(),
# # #                  recovery={'NH3':0.7}, removal={'NH3':0.83, 'K+':0.83, "Na+":0.8}, OPEX_over_CAPEX=0.2):
# # #         SanUnit.__init__(self=self, ID=ID, ins=ins, outs=outs)
# # #         self.recovery = recovery
# # #         self.removal = removal
# # #         self.OPEX_over_CAPEX = OPEX_over_CAPEX


# # #         self.equipment = [
# # #             Electrode('Main_Anode', linked_unit=self, N=1, electrode_type='anode',
# # #                       material='Titanium grid catalyst welded to current collector tab both coated in iridium tantalum mixed metal oxide', surface_area=1, unit_cost=288), #288/unit, 1 unit
# # #             Electrode('Main_Cathode', linked_unit=self, N=1, electrode_type='cathode',
# # #                       material='TIMESETL 3pcs Stainless Steel Woven Wire 20 Mesh - 12"x8"(30x21cm) Metal Mesh Sheet 1mm Hole Great for Air Ventilation - A4', surface_area=30.25, unit_cost=0.847), #in in^2
# # #             Electrode('Current_Collector_Cathode', linked_unit=self, N=1, electrode_type='cathode',
# # #                       material='Stainless Steel 26 gauge 5.5'' x 7''', surface_area=38.5, unit_cost=39.9245), #in unknown units (94/unit, 1 unit)
# # #             Electrode('Reference_Electrode', linked_unit=self, N=1, electrode_type='reference',
# # #                       material='RE-5B Ag/AgCl, 7.5 cm long, with ceramic (MF-2056)', surface_area=1, unit_cost=94), #in unknown units (94/unit, 1 unit)
# # #             Membrane('Cation_Exchange_Membrane', linked_unit=self, N=1,
# # #                      material='CMI-7000S, polystyrene 0.45mm thick [48'' x 20'']',
# # #                      unit_cost=5.5055, surface_area=30.25), # in in^2
# # #             Membrane('Gas_Permeable_Membrane', linked_unit=self, N=1,
# # #                      material='Aquastill 0.3-micron polyethylene membrane', unit_cost=0.968, surface_area=30.25), #in in^2
# # #             ]


# # #     def _run(self):
# # #         influent, catalysts = self.ins
# # #         recovered, removed, residual = self.outs[0], self.outs[1], self.outs[2]

# # #         mixture = WasteStream()
# # #         mixture.mix_from(self.ins)
# # #         residual.copy_like(mixture)

# # #         for chemical, recovery_ratio in self.recovery.items():
# # #             recovered.imass[chemical] = mixture.imass[chemical]*recovery_ratio
            
# # #         for chemical, removal_ratio in self.removal.items():
# # #             recovery_ratio = 0 if recovered.imass[chemical] is None else recovered.imass[chemical]
# # #             removed.imass[chemical] = mixture.imass[chemical]*removal_ratio - mixture.imass[chemical]*recovery_ratio
# # #             residual.imass[chemical] = residual.imass[chemical]-residual.imass[chemical]*removal_ratio

# # #     def _design(self):
# # #         self.add_equipment_design()

# # #     def _cost(self):
# # #         self.add_equipment_cost()
# # #         self.baseline_purchase_costs['Exterior'] = 60.52
# # #         '''
# # #         TOTAL CELL_EXTERIOR_COST = 60.52 USD
# # #         Breakdown:
# # #         Exterior frame (7'' x 7'' x 11/16'')	$21.46
# # #         Interior half-cells (5.5'' x 5.5'' x 11/16'')	$19.87
# # #         Rubber Sheets	$9.196
# # #         Threaded Rods	$2.9696
# # #         Wingnuts	$2.5504
# # #         Flat Washers	$0.728
# # #         Nylon Cable Glands	$3.744
# # #         '''
# # #         self.equip_costs = self.baseline_purchase_costs.values()
# # #         add_OPEX = sum(self.equip_costs)*self.OPEX_over_CAPEX
# # #         recovered, removed = self.outs[0], self.outs[1]

# # #         self.power_utility.rate = recovered.imass['NH3']*0.67577
# # #         # steady state value derived from 17.57 kWh used over 26 hrs
# # #         self._add_OPEX = {'Additional OPEX': add_OPEX}
# # # #%%
# # # class ESAP(SanUnit):
# # #     '''
# # #     Electrochemical stripping, adsorption, and precipitation for nutrient recovery. 
# # #     This unit is able to perform dynamic simulation.

# # #     This unit has the following equipment:
# # #         - :class:`~.equipments.Column`
# # #         - :class:`~.equipments.Machine`
# # #         - :class:`~.equipments.Electrode`
# # #         - :class:`~.equipments.Membrane`

# # <<<<<<< Updated upstream
# =======
# # __all__ = ('ElectrochemicalCell',
# #            'ESAPRecovery',
# #            'ESAPEffluent',
# #            'ESAP',
# #            'ElectrochemicalStrippingAdsorptionPrecipitation',)

# # class ElectrochemicalCell(SanUnit):
# #     '''
# #     Electrochemical cell for nutrient recovery.

# #     This unit has the following equipment:
# #         - :class:`~.equipments.Column`
# #         - :class:`~.equipments.Machine`
# #         - :class:`~.equipments.Electrode`
# #         - :class:`~.equipments.Membrane`

# #     Parameters
# #     ----------
# #     recovery : dict
# #         Keys refer to chemical component IDs. Values refer to recovery fractions (with 1 being 100%) for the respective chemicals.
# #     removal : dict
# #         Keys refer to chemical component IDs. Values refer to removal fractions (with 1 being 100%) for the respective chemicals.
# #     equipment : list(obj)
# #         List of Equipment objects part of the Electrochemical Cell.
# #     OPEX_over_CAPEX : float
# #         Ratio with which operating costs are calculated as a fraction of capital costs

# #     Example
# #     -------
# #     >>> # Set components
# #     >>> import qsdsan as qs
# #     >>> kwargs = dict(particle_size='Soluble',
# #     ...               degradability='Undegradable',
# #     ...               organic=False)
# #     >>> H2O = qs.Component.from_chemical('H2O', phase='l', **kwargs)
# #     >>> NH3 = qs.Component.from_chemical('NH3', phase='g', **kwargs)
# #     >>> NH3.particle_size = 'Dissolved gas'
# #     >>> H2SO4 = qs.Component.from_chemical('H2SO4', phase='l', **kwargs)
# #     >>> AmmoniumSulfate = qs.Component.from_chemical('AmmoniumSulfate', phase='l',
# #     ...                                              **kwargs)
# #     >>> Na_ion = qs.Component.from_chemical('Na+', phase='s', **kwargs)
# #     >>> K_ion = qs.Component.from_chemical('K+', phase='s', **kwargs)
# #     >>> Cl_ion = qs.Component.from_chemical('Cl-', phase='s', **kwargs)
# #     >>> PO4_ion = qs.Component.from_chemical('Phosphate', phase='s', **kwargs)
# #     >>> SO4_ion = qs.Component.from_chemical('Sulfate', phase='s', **kwargs)
# #     >>> C = qs.Component.from_chemical('Carbon', phase='s', **kwargs)
# #     >>> COD = qs.Component.from_chemical('O2', phase='s', **kwargs)
# #     >>> cmps = qs.Components((H2O, NH3, H2SO4, AmmoniumSulfate, Na_ion, K_ion, Cl_ion, PO4_ion, SO4_ion, C, COD))
# #     >>> # Assuming all has the same molar volume as water for demonstration purpose
# #     >>> for cmp in cmps:
# #     ...     cmp.copy_models_from(H2O, names=['V'])
# #     ...     defaulted_properties = cmp.default()
# #     >>> qs.set_thermo(cmps)
# #     >>> # Set waste streams
# #     >>> influent = qs.WasteStream('influent')
# #     >>> influent.set_flow_by_concentration(flow_tot=0.5, concentrations={'NH3':3820,
# #     ...                                     'Na+':1620, 'K+':1470, 'Cl-':3060, 'Phosphate':169,
# #     ...                                     'Sulfate':1680, 'Carbon':1860, 'O2':3460}, units=('mL/min', 'mg/L'))
# #     >>> influent.show()
# #     WasteStream: influent
# #      phase: 'l', T: 298.15 K, P: 101325 Pa
# #      flow (g/hr): H2O        29.5
# #                   NH3        0.115
# #                   Na+        0.0486
# #                   K+         0.0441
# #                   Cl-        0.0918
# #                   Phosphate  0.00507
# #                   Sulfate    0.0504
# #                   Carbon     0.0558
# #                   O2         0.104
# #      WasteStream-specific properties:
# #       pH         : 7.0
# #       Alkalinity : 2.5 mmol/L
# #       TC         : 1860.0 mg/L
# #       TP         : 31.9 mg/L
# #       TK         : 1470.0 mg/L
# #      Component concentrations (mg/L):
# #       H2O              984514.0
# #       NH3              3820.0
# #       Na+              1620.0
# #       K+               1470.0
# #       Cl-              3060.0
# #       Phosphate        169.0
# #       Sulfate          1680.0
# #       Carbon           1860.0
# #       O2               3460.0
# #     >>> catalysts = qs.WasteStream('catalysts', H2SO4=0.00054697, units=('kg/hr'))
# #     >>> # kg/hr imass value derived from 0.145 moles used on average for one, 26 hour experiment in the model cell
# #     >>> catalysts.price = 8.4 #in $/kg
# #     >>> # electricity price for ECS experiments in the default model tested by the Tarpeh Lab
# #     >>> # if want to change the price of electricity,
# #     >>> # use qs.PowerUtility, e.g., qs.PowerUtility.price = 0.1741
    
# #     >>> # Set the unit
# #     >>> U1 = qs.sanunits.ElectrochemicalCell('unit_1', ins=(influent, catalysts), outs=('recovered', 'removed', 'residual'))
# #     >>> # Simulate and look at the results
# #     >>> U1.simulate()
# #     >>> U1.results()
# #     Electrochemical cell                                      Units                              unit_1
# #     Electricity         Power                                    kW                            5.42e-05
# #                         Cost                                 USD/hr                            4.24e-06
# #     Design              Electrode unit_1_Main_Anode - Nu...    None                                   1
# #                         Electrode unit_1_Main_Anode - Ma...    None  titanium grid catalyst welded t...
# #                         Electrode unit_1_Main_Anode - Su...      m2                                   1
# #                         Electrode unit_1_Main_Cathode - ...    None                                   1
# #                         Electrode unit_1_Main_Cathode - ...    None  timesetl 3pcs stainless steel w...
# #                         Electrode unit_1_Main_Cathode - ...      m2                                30.2
# #                         Electrode unit_1_Current_Collect...    None                                   1
# #                         Electrode unit_1_Current_Collect...    None    stainless steel 26 gauge 5.5 x 7
# #                         Electrode unit_1_Current_Collect...      m2                                38.5
# #                         Electrode unit_1_Reference_Elect...    None                                   1
# #                         Electrode unit_1_Reference_Elect...    None  re-5b ag/agcl, 7.5 cm long, wit...
# #                         Electrode unit_1_Reference_Elect...      m2                                   1
# #                         Membrane unit_1_Cation_Exchange_...                                           1
# #                         Membrane unit_1_Cation_Exchange_...          CMI-7000S, polystyrene 0.45mm t...
# #                         Membrane unit_1_Cation_Exchange_...      m2                                30.2
# #                         Membrane unit_1_Gas_Permeable_Me...                                           1
# #                         Membrane unit_1_Gas_Permeable_Me...          Aquastill 0.3-micron polyethyle...
# #                         Membrane unit_1_Gas_Permeable_Me...      m2                                30.2
# #     Purchase cost       unit_1_Main_Anode                       USD                                 288
# #                         unit_1_Main_Cathode                     USD                               0.847
# #                         unit_1_Current_Collector_Cathode        USD                                39.9
# #                         unit_1_Reference_Electrode              USD                                  94
# #                         unit_1_Cation_Exchange_Membrane         USD                                 167
# #                         unit_1_Gas_Permeable_Membrane           USD                                29.3
# #                         Exterior                                USD                                60.5
# #     Total purchase cost                                         USD                                 679
# #     Utility cost                                             USD/hr                            4.24e-06
# #     Additional OPEX                                          USD/hr                                 136
# #     >>> U1.show() # doctest: +ELLIPSIS
# #     ElectrochemicalCell: unit_1
# #     ins...
# #     [0] influent
# #     phase: 'l', T: 298.15 K, P: 101325 Pa
# #     flow (g/hr): H2O        29.5
# #     ...
# #     '''

# #     _N_ins = 2
# #     _N_outs = 3


# #     def __init__(self, ID='', ins=(), outs=(),
# #                  recovery={'NH3':0.7}, removal={'NH3':0.83, 'K+':0.83, "Na+":0.8}, OPEX_over_CAPEX=0.2):
# #         SanUnit.__init__(self=self, ID=ID, ins=ins, outs=outs)
# #         self.recovery = recovery
# #         self.removal = removal
# #         self.OPEX_over_CAPEX = OPEX_over_CAPEX


# #         self.equipment = [
# #             Electrode('Main_Anode', linked_unit=self, N=1, electrode_type='anode',
# #                       material='Titanium grid catalyst welded to current collector tab both coated in iridium tantalum mixed metal oxide', surface_area=1, unit_cost=288), #288/unit, 1 unit
# #             Electrode('Main_Cathode', linked_unit=self, N=1, electrode_type='cathode',
# #                       material='TIMESETL 3pcs Stainless Steel Woven Wire 20 Mesh - 12"x8"(30x21cm) Metal Mesh Sheet 1mm Hole Great for Air Ventilation - A4', surface_area=30.25, unit_cost=0.847), #in in^2
# #             Electrode('Current_Collector_Cathode', linked_unit=self, N=1, electrode_type='cathode',
# #                       material='Stainless Steel 26 gauge 5.5'' x 7''', surface_area=38.5, unit_cost=39.9245), #in unknown units (94/unit, 1 unit)
# #             Electrode('Reference_Electrode', linked_unit=self, N=1, electrode_type='reference',
# #                       material='RE-5B Ag/AgCl, 7.5 cm long, with ceramic (MF-2056)', surface_area=1, unit_cost=94), #in unknown units (94/unit, 1 unit)
# #             Membrane('Cation_Exchange_Membrane', linked_unit=self, N=1,
# #                      material='CMI-7000S, polystyrene 0.45mm thick [48'' x 20'']',
# #                      unit_cost=5.5055, surface_area=30.25), # in in^2
# #             Membrane('Gas_Permeable_Membrane', linked_unit=self, N=1,
# #                      material='Aquastill 0.3-micron polyethylene membrane', unit_cost=0.968, surface_area=30.25), #in in^2
# #             ]


# #     def _run(self):
# #         influent, catalysts = self.ins
# #         recovered, removed, residual = self.outs[0], self.outs[1], self.outs[2]

# #         mixture = WasteStream()
# #         mixture.mix_from(self.ins)
# #         residual.copy_like(mixture)

# #         for chemical, recovery_ratio in self.recovery.items():
# #             recovered.imass[chemical] = mixture.imass[chemical]*recovery_ratio
            
# #         for chemical, removal_ratio in self.removal.items():
# #             recovery_ratio = 0 if recovered.imass[chemical] is None else recovered.imass[chemical]
# #             removed.imass[chemical] = mixture.imass[chemical]*removal_ratio - mixture.imass[chemical]*recovery_ratio
# #             residual.imass[chemical] = residual.imass[chemical]-residual.imass[chemical]*removal_ratio

# #     def _design(self):
# #         self.add_equipment_design()

# #     def _cost(self):
# #         self.add_equipment_cost()
# #         self.baseline_purchase_costs['Exterior'] = 60.52
# #         '''
# #         TOTAL CELL_EXTERIOR_COST = 60.52 USD
# #         Breakdown:
# #         Exterior frame (7'' x 7'' x 11/16'')	$21.46
# #         Interior half-cells (5.5'' x 5.5'' x 11/16'')	$19.87
# #         Rubber Sheets	$9.196
# #         Threaded Rods	$2.9696
# #         Wingnuts	$2.5504
# #         Flat Washers	$0.728
# #         Nylon Cable Glands	$3.744
# #         '''
# #         self.equip_costs = self.baseline_purchase_costs.values()
# #         add_OPEX = sum(self.equip_costs)*self.OPEX_over_CAPEX
# #         recovered, removed = self.outs[0], self.outs[1]

# #         self.power_utility.rate = recovered.imass['NH3']*0.67577
# #         # steady state value derived from 17.57 kWh used over 26 hrs
# #         self._add_OPEX = {'Additional OPEX': add_OPEX}
# # #%%
# # class ESAP(SanUnit):
# #     '''
# #     Electrochemical stripping, adsorption, and precipitation for nutrient recovery. 
# #     This unit is able to perform dynamic simulation.

# #     This unit has the following equipment:
# #         - :class:`~.equipments.Column`
# #         - :class:`~.equipments.Machine`
# #         - :class:`~.equipments.Electrode`
# #         - :class:`~.equipments.Membrane`

# >>>>>>> Stashed changes
# #     Parameters
# #     ----------
# #     ins:
# #         Wastewater stream
# #     outs:
# #         * [0] recovery product
# <<<<<<< Updated upstream
# #         * [1] loss stream (removed by the unit but cannot be recovered)
# =======
# #         * [1] Remainder stream   
# >>>>>>> Stashed changes
# #         * [2] effluent stream
# #     recovery : dict
# #         Keys refer to chemical component IDs. Values refer to recovery fractions (with 1 being 100%) for the respective chemicals.
# #     loss: dict
# #         Keys refer to chemical component IDs. Values refer to loss fractions (with 1 being 100%) for the respective chemicals.
# #     equipment : list(obj)
# #         List of Equipment objects part of the Electrochemical Cell.
# #     N_treatment_capacity: float
# #         kg N/m2/day, designed treatment capacity for a specific CEM area.
# #     OPEX_over_CAPEX : float
# #         Ratio with which operating costs are calculated as a fraction of capital costs
# #     component_ID_NH3: string
# #         The ID for ammonia/ammonium in the influent wastestream
# #     component_ID_P: string
# #         The ID for dissolved phosphate in the influent wastestream
# #     component_ID_Mg: string
# #         The ID for dissolved magnesium in the influent wastestream
# #     specific_energy_consumption: float
# #         kWh/kg N recovered
# <<<<<<< Updated upstream
# #     '''
    
# #     _N_outs = 3
# #     _N_ins = 1
# #     _ins_size_is_fixed = True
# #     _outs_size_is_fixed = True
# #     def __init__(self, ID='', ins=None, outs=(), thermo=None, OPEX_over_CAPEX=0.1, 
# #                  component_ID_NH3 = 'NH3', component_ID_P ='Phosphate', 
# #                  component_ID_Mg = "Mg2+", N_treatment_capacity = 0.362, 
# #                  specific_energy_consumption = 13.9, *, recovery={'NH3':0.7}, 
# #                  loss={'NH3':0.06}, order=None, init_with ='WasteStream', 
# #                  F_BM_default=None, isdynamic=False, **kwargs):
# #         SanUnit.__init__(self, ID=ID, ins=ins, outs=outs, thermo=thermo, 
# #                          init_with = init_with, isdynamic=isdynamic, 
# #                          F_BM_default=F_BM_default, 
# #                          OPEX_over_CAPEX = OPEX_over_CAPEX, 
# #                          # component_ID_NH3 = component_ID_NH3, 
# #                          # component_ID_P = component_ID_P,
# #                          # component_ID_Mg = component_ID_Mg, 
# #                          # N_treatment_capacity = N_treatment_capacity,
# #                          # specific_energy_consumption = specific_energy_consumption,
# #                          **kwargs)
# # =======
# # #     Parameters
# # #     ----------
# # #     ins:
# # #         Wastewater stream
# # #     outs:
# # #         * [0] recovery product
# # #         * [1] Remainder stream   
# # #         * [2] effluent stream
# # #     recovery : dict
# # #         Keys refer to chemical component IDs. Values refer to recovery fractions (with 1 being 100%) for the respective chemicals.
# # #     loss: dict
# # #         Keys refer to chemical component IDs. Values refer to loss fractions (with 1 being 100%) for the respective chemicals.
# # #     equipment : list(obj)
# # #         List of Equipment objects part of the Electrochemical Cell.
# # #     N_treatment_capacity: float
# # #         kg N/m2/day, designed treatment capacity for a specific CEM area.
# # #     OPEX_over_CAPEX : float
# # #         Ratio with which operating costs are calculated as a fraction of capital costs
# # #     component_ID_NH3: string
# # #         The ID for ammonia/ammonium in the influent wastestream
# # #     component_ID_P: string
# # #         The ID for dissolved phosphate in the influent wastestream
# # #     component_ID_Mg: string
# # #         The ID for dissolved magnesium in the influent wastestream
# # #     specific_energy_consumption: float
# # #         kWh/kg N recovered

# # #     '''
# # #     _N_outs = 3
# # #     _N_ins = 1
# # #     _ins_size_is_fixed = True
# # #     _outs_size_is_fixed = True
# # #     def __init__(self, ID='', ins=None, outs=(), thermo=None, OPEX_over_CAPEX=0.1, 
# # #                  component_ID_NH3 = 'NH3', component_ID_P ='Phosphate', 
# # #                  component_ID_Mg = "Mg2+", N_treatment_capacity = 0.362, 
# # #                  specific_energy_consumption = 13.9, *, recovery={'NH3':0.7}, 
# # #                  loss={'NH3':0.06}, order=None, init_with ='WasteStream', 
# # #                  F_BM_default=None, isdynamic=False, **kwargs):
# # #         SanUnit.__init__(self, ID=ID, ins=ins, outs=outs, thermo=thermo, 
# # #                          init_with = init_with, isdynamic=isdynamic, 
# # #                          F_BM_default=F_BM_default, 
# # #                          OPEX_over_CAPEX = OPEX_over_CAPEX, 
# # #                          # component_ID_NH3 = component_ID_NH3, 
# # #                          # component_ID_P = component_ID_P,
# # #                          # component_ID_Mg = component_ID_Mg, 
# # #                          # N_treatment_capacity = N_treatment_capacity,
# # #                          # specific_energy_consumption = specific_energy_consumption,
# # #                          **kwargs)
# # >>>>>>> Stashed changes
        
# # #         self.recovery = recovery
# # #         self.loss = loss
# # #         self._component_ID_NH3 = component_ID_NH3
# # #         self._component_ID_P = component_ID_P
# # #         self._component_ID_Mg = component_ID_Mg
# # #         self._N_treatment_capacity = N_treatment_capacity
# # #         self._specific_energy_consumption = specific_energy_consumption
# # #         self.equipment = [
# # #             Electrode('Anode', linked_unit=self, N=1, electrode_type='anode',
# # #                       material='Ti MMO mesh',surface_area=0.359, unit_cost=2576.19),
# # #             Electrode('Cathode', linked_unit=self, N=1, electrode_type='cathode',
# # #                       material='SS mesh', surface_area=0.359, unit_cost=46.67),
# # #             Membrane('Cation_Exchange_Membrane', linked_unit=self, N=1,
# # #                      material='Selemion CMVN', unit_cost=334.82, lifetime = 744, surface_area= 0.479), #$1016.67/m2
# # #             Membrane('Gas_Permeable_Membrane', linked_unit=self, N=1,
# # #                      material='Omniphobic', unit_cost=279.45, lifetime = 4464, surface_area=0.2394),
# # #             Machine('Pumps',linked_unit=self, N=3, lifetime= 43800, unit_cost =652.36), #pump lifetime 5 years
# # #             ]
    
# # #     @property
# # #     def component_ID_NH3(self):
# # #         '''string 
# # #             The ID for ammonia/ammonium in the influent wastestream '''
# # #         return self._component_ID_NH3
    
# # #     @property
# # #     def component_ID_P(self):
# # #         '''string 
# # #             The ID for dissolved phosphate in the influent wastestream '''
# # #         return self._component_ID_P
    
# # #     @property
# # #     def component_ID_Mg(self):
# # #         '''string 
# # #             The ID for dissolved magnesium in the influent wastestream '''
# # #         return self._component_ID_Mg
        
# # #     @property
# # #     def recovery(self):
# # #         '''dict, keys refer to chemical component IDs. 
# # #         Values refer to recovery fractions (with 1 being 100%) for the respective chemicals. '''
# # #         return self._recovery
    
# # #     @recovery.setter
# # #     def recovery(self, r):
# # #         self._recovery = r
# # #         influent, = self.ins
# # #         components = influent.components
# # #         self._recovery_matrix = components.kwarray(self.recovery)
# # #         self._remaining_matrix = None
        
# # #     @property
# # #     def loss(self):
# # #         '''dict, keys refer to chemical component IDs. 
# # #         Values refer to loss fractions (with 1 being 100%) for the respective chemicals. '''
# # #         return self._loss
    
# # #     @loss.setter
# # #     def loss(self, l):
# # #         self._loss = l
# # #         influent, = self.ins
# # #         components = influent.components
# # #         self._loss_matrix = components.kwarray(self.loss)
# # #         self._remaining_matrix = None
    
# # #     @property
# # #     def N_treatment_capacity(self):
# # #         '''kg N/m2/day, designed treatment capacity for a specific CEM area.'''
# # #         return self._N_treatment_capacity
    
# # #     @N_treatment_capacity.setter
# # #     def N_treatment_capacity(self, c):
# # #         self._N_treatment_capacity = c
    
# # #     @property
# # #     def specific_energy_consumption(self):
# # #         '''kWh/kg N recovered'''
# # #         return self._specific_energy_consumption
    
# # #     @specific_energy_consumption.setter
# # #     def specific_energy_consumption(self, e):
# # #         self._specific_energy_consumption = e
    
# # #     @property
# # #     def recovery_matrix(self):
# # #         '''The mass split [0-1] matrix for components from influent into the recovery_product WasteStream'''
# # #         return self._recovery_matrix
    
# # #     @property
# # #     def loss_matrix(self):
# # #         '''The mass split [0-1] matrix for components from influent into the loss WasteStream'''
# # #         return self._loss_matrix
    
# # #     @property
# # #     def remaining_matrix(self):
# # #         '''The mass split [0-1] matrix for components from influent into the effluent WasteStream'''
# # #         if not hasattr(self,'_remaining_matrix'):
# # #             self._remaining_matrix = None
# # #         if self._remaining_matrix is None:
# # #             influent, = self.ins
# # #             mass = influent.mass.copy()
# # #             components = influent.components
# # #             self._remaining_matrix = np.ones_like(mass) \
# # #                 - components.kwarray(self.recovery) - components.kwarray(self.loss)
# # #         return self._remaining_matrix        
    
# # #     def _run(self):
# # #         influent, = self.ins
# # #         recovery_product, loss, effluent = self.outs
# # #         mass = influent.mass.copy()
# # #         recovery_product.phase = loss.phase = 's'
# # #         # Calculate mass splits
# # #         recovery_product.mass = mass * self.recovery_matrix
# # #         loss.mass = mass * self.loss_matrix
# # #         effluent.mass = mass * self.remaining_matrix
        
# # #     @property
# # #     def state(self):
# # #         '''Component concentrations and total flow rate.'''
# # #         if self._state is None: return None
# # #         else:
# # #             return dict(zip(list(self.components.IDs)+['Q'], self._state))
        
# # #     def _init_state(self):
# # #         influent = self.ins[0]
# # #         # self._state = influent.mass.copy()
# # #         self._state = np.append(influent.mass.copy()*24*1e3, influent.F_vol*24) 
# # #         #convert kg/h to g/d; convert m3/hr to m3/day
# # #         self._dstate = self._state * 0.
        
# # #     def _update_state(self):
# # #         arr = self._state
# # #         for ws in self.outs:
# # #             if ws.state is None: ws.state = np.zeros_like(arr)
# # #         self._outs[0].state = np.zeros_like(arr)
# # #         self._outs[0].state[:-1] = self.recovery_matrix * arr[:-1] # g/d
# # #         self._outs[0].state[-1] = 1
        
# # #         self._outs[1].state[:-1] = self.loss_matrix * arr[:-1] # g/d
# # #         self._outs[1].state[-1] = 1

# # #         self._outs[2].state[:-1] = self.remaining_matrix * arr[:-1]/arr[-1] #mg/L
# # #         self._outs[2].state[-1] = arr[-1] #m3/day

# # #     def _update_dstate(self):
# # #         arr = self._dstate       
# # #         for ws in self.outs:
# # #             if ws.dstate is None: ws.dstate = np.zeros_like(arr)
# # #         self._outs[0].dstate[:-1] = self.recovery_matrix * arr[:-1]
# # #         self._outs[0].dstate[-1] = 0 ##constant derivative = 0
# # #         self._outs[1].dstate[:-1] = self.recovery_matrix * arr[:-1]
# # #         self._outs[1].dstate[-1] = 0 ##constant derivative = 0
# # #         self._outs[2].dstate[:-1] = (self.remaining_matrix * 
# # #                                      ((arr[:-1]*self._outs[2].state[-1]-arr[-1]*self._outs[2].state[:-1])/self._outs[2].state[-1]**2))
# # #         # self._outs[2].dstate[:-1] = (self.remaining_matrix * (arr[:-1]/self._outs[2].state[-1]))
# # #         self._outs[2].dstate[-1] = arr[-1]
    
# # #     @property
# # #     def AE(self):
# # #         if self._AE is None:
# # #             self._compile_AE()
# # #         return self._AE

# # #     def _compile_AE(self):
# # #         _state = self._state
# # #         _dstate = self._dstate
# # #         _update_state = self._update_state
# # #         _update_dstate = self._update_dstate
# # #         def yt(t, QC_ins, dQC_ins):
# # #             _state[-1] = QC_ins[0][-1]
# # #             _state[:-1] = QC_ins[0][:-1]*QC_ins[0][-1]
# # #             _dstate[-1] = dQC_ins[0][-1]
# # #             _dstate[:-1] = dQC_ins[0][:-1]*QC_ins[0][-1] + QC_ins[0][:-1]*dQC_ins[0][-1]
# # #             _update_state()
# # #             _update_dstate()
# # #         self._AE = yt
    
# # #     def _design(self):
# # #         self.add_equipment_design()
# # #         D = self.design_results
# # #         D['Number of ECS towers'] = ceil(self.outs[0].imass[self.component_ID_NH3]/17*14*24\
# # #                                          /self.equipment[2].surface_area/self.N_treatment_capacity)
# # #         D['H2SO4'] = self.outs[0].imol[self.component_ID_NH3] * 0.5 * 98/0.96 *24 #kg/day, assume 1 mol NH3 requires 0.5 mol H2SO4, kg/hr 96% H2SO4
# # #         D['NaOH'] = self.outs[0].imol[self.component_ID_Mg] * 2* 40/0.9 *24 #kg/day, assume 1 mol Mg requires 2 mol NaOH, kg/hr 90% NaOH
# # #         D['KOH'] = self.outs[2].imass['H2O'] * 0.03 * 0.04 *24 #kg/day
# # #         self._units['H2SO4'] = self._units['NaOH'] = self._units['KOH'] = 'kg/day'
# # #         influent = self.ins[0]
# # #         recovery_product, loss, effluent = self.outs
# # #         effluent.imass['H2O'] = influent.imass['H2O'] #assume no water tranport

# # #     def _cost(self):
# # #         self.add_equipment_cost()
# # #         C = self.baseline_purchase_costs
# # #         D = self.design_results
# # #         C['Flanges'] = (47.83+#https://www.mcmaster.com/4881K241 
# # #                         113.79+#https://www.mcmaster.com/95665K324
# # #                         37.5+#https://www.mcmaster.com/4881K239
# # #                         78.98+#https://www.mcmaster.com/95665K322
# # #                         19.15+#https://www.mcmaster.com/4881K236
# # #                         37.65+#https://www.mcmaster.com/95665K217
# # #                         61.53)*2#https://www.mcmaster.com/4881K967
        
# # #         C['TowerWall'] =  (548.55+ #https://www.mcmaster.com/4740K32
# # #                            287.17+#https://www.mcmaster.com/4740K31
# # #                            118.64)#https://www.mcmaster.com/4740K26
# # #         # self.equip_costs = self.baseline_purchase_costs
# # #         for equipment, value in self.baseline_purchase_costs.items():
# # #             self.baseline_purchase_costs[equipment] = value * D['Number of ECS towers']
# # #         add_OPEX = sum(self.baseline_purchase_costs.values())*self.OPEX_over_CAPEX/365/24
# # #         recovered = self.outs[0]

# # #         self.power_utility.rate = (recovered.imass[self.component_ID_NH3]/17*14
# # #                                    *self.specific_energy_consumption) #kW
# # #         # steady state value derived from 17.57 kWh used over 26 hrs
# # #         self._add_OPEX = {'Additional OPEX': add_OPEX}
    
    

# # # #%%
# # # class ESAPRecovery(Splitter):
# # #     '''
# # #     Electrochemical stripping, adsorption, and precipitation for nutrient recovery. 
# # #     This unit is able to perform dynamic simulation.

# # #     This unit has the following equipment:
# # #         - :class:`~.equipments.Column`
# # #         - :class:`~.equipments.Machine`
# # #         - :class:`~.equipments.Electrode`
# # #         - :class:`~.equipments.Membrane`

# # #     Parameters
# # #     ----------
# # #     ins:
# # #         Wastewater stream
# # #     outs:
# # #         * [0] recovery product
# # #         * [1] Remainder stream   
# # #     recovery : dict
# # #         Keys refer to chemical component IDs. Values refer to recovery fractions (with 1 being 100%) for the respective chemicals.
# # #     equipment : list(obj)
# # #         List of Equipment objects part of the Electrochemical Cell.
# # #     N_treatment_capacity: float
# # #         kg N/h designed treatment capacity for 1 tower.
# # #     OPEX_over_CAPEX : float
# # #         Ratio with which operating costs are calculated as a fraction of capital costs
# # #     component_ID_NH3: string
# # #         The ID for ammonia/ammonium in the influent wastestream
# # #     component_ID_P: string
# # #         The ID for dissolved phosphate in the influent wastestream
# # #     component_ID_Mg: string
# # #         The ID for dissolved magnesium in the influent wastestream
# # #     N_prodcut_concentration: float
# # #         molar concentration (mol/L) of N in the final product
# # #     '''
# # #     _ins_size_is_fixed = False
    
# # #     def __init__(self, ID='', ins=None, outs=(), thermo=None, *, 
# # #                  recovery={'NH3':0.8,'Mg2+':0.7,'Phosphate':0.95}, order=None, 
# # #                  init_with ='WasteStream', F_BM_default=None, isdynamic=False, 
# # #                  N_treatment_capacity=0.015, OPEX_over_CAPEX=0.2, component_ID_NH3 = 'NH3',
# # #                  component_ID_P ='Phosphate', component_ID_Mg = "Mg2+", N_prodcut_concentration = 2): #0.015 kg N/h, product 2 mol/L
# # #         Splitter.__init__(self=self, ID=ID, ins=ins, outs=outs,thermo = thermo, split = recovery,
# # #                           order = order, init_with = init_with, F_BM_default = F_BM_default,
# # #                           isdynamic = isdynamic
# # #                           )
# # #         self.recovery = recovery
# # #         self.OPEX_over_CAPEX = OPEX_over_CAPEX
# # #         self.component_ID_NH3 = component_ID_NH3
# # #         self.component_ID_Mg = component_ID_Mg
# # #         self.N_treatment_capacity = N_treatment_capacity

# # #         self.equipment = [
# # #             Electrode('Anode', linked_unit=self, N=1, electrode_type='anode',
# # #                       material='Ti MMO mesh',surface_area=0.359, unit_cost=2576.19),
# # #             Electrode('Cathode', linked_unit=self, N=1, electrode_type='cathode',
# # #                       material='SS mesh', surface_area=0.359, unit_cost=46.67),
# # #             Membrane('Cation_Exchange_Membrane', linked_unit=self, N=1,
# # #                      material='Selemion CMVN', unit_cost=1016.67, surface_area= 0.479), #$1016.67/m2
# # #             Membrane('Gas_Permeable_Membrane', linked_unit=self, N=1,
# # #                      material='Omniphobic', unit_cost=279.45, surface_area=0.2394),
# # #             Machine('Pumps',linked_unit=self, N=3, lifetime= 43800, unit_cost =652.36), #pump lifetime 5 years
# # #             ]

# # #     def _design(self):
# # #         self.add_equipment_design()
# # #         D = self.design_results
# # #         D['Number of ECS towers'] = ceil(self.outs[0].imass['NH3']/17*14/self.N_treatment_capacity)
# # #         D['H2SO4'] = self.outs[0].imol[self.component_ID_NH3] * 0.5 * 98/0.96 #assume 1 mol NH3 requires 0.5 mol H2SO4, kg/hr 96% H2SO4
# # #         D['NaOH'] = self.outs[0].imol[self.component_ID_Mg] * 2* 40/0.9 #assume 1 mol Mg requires 2 mol NaOH, kg/hr 90% NaOH

# # #     def _cost(self):
# # #         self.add_equipment_cost()
# # #         C = self.baseline_purchase_costs
# # #         D = self.design_results
# # #         C['Flanges'] = (47.83+#https://www.mcmaster.com/4881K241 
# # #                         113.79+#https://www.mcmaster.com/95665K324
# # #                         37.5+#https://www.mcmaster.com/4881K239
# # #                         78.98+#https://www.mcmaster.com/95665K322
# # #                         19.15+#https://www.mcmaster.com/4881K236
# # #                         37.65+#https://www.mcmaster.com/95665K217
# # #                         61.53)*2#https://www.mcmaster.com/4881K967
        
# # #         C['TowerWall'] =  (548.55+ #https://www.mcmaster.com/4740K32
# # #                            287.17+#https://www.mcmaster.com/4740K31
# # #                            118.64)#https://www.mcmaster.com/4740K26

# # #         self.equip_costs = self.baseline_purchase_costs.values() * D['Number of ECS towers']
# # #         add_OPEX = sum(self.equip_costs)*self.OPEX_over_CAPEX
# # #         recovered = self.outs[0]

# # #         self.power_utility.rate = recovered.imass['NH3']*0.67577
# # #         # steady state value derived from 17.57 kWh used over 26 hrs
# # #         self._add_OPEX = {'Additional OPEX': add_OPEX}

# # # #%%
# # # class ESAPEffluent(Splitter):
# # #     '''
# # #     Splitting the ESAP effluent. This unit is able to perform dynamic simulation.

# # #     This unit has the following equipment:
# # #         - :class:`~.equipments.Column`
# # #         - :class:`~.equipments.Machine`
# # #         - :class:`~.equipments.Electrode`
# # #         - :class:`~.equipments.Membrane`

# # #     Parameters
# # #     ----------
# # #     ins:
# # #         Inlet fluid to be split
# # #     outs:
# # #         * [0] loss stream during the process
# # #         * [1] effluent stream   
# # #     loss : dict
# # #         Keys refer to chemical component IDs. Values refer to loss fractions (with 1 being 100%) 
# # #         for the respective chemicals. The fraction is respective to the influent of ESAP_effluent unit.
# # #     equipment : list(obj)
# # #         List of Equipment objects part of the Electrochemical Cell.
# # #     OPEX_over_CAPEX : float
# # #         Ratio with which operating costs are calculated as a fraction of capital costs
# # #     '''
# # #     _ins_size_is_fixed = False
    
# # #     def __init__(self, ID='', ins=None, outs=(), thermo=None, *, 
# # #                  loss={'NH3':0.5,'Mg2+':0.3,'Phosphate':0.8}, order=None, 
# # #                  init_with ='WasteStream', F_BM_default=None, isdynamic=False):
# # #         Splitter.__init__(self=self, ID=ID, ins=ins, outs=outs,thermo = thermo, split = loss,
# # #                           order = order, init_with = init_with, F_BM_default = F_BM_default,
# # #                           isdynamic = isdynamic
# # #                           )
# # #         self.loss = loss
        
# # # #%%
# # # class ElectrochemicalStrippingAdsorptionPrecipitation(SanUnit):
# # #     '''
# # #     Electrochemical stripping, adsorption, and precipitation for nutrient recovery.

# # #     This unit has the following equipment:
# # #         - :class:`~.equipments.Column`
# # #         - :class:`~.equipments.Machine`
# # #         - :class:`~.equipments.Electrode`
# # #         - :class:`~.equipments.Membrane`

# # #     Parameters
# # #     ----------
# # #     recovery : dict
# # #         Keys refer to chemical component IDs. Values refer to recovery fractions (with 1 being 100%) for the respective chemicals.
# # #     removal : dict
# # #         Keys refer to chemical component IDs. Values refer to removal fractions (with 1 being 100%) for the respective chemicals.
# # #     equipment : list(obj)
# # #         List of Equipment objects part of the Electrochemical Cell.
# # #     OPEX_over_CAPEX : float
# # #         Ratio with which operating costs are calculated as a fraction of capital costs

# # #     '''
# # #     _ins_size_is_fixed = False
# # #     _N_outs = 3
    
# # #     def __init__(self, ID='', ins=(), outs=(),
# # #                  recovery={'NH3':0.7}, removal={'NH3':0.83, 'K':0.83, "NaCl":0.8}, OPEX_over_CAPEX=0.2):
# # #         SanUnit.__init__(self=self, ID=ID, ins=ins, outs=outs)
# # #         self.recovery = recovery
# # #         self.removal = removal
# # #         self.OPEX_over_CAPEX = OPEX_over_CAPEX


# # #         self.equipment = [
# # #             Electrode('Main_Anode', linked_unit=self, N=1, electrode_type='anode',
# # #                       material='Titanium grid catalyst welded to current collector tab both coated in iridium tantalum mixed metal oxide', surface_area=1, unit_cost=288), #288/unit, 1 unit
# # #             Electrode('Main_Cathode', linked_unit=self, N=1, electrode_type='cathode',
# # #                       material='TIMESETL 3pcs Stainless Steel Woven Wire 20 Mesh - 12"x8"(30x21cm) Metal Mesh Sheet 1mm Hole Great for Air Ventilation - A4', surface_area=30.25, unit_cost=0.847), #in in^2
# # #             Electrode('Current_Collector_Cathode', linked_unit=self, N=1, electrode_type='cathode',
# # #                       material='Stainless Steel 26 gauge 5.5'' x 7''', surface_area=38.5, unit_cost=39.9245), #in unknown units (94/unit, 1 unit)
# # #             Electrode('Reference_Electrode', linked_unit=self, N=1, electrode_type='reference',
# # #                       material='RE-5B Ag/AgCl, 7.5 cm long, with ceramic (MF-2056)', surface_area=1, unit_cost=94), #in unknown units (94/unit, 1 unit)
# # #             Membrane('Cation_Exchange_Membrane', linked_unit=self, N=1,
# # #                      material='CMI-7000S, polystyrene 0.45mm thick [48'' x 20'']',
# # #                      unit_cost=5.5055, surface_area=30.25), # in in^2
# # #             Membrane('Gas_Permeable_Membrane', linked_unit=self, N=1,
# # #                      material='Aquastill 0.3-micron polyethylene membrane', unit_cost=0.968, surface_area=30.25), #in in^2
# # #             ]

# # #     def _run(self):
# # #         influent, = self.ins
# # #         recovered_product, process_loss, effluent = self.outs[0], self.outs[1], self.outs[2]
# # #         effluent.copy_like(self.ins[0])

# # #         for chemical, recovery_ratio in self.recovery.items():
# # #             recovered_product.imass[chemical] = influent.imass[chemical]*recovery_ratio
            
# # #         for chemical, removal_ratio in self.removal.items():
# # #             recovery_ratio = 0 if recovered_product.imass[chemical] is None else recovered_product.imass[chemical]/influent.imass[chemical]
# # #             process_loss.imass[chemical] = influent.imass[chemical]*removal_ratio - influent.imass[chemical]*recovery_ratio
# # #             effluent.imass[chemical] = effluent.imass[chemical]-effluent.imass[chemical]*removal_ratio

# # #     def _design(self):
# # #         self.add_equipment_design()
# # #         D = self.design_results
# # #         D['H2SO4'] = self.outs[0].imol['NH3'] * 98/0.96 #assume 1 mol NH3 requires 1 mol H2SO4, kg/hr 96% H2SO4
# # #         D['NaOH'] = self.outs[0].imol['Mg2+'] * 2* 40/0.9 #assume 1 mol Mg requires 2 mol NaOH, kg/hr 90% NaOH

# # #     def _cost(self):
# # #         self.add_equipment_cost()
# # #         self.baseline_purchase_costs['Exterior'] = 60.52
# # #         '''
# # #         TOTAL CELL_EXTERIOR_COST = 60.52 USD
# # #         Breakdown:
# # #         Exterior frame (7'' x 7'' x 11/16'')	$21.46
# # #         Interior half-cells (5.5'' x 5.5'' x 11/16'')	$19.87
# # #         Rubber Sheets	$9.196
# # #         Threaded Rods	$2.9696
# # #         Wingnuts	$2.5504
# # #         Flat Washers	$0.728
# # #         Nylon Cable Glands	$3.744
# # #         '''
# # #         self.equip_costs = self.baseline_purchase_costs.values()
# # #         add_OPEX = sum(self.equip_costs)*self.OPEX_over_CAPEX
# # #         recovered, removed = self.outs[0], self.outs[1]

# # <<<<<<< Updated upstream
# #         self.power_utility.rate = recovered.imass['NH3']*0.67577
# #         # steady state value derived from 17.57 kWh used over 26 hrs
# #         self._add_OPEX = {'Additional OPEX': add_OPEX}
# # =======
# # #         self.power_utility.rate = recovered.imass['NH3']*0.67577
# # #         # steady state value derived from 17.57 kWh used over 26 hrs
# # #         self._add_OPEX = {'Additional OPEX': add_OPEX}
# # >>>>>>> Stashed changes
# =======

# #     '''
# #     _N_outs = 3
# #     _N_ins = 1
# #     _ins_size_is_fixed = True
# #     _outs_size_is_fixed = True
# #     def __init__(self, ID='', ins=None, outs=(), thermo=None, OPEX_over_CAPEX=0.1, 
# #                  component_ID_NH3 = 'NH3', component_ID_P ='Phosphate', 
# #                  component_ID_Mg = "Mg2+", N_treatment_capacity = 0.362, 
# #                  specific_energy_consumption = 13.9, *, recovery={'NH3':0.7}, 
# #                  loss={'NH3':0.06}, order=None, init_with ='WasteStream', 
# #                  F_BM_default=None, isdynamic=False, **kwargs):
# #         SanUnit.__init__(self, ID=ID, ins=ins, outs=outs, thermo=thermo, 
# #                          init_with = init_with, isdynamic=isdynamic, 
# #                          F_BM_default=F_BM_default, 
# #                          OPEX_over_CAPEX = OPEX_over_CAPEX, 
# #                          # component_ID_NH3 = component_ID_NH3, 
# #                          # component_ID_P = component_ID_P,
# #                          # component_ID_Mg = component_ID_Mg, 
# #                          # N_treatment_capacity = N_treatment_capacity,
# #                          # specific_energy_consumption = specific_energy_consumption,
# #                          **kwargs)
        
# #         self.recovery = recovery
# #         self.loss = loss
# #         self._component_ID_NH3 = component_ID_NH3
# #         self._component_ID_P = component_ID_P
# #         self._component_ID_Mg = component_ID_Mg
# #         self._N_treatment_capacity = N_treatment_capacity
# #         self._specific_energy_consumption = specific_energy_consumption
# #         self.equipment = [
# #             Electrode('Anode', linked_unit=self, N=1, electrode_type='anode',
# #                       material='Ti MMO mesh',surface_area=0.359, unit_cost=2576.19),
# #             Electrode('Cathode', linked_unit=self, N=1, electrode_type='cathode',
# #                       material='SS mesh', surface_area=0.359, unit_cost=46.67),
# #             Membrane('Cation_Exchange_Membrane', linked_unit=self, N=1,
# #                      material='Selemion CMVN', unit_cost=334.82, lifetime = 744, surface_area= 0.479), #$1016.67/m2
# #             Membrane('Gas_Permeable_Membrane', linked_unit=self, N=1,
# #                      material='Omniphobic', unit_cost=279.45, lifetime = 4464, surface_area=0.2394),
# #             Machine('Pumps',linked_unit=self, N=3, lifetime= 43800, unit_cost =652.36), #pump lifetime 5 years
# #             ]
    
# #     @property
# #     def component_ID_NH3(self):
# #         '''string 
# #             The ID for ammonia/ammonium in the influent wastestream '''
# #         return self._component_ID_NH3
    
# #     @property
# #     def component_ID_P(self):
# #         '''string 
# #             The ID for dissolved phosphate in the influent wastestream '''
# #         return self._component_ID_P
    
# #     @property
# #     def component_ID_Mg(self):
# #         '''string 
# #             The ID for dissolved magnesium in the influent wastestream '''
# #         return self._component_ID_Mg
        
# #     @property
# #     def recovery(self):
# #         '''dict, keys refer to chemical component IDs. 
# #         Values refer to recovery fractions (with 1 being 100%) for the respective chemicals. '''
# #         return self._recovery
    
# #     @recovery.setter
# #     def recovery(self, r):
# #         self._recovery = r
# #         influent, = self.ins
# #         components = influent.components
# #         self._recovery_matrix = components.kwarray(self.recovery)
# #         self._remaining_matrix = None
        
# #     @property
# #     def loss(self):
# #         '''dict, keys refer to chemical component IDs. 
# #         Values refer to loss fractions (with 1 being 100%) for the respective chemicals. '''
# #         return self._loss
    
# #     @loss.setter
# #     def loss(self, l):
# #         self._loss = l
# #         influent, = self.ins
# #         components = influent.components
# #         self._loss_matrix = components.kwarray(self.loss)
# #         self._remaining_matrix = None
    
# #     @property
# #     def N_treatment_capacity(self):
# #         '''kg N/m2/day, designed treatment capacity for a specific CEM area.'''
# #         return self._N_treatment_capacity
    
# #     @N_treatment_capacity.setter
# #     def N_treatment_capacity(self, c):
# #         self._N_treatment_capacity = c
    
# #     @property
# #     def specific_energy_consumption(self):
# #         '''kWh/kg N recovered'''
# #         return self._specific_energy_consumption
    
# #     @specific_energy_consumption.setter
# #     def specific_energy_consumption(self, e):
# #         self._specific_energy_consumption = e
    
# #     @property
# #     def recovery_matrix(self):
# #         '''The mass split [0-1] matrix for components from influent into the recovery_product WasteStream'''
# #         return self._recovery_matrix
    
# #     @property
# #     def loss_matrix(self):
# #         '''The mass split [0-1] matrix for components from influent into the loss WasteStream'''
# #         return self._loss_matrix
    
# #     @property
# #     def remaining_matrix(self):
# #         '''The mass split [0-1] matrix for components from influent into the effluent WasteStream'''
# #         if not hasattr(self,'_remaining_matrix'):
# #             self._remaining_matrix = None
# #         if self._remaining_matrix is None:
# #             influent, = self.ins
# #             mass = influent.mass.copy()
# #             components = influent.components
# #             self._remaining_matrix = np.ones_like(mass) \
# #                 - components.kwarray(self.recovery) - components.kwarray(self.loss)
# #         return self._remaining_matrix        
    
# #     def _run(self):
# #         influent, = self.ins
# #         recovery_product, loss, effluent = self.outs
# #         mass = influent.mass.copy()
# #         recovery_product.phase = loss.phase = 's'
# #         # Calculate mass splits
# #         recovery_product.mass = mass * self.recovery_matrix
# #         loss.mass = mass * self.loss_matrix
# #         effluent.mass = mass * self.remaining_matrix
        
# #     @property
# #     def state(self):
# #         '''Component concentrations and total flow rate.'''
# #         if self._state is None: return None
# #         else:
# #             return dict(zip(list(self.components.IDs)+['Q'], self._state))
        
# #     def _init_state(self):
# #         influent = self.ins[0]
# #         # self._state = influent.mass.copy()
# #         self._state = np.append(influent.mass.copy()*24*1e3, influent.F_vol*24) 
# #         #convert kg/h to g/d; convert m3/hr to m3/day
# #         self._dstate = self._state * 0.
        
# #     def _update_state(self):
# #         arr = self._state
# #         for ws in self.outs:
# #             if ws.state is None: ws.state = np.zeros_like(arr)
# #         self._outs[0].state = np.zeros_like(arr)
# #         self._outs[0].state[:-1] = self.recovery_matrix * arr[:-1] # g/d
# #         self._outs[0].state[-1] = 1
        
# #         self._outs[1].state[:-1] = self.loss_matrix * arr[:-1] # g/d
# #         self._outs[1].state[-1] = 1

# #         self._outs[2].state[:-1] = self.remaining_matrix * arr[:-1]/arr[-1] #mg/L
# #         self._outs[2].state[-1] = arr[-1] #m3/day

# #     def _update_dstate(self):
# #         arr = self._dstate       
# #         for ws in self.outs:
# #             if ws.dstate is None: ws.dstate = np.zeros_like(arr)
# #         self._outs[0].dstate[:-1] = self.recovery_matrix * arr[:-1]
# #         self._outs[0].dstate[-1] = 0 ##constant derivative = 0
# #         self._outs[1].dstate[:-1] = self.recovery_matrix * arr[:-1]
# #         self._outs[1].dstate[-1] = 0 ##constant derivative = 0
# #         self._outs[2].dstate[:-1] = (self.remaining_matrix * 
# #                                      ((arr[:-1]*self._outs[2].state[-1]-arr[-1]*self._outs[2].state[:-1])/self._outs[2].state[-1]**2))
# #         # self._outs[2].dstate[:-1] = (self.remaining_matrix * (arr[:-1]/self._outs[2].state[-1]))
# #         self._outs[2].dstate[-1] = arr[-1]
    
# #     @property
# #     def AE(self):
# #         if self._AE is None:
# #             self._compile_AE()
# #         return self._AE

# #     def _compile_AE(self):
# #         _state = self._state
# #         _dstate = self._dstate
# #         _update_state = self._update_state
# #         _update_dstate = self._update_dstate
# #         def yt(t, QC_ins, dQC_ins):
# #             _state[-1] = QC_ins[0][-1]
# #             _state[:-1] = QC_ins[0][:-1]*QC_ins[0][-1]
# #             _dstate[-1] = dQC_ins[0][-1]
# #             _dstate[:-1] = dQC_ins[0][:-1]*QC_ins[0][-1] + QC_ins[0][:-1]*dQC_ins[0][-1]
# #             _update_state()
# #             _update_dstate()
# #         self._AE = yt
    
# #     def _design(self):
# #         self.add_equipment_design()
# #         D = self.design_results
# #         D['Number of ECS towers'] = ceil(self.outs[0].imass[self.component_ID_NH3]/17*14*24\
# #                                          /self.equipment[2].surface_area/self.N_treatment_capacity)
# #         D['H2SO4'] = self.outs[0].imol[self.component_ID_NH3] * 0.5 * 98/0.96 *24 #kg/day, assume 1 mol NH3 requires 0.5 mol H2SO4, kg/hr 96% H2SO4
# #         D['NaOH'] = self.outs[0].imol[self.component_ID_Mg] * 2* 40/0.9 *24 #kg/day, assume 1 mol Mg requires 2 mol NaOH, kg/hr 90% NaOH
# #         D['KOH'] = self.outs[2].imass['H2O'] * 0.03 * 0.04 *24 #kg/day
# #         self._units['H2SO4'] = self._units['NaOH'] = self._units['KOH'] = 'kg/day'
# #         influent = self.ins[0]
# #         recovery_product, loss, effluent = self.outs
# #         effluent.imass['H2O'] = influent.imass['H2O'] #assume no water tranport

# #     def _cost(self):
# #         self.add_equipment_cost()
# #         C = self.baseline_purchase_costs
# #         D = self.design_results
# #         C['Flanges'] = (47.83+#https://www.mcmaster.com/4881K241 
# #                         113.79+#https://www.mcmaster.com/95665K324
# #                         37.5+#https://www.mcmaster.com/4881K239
# #                         78.98+#https://www.mcmaster.com/95665K322
# #                         19.15+#https://www.mcmaster.com/4881K236
# #                         37.65+#https://www.mcmaster.com/95665K217
# #                         61.53)*2#https://www.mcmaster.com/4881K967
        
# #         C['TowerWall'] =  (548.55+ #https://www.mcmaster.com/4740K32
# #                            287.17+#https://www.mcmaster.com/4740K31
# #                            118.64)#https://www.mcmaster.com/4740K26
# #         # self.equip_costs = self.baseline_purchase_costs
# #         for equipment, value in self.baseline_purchase_costs.items():
# #             self.baseline_purchase_costs[equipment] = value * D['Number of ECS towers']
# #         add_OPEX = sum(self.baseline_purchase_costs.values())*self.OPEX_over_CAPEX/365/24
# #         recovered = self.outs[0]

# #         self.power_utility.rate = (recovered.imass[self.component_ID_NH3]/17*14
# #                                    *self.specific_energy_consumption) #kW
# #         # steady state value derived from 17.57 kWh used over 26 hrs
# #         self._add_OPEX = {'Additional OPEX': add_OPEX}
    
    

# # #%%
# # class ESAPRecovery(Splitter):
# #     '''
# #     Electrochemical stripping, adsorption, and precipitation for nutrient recovery. 
# #     This unit is able to perform dynamic simulation.

# #     This unit has the following equipment:
# #         - :class:`~.equipments.Column`
# #         - :class:`~.equipments.Machine`
# #         - :class:`~.equipments.Electrode`
# #         - :class:`~.equipments.Membrane`

# #     Parameters
# #     ----------
# #     ins:
# #         Wastewater stream
# #     outs:
# #         * [0] recovery product
# #         * [1] Remainder stream   
# #     recovery : dict
# #         Keys refer to chemical component IDs. Values refer to recovery fractions (with 1 being 100%) for the respective chemicals.
# #     equipment : list(obj)
# #         List of Equipment objects part of the Electrochemical Cell.
# #     N_treatment_capacity: float
# #         kg N/h designed treatment capacity for 1 tower.
# #     OPEX_over_CAPEX : float
# #         Ratio with which operating costs are calculated as a fraction of capital costs
# #     component_ID_NH3: string
# #         The ID for ammonia/ammonium in the influent wastestream
# #     component_ID_P: string
# #         The ID for dissolved phosphate in the influent wastestream
# #     component_ID_Mg: string
# #         The ID for dissolved magnesium in the influent wastestream
# #     N_prodcut_concentration: float
# #         molar concentration (mol/L) of N in the final product
# #     '''
# #     _ins_size_is_fixed = False
    
# #     def __init__(self, ID='', ins=None, outs=(), thermo=None, *, 
# #                  recovery={'NH3':0.8,'Mg2+':0.7,'Phosphate':0.95}, order=None, 
# #                  init_with ='WasteStream', F_BM_default=None, isdynamic=False, 
# #                  N_treatment_capacity=0.015, OPEX_over_CAPEX=0.2, component_ID_NH3 = 'NH3',
# #                  component_ID_P ='Phosphate', component_ID_Mg = "Mg2+", N_prodcut_concentration = 2): #0.015 kg N/h, product 2 mol/L
# #         Splitter.__init__(self=self, ID=ID, ins=ins, outs=outs,thermo = thermo, split = recovery,
# #                           order = order, init_with = init_with, F_BM_default = F_BM_default,
# #                           isdynamic = isdynamic
# #                           )
# #         self.recovery = recovery
# #         self.OPEX_over_CAPEX = OPEX_over_CAPEX
# #         self.component_ID_NH3 = component_ID_NH3
# #         self.component_ID_Mg = component_ID_Mg
# #         self.N_treatment_capacity = N_treatment_capacity

# #         self.equipment = [
# #             Electrode('Anode', linked_unit=self, N=1, electrode_type='anode',
# #                       material='Ti MMO mesh',surface_area=0.359, unit_cost=2576.19),
# #             Electrode('Cathode', linked_unit=self, N=1, electrode_type='cathode',
# #                       material='SS mesh', surface_area=0.359, unit_cost=46.67),
# #             Membrane('Cation_Exchange_Membrane', linked_unit=self, N=1,
# #                      material='Selemion CMVN', unit_cost=1016.67, surface_area= 0.479), #$1016.67/m2
# #             Membrane('Gas_Permeable_Membrane', linked_unit=self, N=1,
# #                      material='Omniphobic', unit_cost=279.45, surface_area=0.2394),
# #             Machine('Pumps',linked_unit=self, N=3, lifetime= 43800, unit_cost =652.36), #pump lifetime 5 years
# #             ]

# #     def _design(self):
# #         self.add_equipment_design()
# #         D = self.design_results
# #         D['Number of ECS towers'] = ceil(self.outs[0].imass['NH3']/17*14/self.N_treatment_capacity)
# #         D['H2SO4'] = self.outs[0].imol[self.component_ID_NH3] * 0.5 * 98/0.96 #assume 1 mol NH3 requires 0.5 mol H2SO4, kg/hr 96% H2SO4
# #         D['NaOH'] = self.outs[0].imol[self.component_ID_Mg] * 2* 40/0.9 #assume 1 mol Mg requires 2 mol NaOH, kg/hr 90% NaOH

# #     def _cost(self):
# #         self.add_equipment_cost()
# #         C = self.baseline_purchase_costs
# #         D = self.design_results
# #         C['Flanges'] = (47.83+#https://www.mcmaster.com/4881K241 
# #                         113.79+#https://www.mcmaster.com/95665K324
# #                         37.5+#https://www.mcmaster.com/4881K239
# #                         78.98+#https://www.mcmaster.com/95665K322
# #                         19.15+#https://www.mcmaster.com/4881K236
# #                         37.65+#https://www.mcmaster.com/95665K217
# #                         61.53)*2#https://www.mcmaster.com/4881K967
        
# #         C['TowerWall'] =  (548.55+ #https://www.mcmaster.com/4740K32
# #                            287.17+#https://www.mcmaster.com/4740K31
# #                            118.64)#https://www.mcmaster.com/4740K26

# #         self.equip_costs = self.baseline_purchase_costs.values() * D['Number of ECS towers']
# #         add_OPEX = sum(self.equip_costs)*self.OPEX_over_CAPEX
# #         recovered = self.outs[0]

# #         self.power_utility.rate = recovered.imass['NH3']*0.67577
# #         # steady state value derived from 17.57 kWh used over 26 hrs
# #         self._add_OPEX = {'Additional OPEX': add_OPEX}

# # #%%
# # class ESAPEffluent(Splitter):
# #     '''
# #     Splitting the ESAP effluent. This unit is able to perform dynamic simulation.

# #     This unit has the following equipment:
# #         - :class:`~.equipments.Column`
# #         - :class:`~.equipments.Machine`
# #         - :class:`~.equipments.Electrode`
# #         - :class:`~.equipments.Membrane`

# #     Parameters
# #     ----------
# #     ins:
# #         Inlet fluid to be split
# #     outs:
# #         * [0] loss stream during the process
# #         * [1] effluent stream   
# #     loss : dict
# #         Keys refer to chemical component IDs. Values refer to loss fractions (with 1 being 100%) 
# #         for the respective chemicals. The fraction is respective to the influent of ESAP_effluent unit.
# #     equipment : list(obj)
# #         List of Equipment objects part of the Electrochemical Cell.
# #     OPEX_over_CAPEX : float
# #         Ratio with which operating costs are calculated as a fraction of capital costs
# #     '''
# #     _ins_size_is_fixed = False
    
# #     def __init__(self, ID='', ins=None, outs=(), thermo=None, *, 
# #                  loss={'NH3':0.5,'Mg2+':0.3,'Phosphate':0.8}, order=None, 
# #                  init_with ='WasteStream', F_BM_default=None, isdynamic=False):
# #         Splitter.__init__(self=self, ID=ID, ins=ins, outs=outs,thermo = thermo, split = loss,
# #                           order = order, init_with = init_with, F_BM_default = F_BM_default,
# #                           isdynamic = isdynamic
# #                           )
# #         self.loss = loss
        
# # #%%
# # class ElectrochemicalStrippingAdsorptionPrecipitation(SanUnit):
# #     '''
# #     Electrochemical stripping, adsorption, and precipitation for nutrient recovery.

# #     This unit has the following equipment:
# #         - :class:`~.equipments.Column`
# #         - :class:`~.equipments.Machine`
# #         - :class:`~.equipments.Electrode`
# #         - :class:`~.equipments.Membrane`

# #     Parameters
# #     ----------
# #     recovery : dict
# #         Keys refer to chemical component IDs. Values refer to recovery fractions (with 1 being 100%) for the respective chemicals.
# #     removal : dict
# #         Keys refer to chemical component IDs. Values refer to removal fractions (with 1 being 100%) for the respective chemicals.
# #     equipment : list(obj)
# #         List of Equipment objects part of the Electrochemical Cell.
# #     OPEX_over_CAPEX : float
# #         Ratio with which operating costs are calculated as a fraction of capital costs

# #     '''
# #     _ins_size_is_fixed = False
# #     _N_outs = 3
    
# #     def __init__(self, ID='', ins=(), outs=(),
# #                  recovery={'NH3':0.7}, removal={'NH3':0.83, 'K':0.83, "NaCl":0.8}, OPEX_over_CAPEX=0.2):
# #         SanUnit.__init__(self=self, ID=ID, ins=ins, outs=outs)
# #         self.recovery = recovery
# #         self.removal = removal
# #         self.OPEX_over_CAPEX = OPEX_over_CAPEX


# #         self.equipment = [
# #             Electrode('Main_Anode', linked_unit=self, N=1, electrode_type='anode',
# #                       material='Titanium grid catalyst welded to current collector tab both coated in iridium tantalum mixed metal oxide', surface_area=1, unit_cost=288), #288/unit, 1 unit
# #             Electrode('Main_Cathode', linked_unit=self, N=1, electrode_type='cathode',
# #                       material='TIMESETL 3pcs Stainless Steel Woven Wire 20 Mesh - 12"x8"(30x21cm) Metal Mesh Sheet 1mm Hole Great for Air Ventilation - A4', surface_area=30.25, unit_cost=0.847), #in in^2
# #             Electrode('Current_Collector_Cathode', linked_unit=self, N=1, electrode_type='cathode',
# #                       material='Stainless Steel 26 gauge 5.5'' x 7''', surface_area=38.5, unit_cost=39.9245), #in unknown units (94/unit, 1 unit)
# #             Electrode('Reference_Electrode', linked_unit=self, N=1, electrode_type='reference',
# #                       material='RE-5B Ag/AgCl, 7.5 cm long, with ceramic (MF-2056)', surface_area=1, unit_cost=94), #in unknown units (94/unit, 1 unit)
# #             Membrane('Cation_Exchange_Membrane', linked_unit=self, N=1,
# #                      material='CMI-7000S, polystyrene 0.45mm thick [48'' x 20'']',
# #                      unit_cost=5.5055, surface_area=30.25), # in in^2
# #             Membrane('Gas_Permeable_Membrane', linked_unit=self, N=1,
# #                      material='Aquastill 0.3-micron polyethylene membrane', unit_cost=0.968, surface_area=30.25), #in in^2
# #             ]

# #     def _run(self):
# #         influent, = self.ins
# #         recovered_product, process_loss, effluent = self.outs[0], self.outs[1], self.outs[2]
# #         effluent.copy_like(self.ins[0])

# #         for chemical, recovery_ratio in self.recovery.items():
# #             recovered_product.imass[chemical] = influent.imass[chemical]*recovery_ratio
            
# #         for chemical, removal_ratio in self.removal.items():
# #             recovery_ratio = 0 if recovered_product.imass[chemical] is None else recovered_product.imass[chemical]/influent.imass[chemical]
# #             process_loss.imass[chemical] = influent.imass[chemical]*removal_ratio - influent.imass[chemical]*recovery_ratio
# #             effluent.imass[chemical] = effluent.imass[chemical]-effluent.imass[chemical]*removal_ratio

# #     def _design(self):
# #         self.add_equipment_design()
# #         D = self.design_results
# #         D['H2SO4'] = self.outs[0].imol['NH3'] * 98/0.96 #assume 1 mol NH3 requires 1 mol H2SO4, kg/hr 96% H2SO4
# #         D['NaOH'] = self.outs[0].imol['Mg2+'] * 2* 40/0.9 #assume 1 mol Mg requires 2 mol NaOH, kg/hr 90% NaOH

# #     def _cost(self):
# #         self.add_equipment_cost()
# #         self.baseline_purchase_costs['Exterior'] = 60.52
# #         '''
# #         TOTAL CELL_EXTERIOR_COST = 60.52 USD
# #         Breakdown:
# #         Exterior frame (7'' x 7'' x 11/16'')	$21.46
# #         Interior half-cells (5.5'' x 5.5'' x 11/16'')	$19.87
# #         Rubber Sheets	$9.196
# #         Threaded Rods	$2.9696
# #         Wingnuts	$2.5504
# #         Flat Washers	$0.728
# #         Nylon Cable Glands	$3.744
# #         '''
# #         self.equip_costs = self.baseline_purchase_costs.values()
# #         add_OPEX = sum(self.equip_costs)*self.OPEX_over_CAPEX
# #         recovered, removed = self.outs[0], self.outs[1]

# #         self.power_utility.rate = recovered.imass['NH3']*0.67577
# #         # steady state value derived from 17.57 kWh used over 26 hrs
# #         self._add_OPEX = {'Additional OPEX': add_OPEX}
# >>>>>>> Stashed changes
