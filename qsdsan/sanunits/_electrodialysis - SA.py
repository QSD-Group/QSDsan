#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Junhyung Park <junhyungparkenv@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.

'''
from thermosteam.utils import chemicals_user
from thermosteam import settings
from chemicals.elements import molecular_weight as get_mw
from qsdsan import Component, Components, WasteStream, SanUnit, Process, Processes, CompiledProcesses, System
import numpy as np
from qsdsan.utils import ospath, data_path
from scipy.optimize import brenth
from warnings import warn
from math import log10
from qsdsan import Component, Components
from thermosteam import settings
from chemicals.elements import molecular_weight as get_mw
from qsdsan.utils import ospath, data_path

__all__ = ('create_ed_vfa_cmps', 'ED_vfa'
           )

#_path = ospath.join(data_path, 'process_data/_adm1_vfa.tsv')
#_load_components = settings.get_default_chemicals

#%%
# =============================================================================
# ED_vfa-specific components
# =============================================================================

# Define molecular weights for carbon (C) and nitrogen (N)
C_mw = get_mw({'C': 1})
N_mw = get_mw({'N': 1})
Na_mw = get_mw({'Na': 1})
Cl_mw = get_mw({'Cl': 1})
def create_ed_vfa_cmps(set_thermo=True):
    cmps_all = Components.load_default()
    # Acetate (C2)
    S_ac = Component.from_chemical('S_ac', chemical='acetic acid',
                                   description='Acetate',
                                   measured_as='COD',
                                   particle_size='Soluble',
                                   degradability='Readily',
                                   organic=True)
    
    # Ethanol (C2)
    S_et = Component.from_chemical('S_et', chemical='Ethanol',
                                    description='Ethanol',
                                    measured_as='COD',
                                    particle_size='Soluble',
                                    degradability='Readily',
                                    organic=True)
    
    # Lactate (C3)
    S_la = Component.from_chemical('S_la', chemical='lactic acid',
                                description='Lactate',
                                measured_as='COD',
                                particle_size='Soluble',
                                degradability='Readily',
                                organic=True)

    # Propionate (C3)
    S_pro = Component.from_chemical('S_pro', chemical='propionic acid',
                                    description='Propionate',
                                    measured_as='COD',
                                    particle_size='Soluble',
                                    degradability='Readily',
                                    organic=True)

    # Butyrate (C4)
    S_bu = Component.from_chemical('S_bu', chemical='butyric acid',
                                   description='Butyrate',
                                   measured_as='COD',
                                   particle_size='Soluble',
                                   degradability='Readily',
                                   organic=True)

    # Valerate (C5)
    S_va = Component.from_chemical('S_va', chemical='valeric acid',
                                   description='Valerate',
                                   measured_as='COD',
                                   particle_size='Soluble',
                                   degradability='Readily',
                                   organic=True)

    # Glucose (C6)
    S_su = Component.from_chemical('S_su', chemical='glucose',
                                   description='Monosaccharides',
                                   measured_as='COD',
                                   particle_size='Soluble',
                                   degradability='Readily',
                                   organic=True)
    
    # Hexanoate (C6)
    S_he = Component.from_chemical('S_he', chemical='caproic acid',
                                   description='Hexanoate',
                                   measured_as='COD',
                                   particle_size='Soluble',
                                   degradability='Readily',
                                   organic=True)
    
    # How I define molecular mass and unit?
    Na = Component('Na+', formula='Na', i_charge = 1.0,
                       particle_size='Soluble', degradability='Undegradable',
                       organic=False)

    Cl = Component('Cl-', formula='Cl', i_charge = -1.0,
                        particle_size='Soluble', degradability='Undegradable',
                        organic=False)
    
    Fi = Component('Fi', description='Ferricyanide ion', formula='Fe(CN)6', i_charge = -3.0,
                       particle_size='Soluble', degradability='Undegradable',
                       organic=False)
    
    Fo = Component('Fo', description='Ferrocyanide ion', formula='Fe(CN)6', i_charge = -4.0,
                       particle_size='Soluble', degradability='Undegradable',
                       organic=False)

    # Create a Components instance
    cmps_ed_vfa = Components([S_ac, S_et, S_la, S_pro, S_bu, S_va, S_su, S_he, Na, Cl, Fi, Fo, cmps_all.H2O])

    # Compile the components
    cmps_ed_vfa.default_compile()

    # Set thermodynamic settings if specified
    if set_thermo: settings.set_thermo(cmps_ed_vfa)

    return cmps_ed_vfa
cmps = create_ed_vfa_cmps()
# C2 = cmps['S_ac']
# C3 = cmps['S_la', 'S_pro']
# C4 = cmps['S_bu']
# C5 = cmps['S_va']
# C6 = cmps['S_su']
# ed_vfa_cmps = create_ed_vfa_cmps()
# I need to group C2 to C6 later?
#%%
# If we need in industrial application
# feed_stream = qs.WasteStream(ID='feed_stream', T=298.15, P=101325, phase='l',
#                              components={'Water': 965, 'NaCl': 35})
Q_dc = 0.05 # L/hr, Q_dc, Q_ac may be needed
Q_ac = 0.05
HRT_dc = 5 # hr, HRT_dc, HRT_ac may be needed, V_ac should be smaller than V_dc due to upconcentration
HRT_ac = 1

#!!! Below is needed to be revised
inf_dc = WasteStream(ID='inf_dc')
inf_dc.set_flow_by_concentration(concentrations={'S_pro': 5000, 'S_bu': 5000, 'S_he': 5000}, units=('L/hr', 'mg/L'))
eff_dc = WasteStream(ID='eff_dc')
eff_ac = WasteStream(ID='eff_ac')

#WasteStream
#Same as ADM1 Effluent Q
# inf_dc = WasteStream(ID='inf_dc')
# inf_dc.set_flow_by_concentration(flow_tot=5, concentrations={'S_pro': 5000, 'S_bu': 5000, 'S_he': 5000}, units=('L/hr', 'mg/L'))
# inf_ac = WasteStream(ID='inf_ac')
# inf_ac.set_flow_by_concentration(flow_tot=5, concentrations={'Na+': 500, 'Cl-': 500}, units=('L/hr', 'mg/L'))
# eff_dc = WasteStream(ID='eff_dc')               # effluent
# eff_ac = WasteStream(ID='eff_ac')               # effluent
#%%
# SanUnit

# def mass2mol_conversion(cmps):
#     '''conversion factor from kg[measured_as]/m3 to mol[component]/L'''
#     return cmps.i_mass / cmps.chem_MW
# unit_conversion = mass2mol_conversion(cmps)
F=96485.33289  # Faraday's constant in C/mol

class ED_vfa(SanUnit):
    def __init__(self, ID='', ins=None, outs=None, thermo=None, init_with='WasteStream',
                 permselectivity=None,  # Dictionary of permselectivity for each ion pair
                 j=500,  # Current density in A/m^2
                 t=3600,  # Time in seconds
                 A_m=0.5,  # Membrane area in m^2
                 V_dc=Q_dc*HRT_dc,  # Volume of dilute tank in m^3
                 V_ac=Q_ac*HRT_ac,  # Volume of accumulated tank in m^3
                 z_T=1.0,
                 r_m=1.0,  # Membrane resistance in Ohm*m^2
                 r_s=1.0):  # Solution resistance in Ohm*m^2
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.permselectivity = permselectivity or {'S_pro/S_bu': 1.0, 'S_pro/S_he': 1.0, 'S_bu/S_he': 1.0}  # Default permselectivity
        self.j = j
        self.t = t
        self.A_m = A_m
        self.V_dc = V_dc
        self.V_ac = V_ac
        self.z_T = z_T
        self.r_m = r_m
        self.r_s = r_s
        
        # Initialize dictionaries to store transport data
        self.n_T_dict = {}
        self.J_T_dict = {}
        
    _N_ins = 1
    _N_outs = 2

    def _run(self):
        inf_dc = self.ins[0]
        eff_dc, eff_ac = self.outs

        # Calculate total current [A]
        I = self.j * self.A_m
        self.total_current = I
        print(f"Total current (I): {I} A")

        # Obtain the flow rates from the influent streams
        Q_dc = inf_dc.F_vol  # Flow rate from influent dilute stream in m^3/hr
        self.Q_dc = Q_dc / 3600  # Convert to m^3/s

        print(f"Flow rate (Q_dc): {Q_dc} m^3/hr")
        print(f"Converted flow rate (Q_dc): {self.Q_dc} m^3/s")

        # Print volumes of the reactors
        print(f"Dilute tank volume (V_dc): {self.V_dc} m^3")
        print(f"Accumulated tank volume (V_ac): {self.V_ac} m^3")
        
        # eff_dc.copy_like(inf_dc)
        # eff_ac.copy_like(inf_dc)

        initial_concentrations = {
            'S_pro': inf_dc.imol['S_pro'] * 1000 / Q, # = kmole/hr * 1000 * hr/L = mole/L
            'S_bu': inf_dc.imol['S_bu'] * 1000 / Q, # mole/L
            'S_he': inf_dc.imol['S_he'] * 1000 / Q # mole/L
        }
        
        print(f"Initial concentrations: {initial_concentrations} mole/L")

        # Calculate the permselectivity based CE for each ion
        ce_S_pro = 0.2  # Example fixed CE value for S_pro
        ce_S_bu = ce_S_pro / self.permselectivity['S_pro/S_bu'] * (initial_concentrations['S_pro'] / initial_concentrations['S_bu'])
        ce_S_he = ce_S_pro / self.permselectivity['S_pro/S_he'] * (initial_concentrations['S_pro'] / initial_concentrations['S_he'])
        
        # Adjust CE values if their sum exceeds 1
        total_ce = ce_S_pro + ce_S_bu + ce_S_he
        if total_ce > 1:
            scale_factor = 1 / total_ce
            ce_S_pro *= scale_factor
            ce_S_bu *= scale_factor
            ce_S_he *= scale_factor
            
        CE_dict = {'S_pro': ce_S_pro, 'S_bu': ce_S_bu, 'S_he': ce_S_he}
        self.CE_dict = CE_dict
        
        print(f"Current efficiency dictionary: {CE_dict}")
        
        for ion, CE in self.CE_dict.items():
            # Moles of target ion transported [mol]
            n_T = CE * I / (self.z_T * F) * self.t
            self.n_T_dict[ion] = n_T
            print(f"Moles of {ion} transported (n_T): {n_T} mol") # Okay

            # Target ion molar flux [mol/(m^2*s)]
            J_T = CE * I / (self.z_T * F * self.A_m)
            self.J_T_dict[ion] = J_T
            print(f"Target ion molar flux (J_T) for {ion}: {J_T} mol/m^2/s")

            # Initial target ion moles in dilute and accumulated tanks
            # mole = kmole/hr * 1000 * hr/m3 * m3
            n_D_tank_initial = inf_dc.imol[ion] * 1000 / inf_dc.F_vol * self.V_dc
            n_A_tank_initial = 0
            
            print(f"Initial moles in dilute tank (n_D_tank_initial) for {ion}: {n_D_tank_initial} mol")
            print(f"Initial moles in accumulated tank (n_A_tank_initial) for {ion}: {n_A_tank_initial} mol") # Okay
            
            # mol/L below
            C_D_tank = (inf_dc.imol[ion] / 3600 * self.t - J_T * self.A_m * self.t) / (self.V_dc * 1000)
            C_A_tank = J_T * self.A_m * self.t / (self.V_dc * 1000)
            
            # mol/hr below
            eff_dc.imol[ion] = C_D_tank * self.V_dc * 1000 / (self.t / 3600)
            eff_ac.imol[ion] = C_A_tank * self.V_ac * 1000 / (self.t / 3600)
            
            # Ensure non-negative values
            eff_dc.imol[ion] = max(eff_dc.imol[ion], 0)
            eff_ac.imol[ion] = max(eff_ac.imol[ion], 0)
            
            eff_dc.imol[ion] = inf_dc.imol[ion] - eff_ac.imol[ion]
            # # Update moles in dilute and accumulated tanks [mole]
            # n_D_tank = n_D_tank_initial - J_T * self.A_m * self.t
            # n_A_tank = n_A_tank_initial + J_T * self.A_m * self.t

            # # Ensure non-negative moles [mole]
            # n_D_tank = max(n_D_tank, 0)
            # n_A_tank = max(n_A_tank, 0)
            
            # print(f"Updated moles in dilute tank (n_D_tank) for {ion}: {n_D_tank}")
            # print(f"Updated moles in accumulated tank (n_A_tank) for {ion}: {n_A_tank}")

            # # Update effluent streams with moles [kmole/hr]
            # eff_dc.imol[ion] = n_D_tank * 3600 / (1000 * self.t)
            # eff_ac.imol[ion] = n_A_tank * 3600 / (1000 * self.t)

            # # Calculate concentrations in dilute and accumulated tanks
            # C_D_tank = n_D_tank / self.V_dc * 1000 # mol/L
            # C_A_tank = n_A_tank / self.V_ac * 1000 # mol/L

            # # Ensure non-negative concentrations
            # C_D_tank = max(C_D_tank, 0)
            # C_A_tank = max(C_A_tank, 0)
            
            print(f"Concentration in dilute tank (C_D_tank) for {ion}: {C_D_tank}")
            print(f"Concentration in accumulated tank (C_A_tank) for {ion}: {C_A_tank}")
        # # Adjust the mass balance to ensure it matches the influent
        # for comp in inf_dc.chemicals:
        #     if comp.ID not in self.CE_dict:
        #         eff_dc.imass[comp.ID] = inf_dc.imass[comp.ID]

        # for comp in eff_ac.chemicals:
        #     if comp.ID not in self.CE_dict:
        #         eff_ac.imass[comp.ID] = eff_ac.imass[comp.ID]

        # Calculate system resistance [Ohm]
        R_sys = self.A_m * (self.r_m + self.r_s)
        self.R_sys = R_sys

        # Calculate system voltage [V]
        V_sys = R_sys * I
        self.V_sys = V_sys

        # Calculate power consumption [W]
        P_sys = V_sys * I
        self.P_sys = P_sys
        print(f"System resistance (R_sys): {R_sys} Ohm")
        print(f"System voltage (V_sys): {V_sys} V")
        print(f"Power consumption (P_sys): {P_sys} W")

    _units = {
        'Membrane area': 'm2',
        'Total current': 'A',
        'Dilute tank volume': 'm3',
        'Accumulated tank volume': 'm3',
        'System resistance': 'Ohm',
        'System voltage': 'V',
        'Power consumption': 'W'
    }
            
    def _design(self):
        D = self.design_results
        D['Membrane area'] = self.A_m
        D['Total current'] = self.j * self.A_m
        D['Dilute tank volume'] = self.V_dc
        D['Accumulated tank volume'] = self.V_ac
        D['System resistance'] = self.R_sys
        D['System voltage'] = self.V_sys
        D['Power consumption'] = self.P_sys

    def _cost(self):
        self.baseline_purchase_costs['Membrane'] = 100 * self.design_results['Membrane area']  # Assume $100 per m^2 for the membrane
        self.power_utility.consumption = self.design_results['Power consumption'] / 1000  # Assuming kWh consumption based on power
#%%
# Initialize the ED_vfa unit with permselectivity values and different tank volumes
# permselectivity={'S_pro/S_bu': 1.063682, 'S_pro/S_he': 1.555841, 'S_bu/S_he': 1.462693}
ed1 = ED_vfa(
    ID='ED1',
    ins=(inf_dc),
    outs=(eff_dc, eff_ac),
)
#%%
# Simulate the process
ed1.simulate()
print(ed1.results())