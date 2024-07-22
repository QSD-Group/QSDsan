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
import matplotlib.pyplot as plt

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
# I need to group C2 to C6 later?
#%%
# If we need in industrial application
# feed_stream = qs.WasteStream(ID='feed_stream', T=298.15, P=101325, phase='l',
#                              components={'Water': 965, 'NaCl': 35})
Q_dc = 0.3 # L/hr, Q_dc, Q_ac may be needed
Q_ac = 0.3 # 5mL/min flow circulation (NY)
HRT_dc = 3.33 # hr, HRT_dc, HRT_ac may be needed, V_ac should be smaller than V_dc due to upconcentration
HRT_ac = 0.33 # DC: 1L, AC: 0.1L (NY)

#WasteStream
#Same as ADM1 Effluent Q
#S_pro = S_bu = S_he = 25mM (NY&WS)
inf_dc = WasteStream(ID='inf_dc')
inf_dc.set_flow_by_concentration(flow_tot=Q_dc, concentrations={'S_pro': 1852, 'S_bu': 2202, 'S_he': 2904}, units=('L/hr', 'mg/L'))
inf_ac = WasteStream(ID='inf_ac')
inf_ac.set_flow_by_concentration(flow_tot=Q_ac, concentrations={'Na+': 500, 'Cl-': 500}, units=('L/hr', 'mg/L'))
eff_dc = WasteStream(ID='eff_dc')               # effluent
eff_ac = WasteStream(ID='eff_ac')               # effluent
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
                 j=8.23,  # Current density in A/m^2 under 0.8 V (WS)
                 t=3600*6,  # Time in seconds
                 A_m=0.016,  # Membrane area in m^2 (NY&WS)
                 V_dc=Q_dc*HRT_dc,  # Volume of dilute tank in m^3
                 V_ac=Q_ac*HRT_ac,  # Volume of accumulated tank in m^3
                 z_T=1.0,
                 r_m=1.0,  # Membrane resistance in Ohm*m^2
                 r_s=1.0,
                 ce_S_pro=0.1):  # Solution resistance in Ohm*m^2
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.permselectivity = permselectivity or {'S_pro/S_bu': 1.063682, 'S_pro/S_he': 1.555841, 'S_bu/S_he': 1.462693}  # Default permselectivity
        self.j = j
        self.t = t
        self.A_m = A_m
        self.V_dc = V_dc
        self.V_ac = V_ac
        self.z_T = z_T
        self.r_m = r_m
        self.r_s = r_s
        self.ce_S_pro = ce_S_pro
        # Initialize dictionaries to store transport data
        self.n_T_dict = {}
        self.J_T_dict = {}
        # Initializing effluent streams with influent values
        eff_dc.copy_like(inf_dc)
        eff_ac.copy_like(inf_ac)
        
    _N_ins = 2
    _N_outs = 2

    def _run(self):
        inf_dc, inf_ac = self.ins
        eff_dc, eff_ac = self.outs

        # Calculate total current [A]
        I = self.j * self.A_m
        self.total_current = I
        print(f"Total current (I): {I} A")

        # Print volumes of the reactors
        print(f"Dilute tank volume (V_dc): {self.V_dc} m^3")
        print(f"Accumulated tank volume (V_ac): {self.V_ac} m^3")
        
        # eff_dc.copy_like(inf_dc)
        # eff_ac.copy_like(inf_dc)

        initial_concentrations = {
            'S_pro': inf_dc.imol['S_pro'] * 1000 / Q_dc, # = kmole/hr * 1000 * hr/L = mole/L
            'S_bu': inf_dc.imol['S_bu'] * 1000 / Q_dc, # mole/L
            'S_he': inf_dc.imol['S_he'] * 1000 / Q_dc # mole/L
        }
        
        print(f"Initial concentrations: {initial_concentrations} mole/L")

        # Calculate the permselectivity based CE for each ion
        ce_S_pro = self.ce_S_pro  # Example fixed CE value for S_pro
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
        
        print(f"Charge efficiency dictionary: {CE_dict}")
        
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
            C_D_tank = (n_D_tank_initial - J_T * self.A_m * self.t) / (self.V_dc * 1000)
            C_A_tank = J_T * self.A_m * self.t / (self.V_ac * 1000)
            print(f"Concentration in dilute tank (C_D_tank) for {ion}: {C_D_tank} mole/L")
            print(f"Concentration in accumulated tank (C_A_tank) for {ion}: {C_A_tank} mole/L")
            
            # kmole/hr below
            eff_dc.imol[ion] = C_D_tank * self.V_dc / (self.t / 3600)
            eff_ac.imol[ion] = C_A_tank * self.V_ac / (self.t / 3600)
            
            # Ensure non-negative values
            eff_dc.imol[ion] = max(eff_dc.imol[ion], 0)
            eff_ac.imol[ion] = max(eff_ac.imol[ion], 0)

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
    ins=(inf_dc, inf_ac),
    outs=(eff_dc, eff_ac),
)
#%%
# Simulate the process
ed1.simulate()
print(ed1.results())
#%%
# Define the time durations for the simulations (in seconds)
time_durations = [3600 * t for t in range(1, 21)]

# Initialize a dictionary to store the C_A_tank values for each ion
C_A_tank_values = {'S_pro': [], 'S_bu': [], 'S_he': []}

# Loop over the defined time durations and run the simulation for each
for t in time_durations:
    # Create a new instance of the ED_vfa class for each simulation
    ed = ED_vfa(
        ID=f'ED_t_{t}',
        ins=(inf_dc, inf_ac),
        outs=(eff_dc, eff_ac),
        t=t  # Set the time duration for this simulation
    )
    
    # Run the simulation
    ed.simulate()
    
    # Extract the C_A_tank values and store them
    for ion in C_A_tank_values.keys():
        initial_concentration = inf_dc.imol[ion] * 1000 / Q_dc
        CE = ed.CE_dict[ion]
        J_T = CE * ed.total_current / (ed.z_T * F * ed.A_m)
        C_A_tank = J_T * ed.A_m * ed.t / ed.V_ac # M
        C_A_tank_mM = C_A_tank * 1000  # Convert to mM
        C_A_tank_values[ion].append(C_A_tank_mM)
        
# Set font sizes
plt.rcParams.update({
    'font.size': 14,
    'axes.labelsize': 16,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14,
    'legend.fontsize': 14,
    'axes.titlesize': 18,
})

# Plot the results for C_A_tank
plt.figure(figsize=(12, 8))
for ion, values in C_A_tank_values.items():
    plt.plot(range(1, 21), values, label=ion)

plt.xlabel('Time (hours)')
plt.ylabel('C_A_tank (mM)')
plt.title('C_A_tank over time for different ions (in mM)')
plt.legend()
plt.grid(True)
plt.show()
#%%
# Initialize a dictionary to store the C_D_tank values for each ion
C_D_tank_values = {'S_pro': [], 'S_bu': [], 'S_he': []}

# Loop over the defined time durations and run the simulation for each
for t in time_durations:
    # Create a new instance of the ED_vfa class for each simulation
    ed = ED_vfa(
        ID=f'ED_t_{t}',
        ins=(inf_dc, inf_ac),
        outs=(eff_dc, eff_ac),
        t=t  # Set the time duration for this simulation
    )
    
    # Run the simulation
    ed.simulate()
    
    # Extract the C_D_tank values and store them
    for ion in C_D_tank_values.keys():
        initial_concentration = inf_dc.imol[ion] * 1000 / Q_dc
        CE = ed.CE_dict[ion]
        J_T = CE * ed.total_current / (ed.z_T * F * ed.A_m)
        C_D_tank = (initial_concentration - J_T * ed.A_m * ed.t / ed.V_dc) # M
        C_D_tank_mM = C_D_tank *1000
        C_D_tank_values[ion].append(C_D_tank_mM)
        
# Set font sizes
plt.rcParams.update({
    'font.size': 18,
    'axes.labelsize': 18,
    'xtick.labelsize': 16,
    'ytick.labelsize': 16,
    'legend.fontsize': 16,
    'axes.titlesize': 20,
})

# Plot the results for C_D_tank in mM
plt.figure(figsize=(12, 8))
for ion, values in C_D_tank_values.items():
    plt.plot(range(1, 21), values, label=ion)

plt.xlabel('Time (hours)')
plt.ylabel('C_D_tank (mM)')
plt.title('C_D_tank over time for different ions (in mM)')
plt.legend()
plt.grid(True)
plt.show()
#%%
# Define the time duration for the simulation (20 hours in seconds)
time_duration = 3600 * 20

# Define the different membrane area values
A_m_values = [0.016, 0.025, 0.036, 0.049, 0.064, 0.081, 0.1]

# Initialize a dictionary to store the C_A_tank values for each ion and A_m
C_A_tank_results = {ion: [] for ion in ['S_pro', 'S_bu', 'S_he']}

# Loop over the different membrane area values and run the simulation for each
for A_m in A_m_values:
    # Create a new instance of the ED_vfa class for each simulation
    ed = ED_vfa(
        ID=f'ED_A_m_{A_m}',
        ins=(inf_dc, inf_ac),
        outs=(eff_dc, eff_ac),
        t=time_duration,  # Set the time duration for this simulation
        A_m=A_m,  # Set the membrane area for this simulation
        ce_S_pro=0.4  # Set ce_S_pro to 0.1
    )
    
    # Calculate the permselectivity-based CE for each ion
    initial_concentrations = {
        'S_pro': inf_dc.imol['S_pro'] * 1000 / Q_dc,  # mole/L
        'S_bu': inf_dc.imol['S_bu'] * 1000 / Q_dc,  # mole/L
        'S_he': inf_dc.imol['S_he'] * 1000 / Q_dc  # mole/L
    }

    ce_S_pro = ed.ce_S_pro  # Fixed CE value for S_pro
    ce_S_bu = ce_S_pro / ed.permselectivity['S_pro/S_bu'] * (initial_concentrations['S_pro'] / initial_concentrations['S_bu'])
    ce_S_he = ce_S_pro / ed.permselectivity['S_pro/S_he'] * (initial_concentrations['S_pro'] / initial_concentrations['S_he'])
    
    # Adjust CE values if their sum exceeds 1
    total_ce = ce_S_pro + ce_S_bu + ce_S_he
    if (total_ce > 1):
        scale_factor = 1 / total_ce
        ce_S_pro *= scale_factor
        ce_S_bu *= scale_factor
        ce_S_he *= scale_factor

    CE_dict = {'S_pro': ce_S_pro, 'S_bu': ce_S_bu, 'S_he': ce_S_he}
    ed.CE_dict = CE_dict
    
    # Run the simulation
    ed.simulate()
    
    # Extract the C_A_tank values and store them
    for ion in C_A_tank_results.keys():
        CE = ed.CE_dict[ion]
        J_T = CE * ed.total_current / (ed.z_T * F * ed.A_m)
        
        # Calculate C_A_tank and convert to mM
        C_A_tank = J_T * ed.A_m * ed.t / ed.V_ac  # M
        C_A_tank_mM = C_A_tank * 1000  # Convert to mM
        C_A_tank_results[ion].append(C_A_tank_mM)

# Plot the results for C_A_tank in mM as bar graph
fig, ax = plt.subplots(figsize=(12, 8))
bar_width = 0.2
index = np.arange(len(A_m_values))

# Plot C_A_tank
for i, ion in enumerate(C_A_tank_results.keys()):
    ax.bar(index + i * bar_width, C_A_tank_results[ion], bar_width, label=ion)

ax.set_xlabel('Membrane Area (m^2)')
ax.set_ylabel('C_A_tank (mM)')
ax.set_title('C_A_tank over different membrane areas (in mM)')
ax.set_xticks(index + bar_width)
ax.set_xticklabels(A_m_values)
ax.legend()
ax.grid(True)

# Set font sizes
plt.rcParams.update({
    'font.size': 14,
    'axes.labelsize': 16,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14,
    'legend.fontsize': 14,
    'axes.titlesize': 18,
})

plt.tight_layout()
plt.show()
#%%
# Define the time duration for the simulation (20 hours in seconds)
time_duration = 3600 * 20

# Define the different ce_S_pro values
ce_S_pro_values = [0.1, 0.2, 0.3, 0.4, 0.5]

# Define the membrane area
A_m = 0.016  # Example membrane area in m^2

# Initialize a dictionary to store the C_A_tank values for each ion and ce_S_pro
C_A_tank_results = {ion: [] for ion in ['S_pro', 'S_bu', 'S_he']}

# Loop over the different ce_S_pro values and run the simulation for each
for ce_S_pro in ce_S_pro_values:
    # Create a new instance of the ED_vfa class for each simulation
    ed = ED_vfa(
        ID=f'ED_ce_S_pro_{ce_S_pro}',
        ins=(inf_dc, inf_ac),
        outs=(eff_dc, eff_ac),
        t=time_duration,  # Set the time duration for this simulation
        A_m=A_m,  # Set the membrane area for this simulation
        ce_S_pro=ce_S_pro  # Set ce_S_pro
    )
    
    # Calculate the permselectivity-based CE for each ion
    initial_concentrations = {
        'S_pro': inf_dc.imol['S_pro'] * 1000 / Q_dc,  # mole/L
        'S_bu': inf_dc.imol['S_bu'] * 1000 / Q_dc,  # mole/L
        'S_he': inf_dc.imol['S_he'] * 1000 / Q_dc  # mole/L
    }

    ce_S_bu = ce_S_pro / ed.permselectivity['S_pro/S_bu'] * (initial_concentrations['S_pro'] / initial_concentrations['S_bu'])
    ce_S_he = ce_S_pro / ed.permselectivity['S_pro/S_he'] * (initial_concentrations['S_pro'] / initial_concentrations['S_he'])
    
    # Adjust CE values if their sum exceeds 1
    total_ce = ce_S_pro + ce_S_bu + ce_S_he
    if (total_ce > 1):
        scale_factor = 1 / total_ce
        ce_S_pro *= scale_factor
        ce_S_bu *= scale_factor
        ce_S_he *= scale_factor

    CE_dict = {'S_pro': ce_S_pro, 'S_bu': ce_S_bu, 'S_he': ce_S_he}
    ed.CE_dict = CE_dict
    
    # Run the simulation
    ed.simulate()
    
    # Extract the C_A_tank values and store them
    for ion in C_A_tank_results.keys():
        CE = ed.CE_dict[ion]
        J_T = CE * ed.total_current / (ed.z_T * F * ed.A_m)
        
        # Calculate C_A_tank and convert to mM
        C_A_tank = J_T * ed.A_m * ed.t / ed.V_ac  # M
        C_A_tank_mM = C_A_tank * 1000  # Convert to mM
        C_A_tank_results[ion].append(C_A_tank_mM)

# Plot the results for C_A_tank in mM as bar graph
fig, ax = plt.subplots(figsize=(12, 8))
bar_width = 0.15
index = np.arange(len(ce_S_pro_values))

# Plot C_A_tank
for i, ion in enumerate(C_A_tank_results.keys()):
    ax.bar(index + i * bar_width, C_A_tank_results[ion], bar_width, label=ion)

ax.set_xlabel('CE of S_pro')
ax.set_ylabel('C_A_tank (mM)')
ax.set_title('C_A_tank over different CE of S_pro (in mM)')
ax.set_xticks(index + bar_width)
ax.set_xticklabels(ce_S_pro_values)
ax.legend()
ax.grid(True)

# Set font sizes
plt.rcParams.update({
    'font.size': 14,
    'axes.labelsize': 16,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14,
    'legend.fontsize': 14,
    'axes.titlesize': 18,
})

plt.tight_layout()
plt.show()