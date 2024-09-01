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

#_path = ospath.join(data_path, 'process_data/_ed_vfa.tsv')
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
#WasteStream
#Same as ADM1 Effluent Q
inf_dc = WasteStream(ID='inf_dc')
inf_dc.set_flow_by_concentration(flow_tot=5, concentrations={'S_pro': 5000, 'S_bu': 5000, 'S_he': 5000}, units=('L/hr', 'mg/L'))
#fc_eff = WasteStream('FC_Effluent', X_GAO_Gly=.5, H2O=1000, units='kmol/hr')
inf_ac = WasteStream(ID='inf_ac')
inf_ac.set_flow_by_concentration(flow_tot=5, concentrations={'Na+': 500, 'Cl-': 500}, units=('L/hr', 'mg/L'))
#ac_eff = WasteStream('AC_Effluent', X_GAO_Gly=.5, H2O=1000, units='kmol/hr')
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
    def __init__(self, ID='', ins=None, outs=None, thermo=None, isdynamic=True, init_with='WasteStream',
                 CE_dict=None,  # Dictionary of charge efficiencies for each ion
                 j=8.23,   # Current density in A/m^2
                 t=3600*6,  # Time in seconds
                 A_m=0.0016,  # Membrane area in m^2
                 V=0.2/1000,  # Volume of all tanks in m^3
                 z_T=1.0,
                 r_m=3206.875,  # Areal Membrane resistance in Ohm*m^2
                 r_s=4961.875,  # Areal Solution resistance in Ohm*m^2
                 **kwargs
                 ):  # Solution resistance in Ohm*m^2
        SanUnit.__init__(self, ID, ins, outs, thermo, isdynamic=isdynamic)
        self.CE_dict = CE_dict or {'S_pro': 0.4, 'S_bu': 0.37, 'S_he': 0.23}  # Default CE for ions
        self.j = j
        self.t = t
        self.A_m = A_m
        self.V = V
        self.z_T = z_T
        self.r_m = r_m
        self.r_s = r_s
        
    _N_ins = 2
    _N_outs = 2

    @property
    def j(self):
        return self._j
    @j.setter
    def j(self, value):
        if value <= 0:
            raise ValueError("Current density must be positive.")
        self._j = value

    @property
    def t(self):
        return self._t
    @t.setter
    def t(self, value):
        if value <= 0:
            raise ValueError("Time must be positive.")
        self._t = value
        
    @property
    def A_m(self):
        return self._A_m
    @A_m.setter
    def A_m(self, value):
        if value <= 0:
            raise ValueError("Membrane area must be positive.")
        self._A_m = value

    @property
    def V(self):
        return self._V
    @V.setter
    def V(self, value):
        if value <= 0:
            raise ValueError("Volume must be positive.")
        self._V = value
        
    def _run(self):
        inf_dc, inf_ac = self.ins
        eff_dc, eff_ac = self.outs
    
        # Calculate total current [A]
        I = self.j * self.A_m
        self.total_current = I
        print(f"Total current (I): {I} A")
    
        # Obtain the flow rates from the influent streams
        Q_dc = inf_dc.F_vol  # Flow rate from influent dilute stream in m^3/hr
        Q_ac = inf_ac.F_vol  # Flow rate from influent accumulated stream in m^3/hr
        self.Q_dc = Q_dc / 3600  # Convert to m^3/s
        self.Q_ac = Q_ac / 3600  # Convert to m^3/s
    
        self.n_T_dict = {}
        self.J_T_dict = {}
        self.influent_dc_conc = {}
        self.influent_ac_conc = {}
    
        # Print original flow rates [m3/hr]
        print(f"Flow rate (Q_dc): {Q_dc} m^3/hr")
        print(f"Flow rate (Q_ac): {Q_ac} m^3/hr")
    
        # Print the converted flow rates [m3/s]
        print(f"Converted flow rate (Q_dc): {self.Q_dc} m^3/s")
        print(f"Converted flow rate (Q_ac): {self.Q_ac} m^3/s")
    
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
    
    def _design(self):
        D = self.design_results
        D['Membrane area'] = self.A_m
        D['Total current'] = self.j * self.A_m
        D['Tank volume'] = self.V
        D['System resistance'] = self.R_sys
        D['System voltage'] = self.V_sys
        D['Power consumption'] = self.P_sys

    def _cost(self):
        self.baseline_purchase_costs['Membrane'] = 100 * self.design_results['Membrane area']  # Assume $100 per m^2 for the membrane
        self.power_utility.consumption = self.design_results['Power consumption'] / 1000  # Assuming kWh consumption based on power
#%%
# Initialize the ED_vfa unit
ed1 = ED_vfa(
    ID='ED1',
    ins=(inf_dc, inf_ac),
    outs=(eff_dc, eff_ac),
    CE_dict={'S_pro': 0.2, 'S_bu': 0.2, 'S_he': 0.2},  # Separate CE for each ion
    j=5,
    t=288000,
    A_m=0.5,
    V=0.1
)
#%%
# System
'''
# Run the unit
ed_vfa_unit.simulate()
ed_vfa_unit.show()
#%%
'''
# Create the system
sys = System('ED1', path=(ed1,))

# # Simulate the system
# sys.simulate()
# sys._setup()
# sys.converge()
# sys.diagram()
#%%
# Simulation
# Set the dynamic tracker
sys.set_dynamic_tracker(inf_dc, inf_ac, eff_dc, eff_ac, ed1)
# sys
# Simulation settings
t = 288000  # total time for simulation in hours
t_step = 3600  # times at which to store the computed solution in hours
method = 'BDF'  # integration method to use

# Run simulation
sys.simulate(state_reset_hook='reset_cache',
             t_span=(0, t),
             t_eval=np.arange(0, t + t_step, t_step),
             method=method,
             export_state_to='ED_vfa_simulation_results.xlsx')
sys
#%%
# Print the results
print("Effluent dilute stream concentrations (mol/L):")
for ion in ed1.CE_dict.keys():
    eff_dc_conc = eff_dc.imol[ion] / (eff_dc.F_vol * 1000)  # Convert volume from m^3/hr to L/hr for concentration in mol/L
    print(f"{ion}: {eff_dc_conc} mol/L")

print("Effluent accumulated stream concentrations (mol/L):")
for ion in ed1.CE_dict.keys():
    eff_ac_conc = eff_ac.imol[ion] / (eff_ac.F_vol * 1000)  # Convert volume from m^3/hr to L/hr for concentration in mol/L
    print(f"{ion}: {eff_ac_conc} mol/L")
#%%
# Total Mass Calculation
def calculate_total_mass(stream):
    total_mass = 0
    for component in stream.chemicals:
        total_mass += stream.imass[component.ID]
    return total_mass

# Total mass of inf
total_mass_inf_dc = calculate_total_mass(inf_dc)
total_mass_inf_ac = calculate_total_mass(inf_ac)

# Total mass of eff
total_mass_eff_dc = calculate_total_mass(eff_dc)
total_mass_eff_ac = calculate_total_mass(eff_ac)

# Total Mass
print(f"Total mass of inf_dc: {total_mass_inf_dc} kg/hr")
print(f"Total mass of inf_ac: {total_mass_inf_ac} kg/hr")
print(f"Total mass of eff_dc: {total_mass_eff_dc} kg/hr")
print(f"Total mass of eff_ac: {total_mass_eff_ac} kg/hr")

# Check Mass Balance
total_mass_in = total_mass_inf_dc + total_mass_inf_ac
total_mass_out = total_mass_eff_dc + total_mass_eff_ac

print(f"Total mass in (inf_dc + inf_ac): {total_mass_in} kg/hr")
print(f"Total mass out (eff_dc + eff_ac): {total_mass_out} kg/hr")
print(f"Mass balance check: {'Balanced' if total_mass_in == total_mass_out else 'Not balanced'}")