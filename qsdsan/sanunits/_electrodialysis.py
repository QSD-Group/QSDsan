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
                 j=500,   # Current density in A/m^2
                 t=3600,  # Time in seconds
                 A_m=10,  # Membrane area in m^2
                 V=1.0,  # Volume of all tanks in m^3
                 z_T=1.0,
                 r_m=1.0,  # Membrane resistance in Ohm*m^2
                 r_s=1.0,
                 **kwargs
                 ):  # Solution resistance in Ohm*m^2
        SanUnit.__init__(self, ID, ins, outs, thermo, isdynamic=isdynamic)
        self.CE_dict = CE_dict or {'S_pro': 0.1, 'S_bu': 0.1, 'S_he': 0.1}  # Default CE for ions
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
        
    @property
    def n_T(self):
        n_T_dict = {}
        I = self.total_current
        for ion, CE in self.CE_dict.items():
            n_T = CE * I / (self.z_T * F) * self.t
            n_T_dict[ion] = n_T
        return n_T_dict

    @property
    def J_T(self):
        J_T_dict = {}
        I = self.total_current
        for ion, CE in self.CE_dict.items():
            J_T = CE * I / (self.z_T * F * self.A_m)
            J_T_dict[ion] = J_T
        return J_T_dict
    
    def _run(self):
        inf_dc, inf_ac = self.ins
        eff_dc, eff_ac = self.outs
        
        # Calculate total current [A]    
        I = self.j * self.A_m
        self.total_current = I
        print(f"Total current (I): {I} A")
        
        # Obtain the flow rates from the influent streams
        Q_dc = inf_dc.F_vol # Flow rate from influent dilute stream in m^3/hr
        Q_ac = inf_ac.F_vol # Flow rate from influent accumulated stream in m^3/hr
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
        
        # Initializing effluent streams with influent values
        eff_dc.copy_like(inf_dc)
        eff_ac.copy_like(inf_ac)
        print(f"eff_dc: {eff_dc}")
        print(f"eff_ac: {eff_ac}")
        cmps = inf_dc.components
        mass2mol = cmps.i_mass / cmps.chem_MW 
        
        for ion, CE in self.CE_dict.items():        
            # Moles of target ion transported [mol]
            n_T = CE * I / (self.z_T * F) * self.t
            self.n_T_dict[ion] = n_T
            print(f"Moles of {ion} transported (n_T): {n_T} mol")
            
            # Target ion molar flux [mol/(m2*s)]
            J_T = CE * I / (self.z_T * F * self.A_m)
            self.J_T_dict[ion] = J_T
            print(f"Target ion molar flux (J_T) for {ion}: {J_T} mol/m^2/s")

            # Effluent moles for dilute stream [mole]
            idx = cmps.index(ion)
            n_ins_dc = inf_dc.imass[ion] * mass2mol[idx] # kmol/hr
            n_out_dc = (n_ins_dc * self.t / 3600 * 1000) - (self.V * J_T * self.A_m) / self.Q_dc
            

            # Effluent moles for accumulated stream [moles]
            n_out_ac = (self.V * J_T * self.A_m) / self.Q_dc * (1 - np.exp(-self.Q_ac * self.t / self.V))

            # Ensure non-negative effluent moles [mole]
            n_out_dc = max(n_out_dc, 0)
            n_out_ac = max(n_out_ac, 0)

            # Update effluent streams [kmole/hr]
            n_out_dc /= mass2mol[idx]
            eff_dc.imass[ion] = n_out_dc / 1000 / self.t * 3600
            print(f'eff_dc.imol[ion] = {eff_dc.imol[ion]}')
            eff_ac.imol[ion] = n_out_ac / 1000 / self.t * 3600

        # # Ensure non-negative effluent moles for H2O [kmole/hr]
        # eff_dc.imol['H2O'] = max(inf_dc.imol['H2O'], 0)
        # eff_ac.imol['H2O'] = max(inf_ac.imol['H2O'], 0)
        
        # # Adjust the mass balance to ensure it matches the influent
        # Ensure non-negative effluent moles for H2O [kmole/hr]
        # eff_dc.imol['H2O'] = max(inf_dc.imol['H2O'], 0)
        # eff_ac.imol['H2O'] = max(inf_ac.imol['H2O'], 0)
        
        # Adjust the mass balance to ensure it matches the influent
        # for comp in inf_dc.chemicals:
        #     if comp.ID not in self.CE_dict:
        #         eff_dc.imass[comp.ID] = inf_dc.imass[comp.ID]
        
        # for comp in inf_ac.chemicals:
        #     if comp.ID not in self.CE_dict:
        #         eff_ac.imass[comp.ID] = inf_ac.imass[comp.ID]
                        
        # Calculate system resistance [Ohm] -> design function
        R_sys = self.A_m * (self.r_m + self.r_s)
        self.R_sys = R_sys
        
        # Calculate system voltage [V] -> design method
        V_sys = R_sys * I
        self.V_sys = V_sys
        
        # Calculate power consumption [W] -> design
        P_sys = V_sys * I
        self.P_sys = P_sys
        print(f"System resistance (R_sys): {R_sys} Ohm")
        print(f"System voltage (V_sys): {V_sys} V")
        print(f"Power consumption (P_sys): {P_sys} W")

#%%
    def _init_state(self):
        inf_dc, inf_ac = self.ins
        if inf_dc is None or inf_ac is None:
            raise ValueError("Input streams must not be None")
    
        Q_dc = inf_dc.F_vol / 3600
        Q_ac = inf_ac.F_vol / 3600
        # [S_ac, S_et, S_la, S_pro, S_bu, S_va, S_su, S_he, Na, Cl, Fi, Fo, cmps_all.H2O]), total 14
        C_dc = np.array([inf_dc.imol[ion] for ion in inf_dc.chemicals.IDs])
        print(f'C_dc = {C_dc}')
        C_ac = np.array([inf_ac.imol[ion] for ion in inf_ac.chemicals.IDs])
        print(f'C_ac = {C_ac}')
        inf_dc = np.append(C_dc, Q_dc)
        inf_ac = np.append(C_ac, Q_ac)
        # two_Qs = np.array([Q_dc, Q_ac])
        # two_Cs = np.array([C_dc, C_ac])
        # print(f'two_Cs = {two_Cs}')
        # average_Cs = two_Qs @ two_Cs / two_Qs.sum()
        # total_Q = two_Qs.sum()
        self.I = self.j * self.A_m
        print(f"Total current (I): {self.I} A")
    
        CEs = np.array(list(self.CE_dict.values()))
        self.n_T = CEs * self.I / (self.z_T * F) * self.t
        self.J_T = CEs * self.I / (self.z_T * F * self.A_m)
        self.a = self.t / 3600 * 1000
        self.b = self.V * self.J_T * self.A_m
        self.indices = cmps.indices(self.CE_dict.keys())
        self.mass2mol = cmps.i_mass / cmps.chem_MW
        self._state = np.append(inf_dc, inf_ac)
        # self._state = np.append(average_Cs, total_Q)
        self._dstate = np.zeros_like(self._state)  # Ensure _dstate is initialized
    
        print(f"Initial state (average concentrations and total flow rate): {self._state}")
        print(f"Initial state change rate (dstate): {self._dstate}")
#%%
# def _update_state(self):
#     arr = self._state   # retrieving the current state of the SanUnit
#     eff, = self.outs    # assuming this SanUnit has one outlet only
#     eff.state[:] = arr  # assume arr has the same shape as WasteStream.state
# inf_dc, inf_ac = self.ins # this line is wrong because you can't use state of influent in update state
        # # For sludge, the particulate concentrations are multiplied by thickener factor, and
        # # flowrate is multiplied by Qu_factor. The soluble concentrations remains same. 
        # uf, of = self.outs
        # if uf.state is None: uf.state = np.zeros(len(cmps)+1)
        # uf.state[:-1] = self._state[:-1]*cmps.s*1 + self._state[:-1]*cmps.x*thickener_factor
        # uf.state[-1] = self._state[-1]*Qu_factor
        
        # # For effluent, the particulate concentrations are multiplied by thinning factor, and
        # # flowrate is multiplied by Qu_factor. The soluble concentrations remains same. 
        # if of.state is None: of.state = np.zeros(len(cmps)+1)
        # of.state[:-1] = self._state[:-1]*cmps.s*1 + self._state[:-1]*cmps.x*thinning_factor
        # of.state[-1] = self._state[-1]*(1 - Qu_factor)
    def _update_state(self):
        eff_dc, eff_ac = self.outs
    
        # Assuming there are 14 compounds in each stream, update if different
        cmps = eff_dc.components
        n = len(cmps)+1
        n_compounds_dc = len(eff_dc.chemicals.IDs)
        n_compounds_ac = len(eff_ac.chemicals.IDs)
    
        # Print to check the actual number of compounds
        print(f"Number of compounds in dilute stream: {n_compounds_dc}")
        print(f"Number of compounds in accumulated stream: {n_compounds_ac}")
    
        # Unpack the state into effluent dilute and accumulated stream components
        # self._state == [Sac_inf_dc, Set_inf_dc,....., Q_inf_dc,Sac_inf_ac, Set_inf_ac,....., Q_inf_ac]
        eff_dc.state[:] = self._state[:14]
        eff_ac.state[:] = self._state[14:]
        
        # C_dc = self._state[:n_compounds_dc]
        Q_dc = self._state[n_compounds_dc]
        # C_ac = self._state[n_compounds_dc + 1:2 * n_compounds_dc + 1]
        Q_ac = self._state[-1]
        
        # Initialize effluent streams with influent values
        # eff_dc.copy_like(self.ins[0])
        # eff_ac.copy_like(self.ins[1])
    
        # Update effluent flow rates
        # eff_dc.F_vol = Q_dc * 3600  # Convert back to m^3/hr
        # eff_ac.F_vol = Q_ac * 3600  # Convert back to m^3/hr
    
        # Total current calculation (assuming `self.j` and `self.A_m` are defined in the class)
        # I = self.j * self.A_m
        # print(f"Total current (I): {I} A")
    
        # CEs = CE_dict.values()
        # n_T = CEs * I / (self.z_T * F) * self.t
        # J_T = CEs * I / (self.z_T * F * self.A_m)
        # a = self.t / 3600 * 1000
        # b = self.V * J_T * self.A_m
        # indices = cmps.indices(CE_dict.keys())
        
        a, b = self.a, self.b
        indices = self.indices
        mass2mol = self.mass2mol
        n_ins_dc = eff_dc.state[indices] * Q_dc * mass2mol[indices]   # molar flowrate of ions in inf_dc [mol/d]
        n_out_dc = ((n_ins_dc * a) - b / Q_dc) / mass2mol[indices]
        eff_dc.state[indices] = n_out_dc / Q_dc
        
        
        # for ion, CE in self.CE_dict.items():
        #     idx = cmps.index(ion)
        #     # Moles of target ion transported [mol]
        #     n_T = CE * I / (self.z_T * F) * self.t
        #     print(f"Moles of {ion} transported (n_T): {n_T} mol")
    
        #     # Target ion molar flux [mol/(m2*s)]
        #     J_T = CE * I / (self.z_T * F * self.A_m)
        #     print(f"Target ion molar flux (J_T) for {ion}: {J_T} mol/m^2/s")
    
        #     # Effluent moles for dilute stream [mole]
        #     n_ins_dc = eff_dc.state[idx] * Q_dc * mass2mol[idx]  # mol/d
        #     n_out_dc = (n_ins_dc * self.t / 3600 * 1000) - (self.V * J_T * self.A_m) / Q_dc
    
        #     # Effluent moles for accumulated stream [moles]
        #     n_out_ac = (self.V * J_T * self.A_m) / Q_dc * (1 - np.exp(-Q_ac * self.t / self.V))
        #     n_out_ac = b / Q_dc * (1-exp(-c*Q_ac))
    
        #     # Ensure non-negative effluent moles [mole]
        #     n_out_dc = max(n_out_dc, 0)
        #     n_out_ac = max(n_out_ac, 0)
    
        #     # Update effluent streams [kmole/hr]
        #     n_out_dc /= mass2mol # kg/hr
        #     eff_dc.state[idx] = n_out_dc*24 / Q_dc # needs to be in mg/L 
        #     # eff_dc.imol[ion] = n_out_dc / 1000 / self.t * 3600
        #     # eff_ac.imol[ion] = n_out_ac / 1000 / self.t * 3600
        
        # Ensure non-negative effluent moles for H2O [kmole/hr]
        # eff_dc.imol['H2O'] = max(self.ins[0].imol['H2O'], 0)
        # eff_ac.imol['H2O'] = max(self.ins[1].imol['H2O'], 0)
    
        # # Adjust the mass balance to ensure it matches the influent
        # for comp in self.ins[0].chemicals:
        #     if comp.ID not in self.CE_dict:
        #         eff_dc.imass[comp.ID] = self.ins[0].imass[comp.ID]
    
        # for comp in self.ins[1].chemicals:
        #     if comp.ID not in self.CE_dict:
        #         eff_ac.imass[comp.ID] = self.ins[1].imass[comp.ID]
    
        # Print updated effluent concentrations and flow rates for verification
        print(f"Updated effluent dilute stream concentrations: {eff_dc.imol}")
        print(f"Updated effluent dilute stream flow rate: {eff_dc.F_vol} m^3/hr")
        print(f"Updated effluent accumulated stream concentrations: {eff_ac.imol}")
        print(f"Updated effluent accumulated stream flow rate: {eff_ac.F_vol} m^3/hr")

#%%
    def _update_dstate(self):
        pass
        # dy = self._dstate
        # eff_dc, eff_ac = self.outs
    
        # # Assuming there are 14 compounds in each stream, update if different
        # n_compounds_dc = len(eff_dc.chemicals.IDs)
        # n_compounds_ac = len(eff_ac.chemicals.IDs)
    
        # # Print to check the actual number of compounds
        # print(f"Number of compounds in dilute stream: {n_compounds_dc}")
        # print(f"Number of compounds in accumulated stream: {n_compounds_ac}")
    
        # # Unpack the state change into effluent dilute and accumulated stream components
        # dC_dc = dy[:n_compounds_dc]
        # dQ_dc = dy[n_compounds_dc]
        # dC_ac = dy[n_compounds_dc + 1:2 * n_compounds_dc + 1]
        # dQ_ac = dy[-1]
    
        # # Initialize dstate if it's None
        # if eff_dc.dstate is None:
        #     eff_dc.dstate = np.zeros_like(dy[:n_compounds_dc + 1])
        # if eff_ac.dstate is None:
        #     eff_ac.dstate = np.zeros_like(dy[n_compounds_dc + 1:])
    
        # # Update effluent dstate for dilute stream
        # eff_dc.dstate[:-1] = dC_dc
        # eff_dc.dstate[-1] = dQ_dc
    
        # # Update effluent dstate for accumulated stream
        # eff_ac.dstate[:-1] = dC_ac
        # eff_ac.dstate[-1] = dQ_ac
    
        # # Total current calculation (assuming `self.j` and `self.A_m` are defined in the class)
        # I = self.j * self.A_m
        # print(f"Total current (I): {I} A")
    
        # for ion, CE in self.CE_dict.items():
        #     ion_index_dc = eff_dc.chemicals.index(ion)
        #     ion_index_ac = eff_ac.chemicals.index(ion)
    
        #     # Moles of target ion transported [mol]
        #     n_T = CE * I / (self.z_T * F) * self.t
        #     print(f"Moles of {ion} transported (n_T): {n_T} mol")
    
        #     # Target ion molar flux [mol/(m2*s)]
        #     J_T = CE * I / (self.z_T * F * self.A_m)
        #     print(f"Target ion molar flux (J_T) for {ion}: {J_T} mol/m^2/s")
    
        #     # Effluent moles for dilute stream [mole]
        #     n_out_dc = (self.ins[0].imol[ion] * self.t / 3600 * 1000) - (self.V * J_T * self.A_m) / self.Q_dc
    
        #     # Effluent moles for accumulated stream [moles]
        #     n_out_ac = (self.V * J_T * self.A_m) / self.Q_dc * (1 - np.exp(-self.Q_ac * self.t / self.V))
    
        #     # Ensure non-negative effluent moles [mole]
        #     n_out_dc = max(n_out_dc, 0)
        #     n_out_ac = max(n_out_ac, 0)
    
        #     # Update effluent streams dstate [kmole/hr]
        #     eff_dc.dstate[ion_index_dc] = n_out_dc / 1000 / self.t * 3600
        #     eff_ac.dstate[ion_index_ac] = n_out_ac / 1000 / self.t * 3600
    
        # # Ensure non-negative effluent moles for H2O [kmole/hr]
        # eff_dc.dstate[-1] = max(self.ins[0].imol['H2O'], 0)
        # eff_ac.dstate[-1] = max(self.ins[1].imol['H2O'], 0)
    
        # # Adjust the mass balance to ensure it matches the influent
        # for comp in self.ins[0].chemicals:
        #     if comp.ID not in self.CE_dict:
        #         comp_index = eff_dc.chemicals.index(comp.ID)
        #         eff_dc.dstate[comp_index] = self.ins[0].imass[comp.ID]
    
        # for comp in self.ins[1].chemicals:
        #     if comp.ID not in self.CE_dict:
        #         comp_index = eff_ac.chemicals.index(comp.ID)
        #         eff_ac.dstate[comp_index] = self.ins[1].imass[comp.ID]
    
        # # Print updated effluent dstate for verification
        # print(f"Updated dstate for dilute stream: {eff_dc.dstate}")
        # print(f"Updated dstate for accumulated stream: {eff_ac.dstate}")
        # print("dState updated successfully.")
#%%
    def _compile_AE(self):
        _state = self._state
        _dstate = self._dstate
        _update_state = self._update_state
        _update_dstate = self._update_dstate
    
        def yt(t, y_ins, dy_ins):
            # Q_ins = y_ins[:, -1]
            # C_ins = y_ins[:, :-1]
            # dQ_ins = dy_ins[:, -1]
            # dC_ins = dy_ins[:, :-1]
    
            # Q = Q_ins.sum()
            # C = (Q_ins[:, np.newaxis] * C_ins).sum(axis=0) / Q
            # _state[-1] = Q
            # _state[:len(C)] = C
            # Q_dot = dQ_ins.sum()
            # C_dot = ((Q_ins[:, np.newaxis] * dC_ins).sum(axis=0) + (dQ_ins[:, np.newaxis] * C_ins).sum(axis=0) - C * Q_dot) / Q
            # _dstate[-1] = Q_dot
            # _dstate[:len(C_dot)] = C_dot
    
            # y_ins == [[Sac_inf_dc, Set_inf_dc,....., Q_inf_dc],
            #           [Sac_inf_ac, Set_inf_ac,....., Q_inf_ac]]   
            n_rows, n_cols = y_ins.shape
            _state[:n_cols] = y_ins[0]
            _state[n_cols:] = y_ins[1]
            _update_state()
            _update_dstate()
    
            print(f"State and dState compiled at time {t}")
            print(f"State derivative at time {t}: {_dstate}")
            
        self._AE = yt
    
    @property
    def AE(self):
        if self._AE is None:
            self._compile_AE()
        return self._AE
#%%    
    def _design(self):
        D = self.design_results
        D['Membrane area'] = self.A_m
        D['Total current'] = self.j * self.A_m
        D['Tank volume'] = self.V
        D['System resistance'] = self.R_sys
        D['System voltage'] = self.V_sys
        D['Power consumption'] = self.P_sys
        print("Design parameters updated successfully.")
    
    def _cost(self):
        self.baseline_purchase_costs['Membrane'] = 100 * self.design_results['Membrane area']  # Assume $100 per m^2 for the membrane
        self.power_utility.consumption = self.design_results['Power consumption'] / 1000  # Assuming kWh consumption based on power
        print("Cost parameters updated successfully.")

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