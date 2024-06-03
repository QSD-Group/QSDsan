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
from qsdsan import Component, Components, WasteStream, SanUnit, Process, Processes, CompiledProcesses
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
#WasteStream
#Same as ADM1 Effluent Q
inf_dc = WasteStream()
inf_dc.set_flow_by_concentration(flow_tot=7, concentrations={'S_la': .5, 'S_bu': .5, 'S_he': .5}, units=('m3/d', 'mg/L'))
#fc_eff = WasteStream('FC_Effluent', X_GAO_Gly=.5, H2O=1000, units='kmol/hr')
inf_ac = WasteStream()
inf_ac.set_flow_by_concentration(flow_tot=7, concentrations={'Na+': 100, 'Cl-': 100}, units=('m3/d', 'mg/L'))
#ac_eff = WasteStream('AC_Effluent', X_GAO_Gly=.5, H2O=1000, units='kmol/hr')

#%%
# SanUnit
class ED_vfa(SanUnit):
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 CE=0.5,  # Charge efficiency
                 j=500,   # Current density in A/m^2
                 t=3600,  # Time in seconds
                 z_T=1,   # Valence of the target ion
                 F=96485.33289,  # Faraday's constant in C/mol
                 A_m=10,  # Membrane area in m^2
                 V=1.0  # Volume of all tanks in m^3
                ):
        super().__init__(ID, ins, outs, thermo, init_with)
        self.CE = CE
        self.j = j
        self.t = t
        self.z_T = z_T
        self.F = F
        self.A_m = A_m
        self.V = V

    _N_ins = 2
    _N_outs = 2

    def _run(self):
        inf_dc, inf_ac = self.ins
        eff_dc, eff_ac = self.outs

        # Calculate total current
        I = self.j * self.A_m
        
        # Moles of target ion transported
        n_T = self.CE * I / (self.z_T * self.F) * self.t
        
        # Target ion molar flux
        J_T = self.CE * I / (self.z_T * self.F * self.A_m)
        
        # Obtain the flow rate from the influent dilute stream
        Q = inf_dc.F_vol  # Flow rate from influent in m3/hr
        
        # Dilute stream concentration
        C_inf_dc = inf_dc.conc[self.z_T]
        C_out_dc = C_inf_dc - (J_T * self.A_m) / Q
        
        # Accumulated stream concentration
        C_out_ac = (J_T * self.A_m) / Q * (1 - np.exp(-Q * self.t / self.V))
        
        # Update effluent streams
        eff_dc.copy_like(inf_dc)
        eff_dc.conc[self.z_T] = C_out_dc
        
        eff_ac.copy_like(inf_ac)
        eff_ac.conc[self.z_T] = C_out_ac

    _units = {
        'Membrane area': 'm2',
        'Total current': 'A',
        'Accumulated volume': 'm3'
    }

    def _design(self):
        D = self.design_results
        D['Membrane area'] = self.A_m
        D['Total current'] = self.j * self.A_m
        D['Accumulated volume'] = self.V_accum

    def _cost(self):
        self.baseline_purchase_costs['Membrane'] = 100 * self.design_results['Membrane area']  # Assume $100 per m2 for the membrane
        self.power_utility.consumption = self.design_results['Total current'] * self.t / 3600  # Assuming kWh consumption based on current and time

#%%
# Initialize the ED_vfa unit
ed_vfa_unit = ED_vfa(
    ID='ED1',
    ins=(inf_dc, inf_ac),
    outs=(WasteStream(), WasteStream()),
    CE=0.85,
    j=500,
    t=3600,
    z_T=1,
    F=96485.33289,
    A_m=10,
    V=1.0
)
#%%
# Run the unit
ed_vfa_unit.simulate()
ed_vfa_unit.show()
#%%
# Access and display the results
print("Effluent dilute stream concentration:", ed_vfa_unit.outs[0].conc)
print("Effluent accumulated stream concentration:", ed_vfa_unit.outs[1].conc)