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
# ADM1_vfa-specific components
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
# I need to adjust later
fc_inf = WasteStream()
fc_inf.set_flow_by_concentration(flow_tot=100, concentrations={'S_su': .5, 'S_la': .5}, units=('L/hr', 'mg/L'))
#fc_eff = WasteStream('FC_Effluent', X_GAO_Gly=.5, H2O=1000, units='kmol/hr')
ac_inf = WasteStream()
ac_inf.set_flow_by_concentration(flow_tot=100, concentrations={'Na+': 100, 'Cl-': 100}, units=('L/hr', 'mg/L'))
#ac_eff = WasteStream('AC_Effluent', X_GAO_Gly=.5, H2O=1000, units='kmol/hr')

#%%
# SanUnit
# Define the Faraday constant (Coulombs per mole)
F = 96485.33212

class ED_vfa(SanUnit):
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 current_density=1.0,  # default current density in A/m2
                 time=3600,  # default time in seconds
                 z_vfa=1,  # valence of VFAs, typically 1
                 membrane_area=1.0, # default membrane area in m2
                 spacer_thickness=1e-3, # default spacer thickness in m
                 electrode_areal_resistance=1.0, # default in ohm
                 conductivity_dc=1.0, # default in ohm
                 conductivity_ac=1.0, # default in ohm
                ):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.current_density = current_density
        self.time = time
        self.z_vfa = z_vfa
        self.membrane_area = membrane_area
        self.spacer_thickness = spacer_thickness
        self.electrode_areal_resistance = electrode_areal_resistance
        self.conductivity_dc = conductivity_dc
        self.conductivity_ac = conductivity_ac
        
        # Assume a bare module factor of 2
        self.F_BM = {'Membrane': 2}

    _N_ins = 1
    _N_outs = 1

    def _run(self):
        fc_inf, ac_inf = self.ins
        fc_eff, ac_eff = self.outs

        # Calculate the concentrations of target ion (e.g., 'S_su') and counter ion (e.g., 'Na+') in feed and accumulated channels
        C_F_T = fc_inf.imass['S_su'] / fc_inf.F_vol
        C_F_X = fc_inf.imass['Na+'] / fc_inf.F_vol
        C_A_X = ac_inf.imass['Na+'] / ac_inf.F_vol

        # Calculate the ion selectivity based on the known concentrations
        self.ion_selectivity = (C_A_X * C_F_T) / (C_F_X * C_A_T)

        # Calculate the final concentration of the target ion in the accumulated channel using ion selectivity formula
        C_A_T = self.ion_selectivity * (C_F_T / C_F_X) * C_A_X
        ac_inf.imass['S_su'] = C_A_T * ac_inf.F_vol

        # Calculate the total ion flux based on the current density and the membrane area
        membrane_area = self.membrane_area  # m2
        current = self.current_density * membrane_area  # A

        # Calculate the theoretical ion transport
        Q_theoretical = current * self.time / (self.z_vfa * F)  # mol

        # Calculate the actual ion transport based on feed and accumulation concentrations
        actual_transport = 0
        for cmp in fc_inf.components:
            if cmp.ID in ['S_su']:
                transported_mass = fc_inf.imass[cmp.ID] - ac_inf.imass[cmp.ID]
                actual_transport += transported_mass / cmp.MW  # mol

        # Calculate charge efficiency
        self.charge_efficiency = actual_transport / Q_theoretical

        # Calculate total target ion transport rate (Theoretical)
        n_T = self.charge_efficiency * current / (self.z_vfa * F)  # mol/s

        # Calculate actual target ion flux
        J_T = self.ion_selectivity * n_T / membrane_area  # mol/m^2/s

        # Update the concentration in the accumulation channel based on ion transport rate
        for cmp in fc_inf.components:
            if cmp.ID in ['S_su', 'Na+']:
                transported_mass = n_T * cmp.MW * 1000  # g/s
                ac_inf.imass[cmp.ID] += transported_mass * self.time  # update accumulation with ion transport rate
                fc_inf.imass[cmp.ID] -= transported_mass * self.time  # decrease feed concentration accordingly

        ac_eff.copy_like(ac_inf)  # output is the updated accumulation stream
        fc_eff.copy_like(fc_inf)  # output is the updated feed stream

    _units = {
        'Membrane Area': 'm2',
        'Electricity': 'kWh'
    }

    def _design(self):
        D = self.design_results
        total_ion_flux = self.current_density * self.time / (self.z_vfa * F)  # mol/s

        # Assume a membrane area based on ion flux
        membrane_area = total_ion_flux * self.time / (1e-6)  # m2, assuming 1 µm/s flux
        D['Membrane Area'] = membrane_area
        return D

    def _cost(self):
        self.baseline_purchase_costs['Membrane'] = \
            100 * self.design_results['Membrane Area']  # assume $100 per m2 membrane cost

        # Assume the electricity usage is proportional to the current and time
        self.power_utility.consumption = self.current_density * self._design()['Membrane Area'] * self.time * 1e-3  # kWh

    @property
    def ion_selectivity(self):
        '''[float] Selectivity of the target ions over counter ions.'''
        return self._ion_selectivity
    @ion_selectivity.setter
    def ion_selectivity(self, i):
        if i <= 0:
            raise AttributeError('`ion_selectivity` must be positive, '
                                f'the provided value {i} is not valid.')
        self._ion_selectivity = i

    @property
    def charge_efficiency(self):
        '''[float] Charge efficiency of the electrodialysis process.'''
        return self._charge_efficiency

    @property
    def current_density(self):
        '''[float] Electric current density applied to the system (A/m2).'''
        return self._current_density
    @current_density.setter
    def current_density(self, i):
        if i <= 0:
            raise AttributeError('`current_density` must be positive, '
                                f'the provided value {i} is not valid.')
        self._current_density = i

    @property
    def time(self):
        '''[float] Time duration of the electrodialysis process (s).'''
        return self._time
    @time.setter
    def time(self, i):
        if i <= 0:
            raise AttributeError('`time` must be positive, '
                                f'the provided value {i} is not valid.')
        self._time = i

    @property
    def z_vfa(self):
        '''[int] Valence of the VFAs.'''
        return self._z_vfa
    @z_vfa.setter
    def z_vfa(self, i):
        if not isinstance(i, int) or i == 0:
            raise AttributeError('`z_vfa` must be a non-zero integer, '
                                f'the provided value {i} is not valid.')
        self._z_vfa = i

# Usage example:
fc_inf = WasteStream('fc_inf', flow_tot=100, concentrations={'S_su': 0.5, 'S_la': 0.5}, units=('L/hr', 'mg/L'))
ac_inf = WasteStream('ac_inf', flow_tot=100, concentrations={'Na+': 100, 'Cl-': 100}, units=('L/hr', 'mg/L'))
fc_eff = WasteStream('fc_eff')
ac_eff = WasteStream('ac_eff')
ed_unit = ED_vfa(ID='ED_vfa', ins=[fc_inf, ac_inf], outs=[fc_eff, ac_eff], current_density=1.0, time=3600)
ed_unit.simulate()
