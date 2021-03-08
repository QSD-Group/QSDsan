#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/master/LICENSE.txt
for license details.

Ref:
    [1] Trimmer et al., Navigating Multidimensional Social–Ecological System
        Trade-Offs across Sanitation Alternatives in an Urban Informal Settlement.
        Environ. Sci. Technol. 2020, 54 (19), 12641–12653.
        https://doi.org/10.1021/acs.est.0c03296.

'''


# %%

import thermosteam as tmo
from qsdsan import Component, Components

__all__ = ('cmps', )

NH3 = Component.from_chemical('NH3', tmo.Chemical('NH3'), measured_as='N',
                              phase='l', particle_size='Soluble',
                              degradability='Undegradable', organic=False)

NonNH3 = Component.from_chemical('NonNH3', tmo.Chemical('N'), measured_as='N',
                                 phase='l', particle_size='Soluble',
                                 degradability='Undegradable', organic=False,
                                 description='Non-NH3 nitrogen')

P = Component.from_chemical('P', tmo.Chemical('P'),
                            phase='l', particle_size='Soluble',
                            degradability='Undegradable', organic=False)

K = Component.from_chemical('K', tmo.Chemical('K'),
                            phase='l', particle_size='Soluble',
                            degradability='Undegradable', organic=False)

Mg = Component.from_chemical('Mg', tmo.Chemical('Mg'),
                             phase='l', particle_size='Soluble',
                             degradability='Undegradable', organic=False)

Ca = Component.from_chemical('Ca', tmo.Chemical('Ca'),
                             phase='l', particle_size='Soluble',
                             degradability='Undegradable', organic=False)

H2O = Component.from_chemical('H2O', tmo.Chemical('H2O'),
                              phase='l', particle_size='Soluble',
                              degradability='Undegradable', organic=False)

OtherSS = Component('OtherSS', phase='l', particle_size='Soluble',
                    degradability='Undegradable', organic=False,
                    description='Unspecified soluble solids')

N2O = Component.from_chemical('N2O', tmo.Chemical('N2O'),
                              phase='g', particle_size='Dissolved gas',
                              degradability='Undegradable', organic=False)

CH4 = Component.from_chemical('CH4', tmo.Chemical('CH4'),
                              phase='g', particle_size='Dissolved gas',
                              degradability='Slowly', organic=True)

# Below three are for combustion reactions
O2 = Component.from_chemical('O2', tmo.Chemical('O2'),
                             phase='g', particle_size='Dissolved gas',
                             degradability='Undegradable', organic=False)

N2 = Component.from_chemical('N2', tmo.Chemical('N2'),
                             phase='g', particle_size='Dissolved gas',
                             degradability='Undegradable', organic=False)

CO2 = Component.from_chemical('CO2', tmo.Chemical('CO2'),
                              phase='g', particle_size='Dissolved gas',
                              degradability='Undegradable', organic=False)

P4O10 = Component.from_chemical('P4O10', tmo.Chemical('P4O10'),
                                phase='s', particle_size='Particulate',
                                degradability='Undegradable', organic=False)

def add_V_from_rho(cmp, rho):
    V_model = tmo.functional.rho_to_V(rho, cmp.MW)
    try: cmp.V.add_model(V_model)
    except:
        handle = getattr(cmp.V, cmp.locked_state)
        handle.add_model(V_model)

Tissue = Component('Tissue', MW=1, phase='s', particle_size='Particulate',
                    degradability='Undegradable', organic=False,
                    description='Tissue for toilet paper')
# 375 kg/m3 is the average of 250-500 for tissue from
# https://paperonweb.com/density.htm (accessed 2020-11-12)
add_V_from_rho(Tissue, 375)

WoodAsh = Component('WoodAsh', MW=1, phase='s', i_Mg=0.0224, i_Ca=0.3034,
                    particle_size='Particulate', degradability='Undegradable',
                    organic=False, description='Wood ash for desiccant')
add_V_from_rho(WoodAsh, 760)

for i in (Tissue, WoodAsh):
    i.copy_models_from(tmo.Chemical('Glucose'), ('Cn', 'mu'))

Struvite = Component.from_chemical('Struvite',
                                   tmo.Chemical('MagnesiumAmmoniumPhosphate'),
                                   formula='NH4MgPO4·H12O6',
                                   phase='s', particle_size='Particulate',
                                   degradability='Undegradable', organic=False)
# http://www.chemspider.com/Chemical-Structure.8396003.html (accessed 2020-11-19)
add_V_from_rho(Struvite, 1711)
    
HAP = Component.from_chemical('HAP', tmo.Chemical('Hydroxyapatite'),
                              phase='s', particle_size='Particulate',
                              degradability='Undegradable', organic=False)
# Taking the average of 3.1-3.2 g/cm3 from
# https://pubchem.ncbi.nlm.nih.gov/compound/Hydroxyapatite (accessed 2020-11-19)
add_V_from_rho(HAP, 3150)

for cmp in (NonNH3, P, K, Mg, Ca, OtherSS):
    cmp.default()
    cmp.copy_models_from(H2O, ('sigma', 'epsilon', 'kappa', 'Cn', 'mu'))
    add_V_from_rho(cmp, 1e3) # assume the same density as water

cmps = Components((NH3, NonNH3, P, K, Mg, Ca, H2O, OtherSS, N2O, CH4, O2, N2,
                   CO2, P4O10, Tissue, WoodAsh, Struvite, HAP))
for i in cmps:
    if i.HHV is None: i.HHV = 0
    if i.LHV is None: i.LHV = 0
    if i.Hf is None: i.Hf = 0
cmps.compile()

cmps.set_synonym('H2O', 'Water')




