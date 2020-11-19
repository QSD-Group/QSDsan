# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 12:22:58 2020

@author: joy_c
"""

import biosteam as bst
import thermosteam as tmo
import os

os.chdir("C:/Users/joy_c/Dropbox/PhD/Research/QSD/codes_developing/QSD-for-WaSH/sanitation")
from sanitation import Component, Components, WasteStream


path = "C:/Users/joy_c/Dropbox/PhD/Research/QSD/codes_developing/QSD-for-WaSH/sanitation/sanitation/default_data/_components.csv"
cmps = Components.load_from_file(file_type='cvs', path=path)

H2O = Component.from_chemical('H2O', tmo.Chemical('H2O'),
                          i_C=0, i_N=0, i_P=0, i_K=0, i_mass=1,
                          i_charge=0, f_BOD5_COD=0, f_uBOD_COD=0,
                          f_Vmass_Totmass=0,
                          particle_size='Soluble',
                          degradability='Undegradable', organic=False)

cmps.append(H2O)

cmps = Components.load_default(default_compile=True)
ws1 = WasteStream.codstates_inf_model('ws1', flow_tot=1e5)
ws2 = WasteStream.codstates_inf_model('ws2', flow_tot=6.34e5, units=('gal/d', 'g/m3'))
ws3 = WasteStream.codbased_inf_model('ws3', flow_tot=1e5)
