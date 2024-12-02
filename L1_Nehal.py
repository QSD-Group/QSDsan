#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 24 11:14:48 2024

@author: saumitrarai
"""

import qsdsan as qs
from qsdsan import processes as pc, sanunits as su

Q_domestic = 38000 # m3/day
Temp = 273.15+20 # temperature [K]
Q_was = 378.54 # [m3/day]

biomass_IDs = ('X_H', 'X_AUT', 'X_PAO')

domestic_ww = {
   'S_I': 20,
   'X_I': 40,
   'S_F': 45,
   'S_A': 63,
   'X_S': 160,
   'S_NH4': 25,
   'S_PO4': 4.5,
   'X_PP': 0,
   'X_PHA': 10,
   'X_H': 10,
   'X_AUT': 10, 
   'X_PAO': 5, 
   'X_MeOH': 32, 
   'S_ALK':7*12,
    }

default_asm2d_kwargs = dict(iN_SI=0.01, iN_SF=0.03, iN_XI=0.02, iN_XS=0.04, iN_BM=0.07,
            iP_SI=0.0, iP_SF=0.01, iP_XI=0.01, iP_XS=0.01, iP_BM=0.02,
            iTSS_XI=0.75, iTSS_XS=0.75, iTSS_BM=0.9,
            f_SI=0.0, Y_H=0.625, f_XI_H=0.1,
            Y_PAO=0.625, Y_PO4=0.4, Y_PHA=0.2, f_XI_PAO=0.1,
            Y_A=0.24, f_XI_AUT=0.1,
            K_h=3.0, eta_NO3=0.6, eta_fe=0.4, K_O2=0.2, K_NO3=0.5, K_X=0.1,
            mu_H=6.0, q_fe=3.0, eta_NO3_H=0.8, b_H=0.4, K_O2_H=0.2, K_F=4.0,
            K_fe=4.0, K_A_H=4.0, K_NO3_H=0.5, K_NH4_H=0.05, K_P_H=0.01, K_ALK_H=0.1,
            q_PHA=3.0, q_PP=1.5, mu_PAO=1.0, eta_NO3_PAO=0.6, b_PAO=0.2, b_PP=0.2,
            b_PHA=0.2, K_O2_PAO=0.2, K_NO3_PAO=0.5, K_A_PAO=4.0, K_NH4_PAO=0.05,
            K_PS=0.2, K_P_PAO=0.01, K_ALK_PAO=0.1,
            K_PP=0.01, K_MAX=0.34, K_IPP=0.02, K_PHA=0.01,
            mu_AUT=1.0, b_AUT=0.15, K_O2_AUT=0.5, K_NH4_AUT=1.0, K_ALK_AUT=0.5, K_P_AUT=0.01,
            k_PRE=1.0, k_RED=0.6, K_ALK_PRE=0.5, 
            )

cmps = pc.create_asm2d_cmps(False)
cmps.compile()
qs.set_thermo(cmps)
thermo_asm2d = qs.get_thermo()

dom_ww = qs.WasteStream('domestic_wastewater', T=Temp)
dom_ww.set_flow_by_concentration(Q_domestic, 
                                     concentrations=domestic_ww, 
                                     units=('m3/d', 'mg/L'))

effluent = qs.WasteStream('effluent', T=Temp)
sludge = qs.WasteStream('sludge', T=Temp)

asm_kwargs = default_asm2d_kwargs
asm2d = pc.ASM2d()

PC = su.PrimaryClarifier('PC', ins=dom_ww,
                              outs=('effluent', 'sludge'), isdynamic=True,
                              init_with='WasteStream', thermo=thermo_asm2d,
                              sludge_flow_rate=Q_was, 
                              solids_removal_efficiency=0.6)

cmps_adm1 = qs.processes.create_adm1_p_extension_cmps()
cmps_adm1.X_PAO.i_N = 0.07
cmps_adm1.X_PAO.i_P = 0.02
cmps_adm1.refresh_constants()
thermo_adm1 = qs.get_thermo()
adm1 = qs.processes.ADM1_p_extension()

J1 = su.ASM2dtomADM1('J1', upstream=[PC-1], thermo=thermo_adm1, isdynamic=True, 
                         adm1_model=adm1, asm2d_model=asm2d)

DG = su.AnaerobicCSTR(ID='DG', ins = J1.outs[0], outs=['gas', 'sludge_DG'],
                        model=adm1, thermo = thermo_adm1)
DG.algebraic_h2 = True

sys = qs.System('AD_trial', path=[PC, J1, DG])
sys.set_tolerance(rmol=1e-6)
sys.maxiter = 500
sys.simulate(t_span=(0, 0.001))
sys.diagram()

# print("--------------------Influent wastewater properties---------------------")
# dom_ww.show()

# print("--------------------Effluent wastewater properties---------------------")
# effluent.show()
   
# print("---------------------------Untreated Sludge properties------------------------------")
# sludge.show()

# print("---------------------------Digested sludge properties------------------------------")
# DG.outs[1].show()

# print("---------------------------Biogas properties------------------------------")
# DG.outs[0].show()
