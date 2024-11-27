#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 24 11:14:48 2024

@author: saumitrarai
"""

import qsdsan as qs
import numpy as np
from qsdsan import processes as pc, sanunits as su

from qsdsan.utils import (
    time_printer
    )

__all__ = (
    'biomass_IDs',
    'create_system',
    'default_asm2d_kwargs', 
    'steady_state_ad_init_conds',
    'domestic_ww', 'Q_domestic', 'Q_brewery',
    'Q_ras', 'Q_was', 'Temp', 'V_ae',
    )
# %%

# =============================================================================
# Parameters and util functions
# =============================================================================

Q_domestic = 38000      # influent flowrate [m3/d] = 10 MGD (WERF report)

Temp = 273.15+20 # temperature [K]

biomass_IDs = ('X_H', 'X_AUT', 'X_PAO')

domestic_ww = {
   'S_I': 20,
   'X_I': 120,
   'S_F': 45,
   'S_A': 63,
   'X_S': 480,
   'S_NH4': 25,
   'S_PO4': 4.5,
   'X_PP': 100,
   'X_PHA': 100,
   'X_H': 100,
   'X_AUT': 500, 
   'X_PAO': 500, 
   'X_MeOH': 320, 
   'S_ALK':7*12,
    }

cmps = pc.create_asm2d_cmps(False)
cmps.compile()
qs.set_thermo(cmps)
thermo_asm2d = qs.get_thermo()

dom_ww = qs.WasteStream('domestic_wastewater', T=Temp)
dom_ww.set_flow_by_concentration(Q_domestic, 
                                concentrations=domestic_ww, 
                                units=('m3/d', 'mg/L'))
    
effluent = qs.WasteStream('effluent', T=Temp)
WAS = qs.WasteStream('WAS', T=Temp)

C1 = su.PrimaryClarifier('C1', ins=dom_ww, outs=[effluent, WAS], 
                         solids_removal_efficiency=0.6, 
                         sludge_flow_rate=dom_ww.F_vol*24*0.3)
    
sys = qs.System('L1_WERF', path=(C1,))
sys.set_tolerance(rmol=1e-6)
sys.maxiter = 500

sys.set_dynamic_tracker(*sys.products)

if __name__ == '__main__':
    t = 0.001
    method = 'RK23'

# msg = f'Method {method}'
# print(f'\n{msg}\n{"-"*len(msg)}') # long live OCD!
# print(f'Time span 0-{t}d \n')    

sys.simulate()

print("\n")    
print("--------------------Influent wastewater properties---------------------")
dom_ww.show()

print("\n")    
print("--------------------Effluent wastewater properties---------------------")
effluent.show()

print("\n")    
print("---------------------------WAS properties------------------------------")
WAS.show()

sys.diagram()