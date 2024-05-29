# -*- coding: utf-8 -*-
"""
Created on Fri May 24 08:18:43 2024

@author: joy_c
"""
import numpy as np
from numba import njit, cfunc


N = 6
V = np.array([249, 2313, 2366, 2366, 4982, 4982])
denom = np.diag(1/V)
f_in = np.array([[0, 0.8, 0.2, 0,0,0],
                 [1, 0, 0, 0, 0, 0]])
internal_recycles = [(5, 2, 1.5e5)] # from, to, Q
Q_internal = np.zeros(N)
for i_from, i_to, qr in internal_recycles:
    Q_internal[i_to: i_from] += qr
    
DO = [0,0,0,0,2.,2.]
kLa = [0]*6
DOsat = 8.0

#%%
from exposan import bsm1
bsm1.load()
cmps = bsm1.PE.components
asm = bsm1.A1.suspended_growth_model
M_stoi = asm.stoichio_eval()
f_rho = asm.rate_function
DO_id = cmps.index('S_O')

#%%
@njit
def dydt(t, QC_ins, QC, dQC_ins=None):
    y = QC.reshape((N, len(QC)/N))
    Cs = y[:,:-1]
    if any(DO):
        for i in range(N):
            if DO[i] > 0: Cs[i, DO_id] = DO[i]
    elif any(kLa):
        do = Cs[:, DO_id]
        aer = kLa*(DOsat-do)
        aer[aer < 0] = 0.
    Qs = np.dot(QC_ins[:,-1], f_in).cumsum() + Q_internal
    M_ins = f_in.T @ np.diag(QC_ins[:,-1]) @ QC_ins[:,:-1]  # N * n_cmps
    for i_from, i_to, qr in internal_recycles:
        M_ins[i_to] += Cs[i_from] * qr
    # rhos = np.apply_along_axis(f_rho, 1, Cs) # n_zone * n_process
    rhos = np.vstack([f_rho(c) for c in Cs])
    rxn = rhos @ M_stoi
    dy = np.zeros_like(y)
    dy[:,:-1] = denom @ (M_ins - np.diag(Qs) @ Cs) + rxn
    if any(DO): dy[:,DO_id] = 0.
    elif any(kLa): dy[:,DO_id] += aer 
    return dy.flatten()

#%%
import qsdsan.sanunits as su
from qsdsan import System
s = bsm1.sys.flowsheet.stream
u = bsm1.sys.flowsheet.unit
inf = s.wastewater
inf.unlink()
ras = s.RAS
ras.unlink()
s.treated.unlink()
AS = su.PFR('AS', ins=[inf, ras], outs=0-u.C1, V_tanks=[1000]*2+[1333]*3, 
            influent_fractions=[[1.0, 0,0,0,0], [1.0, 0,0,0,0]],
            internal_recycles=[(4,0,55338)], kLa=[0,0,240,240,84], DO_ID='S_O',
            suspended_growth_model=u.A1._model)
AS.set_init_conc(
    S_I=30, S_S=5, X_I=1000, X_S=100, X_BH=500, X_BA=100, X_P=100, S_O=2, S_NO=20,
    S_NH=2, S_ND=1, X_ND=1, S_ALK=84
    )
sys = System('sys', path=(AS, u.C1))

#%%
sys.simulate(t_span=(0,50), method='BDF')
