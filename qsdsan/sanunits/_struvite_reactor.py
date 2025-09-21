

from qsdsan._sanunit import SanUnit
from qsdsan.sanunits._abstract import Splitter
from qsdsan._waste_stream import WasteStream
from qsdsan.equipments import Electrode, Machine, Membrane

from math import ceil
import numpy as np


__all__ = ('SC',)
#%%
class SC(SanUnit):
    """
    Struvite reactor (simplified):
      - Removes phosphate by precipitation with Mg2+ and NH3 to form struvite (MgNH4PO4·6H2O)
      - Recovered outlet [0]: solid struvite
      - Loss outlet [1]: LIQUID (unrecovered removed P)
      - Effluent outlet [2]: LIQUID remainder

    Notes:
      * No crystallization-water bookkeeping; we do NOT adjust free H2O anywhere.
      * Struvite formation uses 1:1:1 molar ratios among removed P, available Mg2+, and NH3.

    Parameters (key)
    ---------------
    removal_P : float in [0,1]
        Target fraction of influent Phosphate removed (not in effluent).
    component_ID_P : str   (default 'S_PO4')
    component_ID_Mg: str   (default 'S_Mg')
    component_ID_NH3: str  (default 'S_NH4')
    component_ID_struvite: str (default 'X_struv')
    """

    _N_ins = 1
    _N_outs = 3
    _ins_size_is_fixed = True
    _outs_size_is_fixed = True

    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                 component_ID_NH3='S_NH4', component_ID_P='S_PO4',
                 component_ID_Mg='S_Mg', component_ID_struvite='X_struv',
                 removal_P=0.90,
                 precip_yield=1.0,
                 order=None, init_with='WasteStream',
                 F_BM_default=None, isdynamic=False, **kwargs):

        super().__init__(ID=ID, ins=ins, outs=outs, thermo=thermo,
                         init_with=init_with, isdynamic=isdynamic,
                         F_BM_default=F_BM_default, **kwargs)

        # IDs / knobs
        self._component_ID_NH3 = component_ID_NH3
        self._component_ID_P   = component_ID_P
        self._component_ID_Mg  = component_ID_Mg
        self._component_ID_struvite = component_ID_struvite

        # phosphate removal target
        self.removal_P = float(removal_P)
        self.precip_yield = float(precip_yield)

        # dynamic caches for split matrices
        self._recovery_matrix  = None
        self._loss_matrix      = None
        self._remaining_matrix = None

    # ---------- convenience ----------
    @property
    def component_ID_NH3(self): return self._component_ID_NH3
    @property
    def component_ID_P(self):   return self._component_ID_P
    @property
    def component_ID_Mg(self):  return self._component_ID_Mg
    @property
    def component_ID_struvite(self): return self._component_ID_struvite

    def _compute_struvite_splits(self, influent):
        """
        Compute per-component split fractions (of INFLUENT) to:
          - recovered (struvite solid),
          - loss (liquid),
          - effluent (liquid).
        Returns (recovery_matrix, loss_matrix, remaining_matrix, struvite_moles)
        aligned to influent.components order.
        """
        cmps = influent.components
        IDs  = cmps.IDs
        n    = len(IDs)
    
        # Start with rec, los = 0; rem = 1 for pass-through by default
        rec = np.zeros(n, dtype=float)
        los = np.zeros(n, dtype=float)
        rem = np.ones(n,  dtype=float)
    
        # Short-hands
        P_id   = self.component_ID_P
        Mg_id  = self.component_ID_Mg
        NH3_id = self.component_ID_NH3
        # STR_id = self.component_ID_struvite  # not needed here, product added later
    
        # Convert to moles for limiting-reagent logic
        imol_in = influent.imol
        mol_P   = float(imol_in[P_id])
        mol_Mg  = float(imol_in[Mg_id])
        mol_NH3 = float(imol_in[NH3_id])
    
        # Target moles of P to remove from liquid
        target_P_removed = max(0.0, min(self.removal_P, 1.0)) * mol_P
    
        # Of that removed P, user wants only a fraction to precipitate
        desired_precip_P = max(0.0, min(self.precip_yield, 1.0)) * target_P_removed
    
        # Actual struvite limited by Mg and NH3 availability (1:1:1 stoichiometry)
        struvite_moles = min(desired_precip_P, mol_Mg, mol_NH3)
    
        # Fractions to RECOVERED (as struvite) for P, Mg, NH3
        rec_frac_P   = (struvite_moles / mol_P)   if mol_P   > 0 else 0.0
        rec_frac_Mg  = (struvite_moles / mol_Mg)  if mol_Mg  > 0 else 0.0
        rec_frac_NH3 = (struvite_moles / mol_NH3) if mol_NH3 > 0 else 0.0
    
        # Any "removed P" not precipitated goes to LOSS (liquid)
        extra_P_removed = max(0.0, target_P_removed - struvite_moles)
        loss_frac_P     = (extra_P_removed / mol_P) if mol_P > 0 else 0.0
    
        # Build per-component splits
        for i, cid in enumerate(IDs):
            if cid == P_id:
                rec[i] = rec_frac_P
                los[i] = loss_frac_P
                rem[i] = 1.0 - rec[i] - los[i]
            elif cid == Mg_id:
                rec[i] = rec_frac_Mg
                los[i] = 0.0
                rem[i] = 1.0 - rec[i]
            elif cid == NH3_id:
                rec[i] = rec_frac_NH3
                los[i] = 0.0
                rem[i] = 1.0 - rec[i]
            else:
                # unaffected: already rem[i] = 1.0 from initialization
                pass
    
        # clip for safety
        rec = np.clip(rec, 0.0, 1.0)
        los = np.clip(los, 0.0, 1.0)
        rem = np.clip(rem, 0.0, 1.0)
    
        return rec, los, rem, struvite_moles



    # ---------- steady run ----------
    def _run(self):
        influent, = self.ins
        recovered, loss, effluent = self.outs

        recovered.phase = 's'
        loss.phase = effluent.phase = 'l'

        recM, losM, remM, struvite_mol = self._compute_struvite_splits(influent)

        mass_in = influent.mass.copy()
        recovered.mass = mass_in * recM
        loss.mass      = mass_in * losM
        effluent.mass  = mass_in * remM

        # Add struvite solid (no free-water adjustments)
        recovered.imol[self.component_ID_struvite] += struvite_mol

        # Cache matrices for dynamics
        self._recovery_matrix  = recM
        self._loss_matrix      = losM
        self._remaining_matrix = remM

    # ---------- matrices for dynamics ----------
    @property
    def recovery_matrix(self):  return self._recovery_matrix
    @property
    def loss_matrix(self):      return self._loss_matrix
    @property
    def remaining_matrix(self): return self._remaining_matrix

    # ---------- dynamics scaffolding ----------
    @property
    def state(self):
        if self._state is None: return None
        return dict(zip(list(self.components.IDs)+['Q'], self._state))

    def _init_state(self):
        influent = self.ins[0]
        self._state  = np.append(influent.mass.copy()*24*1e3, influent.F_vol*24)  # g/d, m3/d
        self._dstate = self._state * 0.

    def _update_state(self):
        arr = self._state
        for ws in self.outs:
            if ws.state is None: ws.state = np.zeros_like(arr)

        self._outs[0].state = np.zeros_like(arr)
        self._outs[0].state[:-1] = (self.recovery_matrix if self.recovery_matrix is not None else 0.0) * arr[:-1]
        self._outs[0].state[-1]  = 1

        self._outs[1].state[:-1] = (self.loss_matrix if self.loss_matrix is not None else 0.0) * arr[:-1]
        self._outs[1].state[-1]  = 1

        remM = self.remaining_matrix if self.remaining_matrix is not None else 1.0
        self._outs[2].state[:-1] = remM * arr[:-1] / arr[-1]   # mg/L
        self._outs[2].state[-1]  = arr[-1]                    # m3/d

    def _update_dstate(self):
        arr = self._dstate
        for ws in self.outs:
            if ws.dstate is None: ws.dstate = np.zeros_like(arr)

        recM = self.recovery_matrix if self.recovery_matrix is not None else 0.0
        losM = self.loss_matrix     if self.loss_matrix     is not None else 0.0
        remM = self.remaining_matrix if self.remaining_matrix is not None else 1.0

        self._outs[0].dstate[:-1] = recM * arr[:-1]
        self._outs[0].dstate[-1]  = 0

        self._outs[1].dstate[:-1] = losM * arr[:-1]
        self._outs[1].dstate[-1]  = 0

        Q  = self._outs[2].state[-1]
        C  = self._outs[2].state[:-1]
        dM = arr[:-1]
        dQ = arr[-1]
        self._outs[2].dstate[:-1] = remM * ((dM*Q - dQ*C) / (Q**2))
        self._outs[2].dstate[-1]  = dQ

    @property
    def AE(self):
        if self._AE is None: self._compile_AE()
        return self._AE

    def _compile_AE(self):
        _state = self._state
        _dstate = self._dstate
        _update_state = self._update_state
        _update_dstate = self._update_dstate

        def yt(t, QC_ins, dQC_ins):
            Q  = QC_ins[0][-1];   C = QC_ins[0][:-1]
            dQ = dQC_ins[0][-1]; dC = dQC_ins[0][:-1]
            _state[-1]  = Q
            _state[:-1] = C * Q
            _dstate[-1]  = dQ
            _dstate[:-1] = dC*Q + C*dQ
            _update_state()
            _update_dstate()
        self._AE = yt