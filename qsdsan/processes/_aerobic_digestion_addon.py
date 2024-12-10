# -*- coding: utf-8 -*-
'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''
import numpy as np
from qsdsan import Process
from thermosteam import settings
_load_components = settings.get_default_chemicals

__all__ = ('ASM_AeDigAddOn',)

class ASM_AeDigAddOn(Process):
    '''
    Creates a `Process` object representing the degradation of particulate
    inert organic materials that typically occur in an aerobic digester. 
    Stoichiometry is determined by rules of element conservation in corresponding
    activated sludge models.

    Parameters
    ----------
    k_dig : float, optional
        The 1st-order degradation rate constant, in d^(-1). The default is 0.04.

    See Also
    --------
    :class:`qsdsan.processes.ASM1`
    :class:`qsdsan.processes.ASM2d`
    :class:`qsdsan.processes.mASM2d`

    Examples
    --------
    >>> import qsdsan.processes as pc
    >>> cmps_asm1 = pc.create_asm1_cmps()
    >>> dig_asm1 = pc.ASM_AeDigAddOn('dig_asm1')
    >>> dig_asm1.show()
    Process: dig_asm1
    [stoichiometry]      X_I: -1
                         X_S: 1
                         S_NH: 0.06
                         S_ALK: 0.0514
    [reference]          X_I
    [rate equation]      X_I*k_dig
    [parameters]         k_dig: 0.04
    [dynamic parameters] 
    
    >>> cmps_masm2d = pc.create_masm2d_cmps(set_thermo=False)
    >>> dig_masm2d = pc.ASM_AeDigAddOn('dig_masm2d', components=cmps_masm2d)
    >>> dig_masm2d.show()
    Process: dig_masm2d
    [stoichiometry]      S_NH4: 0.0265
                         S_PO4: 0.0009
                         S_IC: 0.0433
                         X_I: -1
                         X_S: 1
    [reference]          X_I
    [rate equation]      X_I*k_dig
    [parameters]         k_dig: 0.04
    [dynamic parameters]
    
    '''
    
    def __init__(self, ID, k_dig=0.04, components=None):
        cmps = _load_components(components)
        rxn = 'X_I -> X_S'
        consrv = []
        if 'S_ALK' in cmps.IDs: 
            consrv.append('charge')
            rxn += ' + [?]S_ALK'
        elif 'S_IC' in cmps.IDs:
            consrv.append('C')
            rxn += ' + [?]S_IC'
        
        if 'S_NH' in cmps.IDs:
            consrv.append('N')
            rxn += ' + [?]S_NH'
        elif 'S_NH4' in cmps.IDs:
            consrv.append('N')
            rxn += ' + [?]S_NH4'            
        
        if 'S_PO4' in cmps.IDs:
            consrv.append('P')
            rxn += ' +[?]S_PO4'
            
        super().__init__(ID=ID, reaction=rxn,
                         rate_equation='k_dig*X_I',
                         ref_component='X_I',
                         components=cmps,
                         conserved_for=consrv,
                         parameters=('k_dig',))
        self.k_dig=k_dig
        self._stoichiometry = np.asarray(self._stoichiometry, dtype=float)
    
    @property
    def k_dig(self):
        '''[float] Degradation rate constant, in d^(-1).'''
        return self._k
    @k_dig.setter
    def k_dig(self, k):
        self._k = k
        self.set_parameters(k_dig=k)