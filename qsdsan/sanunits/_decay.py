#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Yalin Li <mailto.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''


# %%

from math import exp
from .. import SanUnit

__all__ = ('Decay',)


class Decay:
    '''For non-steady state degradation.'''

    # Put these as default class attributes
    _t0 = 0
    _tau = None
    _COD_removal = None
    _COD_max_decay = None
    _decay_k_COD = None
    _MCF_decay = None
    _max_CH4_emission = None
    _N_max_decay = None
    _decay_k_N = None
    _N2O_EF_decay = None


    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', F_BM_default=1,
                 degraded_components=('OtherSS',)):
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo,
                         init_with=init_with, F_BM_default=F_BM_default)
        self.degraded_components = degraded_components


    # TODOs: enable COD calc by decay; enable sludge calc with moisture; enable biogas or not
    def _first_order_run(self, degraded_components=('OtherSS',),
                         if_capture_biogas=False, if_N2O_emission=True):
        '''
        A generic :func:`_run` considering first order decay of
        organics and N, including allocation of N to to the ammonia
        and non-ammonia part, as well as the emission of fugitive CH4 and N2O.

        .. note::

            COD degradation will be with `COD_removal` is used if it is given
            (this removal will be treated as the final removal);
            alternatively, `COD_max_decay` is used to calculate the
            final removal with the retention time tau.

        Parameters
        ----------
        degraded_components : tuple
            IDs of components that will degrade (at the same removal as `COD_removal`).
        if_capture_biogas : bool
            If True, self.outs[1] would be biogas;
            if False, the generated biogas will be considered as fugitive emission.
        if_N2O_emission : bool
            Whether to consider N2O emission.
        tau : float
            Retention time of the unit.
        COD_removal : float
            Removal of organics as chemical oxygen demand (COD).
        COD_max_decay : float
            Maximum fraction of COD removed during storage given sufficient time.
        decay_k_COD : float
            Rate constant for COD degradation.
        MCF_decay : float
            Methane correction factor for COD decay,
            [fraction of anaerobic conversion of degraded COD].
        max_CH4_emission : float
            Maximum methane emission, [g CH4/g COD].
        N_max_decay : float
            Maximum degradation of N.
        decay_k_N : float
            First order rate constant for N degradation.
        N2O_EF_decay : float
            N2O emission factor, [fraction of degraded N].
        '''

        waste = self.ins[0]
        if if_capture_biogas:
            treated, biogas, CH4, N2O = self.outs
        else:
            treated, CH4, N2O = self.outs
            biogas = CH4.copy() # this stream doesn't matter
        treated.copy_like(waste)
        biogas.phase = CH4.phase = N2O.phase = 'g'

        # COD removal
        _COD = waste._COD or waste.COD
        COD_deg = _COD*treated.F_vol/1e3*self.COD_removal # kg/hr
        treated._COD *= (1-self.COD_removal)

        degraded_components = degraded_components or self.degraded_components
        treated.imass[degraded_components] *= (1-self.COD_removal)

        CH4_prcd = COD_deg*self.MCF_decay*self.max_CH4_emission
        if if_capture_biogas:
            biogas.imass['CH4'] = CH4_prcd
            CH4.empty()
        else:
            CH4.imass['CH4'] = CH4_prcd
            biogas.empty()

        if if_N2O_emission:
            N_loss = self.first_order_decay(k=self.decay_k_N,
                                            t=self.tau/365,
                                            max_decay=self.N_max_decay)
            N_loss_tot = N_loss*waste.TN/1e3*waste.F_vol
            NH3_rmd, NonNH3_rmd = \
                self.allocate_N_removal(N_loss_tot, waste.imass['NH3'])
            treated.imass ['NH3'] = waste.imass['NH3'] - NH3_rmd
            treated.imass['NonNH3'] = waste.imass['NonNH3'] - NonNH3_rmd
            N2O.imass['N2O'] = N_loss_tot*self.N2O_EF_decay*44/28
        else:
            N2O.empty()



    @staticmethod
    def allocate_N_removal(tot_red, preferred_N):
        '''
        Allocate the total amount of N removal to NH3 and non-NH3 components.

        Parameters
        ----------
        tot_red : float
            Total amount of N to be removed.
        preferred_N : float
            Current content of the N that will be removed first.

        Returns
        -------
        N removal: tuple
            Amount of preferred N to be removed, amount of other N to be removed.
        '''
        if not preferred_N > 0:
            return 0, tot_red
        elif preferred_N > tot_red:
            return tot_red, 0
        else:
            return preferred_N, tot_red-preferred_N


    def first_order_decay(self, k, t, max_decay, t0=None, tot=1):
        r'''
        Calculate first-order degradation loss based on
        `Trimmer et al. <https://doi.org/10.1021/acs.est.0c03296>`_.

        .. math:: C_{deg} = tot * max_{decay}
        .. math:: t_f = t_0 + t
        .. math:: C_{avg} = \frac{C_deg}{k*t} * (e^{-k*t_0}-e^{-k*t_f})
        .. math:: loss = C_{deg} - C_{avg}

        Parameters
        ----------
        k : float
            Degradation rate constant.
        t0 : float
            Degradation time prior to current process.
        t : float
            Degradation time in current process (i.e., tf-t0).
        max_decay : float
            Maximum removal ratio.
        tot : float, optional
            Total degradable amount.
            If set to 1 (default), the return is the relative ratio (i.e., loss/tot).

        Returns
        -------
        loss : float
            Amount lost due to degradation.

        References
        ----------
        [1] Trimmer et al., Navigating Multidimensional Social–Ecological System
        Trade-Offs across Sanitation Alternatives in an Urban Informal Settlement.
        Environ. Sci. Technol. 2020, 54 (19), 12641–12653.
        https://doi.org/10.1021/acs.est.0c03296.
        '''
        t0 = self.t0 if not t0 else t0
        tf = t0 + t
        Cdeg = tot * max_decay
        Cavg = Cdeg/(k*t) * (exp(-k*t0)-exp(-k*tf))
        loss = Cdeg - Cavg
        return loss


    @property
    def t0(self):
        '''[float] Degradation time prior to current process.'''
        return self._t0
    @t0.setter
    def t0(self, i):
        self._t0 = i

    @property
    def tau(self):
        '''[float] Retention time of the unit.'''
        return self._tau
    @tau.setter
    def tau(self, i):
        self._tau = i

    @property
    def COD_removal(self):
        '''[float] Final COD removal.'''
        return self._COD_removal
    @COD_removal.setter
    def COD_removal(self, i):
        self._COD_removal = i

    @property
    def COD_max_decay(self):
        '''[float] Maximum fraction of COD removed during storage given sufficient time.'''
        return self._COD_max_decay
    @COD_max_decay.setter
    def COD_max_decay(self, i):
        self._COD_max_decay = i

    @property
    def decay_k_COD(self):
        '''[float] Rate constant for COD decay.'''
        return self._decay_k_COD
    @decay_k_COD.setter
    def decay_k_COD(self, i):
        self._decay_k_COD = i

    @property
    def MCF_decay(self):
        '''[float] Methane correction factor for COD decay.'''
        return self._MCF_decay
    @MCF_decay.setter
    def MCF_decay(self, i):
        self._MCF_decay = i

    @property
    def max_CH4_emission(self):
        '''[float] Maximum methane emission as a fraction of degraded COD, [g CH4/g COD].'''
        return self._max_CH4_emission
    @max_CH4_emission.setter
    def max_CH4_emission(self, i):
        self._max_CH4_emission = i

    @property
    def N_max_decay(self):
        '''[float] Maximum fraction of N degraded through denitrification during storage given sufficient time.'''
        return self._N_max_decay
    @N_max_decay.setter
    def N_max_decay(self, i):
        self._N_max_decay = i

    @property
    def decay_k_N(self):
        '''[float] Rate constant for N decay.'''
        return self._decay_k_N
    @decay_k_N.setter
    def decay_k_N(self, i):
        self._decay_k_N = i

    @property
    def N2O_EF_decay(self):
        '''[float] N2O emission factor, [fraction of degraded N].'''
        return self._N2O_EF_decay
    @N2O_EF_decay.setter
    def N2O_EF_decay(self, i):
        self._N2O_EF_decay = i