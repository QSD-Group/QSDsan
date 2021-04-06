# -*- coding: utf-8 -*-
'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Joy Cheung <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/master/LICENSE.txt
for license details.
'''

from math import exp
from . import Process

A = 10.19621
B = 1730.63
C = -39.724  # Antoine coefficients for water at 1~100 degree C, converted to SI units
g = 9.81
M_air = 28.964
R = 8.314
P_s = 101325    


__all__ = ('DiffusedAeration',)

#TODO: complete docstring
class DiffusedAeration(Process):
    """
    Create a :class:`DiffusedAeration` object by specifying oxygen mass transfer coefficient 
    at standard condition or air flow rate at field condition. :class:`DiffusedAeration` is a
    subclass of :class:`Process`. A :class:`DiffusedAeration` object is capable of increasing
    the oxygen flow rates of a :class:`WasteStream` object.

    Parameters
    ----------
    ID : str
        A unique identification.
    DO_ID : str
        The component ID of oxygen defined in the :class:`WasteStream` objects.
    KLa_20: float


    
    """
    def __init__(self, ID, DO_ID, KLa_20=None, DOsat_s20=9.09, 
                 alpha=0.8, beta=0.95, F=1.0, theta=1.024, 
                 d_submergence=5, V=None, fine_pore=True, 
                 elevation=0, T_air=293.15, T_water=293.15, 
                 Q_air=None, SOTE=None, CF1=277.6533841):

        super().__init__(self, ID, reaction={DO_ID: 1}, ref_component=DO_ID,
                         conserved_for=(), parameters=('KLa', 'DOsat'),
                         rate_equation='KLa*(DOsat - %s)' % DO_ID)
        self._KLa_20 = KLa_20
        self._DOsat_s20 = DOsat_s20
        self._alpha = alpha
        self._beta = beta
        self._F = F
        self._theta = theta
        self._d_submergence = d_submergence
        self._V = V
        self._fine_pore = fine_pore
        self._elevation = elevation
        self._T_air = T_air
        self._T_water = T_water
        self._Q_air = Q_air
        self._SOTE = SOTE
        self._CF1 = CF1
        self.update()
            
    def _calc_delta(self):
        if self._fine_pore: self._delta = 1 + 0.03858*self._d_submergence
        else: self._delta = 0.99 + 0.0291*self._d_submergence
    
    def _calc_omega(self):
        P_b = exp(-g*M_air*self._elevation/R/self._T_air/1e3) * P_s
        p_v = 10**(A - B/(self._T_water+C))
        p_de = (self._delta-1)*self._d_submergence*(P_s - p_v)
        return (P_b + p_de - p_v)/(P_s + p_de - p_v)
        
    
    def _calc_tau(self):
        x = self._T_water - 273.15
        return -7e-6 * x**3 + 8e-4 * x**2 - 4.35e-2 * x + 1.6052
    
    def _calc_DOsat(self):
        Omega = self._calc_omega()
        tau = self._calc_tau()
        self._DOsat = tau * self._beta * Omega * self._delta * self._DOsat_s20
    
    def _calc_KLa(self, KLa_20=None):
        self._KLa = self._alpha * self._F * self._theta ** (self._T_water - 293.15) * (KLa_20 or self._KLa_20)
    
    def _calc_KLa_20(self):
        Q_air_s = self._Q_air * (293.15/self._T_air) * exp(-g*M_air*self._elevation/R/self._T_air/1e3)
        KLa_20 = Q_air_s * self._SOTE * self._CF1 / (self._V * self._delta * self._DOsat_s20)
        return KLa_20

    def update(self):
        self._calc_delta()
        self._calc_DOsat()
        if self._Q_air and self._SOTE and self._V:
            KLa_20 = self._calc_KLa_20()
            self._calc_KLa(KLa_20)
        elif self._KLa_20:
            self._calc_KLa()
        else: 
            self._KLa = None
        self.set_parameters(KLa=self._KLa, DOsat = self._DOsat)
    
    @property
    def DOsat(self):
        return self._DOsat
    
    @property
    def KLa(self):
        return self._KLa

    @property
    def Q_air(self):
        if self._Q_air: return self._Q_air
        elif self._KLa_20 and self._SOTE and self._V:
            SOTR = self._KLa_20 * self._delta * self._DOsat_s20 * self._V
            Q_air_s = SOTR / self._SOTE /self._CF1
            return Q_air_s * (self._T_air/293.15) / exp(-g*M_air*self._elevation/R/self._T_air/1e3)
        else: return None

    @property
    def KLa_20(self):
        if self._KLa_20: return self._KLa_20
        elif self._Q_air and self._SOTE and self._V:
            return self._calc_KLa_20()
        else: return None
