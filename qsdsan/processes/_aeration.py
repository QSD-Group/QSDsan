# -*- coding: utf-8 -*-
'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Joy Cheung <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

from math import exp
from .. import Process

__all__ = ('DiffusedAeration',)

class DiffusedAeration(Process):
    """
    Create a :class:`DiffusedAeration` object by specifying oxygen mass transfer coefficient 
    at standard condition or air flow rate at field condition. :class:`DiffusedAeration` is a
    subclass of :class:`~.Process`. A :class:`DiffusedAeration` object is capable of increasing
    the oxygen flow rates of a :class:`~.WasteStream` object.

    Parameters
    ----------
    ID : str
        A unique identification.
    DO_ID : str
        The component ID of dissolved oxygen (DO) defined in the :class:`~.WasteStream` objects.
    KLa_20 : float, optional
        Oxygen mass transfer coefficient at standard condition (20 degree C, 
        clean water, new diffuser), [1/d]. The default is None.
    DOsat_s20 : float, optional
        Surface DO saturation concentration at standard condition (20 degree C, 
        1 atm, clean water), [mg/L]. The default is 9.09.
    alpha : float, optional
        Wastewater correction factor for KLa_20, unitless. The default is 0.6.
    beta : float, optional
        Wastewater correction factor for DOsat_s20, unitless. The default is 0.95.
    F : float, optional
        Diffuser fouling factor for KLa_20, unitless. The default is 1.0.
    theta : float, optional
        Temperature correction factor for KLa_20, unitless. The default is 1.024.
    d_submergence : float, optional
        Diffuser submergence depth, [m]. The default is 5.
    V : float, optional
        Reactor volume, [m^3]. The default is None.
    fine_pore : bool, optional
        Whether the diffuser type is fine pore (fine bubble). The default is True.
    elevation : float, optional
        Diffuser elevation, [m]. The default is 0.
    T_air : float, optional
        Air temperature, [K]. The default is 293.15.
    T_water : float, optional
        Water temperature [K]. The default is 293.15.
    Q_air : float, optional
        Airflow rate at field conditions, [m^3/d]. The default is None.
    SOTE : float, optional
        Standard oxygen transfer efficiency (as a fraction), unitless. The default is None.
    CF1 : float, optional
        Oxygen gas content in standard air, [g O2/m3 air]. The default is 277.6533841.    
    """
    
    
    A = 10.19621         # Antoine coefficients for water at 1~100 degree C, 
    B = 1730.63          # corresponding to temperature in K and pressure in Pa
    C = -39.724  
    g = 9.81             # Gravitational acceleration, [m/s^2]
    M_air = 28.964       # Molecular weight of dry air, [g/mol]
    R = 8.314            # Universal gas constant, [J/mol/K]
    P_s = 101325         # Standard sea-level atmospheric pressure, [Pa]

    def __init__(self, ID, DO_ID, KLa_20=None, DOsat_s20=9.09, 
                 alpha=0.6, beta=0.95, F=1.0, theta=1.024, 
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
        P_b = exp(-self.g * self.M_air * self._elevation / self.R / self._T_air/1e3) * self.P_s
        p_v = 10**(self.A - self.B/(self._T_water + self.C))
        p_de = (self._delta-1)*self._d_submergence*(self.P_s - p_v)
        return (P_b + p_de - p_v)/(self.P_s + p_de - p_v)
        
    
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
        Q_air_s = self._Q_air * (293.15/self._T_air) * exp(-self.g * self.M_air * self._elevation / self.R / self._T_air/1e3)
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
        """
        [float] DO saturation concentration at field condition, corrected for water 
        quality, temperature, atmospheric pressure, diffuser submergence depth, [mg/L].
        """
        return self._DOsat
    
    @property
    def KLa(self):
        """
        [float] Oxygen mass transfer coefficient at field condition, corrected for
        water quality, diffuser fouling, temperature, [1/d].
        """
        return self._KLa

    @property
    def Q_air(self):
        """[float] Airflow rate at field conditions, [m^3/d]."""
        if self._Q_air: return self._Q_air
        elif self._KLa_20 and self._SOTE and self._V:
            SOTR = self._KLa_20 * self._delta * self._DOsat_s20 * self._V
            Q_air_s = SOTR / self._SOTE /self._CF1
            return Q_air_s * (self._T_air/293.15) / exp(-self.g * self.M_air * self._elevation / self.R / self._T_air/1e3)
        else: return None

    @property
    def KLa_20(self):
        """
        [float] Oxygen mass transfer coefficient at standard condition 
        (20 degree C, clean water, new diffuser), [1/d].
        """
        if self._KLa_20: return self._KLa_20
        elif self._Q_air and self._SOTE and self._V:
            return self._calc_KLa_20()
        else: return None
