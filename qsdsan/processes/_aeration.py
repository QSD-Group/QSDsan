# -*- coding: utf-8 -*-
'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

from math import exp
from .. import Process

__all__ = ('DiffusedAeration',)

class DiffusedAeration(Process):
    """
    Create a :class:`DiffusedAeration` object by specifying either an oxygen
    mass transfer coefficient (KLa or KLa_20), or air flow rate at field
    condition (Q_air) and reactor volume (V). :class:`DiffusedAeration`
    is a subclass of :class:`~.Process`.

    Parameters
    ----------
    ID : str
        A unique identification.
    DO_ID : str
        The component ID of dissolved oxygen (DO) defined in the :class:`~.WasteStream` objects.
    KLa : float, optional
        Oxygen mass transfer coefficient at field condition, [1/d]. A user-defined
        value of KLa supersedes values calculated with other parameters.
    DOsat : float, optional
        Surface DO saturation concentration at field condition, [mg/L]. A user-defined
        value of DOsat supersedes values calculated with other parameters.
    KLa_20 : float, optional
        Oxygen mass transfer coefficient at standard condition (20 degree C,
        clean water, new diffuser), [1/d]. The default is None.
    Q_air : float, optional
        Airflow rate at field conditions, [m^3/d]. The default is None.
    V : float, optional
        Reactor volume, [m^3]. The default is None.
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
    fine_pore : bool, optional
        Whether the diffuser type is fine pore (fine bubble). The default is True.
    elevation : float, optional
        Diffuser elevation, [m]. The default is 0.
    T_air : float, optional
        Air temperature, [K]. The default is 293.15.
    T_water : float, optional
        Water temperature [K]. The default is 293.15.
    SOTE : float, optional
        Standard oxygen transfer efficiency (as a fraction), unitless. The default is 0.3.
    CF1 : float, optional
        Oxygen gas content in standard air, [g O2/m3 air]. The default is 277.6533841.

    Examples
    --------
    >>> from qsdsan import processes as pc, set_thermo
    >>> cmps = pc.load_asm1_cmps()
    >>> set_thermo(cmps)
    >>> aer = pc.DiffusedAeration('aer', 'S_O', KLa=240, DOsat=8.0, V=1333)
    >>> aer.show()
    Process: aer
    [stoichiometry] S_O: 1
    [reference]     S_O
    [rate equation] KLa*(DOsat - S_O)
    [parameters]    KLa: 240
                    DOsat: 8
    """


    A = 10.19621         # Antoine coefficients for water at 1~100 degree C,
    B = 1730.63          # corresponding to temperature in K and pressure in Pa
    C = -39.724
    g = 9.81             # Gravitational acceleration, [m/s^2]
    M_air = 28.964       # Molecular weight of dry air, [g/mol]
    R = 8.314            # Universal gas constant, [J/mol/K]
    P_s = 101325         # Standard sea-level atmospheric pressure, [Pa]

    def __init__(self, ID, DO_ID, V, KLa=None, DOsat=None,
                 KLa_20=None, Q_air=None, DOsat_s20=9.09,
                 alpha=0.6, beta=0.95, F=1.0, theta=1.024,
                 d_submergence=5, fine_pore=True,
                 elevation=0, T_air=293.15, T_water=293.15,
                 SOTE=0.3, CF1=277.6533841):

        super().__init__(ID, reaction={DO_ID:1}, ref_component=DO_ID,
                         rate_equation='KLa*(DOsat - %s)' % DO_ID,
                         conserved_for=(), parameters=('KLa', 'DOsat'))

        self._V = V
        self._KLa_20 = KLa_20
        self._Q_air = Q_air
        self._DOsat_s20 = DOsat_s20
        self._alpha = alpha
        self._beta = beta
        self._F = F
        self._theta = theta
        self._dsub = d_submergence
        self._fine_pore = fine_pore
        self._elevation = elevation
        self._T_air = T_air
        self._T_water = T_water
        self._SOTE = SOTE
        self._CF1 = CF1
        self.KLa = KLa
        self.DOsat = DOsat

    def _calc_DOsat(self):
        return self.tau * self._beta * self.Omega * self.delta * self._DOsat_s20

    def _calc_KLa(self):
        return self._alpha * self._F * self._theta ** (self._T_water - 293.15) * self.KLa_20

    def _calc_KLa_20(self):
        self._SOTR = SOTR = self._Q_air / self.factor * self._SOTE * self._CF1
        return SOTR / (self._V * self.delta * self._DOsat_s20)

    def _calc_Q_air(self):
        self._SOTR = SOTR = self.KLa_20 * self.delta * self._DOsat_s20 * self._V
        return SOTR / self._SOTE / self._CF1 * self.factor

    @property
    def SOTR(self):
        '''[float] Standard oxygen transfer rate, [g/d]'''
        return self._SOTR

    @property
    def factor(self):
        '''[float] Air flowrate correction factor.'''
        return (self._T_air/293.15) * exp(self.g * self.M_air * self._elevation / self.R / self._T_air/1e3)

    @property
    def delta(self):
        '''[float] Depth correction factor, unitless.'''
        if self._fine_pore: return 1 + 0.03858*self._dsub
        else: return 0.99 + 0.0291*self._dsub

    @property
    def tau(self):
        '''[float] Temperature correction factor, unitless.'''
        x = self._T_water - 293.15
        return -7e-6 * x**3 + 4e-4 * x**2 - 1.98e-2 * x + 1.0

    @property
    def Omega(self):
        '''[float] Pressure correction factor, unitless.'''
        P_b = exp(-self.g * self.M_air * self._elevation / self.R / self._T_air/1e3) * self.P_s
        p_v = 10**(self.A - self.B/(self._T_water + self.C))
        p_de = (self.delta-1)*self._dsub*(self.P_s - p_v)
        return (P_b + p_de - p_v)/(self.P_s + p_de - p_v)

    @property
    def KLa_20(self):
        """
        [float] Oxygen mass transfer coefficient at standard condition
        (20 degree C, clean water, new diffuser), [1/d].
        """
        if self._KLa_20: return self._KLa_20
        elif self._Q_air and self._V:
            return self._calc_KLa_20()
        elif self._KLa:
            return self._KLa / (self._alpha * self._F * self._theta**(self._T_water-293.15))
        else: return None

    @KLa_20.setter
    def KLa_20(self, i):
        self._KLa_20 = i
        self._Q_air = None
        self.KLa = None

    @property
    def Q_air(self):
        """[float] Airflow rate at field conditions, [m^3/d]."""
        return self._Q_air or self._calc_Q_air()

    @Q_air.setter
    def Q_air(self, Q):
        self._Q_air = Q
        self._KLa_20 = None
        self.KLa = None

    @property
    def V(self):
        '''[float] Reactor volume, assuming it is filled, [m^3].'''
        return self._V

    @V.setter
    def V(self, i):
        self._V = i
        # self.KLa = None

    @property
    def DOsat_s20(self):
        '''[float] Surface DO saturation concentration at standard condition
        (20 degree C, 1 atm, clean water), [mg/L].'''
        return self._DOsat_s20

    @DOsat_s20.setter
    def DOsat_s20(self, i):
        self._DOsat_s20 = i
        self.KLa = None
        self.DOsat = None

    @property
    def alpha(self):
        '''[float] Wastewater correction factor for KLa_20, unitless.'''
        return self._alpha

    @alpha.setter
    def alpha(self, i):
        self._alpha = i
        self.KLa = None

    @property
    def beta(self):
        '''[float] Wastewater correction factor for DOsat_s20, unitless.'''
        return self._beta

    @beta.setter
    def beta(self, i):
        self._beta = i
        self.DOsat = None

    @property
    def F(self):
        '''[float] Diffuser fouling factor for KLa_20, unitless.'''
        return self._F

    @F.setter
    def F(self, i):
        self._F = i
        self.KLa = None

    @property
    def theta(self):
        '''[float] Temperature correction factor for KLa_20, unitless.'''
        return self._theta

    @theta.setter
    def theta(self, i):
        self._theta = i
        self.KLa = None

    @property
    def d_submergence(self):
        '''[float] Diffuser submergence depth, [m].'''
        return self._dsub

    @d_submergence.setter
    def d_submergence(self, i):
        self._dsub = i
        self.KLa = None
        self.DOsat = None

    @property
    def fine_pore(self):
        '''[bool] Whether the diffuser type is fine pore (fine bubble). '''
        return self._fine_pore

    @fine_pore.setter
    def fine_pore(self, i):
        self._fine_pore = bool(i)
        self.KLa = None
        self.DOsat = None

    @property
    def elevation(self):
        '''[bool] Diffuser elevation, [m]. '''
        return self._elevation

    @elevation.setter
    def elevation(self, i):
        self._elevation = i
        self.KLa = None
        self.DOsat = None

    @property
    def T_air(self):
        '''[bool] Air temperature, [K]. '''
        return self._T_air

    @T_air.setter
    def T_air(self, i):
        self._T_air = i
        self.KLa = None
        self.DOsat = None

    @property
    def T_water(self):
        '''[bool] Water temperature, [K]. '''
        return self._T_water

    @T_water.setter
    def T_water(self, i):
        self._T_water = i
        self.KLa = None
        self.DOsat = None

    @property
    def SOTE(self):
        '''[bool] Standard oxygen transfer efficiency (as a fraction), unitless.'''
        return self._SOTE

    @SOTE.setter
    def SOTE(self, i):
        self._SOTE = i
        self.KLa = None

    @property
    def CF1(self):
        '''[bool] Oxygen gas content in standard air, [g O2/m3 air].'''
        return self._CF1

    @CF1.setter
    def CF1(self, i):
        self._CF1 = i
        self.KLa = None

    @property
    def KLa(self):
        """
        [float] Oxygen mass transfer coefficient at field condition, corrected for
        water quality, diffuser fouling, temperature, [1/d].
        """
        return self._KLa

    @KLa.setter
    def KLa(self, KLa):
        self._KLa = KLa or self._calc_KLa()
        self.set_parameters(KLa=self._KLa)

    @property
    def DOsat(self):
        """
        [float] DO saturation concentration at field condition, corrected for water
        quality, temperature, atmospheric pressure, diffuser submergence depth, [mg/L].
        """
        return self._DOsat
    @DOsat.setter
    def DOsat(self, DOsat):
        self._DOsat = DOsat or self._calc_DOsat()
        self.set_parameters(DOsat = self._DOsat)