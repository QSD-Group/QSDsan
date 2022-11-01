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

from math import e, log
from sympy import Function, lambdify, Symbol, plot
from thermosteam import Reaction as Rxn

__all__ = ('KineticReaction',)


class KineticReaction(Rxn):
    r'''
    A general class used to to handle reactions with common kinetics
    where the reaction ate is controlled by the molar concentration of a single component
    (i.e., the rate reactant) so the concenration can be analytically solved.
    
    With a rate constant of k and time of t, for nth-order reaction of
    the rate reactant with a starting molar concentration of :math:`[C]_0`,
    the rate law is given as
    
        .. math:: -\frac{d[C]}{dt} = k*[C]^n
    
    and the integrated rate law is (n!=1)
        
        .. math:: [C]^{1-n} = [C]_0^{1-n} - (1-n)*k*t
    
    when n==1, the integrated rate law is
    
        .. math:: [C]= [C]_0*e^{-k*t}
        
    The linear plot that can be used to determine k is :math:`[C]^g` vs. t
    when n!=1 and :math:`ln([C])` vs. t when n==1.
    
    Half-life of the component is (n!=1)
    
        .. math:: t_{\frac{1}{2}} = \frac{[C]_0^g*(1-2^{n-1})}{(1-n)*k}
        
    when n==1, the half-life is
    
        .. math:: t_{\frac{1}{2}} = \frac{ln(2)}{k}
        
    .. note::
        
        Only single-phase reaction is supported, and it only be used to react a stream
        (i.e., property array cannot be used).
    
    Parameters
    ----------
    rate_reactant : str
        ID of the rate-limiting component of the rate equation (None for 0th order),
        for reactions of second order and up, it only applies for the special situation
        where the rate is controlled by the concentration of a single component.
    n : int
        Order of the kinetic reaction, can be a non-negative integer.
    k : float
        Kinetic rate constant, unit should match the one used for C/C0 and t
        (:math:`\frac{{mol}^g}{unit\\ of\\ t}`).
    t : float
        Time of the reaction for conversion calculation.
    reaction : dict or str
        A dictionary of stoichiometric coefficients or a stoichiometric
        equation written as:
        i1 R1 + ... + in Rn -> j1 P1 + ... + jm Pm
    reactant : str
        ID of reactant component, will be set to the `rate_reactant` if not given.
        
    Examples
    --------
    >>> import qsdsan as qs
    >>> from qsdsan.processes import KineticReaction as KRxn
    
    >>> kwargs = dict(phase='g', particle_size='Dissolved gas',
    ...               degradability='Undegradable', organic=False)
    >>> SO2Cl2 = qs.Component('SO2Cl2', **kwargs)
    >>> SO2 = qs.Component('SO2', **kwargs)
    >>> Cl2 = qs.Component('Cl2', **kwargs)
    >>> qs.set_thermo(qs.Components((SO2Cl2, SO2, Cl2)))
    >>> s1 = qs.Stream('s1', SO2Cl2=100, SO2=10, Cl2=5)
    
    Let's look at the decomposition of sulfuryl chloride (SO2Cl2) to SO2 and Cl2,
    which is a first order equation at 320°C [1]_
    
    >>> rxn = KRxn('SO2Cl2', n=1, k=2.2e-5, t=1e5, reaction='SO2Cl2 -> SO2 +  Cl2')
    >>> # The conversion is 0 at this stage (because we don't konw the concentration of SO2Cl2 yet)
    >>> rxn.show()
    KineticReaction (by mol):
    stoichiometry        reactant    X[%]
    SO2Cl2 -> SO2 + Cl2  SO2Cl2      0.00
    >>> # React the stream
    >>> rxn(s1)
    >>> s1.show()
    Stream: s1
     phase: 'l', T: 298.15 K, P: 101325 Pa
     flow (kmol/hr): SO2Cl2  11.1
                     SO2     98.9
                     Cl2     93.9
    >>> rxn.show() # now we know the conversion
    KineticReaction (by mol):
    stoichiometry        reactant    X[%]
    SO2Cl2 -> SO2 + Cl2  SO2Cl2     88.92
    >>> # You can check the rate equation, integrated rate equation, and half-life
    >>> # of the component
    >>> rxn.rate_equation
    -2.2e-5*C(t) - Derivative(C(t), t)
    >>> rxn.integrated_rate_equation.evalf(n=5) # `evalf` is to limit the digits
    0.035543/2.7183**(2.2e-5*t)
    >>> round(rxn.half_life, 2)
    31506.69
    >>> # You can also look at the conversion over time
    >>> fig = rxn.plot_conversion_over_time()
        
    References
    ----------
    .. [1] https://chem.libretexts.org/Bookshelves/General_Chemistry/Map%3A_General_Chemistry_(Petrucci_et_al.)/14%3A_Chemical_Kinetics/14.05%3A_First-Order_Reactions

    See Also
    --------
    `thermosteam.Reaction <https://biosteam.readthedocs.io/en/latest/API/thermosteam/reaction/Reaction.html>`_
    '''
    _t_sym = Symbol('t')
    _C_sym = Function('C')(_t_sym)
    components = Rxn.chemicals
    reaction_components = Rxn.reaction_chemicals
    reset_chemicals = Rxn.reset_chemicals
    
    def __init__(self, rate_reactant, n, k, t, reaction, reactant=None,
                 C0=None, components=None, **kwargs):
        self.n = n
        self.k = k
        self.t = t
        self.C0 = C0
        reactant = reactant or rate_reactant
        basis = kwargs.pop('basis', None) 
        if basis not in ('mol', None):
            raise ValueError(f'Reaction basis can only be "mol", not {basis}.')
        if kwargs.pop('phases', None): raise ValueError('Only single-phase reaction is supported.')
        Rxn.__init__(
            self=self,
            reaction=reaction,
            reactant=reactant,
            X=0, # just to initialize the Reaction, will be updated later
            chemicals=components,
            basis='mol', 
            phases=None,
            **kwargs)
        # Need to be put after Rxn initialization to have the components set
        self.rate_reactant = rate_reactant
        
    def _calculate_X(self):
        '''Calculate the conversion of the rate reactant.'''
        f = lambdify(self._t_sym, self._X_t, 'math')
        return f(self.t)
    
    @property
    def _X_t(self):
        '''Rate reactant conversion over time.'''
        return 1 - self.integrated_rate_equation/self.C0
        
        
    def __call__(self, stream):
        # Update C0
        rate_reactant = self.rate_reactant
        self.C0 = stream.imol[rate_reactant]/stream.F_vol # [kmol/hr]/[m3/hr]=[kmol/m3]=[mole/L]=M
        # Make the converion 0 for an empty stream
        self._X = 0 if stream.F_mass == 0 else self._calculate_X()
        Rxn.__call__(self, stream)
        
    def plot_conversion_over_time(self, **kwargs):
        '''
        Plot concentrations of the reactants and products of this kinetic reaction.
        All keyword arguments will be passed to :func:`sympy.plot`.
        '''
        ylabel = kwargs.pop('ylabel', 'Conversion')
        plot(self._X_t, (self._t_sym, 0, self.t), ylabel=ylabel)
        
    @property
    def rate_reactant(self):
        '''
        [str] ID of the rate-limiting component of the rate equation
        (None for 0th order).
        When setting the 
        '''
        return self._rate_reactant
    @rate_reactant.setter
    def rate_reactant(self, i):
        i = None if not i else i
        if i not in self.components.IDs:
            raise ValueError(f'The `rate_reactant` "{i}" is not in `components` '
                             'of this `KineticReaction` object.')
        self._rate_reactant = i
        
    @property
    def n(self):
        '''[int] Order of the kinetic reaction, can be a non-negative integer.'''
        return self._n
    @n.setter
    def n(self, i):
        int_i = int(i)
        if not int_i == i:
            raise ValueError('`n` must be a non-negative integer, '
                             f'the provided value {i} is not allowed.')
        if hasattr(self, '_n'):
            if self._n != int_i and hasattr(self, '_C0'):
                if self._C0:
                    try: self._X = self._calculate_X()
                    except AttributeError: pass # have not provided values for k/t
        self._n = int_i
        
    @property
    def k(self):
        r'''
        [float] Kinetic rate constant, unit should match the one used for C/C0 and t
        (:math:`\frac{{mol}^g}{unit\\ of\\ t}`).
        '''
        return self._k
    @k.setter
    def k(self, i):
        if hasattr(self, '_k'):
            if self._k != i and hasattr(self, '_C0'):
                if self._C0:
                    try: self._X = self._calculate_X()
                    except AttributeError: pass # have not provided values for n/t
        self._k = i
        
    @property
    def t(self):
        '''[float] Time of the reaction for conversion calculation.'''
        return self._t
    @t.setter
    def t(self, i):
        if hasattr(self, '_t'):
            if self._t != i and hasattr(self, '_C0'):
                if self._C0:
                    try: self._X = self._calculate_X()
                    except AttributeError: pass # have not provided values for n/k
        self._t = i
        
    @property
    def C0(self):
        '''
        [float] Molar concentration of the rate reactant at t=0,
        needed for conversion calculation except for first-order reaction.
        '''
        return self._C0
    @C0.setter
    def C0(self, i):
        if self.n == 1:
            if not i: self._C0 = None
            else: self._C0 = float(i)
            return
        if i <= 0:
            raise ValueError('`C0` can only be a positive value, '
                              f'the provided value of {i} is not allowed.')
        if hasattr(self, '_C0'):
            if self._C0 != i:
                try: self._X = self._calculate_X()
                except AttributeError: pass # have not provided values for n/k/t
        self._C0 = i
        
    @property
    def rate_equation(self):
        '''Differential rate equation of the rate reactant.'''
        C, t = self._C_sym, self._t_sym
        return -C.diff(t)-self.k*C**self.n # Eq(-C.diff(t), k*C**(self.n))

    @property
    def integrated_rate_equation(self):
        '''Integrated rate equation of the rate reactant.'''
        n, k, C0 = self.n, self.k, self.C0
        t = self._t_sym
        if n ==0: return C0 - k*t
        elif n == 1: return C0*(e**(-k*t))
        elif n == 2: return 1/(1/C0+k*t)
        else:
            g = 1 - n
            return log((C0**g-g*k*t), base=g)
        
    @property
    def half_life(self):
        '''[float] Half-life of the rate reactant.'''
        n, k = self.n, self.k
        if n == 1: return log(2)/k
        C0 = self.C0
        if not C0: raise ValueError('`C0` is not provided, cannot calculate half life.')
        g = 1 - n
        return (C0**g)*(1-2**-g)/(g*k)