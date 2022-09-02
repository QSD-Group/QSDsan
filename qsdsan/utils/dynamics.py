# -*- coding: utf-8 -*-
'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

from .loading import load_data
from scipy.interpolate import InterpolatedUnivariateSpline, CubicSpline, interp1d
import matplotlib.pyplot as plt
import numpy as np


__all__ = ('ExogenousDynamicVariable', )

class ExogenousDynamicVariable:
    """
    Creates an exogenously dynamic variable by interpolating time-series data
    or providing a function of time.

    Parameters
    ----------
    ID : str
        Unique identifier of the variable.
    t : array_like
        A 1-D array of time data, must be sorted.
    y : array_like
        A 1-D array of the variable values at corresponding time points, length
        must match the t array.
    function : callable, optional
        A function that returns the variable value at a given time. Ignored if `t` and `y`
        are provided.
    interpolator : str or int or callable, optional
        Interpolation method to use. It can be a string (e.g., 'slinear', 'quadratic', 'cubic')
        or an integer within [1,5] to specify the order of a spline interpolation.  
        Other strings will be passed on to `scipy.interpolate.interp1d` as the 'kind'
        argument. It can also be a class in `scipy.interpolate` that takes time-series 
        data as input upon initiation. Interpolant that is not at least 
        first-order differentiable is not recommended (e.g., linear interpolation).
        The default is `scipy.interpolate.CubicSpline`. 
    derivative_approximator : callable, optional
        A function that returns derivative of the variable at given time. If none specified,
        will use the `.derivative()` function of (if available) of the interpolant.
    intpl_kwargs : dict, optional
        Keyword arguments for initiating the interpolant.
    
    .. note::
        Extrapolation is allowed, assuming the time-series data are periodic and continuous 
        at both bounds.

    See Also
    --------
    `scipy.interpolate.CubicSpline <https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.CubicSpline.html>`
    `scipy.interpolate.InterpolatedUnivariateSpline <https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.InterpolatedUnivariateSpline.html>`
    `scipy.interpolate.interp1d <https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.interp1d.html>`
    """
    
    def __init__(self, ID, t=None, y=None, function=None, interpolator=None, 
                 derivative_approximator=None, intpl_kwargs={}):
        self._ID = ID
        self.t_data = t
        self.y_data = y
        self._intpl_kwargs = intpl_kwargs
        self.interpolator = interpolator
        self.derivative_approximator = derivative_approximator
        
        if t is not None and y is not None:
            if y[-1] != y[0]:
                t = np.append(t, 2*t[-1] - t[-2])
                y = np.append(y, y[0])
            self._t_end = t[-1]
            self._f = self._intpl(t, y, **self._intpl_kwargs)
        elif hasattr(function, '__call__'):
            self._t_end = None
            self._f = function
        
        if self._func_dydt is None and hasattr(self._f, 'derivative'):
            self._df = self._f.derivative()
        elif hasattr(self._func_dydt, '__call__'):
            self._df = self._func_dydt
        else:
            self._df = lambda t: 0
    
    @property
    def ID(self):
        return self._ID
    
    @property
    def interpolator(self):
        return self._intpl

    @interpolator.setter
    def interpolator(self, i):
        if self.t_data is None: self._intpl = None
        else:
            isa = isinstance
            if i is None:
                self._intpl = CubicSpline
                self._intpl_kwargs.update({'bc_type':'periodic'})
            elif isa(i, (str, int)):
                if i in ('slinear', 'quadratic', 'cubic', 1, 2, 3, 4, 5):
                    self._intpl = InterpolatedUnivariateSpline
                    if isa(i, str):
                        if i == 'slinear': i = 1
                        elif i == 'quadratic': i = 2
                        else: i = 3
                    self._intpl_kwargs.update({'k':i})
                else:
                    self._intpl = interp1d
                    self._intpl_kwargs.update({'kind': i})
            elif isa(i, type) and hasattr(i, '__call__'):
                self._intpl = i
            else:
                raise TypeError(f'interpolator must be None, or a str, or an int,'
                                f'or a class that takes in the data points and returns'
                                f'a callable interpolator (e.g., from scipy.interpolate), '
                                f'not {type(i)}.')
            
    @property
    def derivative_approximator(self):
        return self._func_dydt

    @derivative_approximator.setter
    def derivative_approximator(self, f):
        if f is None or hasattr(f, '__call__'):
            self._func_dydt = f
        else:
            raise TypeError(f'derivative approximator must be None or a callable that takes in '
                            f'time and returns derivatives as a 1d numpy.array, not {type(f)}.')
            
    def __call__(self, t):
        '''Evaluates the variable at time t.'''
        if self._t_end: t %= self._t_end
        return float(self._f(t))
    
    def derivative(self, t):
        '''Returns the derivative of the variable at time t'''
        if self._t_end: t %= self._t_end
        return float(self._df(t))
    
    def plot(self, t_start=0, t_end=None):
        fig, ax = plt.subplots(figsize=(8, 4.5))
        t, y = self.t_data, self.y_data
        if t is not None: n_eval = len(t)*5
        else: n_eval = 50
        _t_end = t_end or self._t_end or 10
        t_eval = np.linspace(t_start, _t_end, n_eval)
        y_eval = self._f(t_eval)
        if t is not None:
            ax.plot(t, y, 'o', label='data')
            ax.plot(t_eval, y_eval, '-', label='interpolation')
            ax.legend(loc='best')
        else:
            ax.plot(t_eval, y_eval, '-')
        ax.set(xlabel='Time', ylabel=self.ID)
        return fig, ax
        
    @classmethod
    def batch_init(cls, data, interpolator=None, derivative_approximator=None, 
                   intpl_kwargs={}, **load_data_kwargs):
        '''
        Creates a list of :class:`ExogenousDynamicVariable` objects from time-series data.

        Parameters
        ----------
        data : str or :class:`pandas.DataFrame`
            Time series data of multiple variables. If provided as a file,
            file extension should be one of (".cvs", ".xls", "xlsx").
        load_data_kwargs : optional
            Additional keyword arguments passed to `qsdsan.utils.load_data`.
        '''
        if isinstance(data, str):
            idx_col = load_data_kwargs.pop('index_col', None)
            load_data_kwargs['index_col'] = idx_col
            df = load_data(data, **load_data_kwargs)
        else:
            df = data
        df.sort_values('t', inplace=True)
        dct_y = df.loc[:, df.columns!='t'].to_dict('series')
        t = np.asarray(df.t)
        return [cls(k, t, np.asarray(v), interpolator=interpolator, 
                    derivative_approximator=derivative_approximator,
                    intpl_kwargs=intpl_kwargs) \
                for k, v in dct_y.items()]
        
    def __repr__(self):
        return f"<{type(self).__name__}: {self._ID}>"