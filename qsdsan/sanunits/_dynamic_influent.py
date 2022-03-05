# -*- coding: utf-8 -*-
'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

from .. import SanUnit
from ..utils import ospath, load_data, data_path
from scipy.interpolate import InterpolatedUnivariateSpline, CubicSpline, interp1d
import numpy as np

__all__ = ('DynamicInfluent',)
dynamic_inf_path = ospath.join(data_path, 'sanunit_data/_inf_dry_2006.tsv')

class DynamicInfluent(SanUnit):
    """
    A fake SanUnit to generate a dynamic :class:`WasteStream ` at its outlet 
    by interpolating time-series data.

    Parameters
    ----------
    ID : str, optional
        ID for the dynamic influent generator. The default is ''.
    outs : :class:`WasteStream`
        Dynamic influent.
    data_file : str, optional
        The file path for the time-series data. Acceptable file extensions are
        `.xlsx`, `.xls`, `.csv`, `.tsv`. If none specified, will load the default
        time-series data of dry-weather influent with components from ASM1.
    interpolator : str or int or callable, optional
        Interpolator to use. It can be a string (e.g., 'slinear', 'quadratic', 'cubic')
        or an integer within [1,5] to specify the order of a spline interpolation.  
        Other strings will be passed on to `scipy.interpolate.interp1d` as the 'kind'
        argument. It can also be a class in `scipy.interpolate` that takes time-series 
        data as input upon initiation. Interpolant that is not at least 
        first-order differentiable is not recommended (e.g., linear 
        interpolation).The default is `scipy.interpolate.CubicSpline`. 
    derivative_approximator : callable, optional
        Function that returns derivative of state at given time. If none specified,
        will use the `.derivative` method (if available) of the interpolator.
    load_data_kwargs : dict, optional
        Keyword arguments for loading the data file with the `qsdsan.utils.load_data()`
        function. The default is {}.
    intpl_kwargs : dict, optional
        Keyword arguments for initiating the interpolant. The default is {}.

    See Also
    --------
    `scipy.interpolate.CubicSpline <https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.CubicSpline.html>`
    `scipy.interpolate.InterpolatedUnivariateSpline <https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.InterpolatedUnivariateSpline.html>`
    `scipy.interpolate.interp1d <https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.interp1d.html>`
    """
    _N_ins = 0
    _N_outs = 1
    
    def __init__(self, ID='', ins=None, outs=(), data_file=None, interpolator=None, 
                 derivative_approximator=None, thermo=None, init_with='WasteStream', 
                 isdynamic=True, load_data_kwargs={}, intpl_kwargs={}, **kwargs):
        SanUnit.__init__(self, ID, None, outs, thermo, init_with, isdynamic=isdynamic)
        self._intpl_kwargs = intpl_kwargs
        self.interpolator = interpolator
        self.derivative_approximator = derivative_approximator
        self._init_from_file(data_file, **load_data_kwargs)
        for attr, value in kwargs.items():
            setattr(self, attr, value)

    @property
    def interpolator(self):
        return self._intpl

    @interpolator.setter
    def interpolator(self, i):
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
        if (f is None and hasattr(self._intpl, 'derivative')) or hasattr(f, '__call__'):
            self._func_dydt = f
        else:
            raise TypeError(f'derivative approximator must be None (i.e., inferred from '
                            f'a differentiable interpolant), or a callable that takes in '
                            f'time and returns derivatives as a 1d numpy.array, not {type(f)}.')
            
    def _init_from_file(self, file_path=None, **kwargs):
        path = file_path or dynamic_inf_path
        self._data = df = load_data(path, index_col=None, **kwargs)
        df.sort_values('t', inplace=True)
        df_y = df.loc[:, df.columns!='t']
        if any(df_y.iloc[-1,:] != df_y.iloc[0,:]):
            y_end = df_y.iloc[0,:].to_dict()
            y_end['t'] = 2*df.t.iloc[-1] - df.t.iloc[-2]
            df.append(y_end, ignore_index=True)
        self._t_end = df.t.iloc[-1]
        intpl = self._intpl
        ikwargs = self._intpl_kwargs
        y_IDs = self.components.IDs + ('Q',)
        diff_set = set(df.columns) - set(y_IDs) - {'t'}
        if diff_set: 
            raise RuntimeError(f'The data file contains state variable(s) that are'
                               f'inconsistent with the thermo: {diff_set}')
        self._interpolant = [intpl(df.t, df.loc[:,y], **ikwargs) \
                             if y in df.columns else lambda t: 0 \
                             for y in y_IDs]
        if self._func_dydt is None:
            self._derivative = [i.derivative() if hasattr(i, 'derivative') else lambda t: 0 \
                                for i in self._interpolant]
        else: self._derivative = None

    def interpolant(self, t):
        t %= self._t_end
        return np.array([f(t) for f in self._interpolant])
    
    def derivative(self, t):
        if self._derivative is None: return self.derivative_approximator(t)
        else: 
            t %= self._t_end
            return np.array([f(t) for f in self._derivative])

    def _init_state(self):
        self._state = self.interpolant(0)
        self._dstate = self.derivative(0)

    def _update_state(self):
        self._outs[0].state = self._state

    def _update_dstate(self):
        self._outs[0].dstate = self._dstate

    def _run(self):
        '''Only to converge volumetric flows.'''
        out, = self._outs
        y0 = self._data.iloc[0,:].to_dict()
        Q = y0.pop('Q')
        y0.pop('t', None)
        y0.pop('H2O', None)
        out.set_flow_by_concentration(flow_tot=Q, concentrations=y0, units=('m3/d', 'mg/L'))

    @property
    def AE(self):
        if self._AE is None:
            self._compile_AE()
        return self._AE

    def _compile_AE(self):
        _state = self._state
        _dstate = self._dstate
        f = self.interpolant
        f_dy = self.derivative
        _update_state = self._update_state
        _update_dstate = self._update_dstate      

        def yt(t, QC_ins, dQC_ins):
            _state[:] = f(t)
            _update_state()
            _dstate[:] = f_dy(t)
            _update_dstate()

        self._AE = yt
