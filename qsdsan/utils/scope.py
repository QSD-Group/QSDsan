# -*- coding: utf-8 -*-
'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Joy Zhang <joycheung1994@gmail.com>

Part of this module is based on the Thermosteam package:
https://github.com/BioSTEAMDevelopmentGroup/thermosteam

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

from biosteam.utils import Scope
import matplotlib.pyplot as plt


__all__ = ('WasteStreamScope', 'SanUnitScope')

#!!! TODO: add link to "see also"
class WasteStreamScope(Scope):    
    """
    A tracker of the dynamic component concentrations and volumetric flowrates
    of a :class:`WasteStream`.

    Parameters
    ----------
    wastestream : :class:`WasteStream`
        The `WasteStream` object to keep track of during dynamic simulation. Must
        have a `state` attribute

    See Also
    --------
    `biosteam.utils.Scope <>`
    `WasteStream <https://qsdsan.readthedocs.io/en/latest/streams.html#id1>`
    """
    def __init__(self, wastestream):
        names = [f'{ID} [mg/L]' for ID in wastestream.components.IDs] + ['Q [m3/d]']
        header = list(zip([wastestream.ID]*len(names), names))
        super().__init__(wastestream, ('state',), header=header)
    
    def __repr__(self):
        return f'<WasteStreamScope: {self.subject.ID}>'
    
    def plot_time_series(self, state_var=()):
        """Plot the time series data of specified state variables."""
        state_var = [state_var] if isinstance(state_var, str) else state_var
        cmps_idx = self.subject.components.index
        ids = [-1 if var == 'Q' else cmps_idx(var) for var in state_var]
        t = self.time_series
        ys = self.record.T
        fig, ax = plt.subplots(figsize=(8, 4.5))
        for i, var in zip(ids, state_var):
            ax.plot(t, ys[i], '-o', label=var)
        ax.legend(loc='best')
        ylabel = 'Concentration [mg/L] or Flowrate [m3/d]' if 'Q' in state_var else 'Concentration [mg/L]'
        ax.set(xlabel='Time [d]', ylabel=ylabel)
        return fig, ax
        
class SanUnitScope(Scope):
    """
    A tracker of the dynamic state variables of a :class:`SanUnit`.

    Parameters
    ----------
    unit : :class:`SanUnit`
        The `SanUnit` object to keep track of during dynamic simulation. Must
        have a '_state' attribute.

    See Also
    --------
    `biosteam.utils.Scope <>`
    `WasteStream <https://qsdsan.readthedocs.io/en/latest/streams.html#id1>`
    """    
    def __init__(self, unit):
        names = unit._state_header
        header = list(zip([unit.ID]*len(names), names))
        super().__init__(unit, ('_state',), header=header)

    def __repr__(self):
        return f'<SanUnitScope: {self.subject.ID}>'
    
    def plot_time_series(self, state_var=()):
        state_var = [state_var] if isinstance(state_var, str) else state_var
        idx_label = [[i, header[1]] for i, header in enumerate(self.header) \
                     if header[1].split(' ')[0] in state_var]
        t = self.time_series
        ys = self.record.T
        fig, ax = plt.subplots(figsize=(8, 4.5))
        for i, label in idx_label:
            ax.plot(t, ys[i], '-o', label=label)
        ax.legend(loc='best')
        ax.set(xlabel='Time [d]', ylabel='State Variable')
        return fig, ax