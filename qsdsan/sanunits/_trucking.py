#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

# %%

from warnings import warn
from .. import currency, SanUnit, ImpactItem, Transportation
from ..utils import auom, copy_attr

__all__ = ('Trucking',)


class Trucking(SanUnit):
    '''
    For transportation of materials with considerations on material loss
    based on `Trimmer et al. <https://doi.org/10.1021/acs.est.0c03296>`_

    Parameters
    ----------
    load_type : str
        Either 'mass' or 'volume'.
    load : float
        Transportation load per trip.
    load_unit : str
        Unit of the load.
    distance : float
        Transportation distance per trip.
    distance_unit : float
        Unit of the distance.
    interval : float
        Time interval between trips.
    fee : float
        Transportation fee per trip.
    fee_unit : str
        Currency of the fee.
    if_material_loss : bool
        If material loss occurs during transportation.
    loss_ratio : float or dict
        Fractions of material losses due to transportation.

    Examples
    --------
    `bwaise systems <https://github.com/QSD-Group/EXPOsan/blob/main/exposan/bwaise/systems.py>`_

    References
    ----------
    [1] Trimmer et al., Navigating Multidimensional Social–Ecological System
    Trade-Offs across Sanitation Alternatives in an Urban Informal Settlement.
    Environ. Sci. Technol. 2020, 54 (19), 12641–12653.
    https://doi.org/10.1021/acs.est.0c03296.
    '''

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 load_type='mass', load=1., load_unit='kg',
                 distance=1., distance_unit='km',
                 interval=1., interval_unit='hr',
                 fee=0., fee_unit=currency,
                 if_material_loss=True, loss_ratio=0.02):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.single_truck = single = \
            Transportation('single_truck', item='Trucking',
                           load_type=load_type, load=load, load_unit=load_unit,
                           distance=distance, distance_unit=distance_unit,
                           interval=interval, interval_unit=interval_unit)
        total = single.copy('total_truck')
        self.transportation = (total,)
        self._update_fee(fee, fee_unit)
        self.if_material_loss = if_material_loss
        self.loss_ratio = loss_ratio

    _N_ins = 1
    _N_outs = 2

    def _update_fee(self, fee=0., unit=''):
        if not unit or unit == currency:
            self._fee = float(fee)
        else:
            converted = auom(unit).convert(float(fee), currency)
            self._fee = converted

    def _run(self):
        transported, loss = self.outs
        transported.copy_like(self.ins[0])
        loss.empty()
        if self.if_material_loss:
            if self._loss_ratio_type == 'float':
                loss.copy_like(transported)
                transported.mass *= 1 - self.loss_ratio
                loss.mass = self.ins[0].mass - transported.mass
            else:
                for cmp, ratio in self.loss_ratio.items():
                    transported.imass[cmp] *= 1 - ratio
                    loss.imass[cmp] = self.ins[0].imass[cmp] - transported.imass[cmp]

    def _design(self):
        single = self.single_truck
        single.item = ImpactItem.get_item('Trucking') # in case it has been updated
        if single.load_type == 'volume':
            factor = auom('m3').conversion_factor(single.default_units['load'])
            N = self.F_vol_in*factor*single.interval/single.load
        else:
            factor = auom('kg').conversion_factor(single.default_units['load'])
            N = self.F_mass_in*factor*single.interval/single.load
        self.design_results['Parallel trucks'] = N

        total, = self.transportation
        copy_attr(total, single, skip=('_ID',)) # in case attributes have been updated
        total.load = single.load * N
        self._add_OPEX = {'Total fee': self.fee/total.interval*N}

    @property
    def fee(self):
        '''[float] Transportation fee per trip.'''
        return self._fee
    @fee.setter
    def fee(self, fee, unit=''):
        self._update_fee(fee, unit)

    @property
    def loss_ratio(self):
        '''
        [float] or [dict] Fractions of material losses due to transportation.
        If a single number is provided, then it is assumed that losses of
        all Components in the WasteStream are the same.

        .. note::

            Pre-set state variable values in :class:`~.WasteStream`
            (e.g., COD set through `_COD` rather than calculated)
            will be retained if the loss ratio is a single number
            (treated like the loss stream is split from the original stream),
            but will be overwritten when the ratio is a dict.

        '''
        return self._loss_ratio
    @loss_ratio.setter
    def loss_ratio(self, i):
        if not self.if_material_loss:
            msg = f'`if_material_loss` is False, the set value {i} is ignored.'
            warn(msg, source=self)
        try:
            self._loss_ratio = float(i)
            self._loss_ratio_type = 'float'
        except:
            if isinstance(i, dict):
                self._loss_ratio = i
                self._loss_ratio_type = 'dict'
            else:
                raise TypeError(f'Only float or dict allowed, not {type(i).__name__}.')