#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>
    Joy Zhang <joycheung1994@gmail.com>

Part of this module is based on the BioSTEAM package:
https://github.com/BioSTEAMDevelopmentGroup/biosteam

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''


# %%

import numpy as np
from collections import defaultdict
from collections.abc import Iterable
from biosteam.utils.misc import format_title
from . import currency, Unit, Stream, SanStream, WasteStream, \
    Construction, Transportation

__all__ = ('SanUnit',)

def _update_init_with(init_with, ins_or_outs, size):

    if isinstance(init_with, str):
        init_with = dict.fromkeys([f'{ins_or_outs}{n}' for n in range(size)], init_with)
        return init_with

    if len(init_with) == 1 and (('all' or 'All' or 'else' or 'Else') in init_with.keys()):
        init_with = dict.fromkeys([f'{ins_or_outs}{n}' for n in range(size)], init_with)
        return init_with

    if 'else' in init_with.keys():
        new_init_with = dict.fromkeys([f'{ins_or_outs}{n}' for n in range(size)], init_with['else'])
    elif 'Else' in init_with.keys():
        new_init_with = dict.fromkeys([f'{ins_or_outs}{n}' for n in range(size)], init_with['Else'])
    else:
        new_init_with = {}

    for k, v in init_with.items():
        if k in ('Else', 'else'):
            continue
        v_lower = v.lower()
        if v_lower in ('stream', 's'):
            new_init_with[k] = 's'
        elif v_lower in ('sanstream', 'ss'):
            new_init_with[k] = 'ss'
        elif v_lower in ('wastestream', 'ws'):
            new_init_with[k] = 'ws'
        else:
            raise ValueError(f'Stream type for {k} is invalid, ' \
                             'valid values are "Stream" or "s", "SanStream" or "ss", ' \
                             'and "WasteStream" or "ws".')

    return new_init_with

add2list = lambda lst, item: lst.extend(item) if isinstance(item, Iterable) \
    else lst.append(item)

class SanUnit(Unit, isabstract=True):

    '''
    Subclass of :class:`biosteam.Unit`, can be initialized with
    :class:`thermosteam.Stream`, :class:`~.SanStream`, or :class:`~.WasteStream`.

    Parameters
    ----------
    init_with : str or dict
        Which class of stream the :class:`SanUnit` will be initialized with,
        can be "Stream" (shorthanded as "s"), "SanStream" ("ss"), or "WasteStream" ("ws").
        When provided as a str, all streams will be of the same class;
        when provided as a dict, use "ins" or "outs" followed with the order number
        (i.e., ins0, outs-1) as keys; you can use ":" to denote a range (e.g., ins2:4);
        you can also use "else" to specify the stream class for non-provided ones.
    construction : :class:`~.Construction` or iterable
        Contains construction information.
    transportation : :class:`~.Transportation` or iterable
        Contains construction information.
    add_OPEX : float or dict
        Operating expense per hour in addition to utility cost (assuming 100% uptime).
        Float input will be automatically converted to a dict with the key being
        "Additional OPEX".
    uptime_ratio : float
        Uptime of the unit to adjust `add_OPEX`, should be in [0,1]
        (i.e., a unit that is always operating has an uptime_ratio of 1).

        .. note::

            This will not affect the utility (heating, cooling, power) and
            material costs/environmental impacts.

            To account for less utility and material flows, normalize them to
            the `uptime_ratio` of the system.
            For example, if the system operates 100% of time but a pump only works
            50% of the pump at 50 kW. Set the pump `power_utility` to be 50*50%=25 kW.


    equipment_lifetime : int or dict
        Lifetime of this unit (int) or individual equipment within this unit
        (dict) in year.
        It will be used to adjust cost and emission calculation in TEA and LCA.
        Equipment without provided lifetime will be assumed to have the same
        lifetime as the TEA/LCA.
    F_BM_default : float
        If not None, all bare module factors will be default to the set value.

        .. note::

            Regardless of `F_BM_default`, design (F_D), pressure (F_P),
            and material (F_M) factors are all defaulted to 1.

    See Also
    --------
    `biosteam.Unit <https://biosteam.readthedocs.io/en/latest/Unit.html>`_
    `thermosteam.Stream <https://thermosteam.readthedocs.io/en/latest/Stream.html>`_

    '''

    ticket_name = 'SU'

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 construction=(), transportation=(), equipments=(),
                 add_OPEX={}, uptime_ratio=1., lifetime=None, F_BM_default=None,
                 isdynamic=False, **kwargs):

        self._register(ID)
        self._specification = None
        self._load_thermo(thermo)
        self._init_with = init_with
        self._init_ins(ins, init_with)
        self._init_outs(outs, init_with)
        self._init_utils()
        self._init_results()
        self._init_specification()
        self._assert_compatible_property_package()

        self.construction = construction
        self.transportation = transportation
        for equip in equipments:
            equip._linked_unit = self
        self.equipments = equipments

        self.add_OPEX = add_OPEX
        self.uptime_ratio = 1.
        self.lifetime = lifetime

        if F_BM_default:
            F_BM = self.F_BM
            self.F_BM = defaultdict(lambda: F_BM_default)
            self.F_BM.update(F_BM)

        self._isdynamic = isdynamic
        self._state = None
        self._ODE = None

        for attr, val in kwargs.items():
            setattr(self, attr, val)


    def _convert_stream(self, strm_inputs, streams, init_with, ins_or_outs):
        if not streams:
            return []

        if isinstance(init_with, str):
            init_with = dict.fromkeys([f'{ins_or_outs}{n}'
                                       for n in range(len(streams))], init_with)

        # Do not change pre-defined stream types
        if isinstance(strm_inputs, Stream): # input is a stream
            init_with[f'{ins_or_outs}0'] = type(strm_inputs).__name__
         # input is an iterable of stream
        elif isinstance(strm_inputs, Iterable) and not isinstance(strm_inputs, str):
            for n, s in enumerate(strm_inputs):
                if isinstance(s, Stream):
                    init_with[f'{ins_or_outs}{n}'] = type(s).__name__

        init_with = _update_init_with(init_with, ins_or_outs, len(streams))

        converted = []
        for k, v in init_with.items():
            if not ins_or_outs in k:
                continue

            num = k.split(ins_or_outs)[-1]
            s = streams[int(num)]
            if v == 's':
                converted.append(s)
            elif v == 'ss':
                converted.append(SanStream.from_stream(SanStream, s))
            else:
                converted.append(WasteStream.from_stream(WasteStream, s))

        diff = len(converted) - len(streams)
        if diff != 0:
            raise ValueError(f'Type(s) of {diff} stream(s) has/have not been specified.')

        return converted


    def _init_ins(self, ins, init_with):
        super()._init_ins(ins)
        self._ins = self._convert_stream(ins, self.ins, init_with, 'ins')

    def _init_outs(self, outs, init_with):
        super()._init_outs(outs)
        self._outs = self._convert_stream(outs, self.outs, init_with, 'outs')

    def _init_results(self):
        super()._init_results()

    def __repr__(self):
        return f'<{type(self).__name__}: {self.ID}>'

    def _get_stream_info(self, info, ins_or_outs, _stream_info,
                         T, P, flow, composition, N, IDs):
        info += 'ins...\n' if ins_or_outs=='ins' else 'outs...\n'
        i = 0
        for stream in getattr(self, ins_or_outs):
            if not stream:
                info += f'[{i}] {stream}\n'
                i += 1
                continue
            ws_info = stream._wastestream_info() if isinstance(stream, WasteStream) else ''
            if _stream_info:
                stream_info = stream._info(None, T, P, flow, composition, N, IDs) #+ \
                    # '\n' # this breaks the code block in sphinx
                stream_info += ('\n' + ws_info) if ws_info else ''
            else:
                stream_info = stream._wastestream_info()
            su = stream._source if ins_or_outs=='ins' else stream._sink
            index = stream_info.index('\n')
            from_or_to = 'from' if ins_or_outs=='ins' else 'to'
            link_info = f'  {from_or_to}  {type(su).__name__}-{su}\n' if su else '\n'
            info += f'[{i}] {stream.ID}' + link_info + stream_info[index+1:]
            i += 1
        return info


    def _info(self, T, P, flow, composition, N, IDs, _stream_info):
        '''Information of the unit.'''
        if self.ID:
            info = f'{type(self).__name__}: {self.ID}\n'
        else:
            info = f'{type(self).__name__}\n'

        info = self._get_stream_info(info, 'ins', _stream_info,
                                     T, P, flow, composition, N, IDs)

        info = self._get_stream_info(info, 'outs', _stream_info,
                                     T, P, flow, composition, N, IDs)
        info = info.replace('\n ', '\n    ')
        return info[:-1]


    def show(self, T=None, P=None, flow='g/hr', composition=None, N=15, IDs=None, stream_info=True):
        '''Print information of the unit, including waste stream-specific information.'''
        print(self._info(T, P, flow, composition, N, IDs, stream_info))

    def add_equipment_design(self):
        for equip in self.equipments:
            name = equip.name or format_title(type(equip).__name__)
            equip_design = equip._design()
            equip_design = {} if not equip_design else equip_design
            self.design_results.update(equip_design)

            equip_units = {} if not equip.design_units else equip.design_units
            self._units.update(equip_units)
            self.F_BM[name] = equip.F_BM
            if equip.lifetime:
                self._default_equipment_lifetime[name] = equip.lifetime

    def add_equipment_cost(self):
        for equip in self.equipments:
            name = equip.name or format_title(type(equip).__name__)
            self.baseline_purchase_costs[name] = equip._cost()

    def add_construction(self, add_unit=True, add_design=True, add_cost=True,
                         add_lifetime=True):
        '''Batch-adding construction unit, designs, and costs.'''
        for i in self.construction:
            if add_unit:
                self._units[i.item.ID] = i.item.functional_unit
            if add_design:
                self.design_results[i.item.ID] = i.quantity
            if add_cost:
                self.baseline_purchase_costs[i.item.ID] = i.cost
            if add_lifetime and i.lifetime:
                self._default_equipment_lifetime[i.item.ID] = i.lifetime

    @property
    def components(self):
        '''[Components] The :class:`~.Components` object associated with this unit.'''
        return self.chemicals

    @property
    def isdynamic(self):
        '''[bool] Whether the unit runs dynamically within a system.'''
        return self._isdynamic
    @isdynamic.setter
    def isdynamic(self, i):
        self._isdynamic = bool(i)

    def _state_tracer(self):
        states = []
        for inf in self.ins:
            u = inf._source
            state = u._state_locator(u._state)[inf.ID] if u \
                else np.append(inf.conc, inf.get_total_flow('m3/d'))
            states.append(state)
        return np.array(states)

    @property
    def construction(self):
        '''[:class:`~.Construction` or iterable] Contains construction information.'''
        return self._construction
    @construction.setter
    def construction(self, i):
        if isinstance(i, Construction):
            i = (i,)
        else:
            if not isinstance(i, Iterable):
                raise TypeError(
                    f'Only `Construction` object can be included, not {type(i).__name__}.')
            for j in i:
                if not isinstance(j, Construction):
                    raise TypeError(
                        f'Only `Construction` object can be included, not {type(j).__name__}.')
        self._construction = i

    @property
    def transportation(self):
        '''[:class:`~.Transportation` or iterable] Contains transportation information.'''
        return self._transportation
    @transportation.setter
    def transportation(self, i):
        if isinstance(i, Transportation):
            i = (i,)
        else:
            if not isinstance(i, Iterable):
                raise TypeError(
                    f'Only `Transportation` object  can be included, not {type(i).__name__}.')
            for j in i:
                if not isinstance(j, Transportation):
                    raise TypeError(
                        f'Only `Transportation` can be included, not {type(j).__name__}.')
        self._transportation = i

    @property
    def add_OPEX(self):
        '''
        [dict] Operating expense per hour in addition to utility cost.
        Float input will be automatically converted to a dict with the key being
        "Additional OPEX".
        '''
        return {'Additional OPEX': self._add_OPEX} if isinstance(self._add_OPEX, float) \
            else self._add_OPEX
    @add_OPEX.setter
    def add_OPEX(self, i):
        if isinstance(i, float):
            i = {'Additional OPEX': i}
        if not isinstance(i, dict):
            raise TypeError(
                f'add_OPEX can only be float of dict, not {type(i).__name__}.')
        self._add_OPEX = i

    def results(self, with_units=True, include_utilities=True,
                include_total_cost=True, include_installed_cost=False,
                include_zeros=True, external_utilities=(), key_hook=None):

        results = super().results(with_units, include_utilities,
                                  include_total_cost, include_installed_cost,
                                  include_zeros, external_utilities, key_hook)
        if not self.add_OPEX:
            results.loc[('Additional OPEX', ''), :] = ('USD/hr', 0)
        else:
            for k, v in self.add_OPEX.items():
                if not with_units:
                    results.loc[(k, '')] = v
                else:
                    try: results.loc[(k, ''), :] = ('USD/hr', v)
                    # When `results` is a series instead of dataframe,
                    # might not need this
                    except ValueError:
                        results = results.to_frame(name=self.ID)
                        results.insert(0, 'Units', '')
                        results.loc[(k, ''), :] = ('USD/hr', v)
                        results.columns.name = type(self).__name__
        if with_units:
            results.replace({'USD': f'{currency}', 'USD/hr': f'{currency}/hr'},
                            inplace=True)
        return results

    @property
    def uptime_ratio(self):
        '''
        [float] Uptime of the unit to adjust `add_OPEX`, should be in [0,1]
        (i.e., a unit that is always operating).
        '''
        return self._uptime_ratio
    @uptime_ratio.setter
    def uptime_ratio(self, i):
        if not 0 <=i<= 1:
            raise ValueError(f'Uptime must be between 0 and 1 (100%), not {i}.')
        self._uptime_ratio = float(i)

    @property
    def lifetime(self):
        '''
        [int or dict] Lifetime of this unit (int) or individual equipment
        within this unit (dict) in year.
        It will be used to adjust cost and emission calculation in TEA and LCA.
        Equipment without provided lifetime will be assumed to have the same
        lifetime as the TEA/LCA.
        '''
        return self.equipment_lifetime
    @lifetime.setter
    def lifetime(self, i):
        if not i:
            self.equipment_lifetime = {}
        else:
            self.equipment_lifetime = i if isinstance(i, dict) else int(i)