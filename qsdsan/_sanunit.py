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
from warnings import warn
from biosteam.utils import MissingStream, Inlets, Outlets
from . import currency, Unit, Stream, SanStream, WasteStream, System, \
    Construction, Transportation, Equipment
from .utils import SanUnitScope
from biosteam.utils import Scope

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


def _replace_missing_streams(port, missing):
    for idx, s_type in missing:
        if s_type == 's':
            continue
        s = port[idx]
        stream_class = SanStream if s_type=='ss' else WasteStream
        port.replace(s, stream_class.from_stream(stream_class, s))


def _get_inf_state(inf):
    inf._init_state()
    return inf._state


add2list = lambda lst, item: lst.extend(item) if isinstance(item, Iterable) \
    else lst.append(item)

class SanUnit(Unit, isabstract=True):

    '''
    Subclass of :class:`biosteam.Unit`, can be initialized with
    :class:`thermosteam.Stream`, :class:`~.SanStream`, or :class:`~.WasteStream`.

    Parameters
    ----------
    init_with : str or dict
        Which class of stream the :class:`~.SanUnit` will be initialized with,
        can be "Stream" (shorthanded as "s"), "SanStream" ("ss"), or "WasteStream" ("ws").
        When provided as a str, all streams will be of the same class;
        when provided as a dict, use "ins" or "outs" followed with the order number
        (i.e., ins0, outs-1) as keys; you can use ":" to denote a range (e.g., ins2:4);
        you can also use "else" to specify the stream class for non-provided ones.
    construction : Iterable(obj)
        :class:`~.Construction` with construction information.
    transportation : Iterable(obj)
        :class:`~.Transportation` with transportation information.
    equipments: Iterable(obj)
        :class:`~.Equipment` with equipment information.
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


    lifetime : int or dict
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
        for i in (*construction, *transportation, *equipments):
            i._linked_unit = self
        # Make fresh ones for each unit
        self.construction = () if not construction else construction
        self.transportation = () if not transportation else transportation
        self.equipments = () if not equipments else equipments
        self.add_OPEX = add_OPEX.copy()
        self.uptime_ratio = 1.
        self.lifetime = lifetime
        if F_BM_default:
            F_BM = self.F_BM
            self.F_BM = defaultdict(lambda: F_BM_default)
            self.F_BM.update(F_BM)
        # For units with different state headers, should update it in the unit's ``__init__``
        self.isdynamic = isdynamic
        for attr, val in kwargs.items():
            setattr(self, attr, val)


    def _convert_stream(self, strm_inputs, streams, init_with, ins_or_outs):
        isa = isinstance
        if not streams:
            return [], []
        if isa(init_with, str):
            init_with = dict.fromkeys([f'{ins_or_outs}{n}'
                                       for n in range(len(streams))], init_with)
        # Do not change pre-defined stream types
        if isa(strm_inputs, Stream): # input is a stream
            init_with[f'{ins_or_outs}0'] = type(strm_inputs).__name__
        # `input` is an iterable of stream
        elif isa(strm_inputs, Iterable) and not isa(strm_inputs, str):
            for in_or_out, s in zip(ins_or_outs, strm_inputs):
                if isa(s, Stream):
                    init_with[in_or_out] = type(s).__name__

        init_with = _update_init_with(init_with, ins_or_outs, len(streams))

        converted, missing = [], []
        for k, v in init_with.items():
            if not ins_or_outs in k: # leave out outsX when going through ins and vice versa
                continue
            num = int(k.split(ins_or_outs)[-1])
            s = streams[num]
            if isa(s, MissingStream):
                missing.append((num, v))
                continue
            if v == 's':
                converted.append(s)
            elif v == 'ss':
                converted.append(SanStream.from_stream(SanStream, s))
            else:
                converted.append(WasteStream.from_stream(WasteStream, s))

        diff = len(converted) + len(missing) - len(streams)
        if diff != 0:
            raise ValueError(f'Type(s) of {diff} stream(s) has/have not been specified.')

        return converted, missing


    def _init_dynamic(self):
        self._state = None
        self._dstate = None
        self._ins_QC = np.zeros((len(self._ins), len(self.components)+1))
        self._ins_dQC = self._ins_QC.copy()
        self._ODE = None
        self._AE = None
        if not hasattr(self, '_mock_dyn_sys'):
            self._mock_dyn_sys = System(self.ID+'_dynmock', path=(self,))
        if not hasattr(self, '_state_header'):
            self._state_header = [f'{cmp.ID} [mg/L]' for cmp in self.components] + ['Q [m3/d]']
        # Shouldn't need to re-create the mock system everytime
        # if hasattr(self, '_mock_dyn_sys'):
        #     sys = self._mock_dyn_sys
        #     sys.registry.discard(sys)
        # self._mock_dyn_sys = System(self.ID+'_dynmock', path=(self,))


    def _init_ins(self, ins, init_with):
        super()._init_ins(ins)
        converted, missing = self._convert_stream(ins, self.ins, init_with, 'ins')
        _ins = self._ins = Inlets(self, self._N_ins, converted, self._thermo,
                                  self._ins_size_is_fixed, self._stacklevel)
        # Cannot do it within `_convert_stream` as it creates new streams and
        # cause error in `Inlets/Outlets._redock`
        _replace_missing_streams(_ins, missing)


    def _init_outs(self, outs, init_with):
        super()._init_outs(outs)
        converted, missing = self._convert_stream(outs, self.outs, init_with, 'outs')
        _outs = self._outs = Outlets(self, self._N_outs, converted, self._thermo,
                                     self._outs_size_is_fixed, self._stacklevel)
        _replace_missing_streams(_outs, missing)

    def _init_results(self):
        super()._init_results()
        self.add_OPEX = {}

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

    def set_dynamic_tracker(self, *subjects, **kwargs):
        """
        Set up an :class:`SystemScope` object to track the dynamic data.

        Parameters
        ----------
        *subjects :
            Any subjects of the system to track, which must have an `.scope`
            attribute of type :class:`Scope`.
        """
        sys = self._mock_dyn_sys
        if self.isdynamic:
            sys._scope = {'subjects':subjects, 'kwargs':kwargs}
        else:
            warn(f'{self.ID} is not a dynamic unit, cannot set tracker.')

    def simulate(self, t_span=(0, 0), state_reset_hook=True,
                 solver='', **kwargs):
        '''
        Converge mass and energy flows, design, and cost the unit.

        .. note::

            If this unit is a dynamic unit, ODEs will be run after ``_run``
            and/or ``specification``.

        Parameters
        ----------
        t_span : tuple(float, float)
            Integration time span for dynamic units.
        state_reset_hook: str or callable
            Hook function to reset the cache state between simulations (for dynamic systems).
            Can be "reset_cache" or "clear_state" to call `System.reset_cache` or `System.clear_state`,
            or None to avoiding resetting.
        kwargs : dict
            Other keyword arguments that will be passed to ``scipy.integrate.solve_ivp``

        See Also
        --------
        :func:`~.System.simulate`
        
        `scipy.integrate.solve_ivp <https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html>`_
        '''
        super().simulate()
        if self.isdynamic:
            sys = self._mock_dyn_sys
            sys._feeds = self.ins
            sys._products = self.outs
            sys.simulate(t_span=t_span,
                         state_reset_hook=state_reset_hook,
                         **kwargs)
            self._summary()

    def show(self, T=None, P=None, flow='g/hr', composition=None, N=15, IDs=None, stream_info=True):
        '''Print information of the unit, including waste stream-specific information.'''
        print(self._info(T, P, flow, composition, N, IDs, stream_info))

    def add_equipment_design(self):
        unit_design = self.design_results
        unit_units = self._units
        F_BM, F_D, F_P, F_M, lifetime = \
            self.F_BM, self.F_D, self.F_P, self.F_M, self._default_equipment_lifetime
        isa = isinstance
        get = getattr
        def update_unit_attr(unit_attr, equip_ID, equip_attr):
            if isa(equip_attr, dict):
                unit_attr.update(equip_attr)
            else:
                unit_attr[equip_ID] = equip_attr
        for equip in self.equipments:
            equip_ID = equip.ID
            equip_design = equip._design_results = equip._design()
            equip_design = {} if not equip_design else equip_design
            unit_design.update(equip_design)

            equip_units = {} if not equip.units else equip.units
            unit_units.update(equip_units)
            for unit_attr, equip_attr in zip(
                    (F_BM, F_D, F_P, F_M, lifetime),
                    ('F_BM', 'F_D', 'F_P', 'F_M', 'lifetime'),
                    ):
                update_unit_attr(unit_attr, equip_ID, get(equip, equip_attr))


    def add_equipment_cost(self):
        unit_cost = self.baseline_purchase_costs
        for equip in self.equipments:
            equip_cost = equip._baseline_purchase_costs = equip._cost()
            if isinstance(equip_cost, dict):
                unit_cost.update(equip_cost)
            else:
                unit_cost[equip.ID] = equip_cost

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
        hasfield = hasattr
        if hasfield(self, '_isdynamic'):
            if self._isdynamic == bool(i):
                return
        else: self._isdynamic = bool(i)
        if self.hasode:
            self._init_dynamic()
            if hasattr(self, '_mock_dyn_sys'):
                ID = self.ID+'_dynmock'
                System.registry.discard(ID)
            self._mock_dyn_sys = System(self.ID+'_dynmock', path=(self,))

    @property
    def scope(self):
        """A tracker of the unit's time-series data during dynamic simulation."""
        if not hasattr(self, '_scope'):
            self._scope = SanUnitScope(self)
        return self._scope
    
    @scope.setter
    def scope(self, s):
        if not isinstance(s, Scope): 
            raise TypeError(f'{s} must be an {Scope} not {type(s)}.')
        if self is not s.subject:
            raise ValueError(f'The subject of {s} must be {self} not {s.subject}.')
        self._scope = s

    @property
    def hasode(self):
        """Whether this unit's dynamic states are determined by ordinary differential equations."""
        return hasattr(self, '_compile_ODE')

    def reset_cache(self, dynamic_system=False):
        '''Reset cached states for dynamic units.'''
        super().reset_cache()
        if self.hasode or dynamic_system:
            self._init_dynamic()
            for s in self.outs:
                s.empty()

    def get_retained_mass(self, biomass_IDs):
        warn(f'The retained biomass in {self.ID} is ignored.')

    @property
    def construction(self):
        '''Iterable(obj) :class:`~.Construction` with construction information.'''
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
        '''Iterable(obj) :class:`~.Transportation` with transportation information.'''
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
    def equipments(self):
        '''Iterable(obj) :class:`~.Equipment` with equipments information.'''
        return self._equipments
    @equipments.setter
    def equipments(self, i):
        isa = isinstance
        if isa(i, Equipment):
            i = (i,)
        else:
            if not isa(i, Iterable):
                raise TypeError(
                    f'Only `Equipment` object  can be included, not {type(i).__name__}.')
            for j in i:
                if not isa(j, Equipment):
                    raise TypeError(
                        f'Only `Equipment` can be included, not {type(j).__name__}.')
        self._equipments = i

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
        isa = isinstance
        if isa(i, float):
            i = {'Additional OPEX': i}
        if not isa(i, dict):
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