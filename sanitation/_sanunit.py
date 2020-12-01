#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Sanitation Explorer: Sustainable design of non-sewered sanitation technologies
Copyright (C) 2020, Sanitation Explorer Development Group

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the UIUC open-source license. Please refer to 
https://github.com/QSD-for-WaSH/sanitation/blob/master/LICENSE.txt
for license details.

'''

import biosteam as bst
from biosteam import utils
from . import Construction, Transportation
from .utils.piping import WSIns, WSOuts


__all__ = ('SanUnit',)

@utils.registered(ticket_name='SU')
class SanUnit(bst.Unit, isabstract=True):
    '''Subclass of Unit in biosteam, is initialized with WasteStream rather than Stream.'''

    _stacklevel = 7

    #!!! TODO: write a generic doc for SanUnit and let subclasses inherit it    
    def __init__(self, ID='', ins=None, outs=(), thermo=None):
        self._register(ID)
        self._specification = None
        self._load_thermo(thermo)
        self._init_ins(ins)
        self._init_outs(outs)
        self._init_utils()
        self._init_results()
        self._assert_compatible_property_package()
        self._OPEX = None
        self._construction = ()
        self._transportation = ()

    
    def _init_ins(self, ins):
        self._ins = WSIns(self, self._N_ins, ins, self._thermo,
                          self._ins_size_is_fixed, self._stacklevel)
        
    def _init_outs(self, outs):
        self._outs = WSOuts(self, self._N_outs, outs, self._thermo,
                            self._outs_size_is_fixed, self._stacklevel)

    def _init_results(self):
        super()._init_results()
        self._materials = {}

    def _info(self, T, P, flow, composition, N, _stream_info):
        """Information of the unit."""
        if self.ID:
            info = f'{type(self).__name__}: {self.ID}\n'
        else:
            info = f'{type(self).__name__}\n'
        info += 'ins...\n'
        i = 0
        for stream in self.ins:
            if not stream:
                info += f'[{i}] {stream}\n'
                i += 1
                continue
            if _stream_info:
                stream_info = stream._info(T, P, flow, composition, N) + \
                    '\n' + stream._wastestream_info()
            else:
                stream_info = stream._wastestream_info()
            unit = stream._source
            index = stream_info.index('\n')
            source_info = f'  from  {type(unit).__name__}-{unit}\n' if unit else '\n'
            info += f'[{i}] {stream.ID}' + source_info + stream_info[index+1:] + '\n'
            i += 1
        info += 'outs...\n'
        i = 0
        for stream in self.outs:
            if not stream:
                info += f'[{i}] {stream}\n'
                i += 1
                continue
            if _stream_info:
                stream_info = stream._info(T, P, flow, composition, N) + \
                    '\n' + stream._wastestream_info()
            else:
                stream_info = stream._wastestream_info()
            unit = stream._sink
            index = stream_info.index('\n')
            sink_info = f'  to  {type(unit).__name__}-{unit}\n' if unit else '\n'
            info += f'[{i}] {stream.ID}' + sink_info + stream_info[index+1:] + '\n'
            i += 1
        info = info.replace('\n ', '\n    ')
        return info[:-1]
    
    
    _impact = utils.NotImplementedMethod


    def _summary(self):
        '''After system converges, design unit and calculate cost and environmental impacts'''
        self._design()
        self._cost()
        self._impact()
    
    
    def show(self, T=None, P=None, flow='g/hr', composition=None, N=15, stream_info=True):
        """Print information of the unit, including WasteStream-specific information"""
        print(self._info(T, P, flow, composition, N, stream_info))
        
    def get_ws_inv(self):
        ws_inv = (i for i in self.ins if not (i._source and i.CFs is None))
        ws_inv += (i for i in self.outs if not (i._sink and i.CFs is None))
        return ws_inv
    
    def add_construction(self):
        '''Add construction materials and activities to design and result dict'''
        for i in self.construction:
            self.design_results[i.item.ID] = i.quantity
            self._units[i.item.ID] = i.item.functional_unit
    
    @property
    def construction(self):
        '''[tuple] Contains construction information.'''
        return self._construction
    @construction.setter
    def construction(self, i):
        if isinstance(i, Construction):
            i = (i,)
        else:
            if not iter(i):
                raise TypeError(f'Only <Construction> can be included, not {type(i).__name__}.')
            for j in i:
                if not isinstance(j, Construction):
                    raise TypeError(f'Only <Construction> can be included, not {type(j).__name__}.')
        self._construction = i

    @property
    def transportation(self):
        '''[tuple] Contains transportation information.'''
        return self._transportation
    @transportation.setter
    def transportation(self, i):
        if isinstance(i, Transportation):
            i = (i,)
        else:
            if not iter(i):
                raise TypeError(f'Only <Transportation> can be included, not {type(i).__name__}.')
            for j in i:
                if not isinstance(j, Transportation):
                    raise TypeError(f'Only <Transportation> can be included, not {type(j).__name__}.')
        self._transportation = i

    @property
    def OPEX(self):
        '''[float] Total operating expense per hour.'''
        return self._OPEX or self.utility_cost
    
    #!!! Maybe can just check the _impact results
    @property
    def _if_LCA(self):
        '''
        If this SanUnit is included in LCA. False if:
            - All ins and outs have sink and source (i.e., no chemical inputs or emissions).
            - No HeatUtility and PowerUtility objects.
            - self.construction is empty.

        '''
        inputs = tuple(i for i in self.ins if not i._source)
        if bool(inputs): return True
        emissions = tuple(i for i in self.outs if not i._sink)
        if bool(emissions): return True
        if bool(self.heat_utilities): return True
        if bool(self.power_utility): return True
        if bool(self.construction): return True
        return False


    














