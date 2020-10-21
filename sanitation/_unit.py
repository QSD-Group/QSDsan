#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 21 07:46:34 2020

@author: yalinli_cabbi
"""

import biosteam as bst
from biosteam import utils
from .utils.piping import WSIns, WSOuts


__all__ = ('SanUnit',)

@utils.registered(ticket_name='SU')
class SanUnit(bst.Unit, isabstract=True):
    '''Subclass of Unit in biosteam, is initialized with WasteStream rather than Stream'''

    def __init__(self, ID='', ins=None, outs=(), thermo=None):
        self._specification = None
        self._load_thermo(thermo)
        self._init_ins(ins)
        self._init_outs(outs)
        self._init_utils()
        self._init_results()
        self._assert_compatible_property_package()
        self._register(ID)
    
    def _init_ins(self, ins):
        self._ins = WSIns(self, self._N_ins, ins, self._thermo, self._ins_size_is_fixed)
        
    def _init_outs(self, outs):
        self._outs = WSOuts(self, self._N_outs, outs, self._thermo, self._outs_size_is_fixed)