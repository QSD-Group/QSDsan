#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 21 07:46:34 2020

@author: yalinli_cabbi
"""

import biosteam as bst
from .utils.piping import Ins, Outs


__all__ = ('Unit',)

class Unit(bst.Unit, isabstract=True):
    '''Subclass of Unit in biosteam, is initialized with WasteStream rather than Stream'''
    