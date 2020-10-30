#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 21:47:27 2020

@author: yalinli_cabbi
"""

from biosteam.utils.piping import MissingStream, StreamSequence
from .. import WasteStream as WS

# from biosteam.utils import piping
# MissingStream = piping.MissingStream
# StreamSequence = piping.StreamSequence
# __sub__ = piping.__sub__
# __rsub__ = piping.__rsub__
# pipe_info = piping.pipe_info

__all__ = ('MissingWS', 'WSIns', 'WSOuts', 'as_ws')


# %%

def as_ws(ws):
    if isinstance(ws, WS):
        return ws
    elif isinstance(ws, str):
        return WS(ws)
    elif ws is None:
        return MissingWS(None, None)
    

class MissingWS(MissingStream):
    '''
    A subclass of MissingStream in biosteam, Create a MissingWS (MissingWasteStream)
    object that acts as a dummy in Ins and Outs objects until replaced by an
    actual WasteStream object
    '''
    
    def materialize_connection(self, ID=''):
        '''
        Disconnect this MissingWS from any unit operations and replace
        it with an actual WasteStream
        '''
        source = self._source
        sink = self._sink
        assert source and sink, (
            "both a source and a sink is required to materialize connection")
        material_stream = WS(ID, thermo=source.thermo)
        source._outs.replace(self, material_stream)
        sink._ins.replace(self, material_stream)

    # TODO: add others
    @property
    def COD(self):
        return 0.
    @property
    def BOD(self):
        return 0.
    @property
    def TC(self):
        return 0.
    @property
    def TOC(self):
        return 0.
    @property
    def TN(self):
        return 0.
    @property
    def TKN(self):
        return 0.
    @property
    def TP(self):
        return 0.
    @property
    def TK(self):
        return 0.
    @property
    def charge(self):
        return 0.

    def __repr__(self):
        return '<MissingWS>'
    
    def __str__(self):
        return 'missing waste stream'
    

# %%


def n_missing(ub, N):
    assert ub >= N, f"size of waste streams exceeds {ub}"
    return ub - N

class WSSequence(StreamSequence):
    '''Abstract class for a sequence of waste streams, subclass of StreamSequence in biosteam'''
    
    def __init__(self, size, ws, thermo, fixed_size):
        self._size = size
        self._fixed_size = fixed_size
        dock = self._dock
        redock = self._redock
        if ws == ():
            self._streams = [dock(WS(thermo=thermo)) for i in range(size)]
        else:
            if fixed_size:
                self._initialize_missing_streams()
                if ws:
                    if isinstance(ws, str):
                        self._streams[0] = dock(WS(ws, thermo=thermo))
                    elif isinstance(ws, WS):
                        self._streams[0] = redock(ws)
                    else:
                        N = len(ws)
                        n_missing(size, N) # Assert size is not too big
                        self._streams[:N] = [redock(i) if isinstance(i, WS)
                                             else dock(WS(i, thermo=thermo)) for i in ws]
            else:
                if ws:
                    if isinstance(ws, str):
                        self._streams = [dock(WS(ws, thermo=thermo))]
                    elif isinstance(ws, WS):
                        self._streams = [redock(ws)]
                    else:
                        self._streams = [redock(i) if isinstance(i, WS)
                                         else dock(WS(i, thermo=thermo))
                                         for i in ws]
                else:
                    self._initialize_missing_streams()

    def _create_missing_stream(self):
        MissingWS(None, None).show()
        return MissingWS(None, None)
    
    def _create_N_missing_streams(self, N):
        return [self._create_missing_stream() for i in range(N)]
    
    def _initialize_missing_streams(self):
        self._streams = self._create_N_missing_streams(self._size)

    def __setitem__(self, index, item):
        if isinstance(index, int):
            if item:
                assert isinstance(item, WS), (
                    f"'{type(self).__name__}' object can only contain "
                    f"'WasteStream' objects; not '{type(item).__name__}'")
            elif not isinstance(item, MissingWS):
                item = self._create_missing_stream()
            self._set_stream(index, item)
        elif isinstance(index, slice):
            wastestreams = []
            for ws in item:
                if ws:
                    assert isinstance(ws, WS), (
                        f"'{type(self).__name__}' object can only contain "
                        f"'WasteStream' objects; not '{type(ws).__name__}'")
                elif not isinstance(ws, MissingWS):
                    ws = self._create_missing_stream()
                wastestreams.append(ws)
            self._set_streams(index, item)
        else:
            raise TypeError("Only intergers and slices are valid "
                           f"indices for '{type(self).__name__}' objects")



class WSIns(WSSequence):
    '''Create an Ins object which serves as input streams for a Unit object'''
    __slots__ = ('_sink', '_fixed_size')
    
    def __init__(self, sink, size, ws, thermo, fixed_size=True):
        self._sink = sink
        super().__init__(size, ws, thermo, fixed_size)
    
    @property
    def sink(self):
        return self._sink
    
    def _create_missing_stream(self):
        return MissingWS(None, self._sink)
    
    def _dock(self, ws): 
        ws._sink = self._sink
        return ws

    def _redock(self, ws): 
        sink = ws._sink
        if sink:
            ins = sink._ins
            if ins is not self:
                ins.remove(ws)
                ws._sink = self._sink
        else:
            ws._sink = self._sink
        return ws
    
    def _undock(self, ws): 
        ws._sink = None


class WSOuts(WSSequence):
    '''Create an Outs object which serves as output streams for a Unit object'''
    __slots__ = ('_source',)
    
    def __init__(self, source, size, ws, thermo, fixed_size=True):
        self._source = source
        super().__init__(size, ws, thermo, fixed_size)
    
    @property
    def source(self):
        return self._source
    
    def _create_missing_stream(self):
        return MissingWS(self._source, None)
    
    def _dock(self, ws): 
        ws._source = self._source
        return ws

    def _redock(self, ws): 
        source = ws._source
        if source:
            outs = source._outs
            if outs is not self:
                # Remove from source
                outs.remove(ws)
                ws._source = self._source
        else:
            ws._source = self._source
        return ws
    
    def _undock(self, ws): 
        ws._source = None



# WS.__pow__ = WS.__sub__ = __sub__  # Forward pipping
# WS.__rpow__ = WS.__rsub__ = __rsub__ # Backward pipping    
# WS.sink = property(lambda self: self._sink)
# WS.source = property(lambda self: self._source)
# WS._basic_info = lambda self: (f"{type(self).__name__}: {self.ID or ''}"
#                                   f"{pipe_info(self._source, self._sink)}\n")




























