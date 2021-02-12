#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/master/LICENSE.txt
for license details.
'''

from warnings import warn
from biosteam.utils.piping import MissingStream, StreamSequence
from .. import WasteStream as WS
import biosteam as bst

__all__ = ('MissingWS', 'WSIns', 'WSOuts')


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
    A subclass of :class:`biosteam.MissingStream`, create a special object
    object that acts as a dummy in :class:`WSIns` and :class:`WSOuts` objects
    until replaced by an actual :class:`WasteStream` object.
    '''

    def materialize_connection(self, ID=''):
        '''
        Disconnect this missing waste stream from any unit operations and replace
        it with an actual waste stream.
        '''
        source = self._source
        sink = self._sink
        assert source and sink, (
            'both a source and a sink is required to materialize connection.')
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
    assert ub >= N, f'size of waste streams exceeds {ub}.'
    return ub - N

class WSSequence(StreamSequence):
    '''
    Abstract class for a sequence of waste streams, subclass of
    :class:`biosteam.StreamSequence`.
    '''
    
    def __init__(self, size, ws, thermo, fixed_size, stacklevel):
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
                        self._streams[0] = redock(ws, stacklevel)
                    else:
                        N = len(ws)
                        n_missing(size, N) # Assert size is not too big
                        self._streams[:N] = [redock(i, stacklevel+1) if isinstance(i, WS)
                                             else dock(WS(i, thermo=thermo)) for i in ws]
            else:
                if ws:
                    if isinstance(ws, str):
                        self._streams = [dock(WS(ws, thermo=thermo))]
                    elif isinstance(ws, WS):
                        self._streams = [redock(ws, stacklevel)]
                    else:
                        self._streams = [redock(i, stacklevel+1) if isinstance(i, WS)
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
                    f'`{type(self).__name__}` object can only contain '
                    f'`WasteStream` objects; not {type(item).__name__}.')
            elif not isinstance(item, MissingWS):
                item = self._create_missing_stream()
            self._set_stream(int=index, stream=item, stacklevel=3)
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
            self._set_streams(slice=index, streams=item, stacklevel=3)
        else:
            raise TypeError("Only intergers and slices are valid "
                           f"indices for '{type(self).__name__}' objects")


    # Below shouldn't be needed for newer biosteam
    def _set_streams(self, slice, streams, stacklevel):
        all_streams = self._streams
        for stream in all_streams[slice]: self._undock(stream)
        all_streams[slice] = streams
        for stream in all_streams:
            self._redock(stream, stacklevel)
        if self._fixed_size:
            size = self._size
            N_streams = len(all_streams)
            if N_streams < size:
                N_missing = n_missing(size, N_streams)
                if N_missing:
                    all_streams[N_streams: size] = self._create_N_missing_streams(N_missing)

    def _set_stream(self, int, stream, stacklevel):
        self._undock(self._streams[int])
        self._redock(stream, stacklevel)
        self._streams[int] = stream


class WSIns(WSSequence):
    '''Create an object to serve as an input stream for a :class:`SanUnit` object.'''
    __slots__ = ('_sink', '_fixed_size')
    
    def __init__(self, sink, size, ws, thermo, fixed_size, stacklevel):
        self._sink = sink
        super().__init__(size, ws, thermo, fixed_size, stacklevel)
    
    @property
    def sink(self):
        return self._sink
    
    def _create_missing_stream(self):
        return MissingWS(None, self._sink)
    
    def _dock(self, ws): 
        ws._sink = self._sink
        return ws

    def _redock(self, ws, stacklevel): 
        sink = ws._sink
        if sink:
            ins = sink._ins
            if ins is not self:
                ins.remove(ws)
                ws._sink = new_sink = self._sink
                DOCKING_WARNINGS = bst.utils.piping.DOCKING_WARNINGS
                if (DOCKING_WARNINGS 
                    and sink._ID and new_sink._ID
                    and sink._ID != new_sink._ID):
                    warn(f'Inlet waste stream {ws} is undocked from unit {sink}; '
                         f'{ws} is now docked at {self._sink}', 
                         RuntimeWarning, stacklevel)
        else:
            ws._sink = self._sink
        return ws
    
    def _undock(self, ws): 
        ws._sink = None


class WSOuts(WSSequence):
    '''Create an object to serve as an output stream for a :class:`SanUnit` object.'''
    __slots__ = ('_source',)
    
    def __init__(self, source, size, ws, thermo, fixed_size, stacklevel):
        self._source = source
        super().__init__(size, ws, thermo, fixed_size, stacklevel)
    
    @property
    def source(self):
        return self._source
    
    def _create_missing_stream(self):
        return MissingWS(self._source, None)
    
    def _dock(self, ws): 
        ws._source = self._source
        return ws

    def _redock(self, ws, stacklevel): 
        source = ws._source
        if source:
            outs = source._outs
            if outs is not self:
                # Remove from source
                outs.remove(ws)
                ws._source = new_source = self._source
                DOCKING_WARNINGS = bst.utils.piping.DOCKING_WARNINGS
                if (DOCKING_WARNINGS 
                    and source._ID and new_source._ID
                    and source._ID != new_source._ID):
                    warn(f'Outlet waste stream {ws} is undocked from unit {source}; '
                         f'{ws} is now docked at {self._source}.', 
                         RuntimeWarning, stacklevel)
        else:
            ws._source = self._source
        return ws
    
    def _undock(self, ws): 
        ws._source = None























