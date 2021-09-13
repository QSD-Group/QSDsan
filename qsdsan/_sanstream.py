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
from biosteam.utils import MissingStream
from . import Stream

__all__ = ('SanStream', 'MissingSanStream')

setattr = object.__setattr__

class SanStream(Stream):
    '''
    A subclass of :class:`thermosteam.Stream` with additional attributes
    for environmental impacts.

    .. note::

        Parameters below only include the ones additional to those of :class:`thermosteam.Stream`.


    Parameters
    ----------
    stream_impact_item : :class:`StreamImpactItem`
        The :class:`StreamImpactItem` this stream is linked to.

    Examples
    --------
    `WasteStream <https://qsdsan.readthedocs.io/en/latest/tutorials/WasteStream.html>`_

    See Also
    --------
    `thermosteam.Stream <https://thermosteam.readthedocs.io/en/latest/Stream.html>`_

    '''

    __slots__ = (*Stream.__slots__, '_impact_item', '_stream_impact_item')
    ticket_name = 'ss'

    def __init__(self, ID='', flow=(), phase='l', T=298.15, P=101325.,
                 units='kg/hr', price=0., thermo=None, stream_impact_item=None,
                 **component_flows):
        super().__init__(ID=ID, flow=flow, phase=phase, T=T, P=P,
                         units=units, price=price, thermo=thermo,
                         **component_flows)

        if stream_impact_item:
            stream_impact_item._linked_stream = self
        self._stream_impact_item = stream_impact_item
        self._impact_item = self._stream_impact_item

    def copy(self, new_ID=''):
        '''
        Copy the information of another stream.

        Parameters
        ----------
        new_ID : str
            ID of the new stream, a default ID will be assigned if not provided.


        .. note::

            [1] Price of the original stream is not copied.

            [2] If the original stream has an :class:`~.StreamImpactItem`,
            then a new :class:`~.StreamImpactItem` will be created for the new stream
            and the new impact item will be linked to the original impact item.
        '''

        new = super().copy(ID=new_ID)
        if hasattr(self, '_stream_impact_item'):
            if self.stream_impact_item is not None:
                self.stream_impact_item.copy(stream=new)
        return new

    __copy__ = copy

    @staticmethod
    def from_stream(cls, stream, ID='', **kwargs):
        '''
        Cast a :class:`thermosteam.Stream` or :class:`biosteam.utils.MissingStream`
        to :class:`SanStream` or :class:`MissingSanStream`.

        Parameters
        ----------
        cls : obj
            class of the stream to be created.
        stream : :class:`thermosteam.Stream`
            The original stream.
        ID : str
            If not provided, will use default ID.
        kwargs
            Additional properties of the new stream.

        Examples
        --------
        >>> import qsdsan as qs
        >>> cmps = qs.Components.load_default()
        >>> qs.set_thermo(cmps)
        >>> s = qs.Stream('s', H2O=100, price=5)
        >>> s.show()
        Stream: s
         phase: 'l', T: 298.15 K, P: 101325 Pa
         flow (kmol/hr): H2O  100
        >>> s.price
        5.0
        >>> ss = qs.SanStream.from_stream(qs.SanStream, s, ID='ss', T=350, price=10)
        >>> ss.show()
        SanStream: ss
         phase: 'l', T: 350 K, P: 101325 Pa
         flow (kmol/hr): H2O  100
        >>> ss.price
        10.0
        '''

        if isinstance(stream, MissingStream):
            new = MissingSanStream.__new__(MissingSanStream)
            return new

        elif not isinstance(stream, cls):
            if not ID:
                stream.registry.discard(stream)
                # stream.registry.untrack((stream,))
            new = cls.__new__(cls)

            new_ID = ID if ID else stream.ID
            if new_ID[0]=='s' and new_ID[1:].isnumeric(): # old ID is default
                new_ID = ''
            new.__init__(ID=new_ID)

            new._link = stream._link

            source = new._source = stream._source
            if source:
                source._outs[source._outs.index(stream)] = new

            sink = new._sink = stream._sink
            if sink:
                sink._ins[sink._ins.index(stream)] = new

            new._thermo = stream._thermo
            new._imol = stream._imol.copy()
            new._thermal_condition = stream._thermal_condition.copy()
            new.reset_cache()
            new.price = 0
            new.stream_impact_item = None

            for attr, val in kwargs.items():
                setattr(new, attr, val)

            stream._link = stream._sink = stream._source = None
            return new

        else:
            return stream

    def mix_from(self, others):
        '''
        Update this stream to be a mixture of other streams,
        initial content of this stream will be ignored.

        Parameters
        ----------
        others : iterable
            Can contain :class:`thermosteam.Stream`, :class:`SanStream`,
            or :class:`~.WasteStream`

        .. note::

            Price and impact item are not included.


        Examples
        --------
        >>> import qsdsan as qs
        >>> cmps = qs.Components.load_default()
        >>> qs.set_thermo(cmps)
        >>> s1 = qs.Stream(H2O=100, price=5, units='kg/hr')
        >>> s2 = qs.SanStream(S_O2=100, units='kg/hr')
        >>> s3 = qs.SanStream()
        >>> s3.mix_from((s1, s2))
        >>> s3.show()
        SanStream: ss2
         phase: 'l', T: 298.15 K, P: 101325 Pa
         flow (kmol/hr): S_O2  3.13
                         H2O   5.55
        '''

        others = [s for s in others if not 'Missing' in type(s).__name__]
        Stream.mix_from(self, others)


    @property
    def stream_impact_item(self):
        '''[:class:`StreamImpactItem`] The :class:`StreamImpactItem` this stream is linked to.'''
        return self._stream_impact_item
    @stream_impact_item.setter
    def stream_impact_item(self, i):
        self._stream_impact_item = i
        if i:
            i.linked_stream = self

    @property
    def impact_item(self):
        '''Same as `stream_impact_item`, has been deprecated.'''
        warn('The property `impact_item` has been changed to `stream_impact_item`, '
             'please use `stream_impact_item` instead.', stacklevel=2)
        return self.stream_impact_item

    @property
    def components(self):
        return self.chemicals


# %%

class MissingSanStream(MissingStream):
    '''
    A subclass of :class:`biosteam.MissingStream`, create a special object
    that acts as a dummy until replaced by an actual :class:`SanStream`.

    .. note::

        Users usually do not need to interact with this class.

    '''

    @property
    def stream_impact_item(self):
        return None

    @property
    def impact_item(self):
        '''Same as `stream_impact_item`, has been deprecated.'''
        warn('The property `impact_item` has been changed to `stream_impact_item`, '
             'please use `stream_impact_item` instead.')
        return self.stream_impact_item

    def __repr__(self):
        return '<MissingSanStream>'

    def __str__(self):
        return 'missing sanstream'