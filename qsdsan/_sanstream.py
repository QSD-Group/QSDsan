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
        if 'impact_item' in component_flows.keys():
            raise ValueError('The keyword `impact_item` is deprecated, '
                             'please use `stream_impact_item` instead.')
        super().__init__(ID=ID, flow=flow, phase=phase, T=T, P=P,
                         units=units, price=price, thermo=thermo,
                         **component_flows)
        if stream_impact_item:
            stream_impact_item._linked_stream = self
        self._stream_impact_item = stream_impact_item
        self._impact_item = self._stream_impact_item

    def copy(self, new_ID='', copy_price=False, copy_impact_item=False):
        '''
        Copy the information of another stream.

        There are three functions related to copying: ``copy``, ``copy_like``, and ``copy_flow``,
        and they have slight differences in using.

        Both ``copy`` and ``copy_like`` makes the new stream the same as the original one
        (other than that the new stream does not have the cost),
        but when using ``copy``, you do now need to pre-create the new stream,
        (i.e., you can just do ``new_stream = original_stream.copy('new_ID')``),
        but to use ``copy_like``, you need to firstly create the new stream, then
        ``new_stream.copy_like(original_stream)``.

        For ``copy_flow``, it is similar to ``copy_like`` in that you need to firstly
        creating the new stream, but unlike ``copy_flow`` that copies properties
        such as temperature, pressure, ``copy_flow`` just copies the mass flow information,
        but you can choose which component to copy.

        Parameters
        ----------
        new_ID : str
            ID of the new stream, a default ID will be assigned if not provided.
        copy_price : bool
            If True, price of the new stream will be set to be the same as
            the original stream.
        copy_impact_item : bool
            If True and the original stream has an :class:`~.StreamImpactItem`,
            then a new :class:`~.StreamImpactItem` will be created for the new stream
            and the new impact item will be linked to the original impact item.
        '''

        new = super().copy(ID=new_ID)
        if copy_price:
            new.price = self.price
        if hasattr(self, '_stream_impact_item'):
            if self.stream_impact_item is not None:
                self.stream_impact_item.copy(stream=new)
            else:
                new._stream_impact_item = None
        return new

    __copy__ = copy


    def copy_like(self, other, copy_price=False, copy_impact_item=False):
        '''
        Copy the information of another stream without creating a new stream.

        Parameters
        ----------
        other : obj
            The stream where mass flows and stream properties will be copied from.
        copy_price : bool
            If True, price of the new stream will be set to be the same as
            the original stream.
        copy_impact_item : bool
            If True and the original stream has an :class:`~.StreamImpactItem`,
            then a new :class:`~.StreamImpactItem` will be created for the new stream
            and the new impact item will be linked to the original impact item.

        See Also
        --------
        :func:`copy` for the differences between ``copy``, ``copy_like``, and ``copy_flow``.
        '''

        Stream.copy_like(self, other)
        if not isinstance(other, SanStream):
            return
        if copy_price:
            other.price = self.price
        if copy_impact_item:
            if hasattr(other, '_stream_impact_item'):
                if other.stream_impact_item is not None:
                    self.stream_impact_item.copy(stream=self)


    def copy_flow(self, other, IDs=..., *, remove=False, exclude=False):
        '''
        Copy only the mass flow of another stream without creating a new stream.

        Parameters
        ----------
        other : obj
            The stream where mass flows will be copied from.
        IDs=... : Iterable(str), defaults to all components.
            IDs of the components to be copied from.
        remove=False: bool, optional
            If True, copied components will be removed from the original stream.
        exclude=False: bool, optional
            If True, exclude designated components when copying.


        See Also
        --------
        :func:`copy` for the differences between ``copy``, ``copy_like``, and ``copy_flow``.
        '''
        Stream.copy_flow(self, other=other, IDs=IDs, remove=remove, exclude=exclude)

        if not isinstance(other, SanStream):
            return

        self._stream_impact_item = None


    def flow_proxy(self, ID=None):
        '''
        Return a new stream that shares flow data with this one.

        Parameters
        ----------
        ID : str
            ID of the new proxy stream.

        Examples
        --------
        >>> from qsdsan import set_thermo, SanStream
        >>> from qsdsan.utils import load_example_cmps
        >>> cmps = load_example_cmps()
        >>> set_thermo(cmps)
        >>> ss1 = SanStream('ss1', Water=100, NaCl=1, price=3.18)
        >>> ss2 = ss1.flow_proxy('ss2')
        >>> ss2.mol is ss1.mol
        True
        >>> ss2.thermal_condition is ss1.thermal_condition
        False
        '''
        new = Stream.flow_proxy(self, ID=ID)
        new._stream_impact_item = None
        return new


    def proxy(self, ID=None):
        '''
        Return a new stream that shares all data with this one.

        Note that unlike other properties, the `price` and `stream_impact_item`
        of the two streams are not connected,
        i.e., the price of the new stream will be the same as the
        original one upon creation, but then they can be different.

        Parameters
        ----------
        ID : str
            ID of the new proxy stream.

        Examples
        --------
        >>> from qsdsan import set_thermo, SanStream
        >>> from qsdsan.utils import load_example_cmps
        >>> cmps = load_example_cmps()
        >>> set_thermo(cmps)
        >>> ss1 = SanStream('ss1', Water=100, NaCl=1, price=3.18)
        >>> ss2 = ss1.proxy('ss2')
        >>> ss2.mol is ss1.mol
        True
        >>> ss2.thermal_condition is ss1.thermal_condition
        True
        >>> ss2.price = 5.2335
        >>> ss1.price
        3.18
        '''
        new = Stream.proxy(self, ID=ID)
        new._stream_impact_item = None
        return new


    @staticmethod
    def degassing(original_stream, receiving_stream=None, gas_IDs=()):
        '''
        Remove all the gas components from the original stream,
        if `receiving_stream` is given, then the gas components will be transfered
        to the receiving stream.

        If `gas_IDs` is not provided, then the gas components wiil be those
        either have `locked_state` == "g" or `particle_size` == "Dissolved gas".

        Parameters
        ----------
        original_stream : None or obj
            The stream where the gas components will be removed.
        receiving_stream : None or obj
            The stream to receive the gas components.
        gas_IDs : Iterable(str)
            IDs of the gas components to be removed, will be set according
            to the component properties if not provided.
        '''
        if not gas_IDs:
            gas_IDs = original_stream.gases if isinstance(original_stream, SanStream) \
                else [i.ID for i in original_stream.components if i.locked_state=='g']
        if receiving_stream:
            receiving_stream.imass[gas_IDs] += original_stream.imass[gas_IDs]
        original_stream.imass[gas_IDs] = 0


    @staticmethod
    def filtering(original_stream, receiving_stream=None, solid_IDs=()):
        '''
        Remove all the solid components from the original stream,
        if `receiving_stream` is given, then the solid components will be transfered
        to the receiving stream.

        If `solid_IDs` is not provided, then the gas components wiil be those
        either have `locked_state` == "g" or `particle_size` == "Particulate".

        Parameters
        ----------
        original_stream : None or obj
            The stream where the gas components will be removed.
        receiving_stream : None or obj
            The stream to receive the gas components.
        solid_IDs : Iterable(str)
            IDs of the solid components to be removed, will be set according
            to the component properties if not provided.
        '''
        if not solid_IDs:
            solid_IDs = original_stream.solids if isinstance(original_stream, SanStream) \
                else [i.ID for i in original_stream.components if i.locked_state=='s']
        if receiving_stream:
            receiving_stream.imass[solid_IDs] += original_stream.imass[solid_IDs]
        original_stream.imass[solid_IDs] = 0


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
        # Missing stream, note that if to make updates here,
        # it's likely that `WasteStream.from_stream` should be updated as well.
        if isinstance(stream, MissingStream):
            new = MissingSanStream.__new__(MissingSanStream)
            new.__init__(stream._source, stream._sink)
            return new
        # An actual stream
        if not isinstance(stream, cls):
            if not ID:
                stream.registry.discard(stream)
                # stream.registry.untrack((stream,))
            new = cls.__new__(cls)
            new_ID = ID if ID else stream.ID
            if new_ID[0]=='s' and new_ID[1:].isnumeric(): # old ID is default
                new_ID = ''
            new.__init__(ID=new_ID)

            new._islinked = stream._islinked

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

            stream._islinked = False
            stream._sink = stream._source = None
            return new
        else:
            return stream

    def mix_from(self, others, **kwargs):
        '''
        Update this stream to be a mixture of other streams,
        initial content of this stream will be ignored.

        Parameters
        ----------
        others : Iteralbe(obj)
            Can contain :class:`thermosteam.Stream`, :class:`SanStream`,
            or :class:`~.WasteStream`

        .. note::

            Price and impact item are not included.


        Examples
        --------
        >>> import qsdsan as qs
        >>> cmps = qs.Components.load_default()
        >>> qs.set_thermo(cmps)
        >>> s1 = qs.Stream('s1', H2O=100, price=5, units='kg/hr')
        >>> s2 = qs.SanStream('s2', S_O2=100, units='kg/hr')
        >>> s3 = qs.SanStream('s3')
        >>> s3.mix_from((s1, s2))
        >>> s3.show()
        SanStream: s3
         phase: 'l', T: 298.15 K, P: 101325 Pa
         flow (kmol/hr): S_O2  3.13
                         H2O   5.55
        '''

        others = [s for s in others if not 'Missing' in type(s).__name__]
        Stream.mix_from(self, others, **kwargs)
        if not hasattr(self, '_stream_impact_item'):
            self._stream_impact_item = None


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
             'please use `stream_impact_item` instead.')
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
    line = 'SanStream'

    def materialize_connection(self, ID=''):
        '''
        Disconnect this missing stream from any unit operations and
        replace it with a material stream.
        '''
        source = self._source
        sink = self._sink
        if not (source or sink):
            raise RuntimeError("either a source or a sink is required to "
                               "materialize connection")
        material_stream = SanStream(ID, thermo=(source or sink).thermo)
        if source: source._outs.replace(self, material_stream)
        if sink: sink._ins.replace(self, material_stream)
        return material_stream

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