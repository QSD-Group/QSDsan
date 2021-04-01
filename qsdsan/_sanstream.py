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

from . import Stream
from biosteam.utils import MissingStream

__all__ = ('SanStream', 'MissingSanStream')

setattr = object.__setattr__

class SanStream(Stream):
    '''
    A subclass of :class:`thermosteam.Stream` with additional attributes
    for environmental impacts.
    
    See Also
    --------
    `thermosteam.Stream <https://thermosteam.readthedocs.io/en/latest/Stream.html>`_
    
    '''

    __slots__ = (*Stream.__slots__, '_impact_item')
    ticket_name = 'ss'

    def __init__(self, ID='', flow=(), phase='l', T=298.15, P=101325.,
                 units='kg/hr', price=0., thermo=None, impact_item=None,
                 **chemical_flows):
        super().__init__(ID=ID, flow=flow, phase=phase, T=T, P=P,
                         units=units, price=price, thermo=thermo,
                         **chemical_flows)
        
        if impact_item:
            impact_item._linked_stream = self
        self._impact_item = impact_item

    def copy(self, ID=''):
        '''
        Copy the information of another stream.
        
        Parameters
        ----------
        ID : str
            ID of the new stream, a default ID will be assigned if not provided.
        
        
        .. note::
            
            [1] Price of the original stream is not copied.
            
            [2] If the original stream has an :class:`~.StreamImpactItem`,
            then a new :class:`~.StreamImpactItem` will be created for the new stream
            and the new impact item will be linked to the original impact item.
        '''
        new = super().copy()
        if self.impact_item:
            self.impact_item.copy(stream=new)
        return new
        
    __copy__ = copy

    @staticmethod
    def from_stream(cls, stream, **kwargs):
        '''
        Cast a :class:`thermosteam.Stream` or :class:`biosteam.utils.MissingStream`
        to the designated equivalent.
        '''
        if isinstance(stream, MissingStream):
            new = MissingSanStream.__new__(MissingSanStream)
            return new
            
        elif not isinstance(stream, cls):
            stream.registry.untrack((stream,))
            ID = '' if (stream.ID[0] == 's' and stream.ID[1:].isnumeric()) else stream.ID
            new = cls.__new__(cls, **kwargs)
            new.ID = ID
            new._link = None
            new._sink = stream._sink
            new._source = stream._source
            new._thermo = stream._thermo
            new._imol = stream._imol.copy()
            new._thermal_condition = stream._thermal_condition.copy()
            new.reset_cache()
            new.price = 0
            new.impact_item = None
                
            stream._link = stream._sink = stream._source = None
            return new        

        else:
            return stream


    @property
    def impact_item(self):
        '''[StreamImpactItem] The :class:`StreamImpactItem` this waste stream is linked to.'''
        return self._impact_item
    @impact_item.setter
    def impact_item(self, i):
        self._impact_item = i
        if i:
            i.linked_stream = self

    @property
    def components(self):
        return self.chemicals


# %%

class MissingSanStream(MissingStream):
    '''
    A subclass of :class:`biosteam.MissingStream`, create a special object
    that acts as a dummy until replaced by an actual :class:`SanStream`.
    '''

    @property
    def impact_item(self):
        return None

    def __repr__(self):
        return '<MissingSanStream>'
    
    def __str__(self):
        return 'missing sanstream'





