#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Joy Zhang <joycheung1994@gmail.com>
    Yalin Li <zoe.yalin.li@gmail.com>

Part of this module is based on the Thermosteam package:
https://github.com/BioSTEAMDevelopmentGroup/thermosteam

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''


# %%

import numpy as np
import flexsolve as flx
from free_properties import PropertyFactory, property_array
from thermosteam import settings, indexer
from biosteam.utils import MissingStream
from . import SanStream, MissingSanStream
from .utils import auom, copy_attr, WasteStreamScope
from biosteam.utils import Scope
from warnings import warn

__all__ = ('WasteStream', 'MissingWasteStream',)


_defined_composite_vars = ('COD', 'BOD5', 'BOD', 'uBOD', 'NOD', 'ThOD', 'cnBOD',
                           'C', 'N', 'P', 'K', 'Mg', 'Ca', 'solids', 'charge')

_common_composite_vars = ('_COD', '_BOD', '_uBOD', '_TC', '_TOC', '_TN',
                          '_TKN', '_TP', '_TK', '_TMg', '_TCa',
                          '_dry_mass', '_charge', '_ThOD', '_cnBOD')

_ws_specific_slots = (*_common_composite_vars,
                      '_pH', '_SAlk', '_ratios',
                      # '_stream_impact_item', (pls keep this here, might be useful in debugging)
                      '_state', '_dstate', '_scope')

# Used in the `composite` method
_default_cmp_IDs = {
    'S_H2', 'S_CH4', 'S_CH3OH', 'S_Ac', 'S_Prop', 'S_F', 'S_U_Inf', 'S_U_E',
    'C_B_Subst', 'C_B_BAP', 'C_B_UAP', 'C_U_Inf', 'X_B_Subst', 'X_OHO_PHA',
    'X_GAO_PHA', 'X_PAO_PHA', 'X_GAO_Gly', 'X_PAO_Gly', 'X_OHO', 'X_AOO',
    'X_NOO', 'X_AMO', 'X_PAO', 'X_MEOLO', 'X_FO', 'X_ACO', 'X_HMO', 'X_PRO',
    'X_U_Inf', 'X_U_OHO_E', 'X_U_PAO_E', 'X_Ig_ISS', 'X_MgCO3', 'X_CaCO3',
    'X_MAP', 'X_HAP', 'X_HDP', 'X_FePO4', 'X_AlPO4', 'X_AlOH', 'X_FeOH',
    'X_PAO_PP_Lo', 'X_PAO_PP_Hi', 'S_NH4', 'S_NO2', 'S_NO3', 'S_PO4',
    'S_K', 'S_Ca', 'S_Mg', 'S_CO3', 'S_N2', 'S_O2', 'S_CAT', 'S_AN', 'H2O'
    }
_specific_groups = {'S_VFA': {'S_Ac', 'S_Prop'},
                    'X_Stor': {'X_OHO_PHA', 'X_GAO_PHA', 'X_PAO_PHA',
                              'X_GAO_Gly', 'X_PAO_Gly'},
                    'X_ANO': {'X_AOO', 'X_NOO'},
                    'X_Bio': {'X_OHO', 'X_AOO', 'X_NOO', 'X_AMO', 'X_PAO',
                              'X_MEOLO', 'X_ACO', 'X_HMO', 'X_PRO', 'X_FO'},
                    'S_NOx': {'S_NO2', 'S_NO3'},
                    'X_PAO_PP': {'X_PAO_PP_Lo', 'X_PAO_PP_Hi'}}

_default_ratios = {'iHi_XPAOPP': 0,
                   'iCB_XCB': 0.15,
                   'iBAP_CB': 0.,
                   'iUAP_CB': 0.,
                   'iCUInf_XCUInf': 0.,
                   'iSUInf_SU': 1.,
                   'iXUOHOE_XUE': None,}


vol_unit = auom('L/hr')
conc_unit = auom('mg/L')


# %%

# =============================================================================
# Util functions
# =============================================================================

_load_components = settings.get_chemicals
_load_thermo = settings.get_thermo

def to_float(stream, slot):
    return 0. if getattr(stream, slot) is None else float(getattr(stream, slot))

# functions for calibrations of N, P contents and BOD:COD ratio of certain components
def _calib_SF_iN(components, concentrations, STKN):
    cmp_c = np.asarray([v for v in concentrations.values()])
    SN = (cmp_c * components.i_N * (components.s + components.c)).sum()
    SF_N = concentrations['S_F'] * components.S_F.i_N
    SNOx_N = concentrations['S_NO2'] * components.S_NO2.i_N + concentrations['S_NO3'] * components.S_NO3.i_N
    other_stkn = SN - SF_N - SNOx_N
    SF_N = STKN - other_stkn
    if SF_N < 0:
        raise ValueError("Negative N content {SF_N} for S_F was estimated.")
    return SF_N/concentrations['S_F']

def _calib_XBsub_iN(components, concentrations, XTKN):
    cmp_c = np.asarray([v for v in concentrations.values()])
    other_xtkn = (cmp_c * components.i_N * components.x).sum() - concentrations['X_B_Subst'] * components.X_B_Subst.i_N
    XB_Subst_N = XTKN - other_xtkn
    if XB_Subst_N < 0:
        raise ValueError(f"Negative N content {XB_Subst_N} for X_B_Subst was estimated.")
    return XB_Subst_N/concentrations['X_B_Subst']


def _calib_XBsub_iP(components, concentrations, TP):
    cmp_c = np.asarray([v for v in concentrations.values()])
    other_p = (cmp_c * components.i_P).sum() - concentrations['X_B_Subst'] * components.X_B_Subst.i_P
    XB_Subst_P = TP - other_p
    if XB_Subst_P < 0:
        raise ValueError(f"Negative P content {XB_Subst_P} for X_B_Subst was estimated.")
    return XB_Subst_P/concentrations['X_B_Subst']

def _calib_XBsub_fBODCOD(components, concentrations, substrate_IDs, BOD):
    cmp_c = np.asarray([v for v in concentrations.values()])
    c_sub = np.asarray([v for k,v in concentrations.items() if k in substrate_IDs])
    XB_sub = components.subgroup(substrate_IDs)
    other_BOD = (cmp_c * (components.x + components.c + components.s) * components.f_BOD5_COD).sum() - (c_sub * XB_sub.f_BOD5_COD).sum()
    fbodtocod_sub = (BOD - other_BOD)/c_sub.sum()
    if fbodtocod_sub > 1 or fbodtocod_sub < 0:
        raise ValueError(f"BOD5-to-COD ratio {fbodtocod_sub} for X_B_Subst and X_Stor was estimated out of range [0,1].")
    return fbodtocod_sub


# Indexer for nicer display
@property
def group_conc_compositions(self):
    raise AttributeError('cannot set group by concentratino')

ComponentConcentrationIndexer, ConcentrationIndexer = \
    indexer._new_Indexer('Concentration', 'mg/L', group_conc_compositions)

ChemicalMolarFlowIndexer = indexer.ChemicalMolarFlowIndexer

@PropertyFactory(slots=('name', 'mol', 'index', 'F_vol', 'MW',
                        'phase', 'phase_container'))
def ConcentrationProperty(self):
    '''Concentration flows, in mg/L (g/m3).'''
    f_mass = self.mol[self.index] * self.MW
    phase = self.phase or self.phase_container.phase
    if phase != 'l':
        raise AttributeError('Concentration only valid for liquid phase.')
    V_sum = self.F_vol
    # if V_sum==0:
    #     raise RuntimeError('WasteStream is empty, concentration cannot be calculated.')
    return 1000. * f_mass / V_sum if f_mass else 0.
@ConcentrationProperty.setter
def ConcentrationProperty(self, value):
    raise AttributeError('Cannot set flow rate by concentration.')

def by_conc(self, TP):
    '''
    Return a ComponentConcentrationIndexer that references this object's
    molar and volume data (volume relies on molar).

    Parameters
    ----------
    TP : ThermalCondition

    '''
    cmps = self.chemicals
    mol = self.data
    F_vol = self.by_volume(TP).data.sum()
    conc = np.zeros_like(mol, dtype=object)
    for i, cmp in enumerate(cmps):
        conc[i] = ConcentrationProperty(cmp.ID, mol, i, F_vol, cmp.MW,
                                        None, self._phase)
    self._data_cache['conc', TP] = \
    conc = ComponentConcentrationIndexer.from_data(
        property_array(conc), self._phase, cmps, False)
    return conc
indexer.ChemicalMolarFlowIndexer.by_conc = by_conc

def by_conc(self, TP):
    '''
    Raise an error for attempt multi-phase usage
    as concentration only valid for liquid).
    '''
    raise AttributeError('Concentration only valid for liquid phase.')

indexer.MolarFlowIndexer.by_conc = by_conc
del by_conc, PropertyFactory


# %%

# =============================================================================
# Define the WasteStream class
# =============================================================================

class WasteStream(SanStream):
    '''
    A subclass of :class:`~.SanStream` with additional attributes and methods
    for wastewater.

    .. note::

        [1] Parameters below only include the ones additional to those of
        :class:`thermosteam.Stream`.

        [2] Properties not defined during during initiation will be calculated
        using :func:`WasteStream.composite`.


    Parameters
    ----------
    pH : float
        pH of the stream, unitless.
    SAlk : float
        Alkalinity, in meq/L (or mmol HCO3-/L). Assumed to be mainly bicarbonate.
    COD : float
        Chemical oxygen demand, in mg/L.
    BOD : float
        Biochemical oxygen demand, in mg/L, same as BOD5.
    uBOD : float
        Ultimate biochemical oxygen demand, in mg/L.
    cnBOD : float
        Carbonaceous nitrogenous BOD, in mg/L. Biochemical oxygen demand including nitrification.
    ThOD : float
        Theoretical oxygen demand, in mg/L.
    TC : float
        Total carbon, in mg/L.
    TOC : float
        Total organic carbon, in mg/L.
    TN : float
        Total nitrogen, in mg/L.
    TKN : float
        Total Kjeldahl nitrogen, in mg/L.
    TP : float
        Total phosphorus, in mg/L.
    TK : float
        Total potassium, in mg/L.
    TMg : float
        Total magnesium, in mg/L.
    TCa : float
        Total calcium, in mg/L.
    dry_mass : float
        Total solids, dry mass of dissolved and suspended solids, in mg/L.
    charge : float
        TO BE IMPLEMENTED
    ratios : float
        TO BE IMPLEMENTED
    stream_impact_item : :class:`StreamImpactItem`
        The :class:`StreamImpactItem` this stream is linked to.
    component_flows : kwargs
        Component flow data.


    Examples
    --------
    `WasteStream <https://qsdsan.readthedocs.io/en/latest/tutorials/3_WasteStream.html>`_

    See Also
    --------
    `thermosteam.Stream <https://thermosteam.readthedocs.io/en/latest/Stream.html>`_
    '''

    __slots__ = SanStream.__slots__ + _ws_specific_slots
    _default_ratios = _default_ratios
    ticket_name = 'ws'

    def __init__(self, ID='', flow=(), phase='l', T=298.15, P=101325.,
                 units='kg/hr', price=0., thermo=None,
                 pH=7., SAlk=2.5, COD=None, BOD=None, uBOD=None,
                 ThOD=None, cnBOD=None,
                 TC=None, TOC=None, TN=None, TKN=None, TP=None, TK=None,
                 TMg=None, TCa=None, dry_mass=None, charge=None, ratios=None,
                 stream_impact_item=None, **component_flows):
        SanStream.__init__(self=self, ID=ID, flow=flow, phase=phase, T=T, P=P,
                           units=units, price=price, thermo=thermo,
                           stream_impact_item=stream_impact_item, **component_flows)

        self._init_ws(pH, SAlk, COD, BOD, uBOD, TC, TOC, TN, TKN,
                      TP, TK, TMg, TCa, ThOD, cnBOD, dry_mass, charge, ratios)

    def _init_ws(self, pH=7., SAlk=None, COD=None, BOD=None,
                 uBOD=None, ThOD=None, cnBOD=None,
                 TC=None, TOC=None, TN=None, TKN=None,
                 TP=None, TK=None, TMg=None, TCa=None,
                 dry_mass=None, charge=None, ratios=None):

        self._pH = pH
        self._SAlk = SAlk
        self._COD = COD
        self._BOD = BOD
        self._uBOD = uBOD
        self._ThOD = ThOD
        self._cnBOD = cnBOD
        self._TC = TC
        self._TOC = TOC
        self._TN = TN
        self._TKN = TKN
        self._TP = TP
        self._TK = TK
        self._TMg = TMg
        self._TCa = TCa
        self._dry_mass = dry_mass
        self._charge = charge
        self._ratios = ratios
        self._state = None
        self._dstate = None

    @staticmethod
    def from_stream(cls, stream, ID=None, **kwargs):
        '''
        Cast a :class:`thermosteam.Stream` or :class:`biosteam.utils.MissingStream`
        to :class:`WasteStream` or :class:`MissingWasteStream`.

        Parameters
        ----------
        cls : obj
            class of the stream to be created.
        stream : :class:`thermosteam.Stream`
            The original stream.
        ID : str
            If not provided, will use the ID of the original stream.
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
        >>> ws = qs.WasteStream.from_stream(qs.WasteStream, s, ID='ws',
        ...                                 T=250, price=8)
        >>> ws.show()
        WasteStream: ws
         phase: 'l', T: 250 K, P: 101325 Pa
         flow (g/hr): H2O  1.8e+06
         WasteStream-specific properties:
          pH         : 7.0
         Component concentrations (mg/L):
          H2O          1012485.4
        >>> ws.price
        8.0
        '''
        # Missing stream
        if isinstance(stream, MissingStream):
            new = MissingWasteStream.__new__(MissingWasteStream)
            new.__init__(stream._source, stream._sink)
            return new

        # An actual stream
        new = SanStream.from_stream(cls, stream, ID)
        for attr, val in kwargs.items():
            setattr(new, attr, val)
        new._init_ws()
        return new


    def show(self, T='K', P='Pa', flow='g/hr', composition=False, N=15,
             stream_info=True, details=True, concentrations='mg/L'):
        '''
        Print information related to this :class:`WasteStream`.

        Parameters
        ----------
        T : str, optional
            The unit for temperature. The default is 'K'.
        P : float, optional
            The unit for pressure. The default is 'Pa'.
        flow : str, optional
            The unit for the flow. The default is 'kg/hr'.
        composition : bool, optional
            Whether to show flow information of different :class:`Component` objects in
            this waste stream as a percentage. The default is False.
        N : int, optional
            Number of components to print out, when left as None,
            the number depends on the default of :class:`thermosteam`. The default is 15.
        stream_info : bool, optional
            Whether to print stream-specific information. The default is True.
        details : bool, optional
            Whether to show the all composite variables of this waste stream. The default is True.

        '''

        info = ''
        # Stream-related specifications
        if stream_info:
            super().show(None, T, P, flow, composition, N)
        else:
            info += self._basic_info()
            display_units = self.display_units
            T_units = T or display_units.T
            P_units = P or display_units.P
            info += self._info_phaseTP(self.phase, T_units, P_units)
        info += self._wastestream_info(details=details, concentrations=concentrations, N=N)
        print(info)

    _ipython_display_ = show


    def _wastestream_info(self, details=True, concentrations=None, N=15):
        _ws_info = ' WasteStream-specific properties:'
        # Wastewater-related properties are not relevant for gas or solids
        if self.phase != 'l':
            _ws_info += ' None for non-liquid waste streams\n'
        elif self.F_mass == 0:
            _ws_info += ' None for empty waste streams\n'
        else:
            _ws_info += '\n'
            # Only non-zero properties are shown
            _ws_info += int(bool(self.pH))*f'  pH         : {self.pH:.1f}\n'
            _ws_info += int(bool(self.SAlk))*f'  Alkalinity : {self.SAlk:.1f} mg/L\n'
            if details:
                _ws_info += int(bool(self.COD))   *f'  COD        : {self.COD:.1f} mg/L\n'
                _ws_info += int(bool(self.BOD))   *f'  BOD        : {self.BOD:.1f} mg/L\n'
                _ws_info += int(bool(self.TC))    *f'  TC         : {self.TC:.1f} mg/L\n'
                _ws_info += int(bool(self.TOC))   *f'  TOC        : {self.TOC:.1f} mg/L\n'
                _ws_info += int(bool(self.TN))    *f'  TN         : {self.TN:.1f} mg/L\n'
                # `TKN` not included as the users need to define that to include in TKN calculation
                # _ws_info += int(bool(self.TKN))   *f'  TKN        : {self.TKN:.1f} mg/L\n'
                _ws_info += int(bool(self.TP))    *f'  TP         : {self.TP:.1f} mg/L\n'
                _ws_info += int(bool(self.TK))    *f'  TK         : {self.TK:.1f} mg/L\n'
                # _ws_info += int(bool(self.charge))*f'  charge     : {self.charge:.1f} mmol/L\n'
            else:
                _ws_info += '  ...\n'
            _ws_info += self._concentration_info(unit=concentrations, N=N)
        return _ws_info

    def _concentration_info(self, unit='mg/L', N=15):
        if not isinstance(unit, str): return ''
        else:
            C_arr = self.get_mass_concentration(unit=unit)
            N_ID = sum(C_arr > 0)
            too_many_components = N_ID > N
            N_max =  min(N_ID, N)
            IDs = self.components.IDs
            lengths = [len(ID) for ID in IDs]
            maxlen = max(lengths) + 2
            line = f' Component concentrations ({unit}):\n'
            n = 0
            for i in range(len(IDs)):
                if n == N_max: break
                if C_arr[i] > 0:
                    n += 1
                    spaces = ' ' * (maxlen - lengths[i])
                    line += '  ' + IDs[i] + spaces + f'{C_arr[i]:.1f}\n'
            if too_many_components:
                line += '  ...\n'
            line = line.rstrip('\n')
            return line


    @property
    def ratios(self):
        '''
        The ratios used for estimating waste stream composition based on user input upon initialization.
        Only meaningful for creating a :class:`WasteStream` object from scratch.
        If not used or specified, default as None.
        '''
        return self._ratios
    @ratios.setter
    def ratios(self, ratios):
        r = self._ratios or WasteStream._default_ratios
        for name, ratio in ratios.items():
            if name not in r.keys():
                raise ValueError(f'Cannot identify ratio named "{name}".'
                                 f'Must be one of {r.keys()}.')
            elif isinstance(ratio, (int, float)) and (ratio > 1 or ratio < 0):
                raise ValueError(f"ratio {name}: {ratio} is out of range [0,1].")
            r[name] = ratio
        self._ratios = r

    def composite(self, variable, subgroup=None, particle_size=None,
                  degradability=None, organic=None, volatile=None,
                  specification=None, unit=None):
        """
        Calculate select composite variable with the provided constraints.

        Parameters
        ----------
        variable : str
            The composite variable to calculate. One of the followings:
            ("COD", "BOD5", "BOD", "uBOD", "NOD", "ThOD", "cnBOD",
            "C", "N", "P", "K", "Mg", "Ca", "solids", "charge").
        subgroup : tuple(str or obj), optional
            Iterable of :class:`CompiledComponents` (or their IDs)
            which the composite variable will be calculated for.
        particle_size : "g", "s", "c", or "x", optional
            Dissolved gas ("g"), soluble ("s"), colloidal ("c"), particulate ("x").
            The default is None.
        degradability : "rb", "sb", "b" or "u", optional
            Readily biodegradable ("rb"), slowly biodegradable ("sb"),
            biodegradable ("b"), or undegradable ("u"). The default is None.
        organic : bool, optional
            Organic (True) or inorganic (False). The default is None.
        volatile : bool, optional
            Volatile (True) or involatile (False). The default is None.
        specification : str, optional
            A group of :class:`Component` defined through :func:`Component.define_group`.

            If using the default components (through :func:`Components.load_default`)
            or a subgroup of these default components,
            the specification can be one of
            ("S_VFA", "X_Stor", "X_ANO", "X_Bio", "S_NOx", "X_PAO_PP").

            You can also use :func:`CompileComponents.define_group` to define
            your own specification groups.

            The composite variable will be calculated for the
            intersection of all designated constrains
            (i.e., `subgroup`, `particle_size`, `degradability`, `organic`, `volatile`, and `specification`).
        unit : str
            The unit that the result will be returned in.
            If not provided, result will be in mg/L except for charge (mmol/L).

        Returns
        -------
        value : float
            The estimated value of the composite variable in the desired unit.

        Examples
        --------
        >>> from qsdsan import set_thermo, Components, WasteStream
        >>> cmps = Components.load_default()
        >>> set_thermo(cmps)
        >>> ws = WasteStream.codstates_inf_model('ws', flow_tot=1000, pH=6.8, COD=500, TP=11)
        >>> # To calculate the particulate BOD (i.e., xBOD) of the WasteStream object,
        >>> # you just need to specify the composite variable as "BOD", and particle size as "x"
        >>> ws.composite('BOD', particle_size='x') # doctest: +ELLIPSIS
        152.83...
        >>> # You can also adjust the unit you want the result to be in
        >>> ws.composite('BOD', particle_size='x', unit='g/L') # doctest: +ELLIPSIS
        0.15283...
        >>> # Biomass COD
        >>> ws.composite('COD', specification='X_Bio')
        0.0
        >>> # Nitrogen as nitrate/nitrite
        >>> ws.composite('N', specification='S_NOx')
        0.0
        >>> # Soluble TKN, note that the `TKN` attribute is pre-defined in the default components
        >>> # as all components that are not "S_N2", "S_NO2", and "S_NO3"
        >>> # For your own components, you can pre-define the group using
        >>> # the `define_group` function (e.g., `cmps.define_group`)
        >>> ws.composite('N', subgroup=cmps.TKN, particle_size='s') # doctest: +ELLIPSIS
        27.77...
        >>> ws.composite('C', organic=True) # doctest: +ELLIPSIS
        160.00...
        >>> ws.composite('solids', particle_size='s') # doctest: +ELLIPSIS
        945.69...

        See Also
        --------
        :func:`CompiledComponents.define_group`

        """
        _get = getattr
        if self.F_vol == 0.:
            return 0.

        if variable not in _defined_composite_vars:
            raise KeyError(f"Undefined composite variable {variable},"
                           f"Must be one of {_defined_composite_vars}.")

        # Can only be used for liquid
        if not self.phase == 'l':
            raise RuntimeError('Only liquid streams can use the `composite` method, '
                               f'the current WasteStream {self.ID} is {self.phase}.')

        isa = isinstance
        all_cmps = self.components
        subgroup = subgroup or all_cmps
        subgroup_IDs = {i if isa(i, str) else i.ID for i in subgroup} # use set for easier element comparison
        subgroup_IDs = subgroup_IDs.difference({'Water', 'H2O'}) # remove water

        if specification:
            try:
                specified_IDs = set(_get(all_cmps, specification))
            except AttributeError: # no pre-defined groups
                try:
                    specified_IDs = _specific_groups[specification]
                except KeyError: # specification not in the default ones
                    raise KeyError(f"Undefined specification {specification}. "
                                   f"`specification` must be one of {_specific_groups.keys()}. "
                                   "Use the `subgroup` argument instead or "
                                   "define the specification group using "
                                   "`CompiledComponents.define_group`.")
                # Issue a warning if the subgroup contains components outside of the default ones
                if not subgroup_IDs.issubset(_default_cmp_IDs):
                    warn(f'{specification} is defined with regards to the set of default component IDs. '
                          'Consider using the `subgroup` argument instead of '
                          '`specification`, or define the specification group using '
                          '`CompiledComponents.define_group` '
                          'if different sets of component IDs are used.')
            IDs = tuple(subgroup_IDs.intersection(specified_IDs))
        else:
            IDs = tuple(subgroup_IDs)

        cmps = all_cmps.subgroup(IDs)
        cmp_c = self.iconc[IDs] # [mg/L]
        exclude_gas = cmps.s + cmps.c + cmps.x

        if variable == 'COD':
            var = cmps.i_COD * cmp_c * exclude_gas * (cmps.i_COD >= 0)
        elif variable == 'uBOD':
            var = cmps.i_COD * cmps.f_uBOD_COD * cmp_c * exclude_gas * (cmps.i_COD >= 0)
        elif variable in ('BOD5', 'BOD'):
            var = cmps.i_COD * cmps.f_BOD5_COD * cmp_c * exclude_gas * (cmps.i_COD >= 0)
        elif variable == 'NOD':
            var = cmps.i_NOD * cmp_c * exclude_gas
        elif variable == 'ThOD':
            var = (cmps.i_NOD + cmps.i_COD * (cmps.i_COD >= 0)) * cmp_c
        elif variable == 'cnBOD':
            var = (cmps.i_NOD + cmps.i_COD * cmps.f_BOD5_COD * (cmps.i_COD >= 0)) * cmp_c * exclude_gas
        elif variable == 'C':
            var = cmps.i_C * cmp_c
        elif variable == 'N':
            var = cmps.i_N * cmp_c * exclude_gas
        elif variable == 'P':
            var = cmps.i_P * cmp_c
        elif variable == 'K':
            var = cmps.i_K * cmp_c
        elif variable == 'Mg':
            var = cmps.i_Mg * cmp_c
        elif variable == 'Ca':
            var = cmps.i_Ca * cmp_c
        elif variable == 'solids':
            var = cmps.i_mass * cmp_c * exclude_gas
            if volatile != None:
                if volatile: var *= cmps.f_Vmass_Totmass
                else: var *= 1-cmps.f_Vmass_Totmass
        else:
            var = cmps.i_charge * cmp_c

        dummy = np.ones_like(cmp_c)
        if particle_size:
            if particle_size == 'g':
                dummy *= 1-exclude_gas
            else:
                dummy *= _get(cmps, particle_size)

        if degradability:
            if degradability == 'u': dummy *= 1-cmps.b
            elif degradability == 'b': dummy *= cmps.b
            elif degradability == 'rb': dummy *= cmps.rb
            else: dummy *= cmps.b-cmps.rb

        if organic != None:
            if organic: dummy *= cmps.org
            else: dummy *= 1-cmps.org

        result = (dummy*var).sum()
        if not unit:
            return result

        converted = auom('mg/L').convert(result, unit) if variable != 'charge' \
            else auom('mmol/L').convert(result, unit)
        return converted


    def _liq_sol_properties(self, prop, value):
        if self.phase != 'g':
            return getattr(self, '_'+prop) or value
        else:
            raise AttributeError(f'{self.phase} phase waste stream does not have {prop}.')

    @property
    def pH(self):
        '''[float] pH, unitless.'''
        return self._liq_sol_properties('pH', 7.)

    @property
    def SAlk(self):
        '''[float] Alkalinity, in meq/L (or mmol HCO3-/L). Assumed to be mainly bicarbonate.'''
        return self._liq_sol_properties('SAlk', 0.)

    @property
    def COD(self):
        '''[float] Chemical oxygen demand in mg/L.'''
        return self._liq_sol_properties('COD', self.composite('COD'))

    @property
    def BOD(self):
        '''[float] Biochemical oxygen demand in mg/L. Same as BOD5.'''
        return self._liq_sol_properties('BOD', self.composite('BOD'))

    @property
    def BOD5(self):
        '''[float] 5-day biochemical oxygen demand, in mg/L. Same as BOD.'''
        return self.BOD

    @property
    def uBOD(self):
        '''[float] Ultimate biochemical oxygen demand, in mg/L.'''
        return self._liq_sol_properties('uBOD', self.composite('uBOD'))

    @property
    def cnBOD(self):
        '''[float] Carbonaceous nitrogenous BOD, in mg/L. Biochemical oxygen demand including nitrification.'''
        return self._liq_sol_properties('cnBOD', self.composite('cnBOD'))

    @property
    def ThOD(self):
        '''[float] Theoretical oxygen demand, in mg/L.'''
        return self._liq_sol_properties('ThOD', self.composite('ThOD'))

    #!!! Maybe include C_frac, etc. to calculate C_mass/F_mass - valid for all phases
    # Or a function to calculate it?
    @property
    def TC(self):
        '''[float] Total carbon, in mg/L.'''
        return self._liq_sol_properties('TC', self.composite('C'))

    @property
    def TOC(self):
        '''[float] Total organic carbon, in mg/L.'''
        return self._liq_sol_properties('TOC', self.composite('C', organic=True))

    @property
    def TN(self):
        '''[float] Total nitrogen, in mg/L.'''
        return self._liq_sol_properties('TN', self.composite('N'))

    @property
    def TKN(self):
        '''[float] Total Kjeldahl nitrogen, in mg/L.'''
        try:
            subgroup = self.components.TKN
        except:
            warn('No `TKN` group defined for the current `CompiledComponents`, '
                 '`TKN` will be calculated using `TN` minus concentrations of '
                 '"S_N2", "S_NO2", and "S_NO3". '
                 'Use `define_group` to define the TKN groups if want to calculate otherwise.')
            subgroup = [i.ID for i in self.components if i.ID not in ('S_N2','S_NO2','S_NO3')]
        return self._liq_sol_properties('TKN', self.composite('N', subgroup=subgroup))

    @property
    def TP(self):
        '''[float] Total phosphorus, in mg/L.'''
        return self._liq_sol_properties('TP', self.composite('P'))

    @property
    def TK(self):
        '''[float] Total potassium, in mg/L.'''
        return self._liq_sol_properties('TK', self.composite('K'))

    @property
    def TMg(self):
        '''[float] Total magnesium, in mg/L.'''
        return self._liq_sol_properties('TMg', self.composite('Mg'))

    @property
    def TCa(self):
        '''[float] Total calcium, in mg/L.'''
        return self._liq_sol_properties('TCa', self.composite('Ca'))

    @property
    def dry_mass(self):
        '''[float] Total solids, dry mass of dissolved and suspended solids, in mg/L.'''
        return self._liq_sol_properties('solids', self.composite('solids'))

    # TODO: calibrate Charge when weak acids are involved
    # @property
    # def charge(self):
    #     return self._liq_sol_properties('charge', self.composite('charge'))


    @property
    def iconc(self):
        '''[Indexer] Mass concentrations, in mg/L (g/m3).'''
        return self._imol.by_conc(self._thermal_condition)

    @property
    def conc(self):
        '''[property_array] Mass concentrations, in mg/L (g/m3).'''
        return self.iconc.data

    @property
    def Conc(self):
        '''[property_array] Mass concentrations, in mg/L (g/m3), same as `conc`.'''
        return self.conc

    @property
    def density(self):
        '''[float] Density of the stream, in mg/L (kg/m3).'''
        return self.F_mass/self.F_vol


    def copy(self, new_ID='', copy_price=False, copy_impact_item=False,
             ws_properties=True):
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
        ws_properties : bool
            Whether to copy wastewater-related properties to the new stream.
        '''
        new = SanStream.copy(self, new_ID=new_ID,
                             copy_price=copy_price,
                             copy_impact_item=copy_impact_item)
        if ws_properties:
            new._init_ws()
            # Skip `_scope`, if users want it to be scoped,
            # then calling the `scope` property will automatically make the property
            new = copy_attr(new, self, skip=(*SanStream.__slots__, '_scope'))
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
        SanStream.copy_like(self, other,
                            copy_price=copy_price,
                            copy_impact_item=copy_impact_item)

        if not isinstance(other, WasteStream):
            return

        for slot in _ws_specific_slots:
            if slot == '_scope': continue # see notes in `copy`
            value = getattr(other, slot)
            setattr(self, slot, value)


    def proxy(self, ID=None):
        '''
        Return a new stream that shares all data with this one.

        Note that unlike other properties, the price of the two streams are not
        connected, i.e., the price of the new stream will be the same as the
        original one upon creation, but then they can be different.

        Parameters
        ----------
        ID : str
            ID of the new proxy stream.

        Examples
        --------
        >>> from qsdsan import set_thermo, WasteStream
        >>> from qsdsan.utils import load_example_cmps
        >>> cmps = load_example_cmps()
        >>> set_thermo(cmps)
        >>> ws1 = WasteStream('ws1', Water=100, NaCl=1)
        >>> ws2 = ws1.proxy('ws2')
        >>> ws2.mol is ws1.mol
        True
        >>> # Note that concentration is always calculated (mass/vol) upon request,
        >>> # so `ws2.conc is ws1.conc` will return False
        >>> import numpy as np
        >>> np.all(ws2.conc==ws1.conc)
        True
        >>> ws1.imol['Water'] = 10000
        >>> np.all(ws2.conc==ws1.conc)
        True
        '''
        new = SanStream.proxy(self, ID=ID)
        for slot in _ws_specific_slots:
            if slot == '_scope': continue # see notes in `copy`
            value = getattr(self, slot)
            setattr(new, slot, value)
        return new


    def mix_from(self, others, **kwargs):
        '''
        Update this stream to be a mixture of other streams,
        initial content of this stream will be ignored.

        Parameters
        ----------
        others : iterable
            Can contain :class:`thermosteam.Stream`, :class:`SanStream`,
            or :class:`~.WasteStream`.
        kwargs : dict
            Other keyword arguments that will be passed to :func:`Stream.mix_from`


        .. note::

            Price and impact item are not included.


        Examples
        --------
        >>> import qsdsan as qs
        >>> cmps = qs.Components.load_default()
        >>> qs.set_thermo(cmps)
        >>> s = qs.Stream('s', H2O=100, price=5, units='kg/hr')
        >>> ss = qs.SanStream('ss', S_O2=100, units='kg/hr')
        >>> ws = qs.WasteStream('ws')
        >>> ws.mix_from((s, ss))
        >>> ws.show()
        WasteStream: ws
         phase: 'l', T: 298.15 K, P: 101325 Pa
         flow (g/hr): S_O2  1e+05
                      H2O   1e+05
         WasteStream-specific properties:
          pH         : 7.0
         Component concentrations (mg/L):
          S_O2         303021.4
          H2O          303021.4
        '''
        SanStream.mix_from(self, others, **kwargs)

        for slot in _ws_specific_slots:
            if not hasattr(self, slot):
                continue
            #!!! This needs reviewing, might not be good to calculate some
            # attributes like pH
            # tot = None
            try:
                tot = sum(to_float(i, slot)*i.F_vol for i in others if hasattr(i, slot))
                setattr(self, slot, tot/self.F_vol)
            except:
                setattr(self, slot, None)


    def get_TDS(self, include_colloidal=True):
        '''
        Total dissolved solids (TDS).

        Parameters
        ----------
        include_colloidal : bool, optional
            Whether to include colloidal components as TDS. The default is True.

        Returns
        -------
        TDS : float
            In mg/L.
        '''
        TDS = self.composite('solids', particle_size='s')
        if include_colloidal:
            TDS += self.composite('solids', particle_size='c')
        return TDS

    def get_TSS(self, include_colloidal=False):
        '''
        Total suspended solids (TSS).

        Parameters
        ----------
        include_colloidal : bool, optional
            Whether to include colloidal components as TSS. The default is False.

        Returns
        -------
        TSS : float
            In mg/L.
        '''
        TSS = self.composite('solids', particle_size='x')
        if include_colloidal:
            TSS += self.composite('solids', particle_size='c')
        return TSS

    def get_VSS(self, include_colloidal=False):
        '''[float] Volatile suspended solids, in mg/L.'''
        VSS = self.composite('solids', particle_size='x', volatile=True)
        if include_colloidal:
            VSS += self.composite('solids', particle_size='c', volatile=True)
        return VSS

    def get_ISS(self):
        '''[float] Inorganic/involatile suspended solids, in mg/L.'''
        return self.composite('solids', particle_size='x', volatile=False)


    def get_mass_concentration(self, unit='g/m3', IDs=None):
        '''
        Get mass concentrations in given unit.

        Parameters
        ----------
        unit : str, optional
            Unit of measure. The default is 'g/m3'.
        IDs : Iterable[str], optional
            IDs of components. When not specified, returns mass concentrations of
            all components in thermo.

        '''
        F_vol = self.F_vol
        if not F_vol: raise RuntimeError(f'{repr(self)} is empty')
        if not IDs: IDs = self.components.IDs
        C = self.imass[IDs]/F_vol*1e3      # in mg/L
        return C*conc_unit.conversion_factor(unit)


    def set_flow_by_concentration(self, flow_tot, concentrations, units=('L/hr', 'mg/L'),
                                  bulk_liquid_ID='H2O', atol=1e-5, maxiter=50):
        '''
        Set the mass flows of the WasteStream by specifying total volumetric flow and
        concentrations as well as identifying the component that constitutes the bulk liquid.

        Parameters
        ----------
        flow_tot : float
            Total volumetric flow of the WasteStream.
        concentrations : dict[str, float]
            Concentrations of components.
        units : Iterable[str]
            The first indicates the unit for the input total flow, the second
            indicates the unit for the input concentrations.
        bulk_liquid_ID : str, optional
            ID of the Component that constitutes the bulk liquid, e.g., the solvent.
            The default is 'H2O'.
        atol : float, optional
            The absolute tolerance of error in estimated WasteStream density.
            The default is 1e-5 kg/L.
        maxiter : int, optional
            The maximum number of iterations to estimate the flow of the bulk-liquid
            component and overall density of the WasteStream. The default is 50.


        '''
        if flow_tot == 0: raise RuntimeError(f'{repr(self)} is empty')
        if bulk_liquid_ID in concentrations.keys():
            for i in range(8):
                warn('\n', stacklevel=i)
            C_bulk = concentrations.pop(bulk_liquid_ID)
            warn(f'ignored concentration specified for {bulk_liquid_ID}:{C_bulk}')

        self.empty()
        f = conc_unit.conversion_factor(units[1]) # convert to mg/L
        Q_tot = flow_tot / vol_unit.conversion_factor(units[0])   # converted to L/hr
        C_arr = np.array(list(concentrations.values()))
        IDs = concentrations.keys()
        M_arr = C_arr/f*Q_tot*1e-6       # mg/L * L/hr /1e6 = kg/hr
        self.set_flow(M_arr, 'kg/hr', tuple(IDs))

        # Density of the mixture should be larger than the pure bulk liquid,
        # give it a factor of 100 for safety
        den = self.components[bulk_liquid_ID].rho(phase=self.phase, T=self.T, P=self.P)*1e-3  # bulk liquid density in [kg/L]
        bulk_ref_phase = self.components[bulk_liquid_ID].phase_ref
        if bulk_ref_phase != 'l':
            warn(f'Reference phase of liquid is "{bulk_ref_phase}", not "l", '
                 'flow might not be set correctly.')

        M_bulk_max = Q_tot * den * 100 # L/hr * kg/L = kg/hr

        M_bulk = flx.IQ_interpolation(f=self._Q_obj_f,
                                      x0=0, x1=M_bulk_max,
                                      xtol=0.01, ytol=atol, # ytol here is actually Q, but should be fine
                                      args=(bulk_liquid_ID, Q_tot),
                                      maxiter=maxiter, checkbounds=False)
        self.set_flow(M_bulk, 'kg/hr', bulk_liquid_ID)


    def _Q_obj_f(self, M_bulk, bulk_liquid_ID, target_Q):
        self.set_flow(M_bulk, 'kg/hr', bulk_liquid_ID)
        return self.F_vol*1e3 - target_Q

    @property
    def state(self):
        if self.isproduct(): return self._state
        else:
            i = self._sink._ins.index(self)
            return self._sink._ins_QC[i]

    @state.setter
    def state(self, y):
        if self.isproduct():
            try: self._state[:] = y
            except TypeError: self._state = y
        else:
            i = self._sink._ins.index(self)
            self._sink._ins_QC[i,:] = y

    @property
    def dstate(self):
        if self.isproduct(): return self._dstate
        else:
            i = self._sink._ins.index(self)
            return self._sink._ins_dQC[i]

    @dstate.setter
    def dstate(self, dy):
        if self.isproduct():
            try: self._dstate[:] = dy
            except TypeError: self._dstate = dy
            # decide whether to record state over time
        else:
            i = self._sink._ins.index(self)
            self._sink._ins_dQC[i,:] = dy

    @property
    def scope(self):
        """A tracker of the wastestream's time-series data during dynamic simulation."""
        if not hasattr(self, '_scope'):
            self._scope = WasteStreamScope(self)
        return self._scope

    @scope.setter
    def scope(self, s):
        if not isinstance(s, Scope):
            raise TypeError(f'{s} must be a `Scope`, not {type(s)}.')
        if self is not s.subject:
            if isinstance(s, WasteStreamScope):
                s = WasteStreamScope(self)
                warn(f'{self} is not {s}.subject, so {s} has been re-initialized '
                     f'as a {WasteStreamScope} object with {self} as its subject.')
            else:
                raise ValueError(f'The subject of {s} must be {self} not {s.subject}.')
        self._scope = s

    def _init_state(self):
        self.state = np.append(self.conc.astype('float64'), self.get_total_flow('m3/d'))
        self.dstate = np.zeros_like(self.state)

    def _state2flows(self):
        Q = self.state[-1]
        Cs = dict(zip(self.components.IDs, self.state[:-1]))
        Cs.pop('H2O', None)
        self.set_flow_by_concentration(Q, Cs, units=('m3/d', 'mg/L'))

    @classmethod
    def codstates_inf_model(cls, ID='', flow_tot=0., units = ('L/hr', 'mg/L'),
                            phase='l', T=298.15, P=101325., price=0., thermo=None,
                            pH=7., SAlk=10., ratios=None,
                            COD=430., TKN=40., TP=10., iVSS_TSS=0.75, iSNH_STKN=0.9,
                            S_NH4=25., S_NO2=0., S_NO3=0., S_PO4=8.,
                            S_Ca=140., S_Mg=50., S_K=28., S_CAT=3., S_AN=12., S_N2=18.,
                            frSUInf=0.05, frSF=0.2, frXCUInf=0.13,
                            frSUE=0., frSCH3OH=0., frSAc=0., frSProp=0.,
                            frXOHO=0., frXAOO=0., frXNOO=0., frXAMO=0., frXPAO=0.,
                            frXPRO=0., frXACO=0., frXHMO=0., frXMEOLO=0., frXFO=0.,
                            frXOHO_PHA=0., frXGAO_PHA=0., frXPAO_PHA=0.,
                            frXGAO_Gly=0., frXPAO_Gly=0., frXU_OHO_E=0., frXU_PAO_E=0.,
                            X_FePO4=0., X_AlPO4=0., X_FeOH=0., X_AlOH=0.,
                            X_MAP=0., X_HAP=0., X_HDP=0., X_PAO_PP=0.,
                            X_MgCO3=0., X_CaCO3=0., DO=0., S_H2=0., S_CH4=0.):
        '''
        Create a waste stream using the COD states influent characterization model as in [1]_,
        requires the default component set (can be created by :func:`~.Components.load_default`)

        Normal attributes of the waste streams (e.g., phase, T, Pa, price) can be passed as usual.

        Parameters
        ----------
        ID : str
            ID of the stream.
        flow_tot : float
            Total flowrates.
        units : tuple(str, str)
            Unit of the flowrate and component concentrations.

        Examples
        --------
        >>> from qsdsan import Components, WasteStream, set_thermo
        >>> cmps = Components.load_default()
        >>> set_thermo(cmps)
        >>> ws = WasteStream.codstates_inf_model('ws_codstates', flow_tot=100)
        >>> ws.show()
        WasteStream: ws_codstates
         phase: 'l', T: 298.15 K, P: 101325 Pa
         flow (g/hr): S_F        8.6
                      S_U_Inf    2.15
                      C_B_Subst  4
                      X_B_Subst  22.7
                      X_U_Inf    5.59
                      X_Ig_ISS   5.23
                      S_NH4      2.5
                      S_PO4      0.8
                      S_K        2.8
                      S_Ca       14
                      S_Mg       5
                      S_CO3      12
                      S_N2       1.8
                      S_CAT      0.3
                      S_AN       1.2
                      ...
         WasteStream-specific properties:
          pH         : 7.0
          Alkalinity : 10.0 mg/L
          COD        : 430.0 mg/L
          BOD        : 221.8 mg/L
          TC         : 265.0 mg/L
          TOC        : 137.6 mg/L
          TN         : 40.0 mg/L
          TP         : 10.0 mg/L
          TK         : 28.0 mg/L
         Component concentrations (mg/L):
          S_F          86.0
          S_U_Inf      21.5
          C_B_Subst    40.0
          X_B_Subst    226.6
          X_U_Inf      55.9
          X_Ig_ISS     52.3
          S_NH4        25.0
          S_PO4        8.0
          S_K          28.0
          S_Ca         140.0
          S_Mg         50.0
          S_CO3        120.0
          S_N2         18.0
          S_CAT        3.0
          S_AN         12.0
          ...

        References
        ----------
        [1] Hydromantis. GPS-X Technical Reference - v8.0.
        Hydromantis Environmental Software Solutions, Inc. 2019.
        https://www.hydromantis.com/help/GPS-X/docs/8.0/Technical/index.html

        '''
        if thermo: cmps = thermo.chemicals
        else:
            cmps = _load_components()
            thermo = _load_thermo()

        cmp_dct = dict.fromkeys(cmps.IDs, 0.)

        new = cls(ID=ID, phase=phase, T=T, P=P, units='kg/hr', price=price,
                  thermo=thermo, pH=pH, SAlk=SAlk)

        if ratios: new.ratios = ratios
        else: new.ratios = WasteStream._default_ratios
        r = new._ratios

        #************ user-defined states **************
        cmp_dct['S_H2'] = S_H2
        cmp_dct['S_CH4'] = S_CH4
        cmp_dct['S_N2'] = S_N2
        cmp_dct['S_O2'] = DO
        cmp_dct['S_NH4'] = S_NH4
        cmp_dct['S_NO2'] = S_NO2
        cmp_dct['S_NO3'] = S_NO3
        cmp_dct['S_PO4'] = S_PO4
        cmp_dct['S_CO3'] = SAlk * 12 * conc_unit.conversion_factor(units[1])             # 1 meq/L SAlk ~ 1 mmol/L HCO3- ~ 12 mg C/L (12 mg C/mmol HCO3-)
        cmp_dct['S_Ca'] = S_Ca
        cmp_dct['S_Mg'] = S_Mg
        cmp_dct['S_K'] = S_K
        cmp_dct['X_MAP'] = X_MAP
        cmp_dct['X_HAP'] = X_HAP
        cmp_dct['X_HDP'] = X_HDP
        cmp_dct['X_FePO4'] = X_FePO4
        cmp_dct['X_AlPO4'] = X_AlPO4
        cmp_dct['X_FeOH'] = X_FeOH
        cmp_dct['X_AlOH'] = X_AlOH
        cmp_dct['X_MgCO3'] = X_MgCO3
        cmp_dct['X_CaCO3'] = X_CaCO3
        cmp_dct['S_CAT'] = S_CAT
        cmp_dct['S_AN'] = S_AN

        #************ organic components **************
        cmp_dct['S_CH3OH'] = COD * frSCH3OH
        cmp_dct['S_Ac'] = COD * frSAc
        cmp_dct['S_Prop'] = COD * frSProp
        cmp_dct['S_F'] = COD * frSF
        cmp_dct['S_U_Inf'] = COD * frSUInf
        cmp_dct['S_U_E'] = COD * frSUE

        S_Org = sum([v for k,v in cmp_dct.items() if k in ('S_CH3OH','S_Ac','S_Prop','S_F','S_U_Inf','S_U_E')])

        XC_U_Inf = COD * frXCUInf
        cmp_dct['C_U_Inf'] = XC_U_Inf * r['iCUInf_XCUInf']
        cmp_dct['X_U_Inf'] = XC_U_Inf * (1 - r['iCUInf_XCUInf'])

        cmp_dct['X_OHO'] = COD * frXOHO
        cmp_dct['X_AOO'] = COD * frXAOO
        cmp_dct['X_NOO'] = COD * frXNOO
        cmp_dct['X_AMO'] = COD * frXAMO
        cmp_dct['X_PAO'] = COD * frXPAO
        cmp_dct['X_ACO'] = COD * frXACO
        cmp_dct['X_HMO'] = COD * frXHMO
        cmp_dct['X_PRO'] = COD * frXPRO
        cmp_dct['X_MEOLO'] = COD * frXMEOLO
        cmp_dct['X_FO'] = COD * frXFO

        X_Bio = sum([v for k,v in cmp_dct.items() if k.startswith('X_') and k.endswith('O')])

        cmp_dct['X_OHO_PHA'] = COD * frXOHO_PHA
        cmp_dct['X_GAO_PHA'] = COD * frXGAO_PHA
        cmp_dct['X_PAO_PHA'] = COD * frXPAO_PHA
        cmp_dct['X_GAO_Gly'] = COD * frXGAO_Gly
        cmp_dct['X_PAO_Gly'] = COD * frXPAO_Gly

        X_Stor = sum([v for k,v in cmp_dct.items() if k.endswith(('PHA','Gly'))])

        cmp_dct['X_U_OHO_E'] = COD * frXU_OHO_E
        cmp_dct['X_U_PAO_E'] = COD * frXU_PAO_E

        X_U_E = cmp_dct['X_U_OHO_E'] + cmp_dct['X_U_PAO_E']

        XC_B = COD - S_Org - XC_U_Inf - X_U_E
        C_B = XC_B * r['iCB_XCB']
        cmp_dct['C_B_BAP'] = C_B * r['iBAP_CB']
        cmp_dct['C_B_UAP'] = C_B * r['iUAP_CB']
        cmp_dct['C_B_Subst'] = C_B - cmp_dct['C_B_BAP'] - cmp_dct['C_B_UAP']

        cmp_dct['X_B_Subst'] = XC_B - C_B - X_Bio - X_Stor

        cmp_c = np.asarray([v for v in cmp_dct.values()])
        VSS = (cmp_c * cmps.i_mass * cmps.f_Vmass_Totmass * cmps.x * cmps.org).sum()
        TSS = VSS/iVSS_TSS
        X_Org_ISS = (cmp_c * cmps.i_mass * (1-cmps.f_Vmass_Totmass) * cmps.x * cmps.org).sum()

        del S_Org, XC_U_Inf, X_Bio, X_Stor, X_U_E, XC_B, C_B

        #************ inorganic components **************
        cmp_dct['X_PAO_PP_Hi'] = X_PAO_PP * r['iHi_XPAOPP']
        cmp_dct['X_PAO_PP_Lo'] = X_PAO_PP * (1 - r['iHi_XPAOPP'])

        ISS = TSS - VSS
        cmp_c = np.asarray([v for v in cmp_dct.values()])
        other_ig_iss = (cmp_c * cmps.i_mass * cmps.x * (1-cmps.org)).sum()
        cmp_dct['X_Ig_ISS'] = ISS - X_Org_ISS - other_ig_iss

        del ISS, VSS, TSS, X_Org_ISS, other_ig_iss, cmp_c

        # TODO: calibrate pH, SAlk, SCAT, SAN
        bad_vars = {k:v for k,v in cmp_dct.items() if v<0}
        if len(bad_vars) > 0:
            raise ValueError(f"The following state variable(s) was found negative: {bad_vars}.")

        del bad_vars

        #************ calibrate XB_subst, SF's N, P content *************
        if S_NH4 > 0 and cmp_dct['S_F'] > 0:
            SFi_N = _calib_SF_iN(cmps, cmp_dct, S_NH4/iSNH_STKN)
        XB_Substi_N = _calib_XBsub_iN(cmps, cmp_dct, TKN - S_NH4/iSNH_STKN)
        XB_Substi_P = _calib_XBsub_iP(cmps, cmp_dct, TP)

        if S_NH4 > 0 and cmp_dct['S_F'] > 0: cmps.S_F.i_N = SFi_N
        cmps.X_B_Subst.i_N = XB_Substi_N
        cmps.X_B_Subst.i_P = XB_Substi_P
        cmps.refresh_constants()

        #************ convert concentrations to flow rates *************
        cmp_dct.pop('H2O', None) # avoid warning related to H2O as the bulk liquid
        new.set_flow_by_concentration(flow_tot, cmp_dct, units)
        new.ratios = r

        return new


    @classmethod
    def codbased_inf_model(cls, ID='', flow_tot=0., units = ('L/hr', 'mg/L'),
                           phase='l', T=298.15, P=101325., price=0., thermo=None,
                           pH=7., SAlk=10., ratios=None,
                           COD=430., TKN=40., TP=10., iVSS_TSS=0.75, iSNH_STKN=0.9,
                           iSCOD_COD=0.25, iSBOD_SCOD=0.50, iBOD_COD=0.58,
                           S_NH4=25., S_NO2=0., S_NO3=0., S_PO4=8.,
                           S_Ca=140., S_Mg=50., S_K=28., S_CAT=3., S_AN=12., S_N2=18.,
                           C_B=40., S_CH3OH=0., S_Ac=0., S_Prop=0.,
                           X_OHO_PHA=0., X_GAO_PHA=0., X_PAO_PHA=0.,
                           X_GAO_Gly=0., X_PAO_Gly=0., X_U_OHO_E=0., X_U_PAO_E=0.,
                           X_OHO=0., X_AOO=0., X_NOO=0., X_AMO=0., X_PAO=0.,
                           X_PRO=0., X_ACO=0., X_HMO=0., X_MEOLO=0., X_FO=0.,
                           X_FePO4=0., X_AlPO4=0., X_FeOH=0., X_AlOH=0.,
                           X_MAP=0., X_HAP=0., X_HDP=0., X_PAO_PP=0.,
                           X_MgCO3=0., X_CaCO3=0., DO=0., S_H2=0., S_CH4=0.):
        '''
        Create a waste stream using the COD-based influent characterization model as in [1]_,
        requires the default component set (can be created by :func:`~.Components.load_default`)

        Normal attributes of the waste streams (e.g., phase, T, Pa, price) can be passed as usual.

        Parameters
        ----------
        ID : str
            ID of the stream.
        flow_tot : float
            Total flowrates.
        units : tuple(str, str)
            Unit of the flowrate and component concentrations.

        Examples
        --------
        >>> from qsdsan import Components, WasteStream, set_thermo
        >>> cmps = Components.load_default()
        >>> set_thermo(cmps)
        >>> ws = WasteStream.codbased_inf_model('ws_codinf', flow_tot=100)
        >>> ws.show()
        WasteStream: ws_codinf
         phase: 'l', T: 298.15 K, P: 101325 Pa
         flow (g/hr): S_F        7.5
                      S_U_Inf    3.25
                      C_B_Subst  4
                      X_B_Subst  22.7
                      X_U_Inf    5.58
                      X_Ig_ISS   5.23
                      S_NH4      2.5
                      S_PO4      0.8
                      S_K        2.8
                      S_Ca       14
                      S_Mg       5
                      S_CO3      12
                      S_N2       1.8
                      S_CAT      0.3
                      S_AN       1.2
                      ...
         WasteStream-specific properties:
          pH         : 7.0
          Alkalinity : 10.0 mg/L
          COD        : 430.0 mg/L
          BOD        : 249.4 mg/L
          TC         : 265.0 mg/L
          TOC        : 137.6 mg/L
          TN         : 40.0 mg/L
          TP         : 10.0 mg/L
          TK         : 28.0 mg/L
         Component concentrations (mg/L):
          S_F          75.0
          S_U_Inf      32.5
          C_B_Subst    40.0
          X_B_Subst    226.7
          X_U_Inf      55.8
          X_Ig_ISS     52.3
          S_NH4        25.0
          S_PO4        8.0
          S_K          28.0
          S_Ca         140.0
          S_Mg         50.0
          S_CO3        120.0
          S_N2         18.0
          S_CAT        3.0
          S_AN         12.0
          ...

        References
        ----------
        [1] Hydromantis. GPS-X Technical Reference - v8.0.
        Hydromantis Environmental Software Solutions, Inc. 2019.
        https://www.hydromantis.com/help/GPS-X/docs/8.0/Technical/index.html

        '''

        if thermo: cmps = thermo.chemicals
        else:
            cmps = _load_components()
            thermo = _load_thermo()

        cmp_dct = dict.fromkeys(cmps.IDs, 0.)

        new = cls(ID=ID, phase=phase, T=T, P=P, units='kg/hr', price=price,
                  thermo=thermo, pH=pH, SAlk=SAlk)

        if ratios: new.ratios = ratios
        else: new.ratios = WasteStream._default_ratios
        r = new._ratios

        #************ user-defined inorganic states **************
        cmp_dct['S_H2'] = S_H2
        cmp_dct['S_CH4'] = S_CH4
        cmp_dct['S_N2'] = S_N2
        cmp_dct['S_O2'] = DO
        cmp_dct['S_NH4'] = S_NH4
        cmp_dct['S_NO2'] = S_NO2
        cmp_dct['S_NO3'] = S_NO3
        cmp_dct['S_PO4'] = S_PO4
        cmp_dct['S_CO3'] = SAlk * 12 * conc_unit.conversion_factor(units[1])       # 1 meq/L SAlk ~ 1 mmol/L HCO3- ~ 12 mg C/L (12 mg C/mmol HCO3-)
        cmp_dct['S_Ca'] = S_Ca
        cmp_dct['S_Mg'] = S_Mg
        cmp_dct['S_K'] = S_K
        cmp_dct['X_MAP'] = X_MAP
        cmp_dct['X_HAP'] = X_HAP
        cmp_dct['X_HDP'] = X_HDP
        cmp_dct['X_FePO4'] = X_FePO4
        cmp_dct['X_AlPO4'] = X_AlPO4
        cmp_dct['X_FeOH'] = X_FeOH
        cmp_dct['X_AlOH'] = X_AlOH
        cmp_dct['X_MgCO3'] = X_MgCO3
        cmp_dct['X_CaCO3'] = X_CaCO3
        cmp_dct['S_CAT'] = S_CAT
        cmp_dct['S_AN'] = S_AN

        #************ organic components **************
        sCOD = COD * iSCOD_COD
        sBOD = sCOD * iSBOD_SCOD
        cmp_dct['S_CH3OH'] = S_CH3OH
        cmp_dct['S_Ac'] = S_Ac
        cmp_dct['S_Prop'] = S_Prop

        cmp_c = np.asarray([v for v in cmp_dct.values()])
        other_sBOD = (cmp_c * cmps.s * cmps.org * cmps.f_BOD5_COD).sum()
        cmp_dct['S_F'] = (sBOD - other_sBOD)/cmps.S_F.f_BOD5_COD

        cmp_c = np.asarray([v for v in cmp_dct.values()])
        S_U = sCOD - (cmp_c * cmps.s * cmps.org).sum()
        cmp_dct['S_U_Inf'] = S_U * r['iSUInf_SU']
        cmp_dct['S_U_E'] = S_U - cmp_dct['S_U_Inf']

        cmp_dct['C_B_BAP'] = C_B * r['iBAP_CB']
        cmp_dct['C_B_UAP'] = C_B * r['iUAP_CB']
        cmp_dct['C_B_Subst'] = C_B - cmp_dct['C_B_BAP'] - cmp_dct['C_B_UAP']

        XC_B = C_B/r['iCB_XCB']
        XC_U = COD - sCOD - XC_B
        cmp_dct['X_U_OHO_E'] = X_U_OHO_E
        cmp_dct['X_U_PAO_E'] = X_U_PAO_E

        XC_U_Inf = XC_U - X_U_OHO_E - X_U_PAO_E
        cmp_dct['C_U_Inf'] = XC_U_Inf * r['iCUInf_XCUInf']
        cmp_dct['X_U_Inf'] = XC_U_Inf * (1 - r['iCUInf_XCUInf'])

        cmp_dct['X_OHO'] = X_OHO
        cmp_dct['X_AOO'] = X_AOO
        cmp_dct['X_NOO'] = X_NOO
        cmp_dct['X_AMO'] = X_AMO
        cmp_dct['X_PAO'] = X_PAO
        cmp_dct['X_ACO'] = X_ACO
        cmp_dct['X_HMO'] = X_HMO
        cmp_dct['X_PRO'] = X_PRO
        cmp_dct['X_MEOLO'] = X_MEOLO
        cmp_dct['X_FO'] = X_FO
        cmp_dct['X_OHO_PHA'] = X_OHO_PHA
        cmp_dct['X_GAO_PHA'] = X_GAO_PHA
        cmp_dct['X_PAO_PHA'] = X_PAO_PHA
        cmp_dct['X_GAO_Gly'] = X_GAO_Gly
        cmp_dct['X_PAO_Gly'] = X_PAO_Gly

        cmp_c = np.asarray([v for v in cmp_dct.values()])
        cmp_dct['X_B_Subst'] = COD - (cmp_c * cmps.org * (cmps.s + cmps.c + cmps.x)).sum()

        cmp_c = np.asarray([v for v in cmp_dct.values()])
        VSS = (cmp_c * cmps.i_mass * cmps.f_Vmass_Totmass * cmps.x * cmps.org).sum()
        TSS = VSS/iVSS_TSS
        X_Org_ISS = (cmp_c * cmps.i_mass * (1-cmps.f_Vmass_Totmass) * cmps.x * cmps.org).sum()

        del sCOD, sBOD, other_sBOD, S_U, XC_U, XC_B, XC_U_Inf

        #************ inorganic components **************
        cmp_dct['X_PAO_PP_Hi'] = X_PAO_PP * r['iHi_XPAOPP']
        cmp_dct['X_PAO_PP_Lo'] = X_PAO_PP * (1 - r['iHi_XPAOPP'])

        ISS = TSS - VSS
        cmp_c = np.asarray([v for v in cmp_dct.values()])
        other_ig_iss = (cmp_c * cmps.i_mass * cmps.x * (1-cmps.org)).sum()
        cmp_dct['X_Ig_ISS'] = ISS - X_Org_ISS - other_ig_iss

        del ISS, VSS, TSS, X_Org_ISS, other_ig_iss, cmp_c

        # TODO: calibrate pH, SAlk, SCAT, SAN
        bad_vars = {k:v for k,v in cmp_dct.items() if v<0}
        if len(bad_vars) > 0:
            raise ValueError(f"The following state variable(s) was found negative: {bad_vars}.")

        del bad_vars

        #************ calibrate XB_subst, SF's N, P content *************
        if S_NH4 > 0 and cmp_dct['S_F'] > 0:
            SFi_N = _calib_SF_iN(cmps, cmp_dct, S_NH4/iSNH_STKN)
        XB_Substi_N = _calib_XBsub_iN(cmps, cmp_dct, TKN - S_NH4/iSNH_STKN)
        XB_Substi_P = _calib_XBsub_iP(cmps, cmp_dct, TP)

        BOD = COD * iBOD_COD
        sub_IDs = ('X_B_Subst', 'X_OHO_PHA', 'X_GAO_PHA', 'X_PAO_PHA', 'X_GAO_Gly', 'X_PAO_Gly')
        fbodtocod_sub = _calib_XBsub_fBODCOD(cmps, cmp_dct, sub_IDs, BOD)

        if S_NH4 > 0 and cmp_dct['S_F'] > 0: cmps.S_F.i_N = SFi_N
        cmps.X_B_Subst.i_N = XB_Substi_N
        cmps.X_B_Subst.i_P = XB_Substi_P
        for i in sub_IDs: cmps[i].f_BOD5_COD = fbodtocod_sub
        cmps.refresh_constants()

        #************ convert concentrations to flow rates *************
        cmp_dct.pop('H2O', None) # avoid warning related to H2O as the bulk liquid
        new.set_flow_by_concentration(flow_tot, cmp_dct, units)
        new.ratios = r

        return new


    @classmethod
    def bodbased_inf_model(cls, ID='', flow_tot=0., units = ('L/hr', 'mg/L'),
                           phase='l', T=298.15, P=101325., price=0., thermo=None,
                           pH=7., SAlk=10., ratios=None,
                           BOD=250., TKN=40., TP=10., iVSS_TSS=0.75, iSNH_STKN=0.9,
                           iSBOD_BOD=0.25, iSBOD_SCOD=0.50, iBOD_COD=0.58,
                           S_N2=18., S_NH4=25., S_NO2=0., S_NO3=0., S_PO4=8.,
                           S_Ca=140., S_Mg=50., S_K=28., S_CAT=3., S_AN=12.,
                           C_B=40., S_CH3OH=0., S_Ac=0., S_Prop=0.,
                           X_OHO_PHA=0., X_GAO_PHA=0., X_PAO_PHA=0.,
                           X_GAO_Gly=0., X_PAO_Gly=0., X_U_OHO_E=0., X_U_PAO_E=0.,
                           X_OHO=0., X_AOO=0., X_NOO=0., X_AMO=0., X_PAO=0.,
                           X_PRO=0., X_ACO=0., X_HMO=0., X_MEOLO=0., X_FO=0.,
                           X_FePO4=0., X_AlPO4=0., X_FeOH=0., X_AlOH=0.,
                           X_MAP=0., X_HAP=0., X_HDP=0., X_PAO_PP=0.,
                           X_MgCO3=0., X_CaCO3=0., DO=0., S_H2=0., S_CH4=0.):
        '''
        Create a waste stream using the BOD influent characterization model as in [1]_,
        requires the default component set (can be created by :func:`~.Components.load_default`)

        Normal attributes of the waste streams (e.g., phase, T, Pa, price) can be passed as usual.

        Parameters
        ----------
        ID : str
            ID of the stream.
        flow_tot : float
            Total flowrates.
        units : tuple(str, str)
            Unit of the flowrate and component concentrations.

        Examples
        --------
        >>> from qsdsan import Components, WasteStream, set_thermo
        >>> cmps = Components.load_default()
        >>> set_thermo(cmps)
        >>> ws = WasteStream.bodbased_inf_model('ws_bodinf', flow_tot=100)
        >>> ws.show()
        WasteStream: ws_bodinf
         phase: 'l', T: 298.15 K, P: 101325 Pa
         flow (g/hr): S_F        8.72
                      S_U_Inf    3.78
                      C_B_Subst  4
                      X_B_Subst  22.7
                      X_U_Inf    3.94
                      X_Ig_ISS   4.93
                      S_NH4      2.5
                      S_PO4      0.8
                      S_K        2.8
                      S_Ca       14
                      S_Mg       5
                      S_CO3      12
                      S_N2       1.8
                      S_CAT      0.3
                      S_AN       1.2
                      ...
         WasteStream-specific properties:
          pH         : 7.0
          Alkalinity : 10.0 mg/L
          COD        : 431.0 mg/L
          BOD        : 250.0 mg/L
          TC         : 264.9 mg/L
          TOC        : 137.9 mg/L
          TN         : 40.0 mg/L
          TP         : 10.0 mg/L
          TK         : 28.0 mg/L
         Component concentrations (mg/L):
          S_F          87.2
          S_U_Inf      37.8
          C_B_Subst    40.0
          X_B_Subst    226.7
          X_U_Inf      39.4
          X_Ig_ISS     49.3
          S_NH4        25.0
          S_PO4        8.0
          S_K          28.0
          S_Ca         140.0
          S_Mg         50.0
          S_CO3        120.0
          S_N2         18.0
          S_CAT        3.0
          S_AN         12.0
          ...

        References
        ----------
        [1] Hydromantis. GPS-X Technical Reference - v8.0.
        Hydromantis Environmental Software Solutions, Inc. 2019.
        https://www.hydromantis.com/help/GPS-X/docs/8.0/Technical/index.html

        '''

        if thermo: cmps = thermo.chemicals
        else:
            cmps = _load_components()
            thermo = _load_thermo()

        cmp_dct = dict.fromkeys(cmps.IDs, 0.)

        new = cls(ID=ID, phase=phase, T=T, P=P, units='kg/hr', price=price,
                  thermo=thermo, pH=pH, SAlk=SAlk)

        if ratios: new.ratios = ratios
        else: new.ratios = WasteStream._default_ratios
        r = new._ratios

        #************ user-defined inorganic states **************
        cmp_dct['S_H2'] = S_H2
        cmp_dct['S_CH4'] = S_CH4
        cmp_dct['S_N2'] = S_N2
        cmp_dct['S_O2'] = DO
        cmp_dct['S_NH4'] = S_NH4
        cmp_dct['S_NO2'] = S_NO2
        cmp_dct['S_NO3'] = S_NO3
        cmp_dct['S_PO4'] = S_PO4
        cmp_dct['S_CO3'] = SAlk * 12 * conc_unit.conversion_factor(units[1])       # 1 meq/L SAlk ~ 1 mmol/L HCO3- ~ 12 mg C/L (12 mg C/mmol HCO3-)
        cmp_dct['S_Ca'] = S_Ca
        cmp_dct['S_Mg'] = S_Mg
        cmp_dct['S_K'] = S_K
        cmp_dct['X_MAP'] = X_MAP
        cmp_dct['X_HAP'] = X_HAP
        cmp_dct['X_HDP'] = X_HDP
        cmp_dct['X_FePO4'] = X_FePO4
        cmp_dct['X_AlPO4'] = X_AlPO4
        cmp_dct['X_FeOH'] = X_FeOH
        cmp_dct['X_AlOH'] = X_AlOH
        cmp_dct['X_MgCO3'] = X_MgCO3
        cmp_dct['X_CaCO3'] = X_CaCO3
        cmp_dct['S_CAT'] = S_CAT
        cmp_dct['S_AN'] = S_AN

        #************ organic components **************
        COD = BOD / iBOD_COD
        sBOD = BOD * iSBOD_BOD
        cmp_dct['S_CH3OH'] = S_CH3OH
        cmp_dct['S_Ac'] = S_Ac
        cmp_dct['S_Prop'] = S_Prop

        cmp_c = np.asarray([v for v in cmp_dct.values()])
        other_sBOD = (cmp_c * cmps.s * cmps.org * cmps.f_BOD5_COD).sum()
        cmp_dct['S_F'] = (sBOD - other_sBOD)/cmps.S_F.f_BOD5_COD

        sCOD = sBOD / iSBOD_SCOD
        cmp_c = np.asarray([v for v in cmp_dct.values()])
        S_U = sCOD - (cmp_c * cmps.s * cmps.org).sum()
        cmp_dct['S_U_Inf'] = S_U * r['iSUInf_SU']
        cmp_dct['S_U_E'] = S_U - cmp_dct['S_U_Inf']

        cmp_dct['C_B_BAP'] = C_B * r['iBAP_CB']
        cmp_dct['C_B_UAP'] = C_B * r['iUAP_CB']
        cmp_dct['C_B_Subst'] = C_B - cmp_dct['C_B_BAP'] - cmp_dct['C_B_UAP']

        XC_B = C_B/r['iCB_XCB']
        XC_U = COD - sCOD - XC_B
        cmp_dct['X_U_OHO_E'] = X_U_OHO_E
        cmp_dct['X_U_PAO_E'] = X_U_PAO_E

        XC_U_Inf = XC_U - X_U_OHO_E - X_U_PAO_E
        cmp_dct['C_U_Inf'] = XC_U_Inf * r['iCUInf_XCUInf']
        cmp_dct['X_U_Inf'] = XC_U_Inf * (1 - r['iCUInf_XCUInf'])

        cmp_dct['X_OHO'] = X_OHO
        cmp_dct['X_AOO'] = X_AOO
        cmp_dct['X_NOO'] = X_NOO
        cmp_dct['X_AMO'] = X_AMO
        cmp_dct['X_PAO'] = X_PAO
        cmp_dct['X_ACO'] = X_ACO
        cmp_dct['X_HMO'] = X_HMO
        cmp_dct['X_PRO'] = X_PRO
        cmp_dct['X_MEOLO'] = X_MEOLO
        cmp_dct['X_FO'] = X_FO
        cmp_dct['X_OHO_PHA'] = X_OHO_PHA
        cmp_dct['X_GAO_PHA'] = X_GAO_PHA
        cmp_dct['X_PAO_PHA'] = X_PAO_PHA
        cmp_dct['X_GAO_Gly'] = X_GAO_Gly
        cmp_dct['X_PAO_Gly'] = X_PAO_Gly

        cmp_c = np.asarray([v for v in cmp_dct.values()])
        cmp_dct['X_B_Subst'] = COD - (cmp_c * cmps.org * (cmps.s + cmps.c + cmps.x)).sum()

        cmp_c = np.asarray([v for v in cmp_dct.values()])
        VSS = (cmp_c * cmps.i_mass * cmps.f_Vmass_Totmass * cmps.x * cmps.org).sum()
        TSS = VSS/iVSS_TSS
        X_Org_ISS = (cmp_c * cmps.i_mass * (1-cmps.f_Vmass_Totmass) * cmps.x * cmps.org).sum()

        del sCOD, sBOD, other_sBOD, S_U, XC_U, XC_B, XC_U_Inf

        #************ inorganic components **************
        cmp_dct['X_PAO_PP_Hi'] = X_PAO_PP * r['iHi_XPAOPP']
        cmp_dct['X_PAO_PP_Lo'] = X_PAO_PP * (1 - r['iHi_XPAOPP'])

        ISS = TSS - VSS
        cmp_c = np.asarray([v for v in cmp_dct.values()])
        other_ig_iss = (cmp_c * cmps.i_mass * cmps.x * (1-cmps.org)).sum()
        cmp_dct['X_Ig_ISS'] = ISS - X_Org_ISS - other_ig_iss

        del ISS, VSS, TSS, X_Org_ISS, other_ig_iss, cmp_c

        # TODO: calibrate pH, SAlk, SCAT, SAN
        bad_vars = {k:v for k,v in cmp_dct.items() if v<0}
        if len(bad_vars) > 0:
            raise ValueError(f"The following state variable(s) was found negative: {bad_vars}.")

        del bad_vars

        #************ calibrate XB_subst, SF's N, P content *************
        if S_NH4 > 0 and cmp_dct['S_F'] > 0:
            SFi_N = _calib_SF_iN(cmps, cmp_dct, S_NH4/iSNH_STKN)
        XB_Substi_N = _calib_XBsub_iN(cmps, cmp_dct, TKN - S_NH4/iSNH_STKN)
        XB_Substi_P = _calib_XBsub_iP(cmps, cmp_dct, TP)

        BOD = COD * iBOD_COD
        sub_IDs = ('X_B_Subst', 'X_OHO_PHA', 'X_GAO_PHA', 'X_PAO_PHA', 'X_GAO_Gly', 'X_PAO_Gly')
        fbodtocod_sub = _calib_XBsub_fBODCOD(cmps, cmp_dct, sub_IDs, BOD)

        if S_NH4 > 0 and cmp_dct['S_F'] > 0: cmps.S_F.i_N = SFi_N
        cmps.X_B_Subst.i_N = XB_Substi_N
        cmps.X_B_Subst.i_P = XB_Substi_P
        for i in sub_IDs: cmps[i].f_BOD5_COD = fbodtocod_sub
        cmps.refresh_constants()

        #************ convert concentrations to flow rates *************
        cmp_dct.pop('H2O', None) # avoid warning related to H2O as the bulk liquid
        new.set_flow_by_concentration(flow_tot, cmp_dct, units)
        new.ratios = r

        return new


    @classmethod
    def sludge_inf_model(cls, ID='', flow_tot=0., units = ('L/hr', 'mg/L'),
                         phase='l', T=298.15, P=101325., price=0., thermo=None,
                         pH=7., SAlk=10., ratios=None,
                         TSS=1e4, TKN=750., TP=250., S_NH4=100., S_PO4=50.,
                         iVSS_TSS=0.65, iscCOD_COD=0.01, iSNH_STKN=0.9,
                         frXUInf_VSS=0.4, frXUE_VSS=0.3, frXOHO_VSS=0.2,
                         frXAOO_VSS=0.01, frXNOO_VSS=0.01,frXPAO_VSS=0.01,
                         frCB_scCOD=0.1, frSU_scCOD=0.8,
                         S_Ca=140., S_Mg=50., S_K=28., S_CAT=3., S_AN=12., S_N2=18.,
                         frXACO_VSS=0., frXHMO_VSS=0., frXPRO_VSS=0., frXFO_VSS=0.,
                         frXMEOLO_VSS=0., frXAMO_VSS=0.,
                         frXOHO_PHA_VSS=0., frXGAO_PHA_VSS=0., frXPAO_PHA_VSS=0.,
                         frXGAO_Gly_VSS=0., frXPAO_Gly_VSS=0.,
                         frSCH3OH_scCOD=0., frSAc_scCOD=0., frSProp_scCOD=0.,
                         S_NO2=0., S_NO3=0., X_PAO_PP=0., X_FeOH=0., X_AlOH=0.,
                         X_FePO4=0., X_AlPO4=0., X_MAP=0., X_HAP=0., X_HDP=0.,
                         X_MgCO3=0., X_CaCO3=0., DO=0., S_H2=0., S_CH4=0.):
        '''
        Create a waste stream using the sludge influent characterization model as in [1]_,
        requires the default component set (can be created by :func:`~.Components.load_default`)

        Normal attributes of the waste streams (e.g., phase, T, Pa, price) can be passed as usual.

        Parameters
        ----------
        ID : str
            ID of the stream.
        flow_tot : float
            Total flowrates.
        units : tuple(str, str)
            Unit of the flowrate and component concentrations.

        Examples
        --------
        >>> from qsdsan import Components, WasteStream, set_thermo
        >>> cmps = Components.load_default()
        >>> set_thermo(cmps)
        >>> ws = WasteStream.sludge_inf_model('ws_sludgeinf', flow_tot=100)
        >>> ws.show()
        WasteStream: ws_sludgeinf
         phase: 'l', T: 298.15 K, P: 101325 Pa
         flow (g/hr): S_F        1.08
                      S_U_Inf    8.65
                      C_B_Subst  1.08
                      X_B_Subst  81.9
                      X_OHO      192
                      X_AOO      9.62
                      X_NOO      9.62
                      X_PAO      9.62
                      X_U_Inf    468
                      X_U_OHO_E  285
                      X_U_PAO_E  14.3
                      X_Ig_ISS   298
                      S_NH4      10
                      S_PO4      5
                      S_K        2.8
                      ...
         WasteStream-specific properties:
          pH         : 7.0
          Alkalinity : 10.0 mg/L
          COD        : 10814.4 mg/L
          BOD        : 1744.3 mg/L
          TC         : 4246.5 mg/L
          TOC        : 3702.2 mg/L
          TN         : 750.0 mg/L
          TP         : 250.0 mg/L
          TK         : 52.9 mg/L
         Component concentrations (mg/L):
          S_F          10.8
          S_U_Inf      86.5
          C_B_Subst    10.8
          X_B_Subst    819.0
          X_OHO        1924.0
          X_AOO        96.2
          X_NOO        96.2
          X_PAO        96.2
          X_U_Inf      4680.0
          X_U_OHO_E    2852.0
          X_U_PAO_E    142.6
          X_Ig_ISS     2980.6
          S_NH4        100.0
          S_PO4        50.0
          S_K          28.0
          ...

        References
        ----------
        [1] Hydromantis. GPS-X Technical Reference - v8.0.
        Hydromantis Environmental Software Solutions, Inc. 2019.
        https://www.hydromantis.com/help/GPS-X/docs/8.0/Technical/index.html

        '''

        if thermo: cmps = thermo.chemicals
        else:
            cmps = _load_components()
            thermo = _load_thermo()

        cmp_dct = dict.fromkeys(cmps.IDs, 0.)

        new = cls(ID=ID, phase=phase, T=T, P=P, units='kg/hr', price=price,
                  thermo=thermo, pH=pH, SAlk=SAlk)

        if ratios: new.ratios = ratios
        else: new.ratios = WasteStream._default_ratios
        r = new._ratios

        #************ user-defined inorganic states **************
        cmp_dct['S_H2'] = S_H2
        cmp_dct['S_CH4'] = S_CH4
        cmp_dct['S_N2'] = S_N2
        cmp_dct['S_O2'] = DO
        cmp_dct['S_NH4'] = S_NH4
        cmp_dct['S_NO2'] = S_NO2
        cmp_dct['S_NO3'] = S_NO3
        cmp_dct['S_PO4'] = S_PO4
        cmp_dct['S_CO3'] = SAlk * 12 * conc_unit.conversion_factor(units[1])       # 1 meq/L SAlk ~ 1 mmol/L HCO3- ~ 12 mg C/L (12 mg C/mmol HCO3-)
        cmp_dct['S_Ca'] = S_Ca
        cmp_dct['S_Mg'] = S_Mg
        cmp_dct['S_K'] = S_K
        cmp_dct['X_MAP'] = X_MAP
        cmp_dct['X_HAP'] = X_HAP
        cmp_dct['X_HDP'] = X_HDP
        cmp_dct['X_FePO4'] = X_FePO4
        cmp_dct['X_AlPO4'] = X_AlPO4
        cmp_dct['X_FeOH'] = X_FeOH
        cmp_dct['X_AlOH'] = X_AlOH
        cmp_dct['X_MgCO3'] = X_MgCO3
        cmp_dct['X_CaCO3'] = X_CaCO3
        cmp_dct['S_CAT'] = S_CAT
        cmp_dct['S_AN'] = S_AN

        #************ particulate components **************
        VSS = TSS * iVSS_TSS
        cmp_dct['X_U_Inf'] = frXUInf_VSS * VSS

        if r['iXUOHOE_XUE']: frOHO = r['iXUOHOE_XUE']
        else:
            try: frOHO = frXOHO_VSS/(frXOHO_VSS + frXPAO_VSS)
            except ZeroDivisionError: frOHO = 0.5
        cmp_dct['X_U_OHO_E'] = frXUE_VSS * VSS * frOHO
        cmp_dct['X_U_PAO_E'] = frXUE_VSS * VSS * (1-frOHO)

        cmp_dct['X_OHO'] = frXOHO_VSS * VSS
        cmp_dct['X_AOO'] = frXAOO_VSS * VSS
        cmp_dct['X_NOO'] = frXNOO_VSS * VSS
        cmp_dct['X_AMO'] = frXAMO_VSS * VSS
        cmp_dct['X_PAO'] = frXPAO_VSS * VSS
        cmp_dct['X_ACO'] = frXACO_VSS * VSS
        cmp_dct['X_HMO'] = frXHMO_VSS * VSS
        cmp_dct['X_PRO'] = frXPRO_VSS * VSS
        cmp_dct['X_MEOLO'] = frXMEOLO_VSS * VSS
        cmp_dct['X_FO'] = frXFO_VSS * VSS
        cmp_dct['X_OHO_PHA'] = frXOHO_PHA_VSS * VSS
        cmp_dct['X_GAO_PHA'] = frXGAO_PHA_VSS * VSS
        cmp_dct['X_PAO_PHA'] = frXPAO_PHA_VSS * VSS
        cmp_dct['X_GAO_Gly'] = frXGAO_Gly_VSS * VSS
        cmp_dct['X_PAO_Gly'] = frXPAO_Gly_VSS * VSS

        cmp_c = np.asarray([v for v in cmp_dct.values()])
        cmp_dct['X_B_Subst'] = VSS - (cmp_c * cmps.org * cmps.x).sum()

        # convert gVSS to gCOD
        for cmp in cmps:
            if cmp.organic and cmp.particle_size == 'Particulate':
                cmp_dct[cmp.ID] /= cmp.i_mass * cmp.f_Vmass_Totmass

        cmp_dct['X_PAO_PP_Hi'] = X_PAO_PP * r['iHi_XPAOPP']
        cmp_dct['X_PAO_PP_Lo'] = X_PAO_PP * (1 - r['iHi_XPAOPP'])

        cmp_c = np.asarray([v for v in cmp_dct.values()])
        ig_ISS = TSS - (cmp_c * cmps.i_mass * cmps.x * cmps.org).sum()
        other_ig_iss = (cmp_c * cmps.i_mass * cmps.x * (1-cmps.org)).sum()
        cmp_dct['X_Ig_ISS'] = ig_ISS - other_ig_iss

        del other_ig_iss, cmp_c

        #*********** soluble and colloidal components *************
        cmp_c = np.asarray([v for v in cmp_dct.values()])
        xCOD = (cmp_c * cmps.x * cmps.org).sum()
        scCOD = xCOD * iscCOD_COD / (1-iscCOD_COD)

        cmp_dct['S_CH3OH'] = frSCH3OH_scCOD * scCOD
        cmp_dct['S_Ac'] = frSAc_scCOD * scCOD
        cmp_dct['S_Prop'] = frSProp_scCOD * scCOD
        cmp_dct['S_U_Inf'] = frSU_scCOD * scCOD * r['iSUInf_SU']
        cmp_dct['S_U_E'] = frSU_scCOD * scCOD - cmp_dct['S_U_Inf']

        C_B = frCB_scCOD * scCOD
        cmp_dct['C_B_BAP'] = C_B * r['iBAP_CB']
        cmp_dct['C_B_UAP'] = C_B * r['iUAP_CB']
        cmp_dct['C_B_Subst'] = C_B - cmp_dct['C_B_BAP'] - cmp_dct['C_B_UAP']
        cmp_dct['C_U_Inf'] = cmp_dct['X_U_Inf'] * r['iCUInf_XCUInf'] / (1-r['iCUInf_XCUInf'])

        cmp_c = np.asarray([v for v in cmp_dct.values()])
        cmp_dct['S_F'] = scCOD - (cmp_c * (cmps.s + cmps.c) * cmps.org).sum()

        # TODO: calibrate pH, SAlk, SCAT, SAN
        bad_vars = {k:v for k,v in cmp_dct.items() if v<0}
        if len(bad_vars) > 0:
            raise ValueError(f"The following state variable(s) was found negative: {bad_vars}.")

        del bad_vars

        #************ calibrate XB_subst, SF's N, P content *************
        if S_NH4 > 0 and cmp_dct['S_F'] > 0:
            SFi_N = _calib_SF_iN(cmps, cmp_dct, S_NH4/iSNH_STKN)
        XB_Substi_N = _calib_XBsub_iN(cmps, cmp_dct, TKN - S_NH4/iSNH_STKN)
        XB_Substi_P = _calib_XBsub_iP(cmps, cmp_dct, TP)

        if S_NH4 > 0 and cmp_dct['S_F'] > 0: cmps.S_F.i_N = SFi_N
        cmps.X_B_Subst.i_N = XB_Substi_N
        cmps.X_B_Subst.i_P = XB_Substi_P
        cmps.refresh_constants()

        #************ convert concentrations to flow rates *************
        cmp_dct.pop('H2O', None) # avoid warning related to H2O as the bulk liquid
        new.set_flow_by_concentration(flow_tot, cmp_dct, units)
        new.ratios = r

        return new



# %%

class MissingWasteStream(MissingSanStream):
    '''
    A subclass of :class:`MissingSanStream`, create a special object
    that acts as a dummy until replaced by an actual :class:`WasteStream`.

    .. note::

        Users usually do not need to interact with this class.
    '''
    line = 'WasteStream'

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
        material_stream = WasteStream(ID, thermo=(source or sink).thermo)
        if source: source._outs.replace(self, material_stream)
        if sink: sink._ins.replace(self, material_stream)
        return material_stream

    @property
    def pH(self):
        '''[float] pH, unitless.'''
        return 7.

    @property
    def SAlk(self):
        '''[float] Alkalinity.'''
        return 0.

    @property
    def COD(self):
        '''[float] Chemical oxygen demand.'''
        return 0.

    @property
    def BOD(self):
        '''[float] Biochemical oxygen demand, same as BOD5.'''
        return 0.

    @property
    def BOD5(self):
        '''[float] 5-day biochemical oxygen demand, same as BOD.'''
        return self.BOD

    @property
    def uBOD(self):
        '''[float] Ultimate biochemical oxygen demand.'''
        return 0.

    @property
    def cnBOD(self):
        '''[float] Carbonaceous nitrogenous BOD.'''
        return 0.

    @property
    def ThOD(self):
        '''[float] Theoretical oxygen demand.'''
        return 0.

    @property
    def TC(self):
        '''[float] Total carbon.'''
        return 0.

    @property
    def TOC(self):
        '''[float] Total organic carbon.'''
        return 0.

    @property
    def TN(self):
        '''[float] Total nitrogen.'''
        return 0.

    @property
    def TKN(self):
        '''[float] Total Kjeldahl nitrogen.'''
        return 0.

    @property
    def TP(self):
        '''[float] Total phosphorus.'''
        return 0.

    @property
    def TK(self):
        '''[float] Total potassium.'''
        return 0.

    @property
    def TMg(self):
        '''[float] Total magnesium.'''
        return 0.

    @property
    def TCa(self):
        '''[float] Total calcium.'''
        return 0.

    @property
    def dry_mass(self):
        '''[float] Total solids.'''
        return 0.

    @property
    def density(self):
        '''[float] Density of the stream, in mg/L (kg/m3).'''
        return 0.

    #!!! Keep this up-to-date with WasteStream
    # @property
    # def charge(self):
    #     return 0.

    def __repr__(self):
        return '<MissingWasteStream>'

    def __str__(self):
        return 'missing waste stream'