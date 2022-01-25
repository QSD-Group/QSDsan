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

import numpy as np, flexsolve as flx
from warnings import warn
from math import ceil
from biosteam import Splitter, SolidsCentrifuge
from . import Decay
from .. import SanUnit
from ..sanunits import Pump
from ..utils import ospath, load_data, data_path, dct_from_str

__all__ = (
    'SludgeThickening',
    'BeltThickener',
    'SludgeCentrifuge',
    'SludgeSeparator',
    )


# %%

class SludgeThickening(SanUnit, Splitter):
    '''
    A generic class for concentrating (i.e., thickening) of sludge
    from wastewater treatment processes based on
    `Shoener et al. <https://doi.org/10.1039/C5EE03715H>`_

    The 0th outs is the water-rich supernatant (effluent) and
    the 1st outs is the solid-rich sludge.

    Two pumps (one for the supernatant and one for sludge) are included.

    Separation split is determined by the moisture (i.e., water)
    content of the sludge, soluble components will have the same split as water,
    insolubles components will all go to the retentate.
    
    Note that if the moisture content of the incoming feeds are smaller than
    the target moisture content, the target moisture content will be ignored.

    Parameters
    ----------
    ins : Iterable(obj)
        Dilute sludge stream.
    outs : Iterable(obj)
        Water/bulk-liquid-rich stream, concentrated sludge.
    sludge_moisture : float
        Moisture content of the sludge, [wt% water].
    solids : Iterable(str)
        IDs of the solid components.
        If not provided, will be set to the default `solids` attribute of the components.
    disposal_cost : float
        Disposal cost of the dewatered solids. [$/kg].

    References
    ----------
    [1] Shoener et al., Design of Anaerobic Membrane Bioreactors for the
    Valorization of Dilute Organic Carbon Waste Streams.
    Energy Environ. Sci. 2016, 9 (3), 1102–1112.
    https://doi.org/10.1039/C5EE03715H.
    '''

    _graphics = Splitter._graphics
    _ins_size_is_fixed = False
    _N_outs = 2
    auxiliary_unit_names = ('effluent_pump', 'sludge_pump')

    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream',
                 sludge_moisture=0.96, solids=(), 
                 disposal_cost=125/907.18474): # from $/U.S. ton
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with=init_with)
        self.sludge_moisture = sludge_moisture
        cmps = self.components
        self.solids = solids or tuple((cmp.ID for cmp in cmps.solids))
        self.solubles = tuple([i.ID for i in cmps if i.ID not in self.solids])
        self.disposal_cost = disposal_cost
        eff = self._eff = self.outs[0].proxy(f'{ID}_eff')
        sludge = self._sludge = self.outs[1].proxy(f'{ID}_sludge')
        self.effluent_pump = Pump(f'{self.ID}_eff', ins=eff, init_with=init_with)
        self.sludge_pump = Pump(f'{self.ID}_sludge', ins=sludge, init_with=init_with)
        self._mixed = self.ins[0].copy()
        self._set_split_at_mc()


    def _set_split_at_mc(self):
        mixed = self._mixed
        eff, sludge = self.outs
        mixed.mix_from(self.ins)
        mc = self.sludge_moisture
        mixed_F_mass = mixed.F_mass
        if mixed_F_mass == 0: # empty streams
            split = 0
        else:
            mixed_mc = mixed.imass['Water']/mixed_F_mass
            if mixed_mc < mc: # not enough water in the feeds
                warn(f'The set `sludge_moisture` {mc} is smaller than the '
                     f'moisture content in the influent ({mixed_mc:.3f}) and is ignored.')
                sludge.copy_like(mixed)
                eff.empty()
                split = 0
            else:
                solubles, solids = self.solubles, self.solids
                sludge.copy_flow(mixed)
                solids_mass = sludge.imass[solids].sum()
                sludge.imass[('Water', *solubles)] *= \
                    (solids_mass/(1-mc)-solids_mass)/(mixed_F_mass-solids_mass)
                eff.mass = mixed.mass - sludge.mass
                split = mixed.mass.value.copy()
                idx = np.where(mixed.mass!=0)
                split[idx] = eff.mass[idx]/mixed.mass[idx]
        self._isplit = self.thermo.chemicals.isplit(split)
        

    @staticmethod
    def _mc_at_split(split, solubles, mixed, eff, sludge, target_mc):
        eff.imass[solubles] = mixed.imass[solubles] * split
        sludge.imass[solubles] = mixed.imass[solubles] - eff.imass[solubles]
        mc = sludge.imass['Water'] / sludge.F_mass
        return mc-target_mc


    def _run(self):
        eff, sludge = self.outs
        solubles, solids = self.solubles, self.solids

        mixed = self._mixed
        mixed.mix_from(self.ins)
        eff.T = sludge.T = mixed.T
        sludge.copy_flow(mixed, solids, remove=True) # all solids go to sludge
        eff.copy_flow(mixed, solubles)

        flx.IQ_interpolation(
            f=self._mc_at_split, x0=1e-3, x1=1.-1e-3,
            args=(solubles, mixed, eff, sludge, self.sludge_moisture),
            checkbounds=False)
        self._set_split_at_mc()


    def _cost(self):
        m_solids = self.outs[-1].F_mass
        self.add_OPEX = {'Sludge disposal': m_solids*self.disposal_cost}
        power = 0
        for p in (self.effluent_pump, self.sludge_pump):
            p.simulate()
            power += p.power_utility.rate
        self.power_utility.rate = power


class BeltThickener(SludgeThickening):
    '''
    Gravity belt thickener (GBT) designed based on the manufacture specification
    data sheet. [1]_

    The 0th outs is the water-rich supernatant (effluent) and
    the 1st outs is the solid-rich sludge.

    Key parameters include:

        - Capacity: 80-100 m3/h.
        - Influent solids concentration: 0.2-1%.
        - Sludge cake moisture content: 90-96%.
        - Motor power: 3 (driving motor) and 1.1 (agitator motor) kW.
        - Belt width: 2.5 m.
        - Weight: 2350 kg.
        - Quote price: $3680 ea for three or more sets.

    The bare module (installation) factor is from Table 25 in Humbird et al. [2]_
    (solids handling equipment).

    Parameters
    ----------
    sludge_moisture : float
        Moisture content of the thickened sludge, [wt% water].
    solids : Iterable(str)
        IDs of the solid components.
        If not provided, will be set to the default `solids` attribute of the components.
    max_capacity : float
        Maximum hydraulic loading per belt thickener, [m3/h].
    power_demand : float
        Total power demand of each belt thickener, [kW].

    References
    ----------
    .. [1] `Industrial filtering equipment gravity thickener rotary thickening belt filter press. \
        <https://www.alibaba.com/product-detail/Industrial-filtering-equipment-gravity-thickener-rotary_60757627922.html?spm=a2700.galleryofferlist.normal_offer.d_title.78556be9t8szku>`_
        Data obtained on 7/21/2021.
    .. [2] Humbird et al., Process Design and Economics for Biochemical Conversion of
        Lignocellulosic Biomass to Ethanol: Dilute-Acid Pretreatment and Enzymatic
        Hydrolysis of Corn Stover; Technical Report NREL/TP-5100-47764;
        National Renewable Energy Lab (NREL), 2011.
        https://www.nrel.gov/docs/fy11osti/47764.pdf
    '''

    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 sludge_moisture=0.96, solids=(),
                 max_capacity=100, power_demand=4.1):
        SludgeThickening.__init__(self, ID, ins, outs, thermo,
                                sludge_moisture=sludge_moisture,
                                solids=solids)
        self.max_capacity = max_capacity
        self.power_demand = power_demand


    def _design(self):
        self._N_thickener = N = ceil(self._mixed.F_vol/self.max_capacity)
        self.design_results['Number of thickeners'] = N
        self.F_BM['Thickeners'] = 1.7 # ref [2]
        self.baseline_purchase_costs['Thickeners'] = 4000 * N
        
    def _cost(self):
        super()._cost()
        self.power_utility.rate += self.power_demand * self.N_thickener 


    @property
    def N_thickener(self):
        '''[int] Number of required belt thickeners.'''
        return self._N_thickener


class SludgeCentrifuge(SludgeThickening, SolidsCentrifuge):
    '''
    Solid centrifuge for sludge dewatering.

    `_run` and `_cost` are based on `SludgeThickening` and `_design`
    is based on `SolidsCentrifuge`.

    The 0th outs is the water-rich supernatant (effluent) and
    the 1st outs is the solid-rich sludge.

    Parameters
    ----------
    sludge_moisture : float
        Moisture content of the thickened sludge, [wt% water].
    solids : Iterable(str)
        IDs of the solid components.
        If not provided, will be set to the default `solids` attribute of the components.
    '''

    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 sludge_moisture=0.8, solids=(),
                 centrifuge_type='scroll_solid_bowl'):
        SludgeThickening.__init__(self, ID, ins, outs, thermo,
                                sludge_moisture=sludge_moisture,
                                solids=solids)
        self.centrifuge_type = centrifuge_type

    _run = SludgeThickening._run

    _design = SolidsCentrifuge._design

    _cost = SludgeThickening._cost


# %%

separator_path = ospath.join(data_path, 'sanunit_data/_sludge_separator.tsv')
allocate_N_removal = Decay.allocate_N_removal

class SludgeSeparator(SanUnit):
    '''
    For sludge separation based on
    `Trimmer et al. <https://doi.org/10.1021/acs.est.0c03296>`_,
    note that no default cost or environmental impacts are included.

    Parameters
    ----------
    ins : WasteStream
        Waste for treatment.
    outs : WasteStream
        Liquid, settled solids.
    split : float or dict
        Fractions of material retention in the settled solids.
        Default values will be used if not given.
    settled_frac : float
        Fraction of influent that settles as solids.
        The default value will be used if not given.

    Examples
    --------
    `bwaise systems <https://github.com/QSD-Group/EXPOsan/blob/main/exposan/bwaise/systems.py>`_

    References
    ----------
    [1] Trimmer et al., Navigating Multidimensional Social–Ecological System
    Trade-Offs across Sanitation Alternatives in an Urban Informal Settlement.
    Environ. Sci. Technol. 2020, 54 (19), 12641–12653.
    https://doi.org/10.1021/acs.est.0c03296.

    '''

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 split=None, settled_frac=None, **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with, **kwargs)

        data = load_data(path=separator_path)
        self.split = split or dct_from_str(data.loc['split']['expected'])
        self.settled_frac = settled_frac or float(data.loc['settled_frac']['expected'])
        del data


    _N_ins = 1
    _outs_size_is_fixed = False

    def _adjust_solid_water(self, influent, liq, sol):
        sol.imass['H2O'] = 0
        sol.imass['H2O'] = influent.F_mass * self.settled_frac - sol.F_mass
        if sol.imass['H2O'] < 0:
            sol.imass['H2O'] = 0
            msg = 'Negative water content calculated for settled solids, ' \
                'try smaller split or larger settled_frac.'
            warn(msg)
        liq.imass['H2O'] = influent.imass['H2O'] - sol.imass['H2O']
        return liq, sol

    def _run(self):
        waste = self.ins[0]
        liq, sol = self.outs[0], self.outs[1]

        # Retention in the settled solids
        sol_COD = liq_COD = None
        split = self.split
        if self._split_type == 'float':
            liq.copy_like(waste)
            sol.copy_like(waste)
            sol.mass *= self.split
            liq.mass -= sol.mass
        else:
            for var in self.split.keys():
                if var == 'TS':
                    sol.imass['OtherSS'] = split[var] * waste.imass['OtherSS']
                elif var == 'COD':
                    _COD = waste._COD or waste.COD
                    sol_COD = split[var] * _COD * waste.F_vol
                    liq_COD = _COD * waste.F_vol - sol_COD
                elif var == 'N':
                    N_sol = split[var]*(waste.imass['NH3']+waste.imass['NonNH3'])
                    NonNH3_rmd, NH3_rmd = \
                        allocate_N_removal(N_sol, waste.imass['NonNH3'])
                    sol.imass['NonNH3'] = NonNH3_rmd
                    sol.imass['NH3'] = NH3_rmd
                else:
                    sol.imass[var] = split[var] * waste.imass[var]
            liq.mass = waste.mass - sol.mass

        # Adjust total mass of of the settled solids by changing water content.
        liq, sol = self._adjust_solid_water(waste, liq, sol)
        sol._COD = sol._COD if not sol_COD else sol_COD / sol.F_vol
        liq._COD = liq._COD if not liq_COD else liq_COD / liq.F_vol


    @property
    def split(self):
        '''
        [float] or [dict] Fractions of material retention in the settled solids
        before degradation. If a single number is provided, then it is assumed
        that retentions of all Components in the WasteStream are the same.

        .. note::

            Set state variable values (e.g., COD) will be retained if the retention
            ratio is a single number (treated like the loss stream is split
            from the original stream), but not when the ratio is a dict.

        '''
        return self._split
    @split.setter
    def split(self, i):
        try:
            self._split = float(i)
            self._split_type = 'float'
        except:
            if isinstance(i, dict):
                self._split = i
                self._split_type = 'dict'
            else:
                raise TypeError(f'Only float or dict allowed, not {type(i).__name__}.')

    @property
    def settled_frac(self):
        '''[float] Fraction of influent that settles as solids.'''
        return self._settled_frac
    @settled_frac.setter
    def settled_frac(self, i):
        self._settled_frac = i