#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:

    Yalin Li <mailto.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''


# %%

import biosteam as bst, qsdsan as qs
from datetime import date
from biosteam import TEA as BSTTEA

__all__ = ('TEA', 'SimpleTEA',)

conflict_slots = ('lang_factor', 'system', 'units', 'feeds', 'products')

default_kwargs = dict(
    startup_months=0,
    startup_FOCfrac=1,
    startup_VOCfrac=1,
    startup_salesfrac=1,
    WC_over_FCI=0,
    finance_interest=0,
    finance_years=0,
    finance_fraction=0,
    )


class TEA(BSTTEA):
    '''
    Techno-economic analysis (TEA) with a simplified
    capital cost structure, unit-level operating cost components, and
    annualized cost metrics. Discounted cash flow is also included for 
    net present value (NPV) and internal rate of return (IRR) calculations.

    Key design choices:

    - Uses ``start_year`` + ``lifetime`` to indicate project duration.
    - Uses ``uptime_ratio`` (fraction in [0, 1]) to indicate operating time.
    - Collapses the capital cost hierarchy so DPI (direct permanent investment) 
      = TDC (total depreciable capital) = FCI (fixed capital investment) = installed equipment
      cost by default (no indirect cost adders applied unless ``_DPI``/``_TDC``/``_FCI``
      are overridden in a subclass).
    - Exposes ``CAPEX`` as a direct override for installed equipment cost, and
      ``lang_factor`` as an alternative to bare-module factors.
    - Decomposes fixed operating cost (FOC) into ``annual_maintenance`` (fraction of FCI),
      ``annual_labor``, and additional operating expenditures from individual units
      (``unit_add_OPEX``) and at the system level (``system_add_OPEX``).
    - Adds annualized cost properties: ``annualized_NPV``, ``annualized_CAPEX``,
      ``annualized_equipment_cost``, and ``EAC`` (equivalent annual cost).
    - Defaults to 100% equity financing, though loan financing is still available
      via ``finance_interest``, ``finance_years``, and ``finance_fraction``.

    Parameters
    ----------
    system : obj
        The system this TEA is conducted for.
    discount_rate : float
        Discount rate used in the discounted cash flow analysis.

        .. note::

            Herein ``discount_rate`` equals ``IRR`` (internal rate of return).
            Technically, IRR equals the discount rate only when NPV is 0.

    income_tax : float
        Combined tax rate (e.g., sum of national, state, and local levels)
        applied to net earnings.
    start_year : int
        Calendar year in which the system begins operation.
    lifetime : int
        Total operating lifetime of the system, [yr].

        .. note::

            The depreciation schedule must fit within the lifetime (its length
            must be <= ``lifetime``). The default ``'SL'`` (straight line) spans the
            whole lifetime and always fits, so there is no minimum. MACRS schedules
            run one year longer than their name (IRS half-year convention), e.g.
            ``'MACRS5'`` is a 6-year schedule (needs ``lifetime >= 6``) and
            ``'MACRS7'`` needs ``lifetime >= 8``. See ``depreciation``.

    uptime_ratio : float
        Fraction of time the system is operating, in [0, 1].
        A continuously operating system has ``uptime_ratio = 1``.

        .. note::

            If a unit has a different ``uptime_ratio`` than the system, the unit's
            value is used only when scaling its ``add_OPEX``. Utility and material
            costs are always scaled to the system's operating hours. Flows that
            do not match the system ``uptime_ratio`` should be normalized before
            being assigned to the unit (e.g., a pump that runs 50% of the time at
            50 kW should have ``power_utility`` set to 25 kW).

    CEPCI : float, optional
        Chemical Engineering Plant Cost Index used for equipment cost scaling.
        If None (default), the current ``qsdsan.CEPCI`` (i.e., ``biosteam.CE``) is
        left unchanged; pass a value (e.g., ``qsdsan.CEPCI_by_year[2023]``) to set it.
    CAPEX : float
        Total capital expenditure. When provided, overrides ``installed_equipment_cost``.
    lang_factor : float or None
        Multiplier applied to total equipment purchase cost to estimate installed
        cost. Mutually exclusive with ``CAPEX``; leave as ``None`` when ``CAPEX``
        is provided. If neither is given, installed cost is summed from each unit's
        bare-module factors.
    annual_maintenance : float
        Annual maintenance cost as a fraction of fixed capital investment (FCI).
    annual_labor : float
        Annual labor cost [USD/yr].
    system_add_OPEX : float or dict
        Additional annual operating expenditure at the system level, on top of
        the ``add_OPEX`` of individual units. A float is automatically converted
        to a dict keyed ``"System additional OPEX"``.
    depreciation : str
        Depreciation schedule: ``'SL'`` (straight line, default), ``'DDB'``
        (double-declining balance), ``'SYD'`` (sum-of-years-digits), or a MACRS
        schedule (``'MACRS3'``, ``'MACRS5'``, ``'MACRS7'``, ``'MACRS10'``, ...).
        The schedule length must be <= ``lifetime``. Depreciation only affects
        results when there is taxable income to shield (i.e. ``income_tax`` > 0
        and positive net earnings).
    construction_schedule : tuple
        Fraction of total capital invested in each year prior to start-up; must
        sum to 1. Use the default ``(0, 1)`` if no staged construction is needed.
    accumulate_interest_during_construction : bool
        If ``False`` (default), loan interest accrued during construction is paid
        from equity and not rolled into the loan principal.
        If ``True``, accrued interest is capitalized onto the loan principal.
        See https://github.com/BioSTEAMDevelopmentGroup/biosteam/issues/180
        for details.
    simulate_system : bool
        Whether to simulate the system before creating the TEA object.
    simulate_kwargs : dict
        Keyword arguments passed to ``system.simulate()`` when ``simulate_system``
        is ``True``.
    tea_kwargs
        Additional keyword arguments for the underlying cash flow model.
        Defaults (in parentheses): ``startup_months`` (0), ``startup_FOCfrac`` (1),
        ``startup_VOCfrac`` (1), ``startup_salesfrac`` (1), ``WC_over_FCI`` (0),
        ``finance_interest`` (0), ``finance_years`` (0), ``finance_fraction`` (0).

    Examples
    --------
    A system should be constructed prior to TEA, here we import a pre-constructed one.

    >>> import qsdsan as qs
    >>> from qsdsan.utils import create_example_system
    >>> sys = create_example_system()
    >>> # Uncomment the line below to see the system diagram
    >>> # sys.diagram()
    >>> sys.simulate()
    >>> sys.show() # doctest: +ELLIPSIS
    System: sys
    ins...
    [0] salt_water...
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): H2O   111
                        NaCl  0.856
    [1] methanol...
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow...Methanol
    [2] ethanol...
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow...Ethanol
    outs...
    [0] alcohols...
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Methanol  0.624
                        Ethanol   0.217
    [1] waste_brine...
        phase: 'l', T: 350 K, P: 101325 Pa
        flow (kmol/hr): H2O   88.8
                        NaCl  0.684
    >>> tea = qs.TEA(system=sys, discount_rate=0.05, start_year=2021,
    ...              lifetime=10, uptime_ratio=0.9,
    ...              system_add_OPEX=0.03)
    >>> # Your results maybe slightly different depending on the version of
    >>> # QSDsan's dependent packages (e.g., thermo)
    >>> tea.show() # doctest: +ELLIPSIS
    TEA: sys
    NPV  : -259,...

    See Also
    --------
    `TEA <https://qsdsan.readthedocs.io/en/latest/tutorials/7_TEA.html>`_
    '''

    __slots__ = (*(i for i in BSTTEA.__slots__ if i not in conflict_slots),
                 '_system', '_units', '_feeds', '_products',
                 '_discount_rate', '_start_year', '_lifetime',
                 '_uptime_ratio', '_operating_hours', '_CAPEX', '_lang_factor',
                 '_annual_maintenance', '_annual_labor', '_system_add_OPEX')

    def __init__(self, system, discount_rate=0.05, income_tax=0.,
                 CEPCI=None, start_year=date.today().year,
                 lifetime=10, uptime_ratio=1.,
                 CAPEX=0., lang_factor=None,
                 annual_maintenance=0., annual_labor=0., system_add_OPEX={},
                 depreciation='SL',
                 construction_schedule=(0, 1), accumulate_interest_during_construction=False,
                 simulate_system=True, simulate_kwargs={},
                 **tea_kwargs):
        # Set the cost index (CEPCI) before simulation so it applies to costing; if not
        # provided, leave the current `qsdsan.CEPCI` (i.e., `biosteam.CE`) untouched
        # rather than resetting it (previously the default froze `bst.CE` at import time).
        if CEPCI is not None: self.CEPCI = CEPCI
        if simulate_system: system.simulate(**simulate_kwargs)
        self.system = system
        system._TEA = self
        # IRR (internal rate of return) is the discount rate when net present value is 0
        self.IRR = discount_rate
        self._IRR = discount_rate # guess IRR for solve_IRR method
        self.income_tax = income_tax
        self._sales = 0 # guess cost for solve_price method
        self._depreciation = None # initialize this attribute
        self.start_year = start_year
        self.lifetime = lifetime
        self.uptime_ratio = 1.
        self.annual_maintenance = annual_maintenance
        self.annual_labor = annual_labor
        self.system_add_OPEX = {}.copy() if not system_add_OPEX else system_add_OPEX
        self.depreciation = depreciation
        self.construction_schedule = construction_schedule
        self.accumulate_interest_during_construction = accumulate_interest_during_construction 
        default_kwargs.update(tea_kwargs)
        for k, v in default_kwargs.items():
            setattr(self, k, v)
        self._lang_factor = None
        self._CAPEX = CAPEX
        self.lang_factor = lang_factor

    def __repr__(self):
        return f'<{type(self).__name__}: {self.system.ID}>'

    def show(self):
        c = self.currency
        info = f'{type(self).__name__}: {self.system.ID}'
        info += f'\nNPV  : {self.NPV:,.0f} {c} at {self.discount_rate:.1%} discount rate'
        print(info)

    _ipython_display_ = show

    def _add_first_replacement_costs(self, nontaxable_cashflow):
        system = self.system
        units = system.unit_capital_costs.values() if isinstance(system, bst.AgileSystem) else system.cost_units
        lang_factor = self.lang_factor
        start = self._start
        end = start + self._years
        for unit in units:
            lifetime = unit.lifetime
            if not lifetime:
                continue
            if lang_factor:
                installed_costs = {i: j*lang_factor for i, j in unit.purchase_costs.items()}
            else:
                installed_costs = unit.installed_costs
            if isinstance(lifetime, int):
                replacement_index = start + lifetime
                if replacement_index < end:
                    nontaxable_cashflow[replacement_index] -= sum(installed_costs.values())
            elif isinstance(lifetime, dict):
                for name, installed_cost in installed_costs.items():
                    equipment_lifetime = lifetime.get(name)
                    if equipment_lifetime:
                        replacement_index = start + equipment_lifetime
                        if replacement_index < end:
                            nontaxable_cashflow[replacement_index] -= installed_cost

    def _taxable_nontaxable_depreciation_cashflows(self):
        taxable_cashflow, nontaxable_cashflow, depreciation = super()._taxable_nontaxable_depreciation_cashflows()
        self._add_first_replacement_costs(nontaxable_cashflow)
        return taxable_cashflow, nontaxable_cashflow, depreciation


    def _DPI(self, installed_equipment_cost):
        return installed_equipment_cost

    def _TDC(self, DPI):
        return DPI

    def _FCI(self, TDC):
        return TDC

    def _FOC(self, FCI):
        return FCI*self.annual_maintenance+self.annual_labor+self.total_add_OPEX

    @property
    def system(self):
        '''[:class:`biosteam.System`] The system this TEA is conducted for.'''
        return self._system
    @system.setter
    def system(self, i):
        if i:
            self._system = i
            self._units = sorted([j for j in i.units if j._design or j._cost],
                                key=lambda x: x.line)
            self._feeds = i.feeds
            self._products = i.products
            try:
                self.system._TEA = self
            except AttributeError:
                pass

    @property
    def units(self):
        '''[:class:`qsdsan.SanUnit`] Units in the system.'''
        return self._units

    @property
    def feeds(self):
        '''[:class:`qsdsan.WasteStream`] System feed streams.'''
        return self._feeds

    @property
    def products(self):
        '''[:class:`qsdsan.WasteStream`] System product streams.'''
        return self._products

    @property
    def discount_rate(self):
        '''[float] Interest rate used in discounting, same as `IRR` in :class:`biosteam.TEA`.'''
        return self.IRR
    @discount_rate.setter
    def discount_rate(self, i):
        self.IRR = i

    @property
    def CEPCI(self):
        '''[float] Chemical Engineering Plant Cost Index.'''
        return bst.CE
    @CEPCI.setter
    def CEPCI(self, i):
        bst.CE = i

    @property
    def CEPCI_by_year(self):
        '''
        [dict] Chemical Engineering Plant Cost Index with key being the year
        and values being the index. Same as ``qsdsan.CEPCI_by_year``.
        '''
        return qs.CEPCI_by_year

    @property
    def start_year(self):
        '''[int] Start year of the system.'''
        return self._start_year
    @start_year.setter
    def start_year(self, i):
        self._start_year = i

    @property
    def lifetime(self):
        '''[int] Total lifetime of the system, [yr]. Currently `biosteam` only supports int.'''
        return self._lifetime
    @lifetime.setter
    def lifetime(self, i):
        self._lifetime = self._years = int(i)
        self._duration = (int(self.start_year), int(self.start_year+self.lifetime))
        self.depreciation = self.depreciation

    @property
    def duration(self):
        '''[int] Duration of the system based on start_year and lifetime.'''
        return self._duration

    @property
    def uptime_ratio(self):
        '''[float] Fraction of time that the system is operating.'''
        return self._uptime_ratio
    @uptime_ratio.setter
    def uptime_ratio(self, i):
        if 0 <=i<= 1:
            self._uptime_ratio = float(i)
            self._operating_days = 365*float(i)
            self._operating_hours = self.system.operating_hours = self._operating_days * 24
        else:
            raise ValueError('`uptime_ratio` must be in [0,1].')

    @property
    def operating_days(self):
        '''[float] Equivalent operating days calculated based on `uptime_ratio`.'''
        return self._operating_days
    @operating_days.setter
    def operating_days(self, i):
        raise AttributeError('Set `uptime_ratio` instead.')

    @property
    def operating_hours(self):
        '''[float] Equivalent operating hours calculated based on `uptime_ratio`.'''
        return self._operating_hours
    @operating_hours.setter
    def operating_hours(self, i):
        raise AttributeError('Set `uptime_ratio` instead.')

    @property
    def lang_factor(self):
        '''[float] A factor to estimate the total installation cost based on equipment purchase cost.'''
        return self._lang_factor or None
    @lang_factor.setter
    def lang_factor(self, i):
        if self.CAPEX is not None:
            if i is not None:
                raise AttributeError('`CAPEX` provided, `lang_factor` cannot be set. '
                                     'The calculated `lang_factor is` '
                                     f'{self.installed_equipment_cost/self.purchase_cost:.1f}.')
            else:
                self._lang_factor = None
        elif i >=1:
            self._lang_factor = float(i)
        else:
            raise ValueError('`lang_factor` must >= 1.')

    @property
    def currency(self):
        '''[str] TEA currency, same with `qsdsan.currency`.'''
        return qs.currency
    @currency.setter
    def currency(self, i):
        raise AttributeError('Currency can only be changed through `qsdsan.currency`.')

    @property
    def installed_equipment_cost(self):
        '''[float] Sum of installed cost of all units in the system, is the same as `CAPEX` if `CAPEX` is provided.'''
        if self._CAPEX:
            return self._CAPEX
        if self.lang_factor:
            return self.purchase_cost*self.lang_factor
        return sum([u.installed_cost for u in self.units])

    @property
    def DPI(self):
        '''[float] Direct permanent investment, calculated using `self._DPI` as a function of `installed_equipment_cost`.'''
        return self._DPI(self.installed_equipment_cost)

    @property
    def TDC(self):
        '''[float] Total depreciable capital, calculated using `self._TDC` as a function of `self.DPI`.'''
        return self._TDC(self.DPI)

    @property
    def FCI(self):
        '''[float] Fixed capital investment, calculated using `self._FCI` as a function of `self.TDC`.'''
        return self._FCI(self.TDC)

    @property
    def TCI(self):
        '''[float] Total capital investment, calculated as `self._FCI*(1+self.WC_over_FCI)` (WC for working capital).'''
        return self.FCI*(1+self.WC_over_FCI)

    @property
    def CAPEX(self):
        '''[float] Capital expenditure, if not provided, is set to be the same as `self.TCI`.'''
        return self.TCI

    @property
    def annual_maintenance(self):
        '''[float] Annual maintenance cost as a fraction of fixed capital investment.'''
        return self._annual_maintenance
    @annual_maintenance.setter
    def annual_maintenance(self, i):
        if 0 <= i <= 1:
            self._annual_maintenance = float(i)
        else:
            raise ValueError('`annual_maintenance` must be in [0,1].')

    @property
    def annual_labor(self):
        '''[float] Annual labor cost.'''
        return self._annual_labor
    @annual_labor.setter
    def annual_labor(self, i):
        self._annual_labor = float(i)

    @property
    def unit_add_OPEX(self):
        '''[float] Sum of `add_OPEX` for all units in the system.'''
        tot = 0
        for u in self.units:
            add_OPEX = sum(v for v in u.add_OPEX.values())
            tot += add_OPEX*u.uptime_ratio/self.uptime_ratio*self._operating_hours
        return tot

    @property
    def system_add_OPEX(self):
        '''
        [dict] Annual additional system-wise operating expenditure
        (on top of the `add_OPEX` of each unit).
        Float input will be automatically converted to a dict with the key being
        "System additional OPEX".
        '''
        return {'System additional OPEX': self._system_add_OPEX} \
               if isinstance(self._system_add_OPEX, float) else self._system_add_OPEX
        return self._system_add_OPEX
    @system_add_OPEX.setter
    def system_add_OPEX(self, i):
        if isinstance(i, float):
            i = {'System additional OPEX': i}
        if not isinstance(i, dict):
            raise TypeError('system_add_OPEX can only be float of dict, ' \
                            f'not {type(i).__name__}.')
        self._system_add_OPEX = i

    @property
    def total_add_OPEX(self):
        '''[float] Sum of `unit_add_OPEX` and `system_add_OPEX`.'''
        return self.unit_add_OPEX+sum(v for v in self.system_add_OPEX.values())

    @property
    def FOC(self):
        '''
        [float] Fixed operating cost, including maintenance, labor, and any additional
        operating expenditure other than chemical inputs and utilities.
        '''
        return self._FOC(self.FCI)

    def _get_annuity_factor(self, yrs=None):
        yrs = yrs or self._years
        r = self.discount_rate
        return (1-(1+r)**(-yrs))/r

    @property
    def annualized_NPV(self):
        r'''
        [float] Annualized NPV calculated as:

        .. math::

            annualized\ NPV = \frac{NPV*r}{(1-(1+r)^{-lifetime})}
        '''
        return self.NPV/self._get_annuity_factor()

    @property
    def annualized_CAPEX(self):
        r'''
        [float] Annualized capital expenditure calculated through annualized NPV as:

        .. math::

            annualized\ capital\ cost = annual\ net\ earning - annualized\ NPV

        .. note::

            Read the `tutorial <https://github.com/QSD-Group/QSDsan/blob/main/docs/source/tutorials/7_TEA.ipynb>`_
            about the difference between `annualized_CAPEX` and `annualized_equipment_cost`.
        '''
        return self.net_earnings-self.annualized_NPV


    def get_unit_annualized_equipment_cost(self, units=None):
        r'''
        Annualized equipment cost representing the sum of the annualized
        cost of each equipment, which is calculated as:

        .. math::

            annualized\ equipment\ cost = \frac{equipment\ installed\ cost}{(1-(1+r)^{-lifetime})}

        .. note::

            Read the `tutorial <https://github.com/QSD-Group/QSDsan/blob/main/docs/source/tutorials/7_TEA.ipynb>`_
            about the difference between `annualized_CAPEX` and `annualized_equipment_cost`.
        '''
        units = units or self.units
        try: iter(units)
        except: units = (units,)
        cost = 0
        get_A = self._get_annuity_factor

        for unit in units:
            lifetime = unit.lifetime or self.lifetime
            # no unit equipment lifetime or unit equipment lifetime given as a number
            # (i.e., no individual equipment lifetime)
            if not isinstance(lifetime, dict):
                cost += unit.installed_cost/get_A(lifetime)
            else:
                lifetime_dct = dict.fromkeys(unit.purchase_costs.keys())
                lifetime_dct.update(lifetime)
                for equip, cost in unit.purchase_costs.items():
                    factor = unit.F_BM[equip]*\
                        unit.F_D.get(equip, 1.)*unit.F_P.get(equip, 1.)*unit.F_M.get(equip, 1.)
                    # for equipment that does not have individual lifetime
                    # use the unit lifetime or TEA lifetime
                    equip_lifetime = lifetime_dct[equip] or self.lifetime
                    cost += factor*cost/get_A(equip_lifetime)
        return cost


    @property
    def annualized_equipment_cost(self):
        '''
        [float] Annualized equipment cost representing the sum of the annualized
        cost of each equipment, calculated using ``get_unit_annualized_equipment_cost``.
        '''
        return self.get_unit_annualized_equipment_cost()


    @property
    def EAC(self):
        '''
        [float] Equivalent annual cost calculated as the sum of `annualized_CAPEX` and
        `AOC` (annual operating cost).

        .. note::

            Sales are not included.
        '''
        return self.annualized_CAPEX+self.AOC
    

class SimpleTEA(TEA):
    '''
    .. deprecated::
        Use :class:`TEA` instead. ``SimpleTEA`` is an alias kept for backward
        compatibility and will be removed in a future version.
    '''
    def __init__(self, *args, **kwargs):
        import warnings
        warnings.warn(
            'SimpleTEA is deprecated and will be removed in a future version; '
            'use TEA instead.',
            DeprecationWarning,
            stacklevel=2,
        )
        super().__init__(*args, **kwargs)
