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

import qsdsan as qs
from datetime import date
from biosteam import TEA

__all__ = ('SimpleTEA',)

conflict_slots = ('lang_factor', 'system', 'units', 'feeds', 'products')

class SimpleTEA(TEA):
    '''
    Calculate an annualized cost for simple economic analysis that does not
    include loan payment (i.e., 100% equity).

    Parameters
    ----------
    system : :class:`biosteam.System`
        The system this TEA is conducted for.
    discount_rate : float
        Interest rate used in discounted cash flow analysis.

        .. note::

            Herein `discount_rate` equals to `IRR` (internal rate of return).
            Although theoretically, IRR is the discount rate only when the
            net present value (NPV) is 0.

    income_tax : float
        Combined tax (e.g., sum of national, state, local levels) for net earnings.
    start_year : int
        Start year of the system.
    lifetime : int
        Total lifetime of the system, [yr]. Currently `biosteam` only supports int.

        .. note::

            As :class:`SimpleTEA` is a subclass of :class:`biosteam.TEA`,
            and :class:`biosteam.TEA` currently only supports certain
            depreciation schedules, lifetime must be larger than or equal to 6.


    uptime_ratio : float
        Fraction of time that the system is operating, should be in [0,1]
        (i.e., a system that is always operating has an uptime_ratio of 1).

        .. note::

            If a unit has an `uptime_ratio` that is different from the `uptime_ratio`
            of the system, the `uptime_ratio` of the unit will be used in calculating
            the additional operation expenses (provided in `unit.add_OPEX`).

            However, `uptime_ratio` of the unit will not affect the utility
            (heating, cooling, power) and material costs/environmental impacts.

            For example, if the `uptime_ratio` of the system and the unit are
            1 and 0.5, respectively, then in calculating operating expenses
            associated with the unit:

                - Utility and material costs/environmental impacts will be calculated for 1*24*365 hours per year.
                - Additional operaitng expenses will be calculated for 0.5*24*365 hours per year.

            If utility and material flows are not used at the same `uptime_ratio`
            as the system, they should be normalized to be the same.
            For example, if the system operates 100% of time but a pump only works
            50% of the pump at 50 kW. Set the pump `power_utility` to be 50*50%=25 kW.


    CAPEX : float
        Capital expenditure, if not provided, is set to be the same as `installed_equipment_cost`.
    lang_factor : float or None
        A factor to estimate the total installation cost based on equipment purchase cost,
        leave as ``None`` if providing ``CAPEX``.
        If neither ``CAPEX`` nor ``lang_factor`` is provided,
        ``installed_equipment_cost`` will be calculated as the sum of purchase costs
        of all units within the system.
    annual_maintenance : float
        Annual maintenance cost as a fraction of fixed capital investment.
    annual_labor : float
        Annual labor cost.
    system_add_OPEX : float or dict
        Annual additional system-wise operating expenditure (on top of the `add_OPEX` of each unit).
        Float input will be automatically converted to a dict with the key being
        "System additional OPEX".
    construction_schedule : tuple or None
        Construction progress, must sum up to 1, leave as `None` will assume the system finishes within one year.

    Examples
    --------
    A system should be constructed prior to TEA, here we import a pre-constructed one.

    >>> import qsdsan as qs
    >>> from qsdsan.utils import load_example_cmps, load_example_sys
    >>> cmps = load_example_cmps()
    >>> sys = load_example_sys(cmps)
    >>> # Uncomment the line below to see the system diagram
    >>> # sys.diagram()
    >>> sys.simulate()
    >>> sys.show()
    System: sys
    ins...
    [0] salt_water
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): H2O   111
                        NaCl  0.856
    [1] methanol
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Methanol  0.624
    [2] ethanol
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Ethanol  0.217
    outs...
    [0] alcohols
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Methanol  0.624
                        Ethanol   0.217
    [1] waste_brine
        phase: 'l', T: 350 K, P: 101325 Pa
        flow (kmol/hr): H2O   88.8
                        NaCl  0.684
    >>> tea = qs.SimpleTEA(system=sys, discount_rate=0.05, start_year=2021,
    ...                     lifetime=10, uptime_ratio=0.9,
    ...                     system_add_OPEX=0.03)
    >>> tea.show()
    SimpleTEA: sys
    NPV  : -258,730 USD at 5.0% discount rate

    See Also
    --------
    `SanUnit and System <https://qsdsan.readthedocs.io/en/latest/tutorials/SanUnit_and_System.html>`_

    `TEA and LCA <https://qsdsan.readthedocs.io/en/latest/tutorials/TEA_and_LCA.html>`_
    '''

    __slots__ = (*(i for i in TEA.__slots__ if i not in conflict_slots),
                 '_system', '_units', '_feeds', '_products',
                 '_discount_rate', '_start_year', '_lifetime',
                 '_uptime_ratio', '_operating_hours', '_CAPEX', '_lang_factor',
                 '_annual_maintenance', '_annual_labor', '_system_add_OPEX')

    def __init__(self, system, discount_rate=0.05, income_tax=0.,
                 start_year=date.today().year, lifetime=10, uptime_ratio=1.,
                 CAPEX=0., lang_factor=None,
                 annual_maintenance=0., annual_labor=0., system_add_OPEX={},
                 depreciation='SL', construction_schedule=(1,)):
        system.simulate()
        self.system = system
        system._TEA = self
        self.income_tax = income_tax
        # IRR (internal rate of return) is the discount rate when net present value is 0
        self.IRR = discount_rate
        self._IRR = discount_rate # guess IRR for solve_IRR method
        self._sales = 0 # guess cost for solve_price method
        self._depreciation = None # initialize this attribute
        self.start_year = start_year
        self.lifetime = lifetime
        self.uptime_ratio = 1.
        self._lang_factor = None
        self._CAPEX = CAPEX
        self.lang_factor = lang_factor
        self.annual_maintenance = annual_maintenance
        self.annual_labor = annual_labor
        self.system_add_OPEX = system_add_OPEX
        self.depreciation = depreciation
        self.construction_schedule = construction_schedule

        ########## Not relevant to SimpleTEA but required by TEA ##########
        self.startup_months = 0.
        self.startup_FOCfrac = 0.
        self.startup_VOCfrac = 0.
        self.startup_salesfrac = 0.
        self.WC_over_FCI = 0.
        self.finance_interest = 0.
        self.finance_years = 0
        self.finance_fraction = 0


    def __repr__(self):
        return f'<{type(self).__name__}: {self.system.ID}>'

    def show(self):
        c = self.currency
        info = f'{type(self).__name__}: {self.system.ID}'
        info += f'\nNPV  : {self.NPV:,.0f} {c} at {self.discount_rate:.1%} discount rate'
        print(info)

    _ipython_display_ = show


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
        '''[float] Direct permanent investment, same as `installed_equipment_cost`.'''
        return self._DPI(self.installed_equipment_cost)

    @property
    def TDC(self):
        '''[float] Total depreciable capital, same as `installed_equipment_cost`.'''
        return self._TDC(self.DPI)

    @property
    def FCI(self):
        '''[float] Fixed capital investment, same as `installed_equipment_cost`.'''
        return self._FCI(self.TDC)

    @property
    def TCI(self):
        '''[float] Total capital investment, same as `installed_equipment_cost`.'''
        return self.FCI

    @property
    def CAPEX(self):
        '''[float] Capital expenditure, if not provided, is set to be the same as `installed_equipment_cost`.'''
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
        return {'System dditional OPEX': self._system_add_OPEX} \
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
        operaitng expenditure other than chemical inputs and utilities.
        '''
        return self._FOC(self.FCI)

    def _get_annuity_factor(self, yrs):
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

    @property
    def annualized_equipment_cost(self):
        r'''
        [float] Annualized equipment cost representing the sum of the annualized
        cost of each equipment, which is as:

        .. math::

            annualized\ equipment\ cost = \frac{equipment\ installed\ cost}{(1-(1+r)^{-lifetime})}

        .. note::

            Read the `tutorial <https://github.com/QSD-Group/QSDsan/blob/main/docs/source/tutorials/7_TEA.ipynb>`_
            about the difference between `annualized_CAPEX` and `annualized_equipment_cost`.
        '''
        units = self.units
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
    def EAC(self):
        '''
        [float] Equvalent annual cost calculated as the sum of `annualized_CAPEX` and
        `AOC` (annual operating cost).
        '''
        return self.annualized_CAPEX+self.AOC