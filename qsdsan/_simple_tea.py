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

from datetime import date
import qsdsan as qs
from biosteam import TEA

__all__ = ('SimpleTEA',)

conflict_slots = ('lang_factor', 'system', 'units', 'feeds', 'products')

class SimpleTEA(TEA):    
    '''
    Calculate an annualized cost for simple economic analysis that does not
    include loan payment (i.e., 100% equity) and taxes.

    Parameters
    ----------
    system : :class:`biosteam.System`
        The system this TEA is conducted for.
    discount_rate : float
        Interest rate used in discounted cash flow analysis.
    start_year : int
        Start year of the system.
    lifetime : int
        Total lifetime of the system, [yr]. Currently `biosteam` only supports int.
    uptime_ratio : float
        Fraction of time that the system is operating.
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

    '''
    
    __slots__ = (*(i for i in TEA.__slots__ if i not in conflict_slots),
                 '_system', '_units', '_feeds', '_products',
                 '_discount_rate', '_start_year', '_lifetime',
                 '_uptime_ratio', '_operating_hours', '_CAPEX', '_lang_factor',
                 '_annual_maintenance', '_annual_labor', '_system_add_OPEX')
    
    def __init__(self, system, discount_rate=0.05,
                 start_year=date.today().year, lifetime=10, uptime_ratio=1., 
                 CAPEX=0., lang_factor=None,
                 annual_maintenance=0., annual_labor=0., system_add_OPEX={},
                 construction_schedule=None):
        system.simulate()
        self.system = system
        system._TEA = self
        self.discount_rate = discount_rate
        # IRR (internal rate of return) is the discount rate when net present value is 0
        self.IRR = discount_rate
        self._IRR = discount_rate # guess IRR for solve_IRR method
        self._sales = 0 # guess cost for solve_price method
        self.start_year = start_year
        self.lifetime = lifetime
        self.uptime_ratio = 1.
        self._lang_factor = None
        self._CAPEX = CAPEX
        self.lang_factor = lang_factor
        self.annual_maintenance = annual_maintenance
        self.annual_labor = annual_labor
        self.system_add_OPEX = system_add_OPEX
        if not construction_schedule:
            construction_schedule = (1,)
        self.construction_schedule = construction_schedule
        
        ########## Not relevant to SimpleTEA but required by TEA ##########
        # From U.S. IRS for tax purpose, won't matter when tax set to 0
        # Based on IRS Publication 946 (2019), MACRS15 should be used for
        # municipal wastewater treatment plant, but the system lifetime is
        # just 10 yrs or shorter, so changed to a shorter one
        self.depreciation = 'MACRS7'
        self.income_tax = 0.
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
        info += f'\nEAC  : {self.EAC:,.0f} {c}/yr'
        info += f'\nCAPEX: {self.CAPEX:,.0f} {c} (annualized to {self.annualized_CAPEX:,.0f} {c}/yr)'
        info += f'\nAOC  : {self.AOC:,.0f} {c}/yr'
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

    def get_unit_annualized_CAPEX(self, units):
        try: iter(units)
        except: units = (units, )
        CAPEX = 0
        r = self.discount_rate
        for unit in units:
            lifetime = unit.lifetime or self.lifetime
            if not unit._default_equipment_lifetime:
                CAPEX += unit.installed_cost*r/(1-(1+r)**(-lifetime))
            else:
                lifetime_dct = dict.fromkeys(unit.purchase_costs.keys(), lifetime)
                lifetime_dct.update(unit._default_equipment_lifetime)
                for equip, cost in unit.purchase_costs.items():
                    factor = unit.F_BM[equip]*\
                        unit.F_D.get(equip, 1.)*unit.F_P.get(equip, 1.)*unit.F_M.get(equip, 1.)
                    CAPEX += factor*cost*r/(1-(1+r)**(-lifetime_dct[equip]))
        return CAPEX

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
        '''[float] Interest rate used in discounted cash flow analysis.'''
        return self._discount_rate
    @discount_rate.setter
    def discount_rate(self, i):
        if 0 <= i <= 1:
            self._discount_rate = float(i)
        else:
            raise ValueError('`discount_rate` must be in [0,1].')
    
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
        return int(self._lifetime)
    @lifetime.setter
    def lifetime(self, i):
        self._lifetime = self._years = int(i)
        self._duration = (int(self.start_year), int(self.start_year+self.lifetime))

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
        if 0 <= i <= 1:
            self._uptime_ratio = float(i)
            self._operating_days = 365*float(i)
            self._operating_hours = self._operating_days * 24
        else:
            raise ValueError('`uptime_ratio` must be in [0,1].')

    @property
    def operating_days(self):
        '''[float] Equivalent operating days calculated based on `uptime_ratio`.'''
        return self._operating_days

    @property
    def operating_hours(self):
        '''[float] Equivalent operating hours calculated based on `uptime_ratio`.'''
        return self._operating_hours
    
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
    
    @property
    def annualized_CAPEX(self):
        '''[float] Annualized capital expenditure.'''
        return self.get_unit_annualized_CAPEX(self.units)

    
    @property
    def EAC(self):
        '''[float] Equivalent annual cost of the system.'''
        return self.annualized_CAPEX+self.AOC

    












