#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Lewis Rowles <stetsonsc@gmail.com>
    Yalin Li <zoe.yalin.li@gmail.com>
    Lane To <lane20@illinois.edu>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

# %%

from . import SludgeSeparator
from .. import Construction
from ..utils import load_data, data_path, ospath, price_ratio

__all__ = ('ScrewPress',)

screw_path = ospath.join(data_path,'sanunit_data/_screw_press.tsv')

#!!! See GitHub issue here for question about the data references:
# https://github.com/QSD-Group/QSDsan/issues/58

@price_ratio(default_price_ratio=1)
class ScrewPress(SludgeSeparator):
    '''
    Screw Press is used for dewatering where sludge, conditioned with cationic
    polymer, is fed into the unit. Sludge is continuously dewatered as it travels
    along the screw.

    The following components should be included in system thermo object for simulation:
    'H2O', `Polyacrylamide`.

    The following impact items should be pre-constructed for life cycle assessment:
    `Steel`.

    Parameters
    ----------
    ins : Iterable(stream)
        Waste for treatment (e.g., wastewater or latrine sludge) and polymer added for treatment.
    outs : Iterable(stream)
        Treated liquids and solids.

    References
    ----------
    [1] Tchobanoglous, G.; Stensel, H. D.; Tsuchihashi, R.; Burton, F.; Abu-Orf,
    M.; Bowden, G.; Pfrang, W. Wastewater Engineering: Treatment and Resource
    Recovery, 5th ed.; Metcalf & Eddy, Inc., AECOM, McGraw-Hill: New York, 2014.

    '''

    def __init__(self, ID='', ins=None, outs=(),thermo=None, init_with='WasteStream',
                 split=None, settled_frac=None, if_N2O_emission=False, **kwargs):
        SludgeSeparator.__init__(self, ID, ins, outs, thermo, init_with,
                                 split, settled_frac, F_BM_default=1)

        self.construction = (
            Construction('steel', linked_unit=self, item='Steel', quantity_unit='kg'))

        data = load_data(path=screw_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)

        #!!! What is this doing?
        #**** how can I use this distribution for solids capture?
        #change name in csv to settled_frac
        # self.dewatering_solids_capture = self.settled_frac

    _N_ins = 2
    _N_outs = 2

    def _run(self):
        waste, polymer = self.ins
        liq, cake_sol = self.outs
        SludgeSeparator._run(self)

        #*** is this needed in both places???? Will it be pulled into the function? Isnt cake solids TS needed too?
        # self.dewatering_solids_capture = self.settled_frac
        liq, cake_sol = self._adjust_solid_water(waste, liq, cake_sol)

        #!!! Notes still useful?
        # need the flowrate of the liquid stream or use equations below with dewatering_solids_flowrate
        # self.dewatering_solids_flowrate = liq.F_vol

        # note to self change these once confirm outputs from above function...

        # Does the funtion above for self._adjust_solid_water take care of these equations?
        # self.waste_sol_flow = waste.F_mass # kg TS/hr
        # dewatering_solids_conc = (self.waste_sol_flow) * self.dewatering_solids_capture    # kg TS/hr
        # dewatering_solids_flowrate = dewatering_solids_conc / (1.06 * self.cake_solids_TS * 1000) # m3 / hr
        # liq.F_vol = self.dewatering_centrate_flowrate = (waste.F_vol - dewatering_solids_flowrate)   # m3 / hr
        # dewatering_centrate_conc = (self.waste_sol_flow) * (1 - self.dewatering_solids_capture) # kg / d
        # dewatering_total_solids = dewatering_solids_conc / dewatering_solids_flowrate # kg / m3

        waste_TS = waste.F_mass - waste.imass['H2O']
        polymer.imass['Polyacrylamide'] = self.dewatering_polymer_dose * waste_TS / 1000 # kg polymer / hr


    def _design(self):
        self.construction[0].quantity = self.dewatering_screw_press_steel
        self.add_construction()

    def _cost(self):
        self.baseline_purchase_costs['Screw Press'] = self.dewatering_screw_press_cost*self.price_ratio
        self.power_utility(self.dewatering_energy_demand*self.outs[1].F_mass) # kWh/hr