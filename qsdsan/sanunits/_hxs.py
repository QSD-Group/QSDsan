#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>
    Lewis Rowles <stetsonsc@gmail.com>
    Joy Zhang <joycheung1994@gmail.com>
    Lane To <lane20@illinois.edu>

Part of this module is based on the biosteam package:
https://github.com/BioSTEAMDevelopmentGroup/biosteam

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

from biosteam.units import HXprocess, HXutility
from .. import SanUnit, Construction
from ..utils import ospath, load_data, data_path, price_ratio

__all__ = ('HXprocess', 'HXutility', 'OilHeatExchanger', )


# %%

# =============================================================================
# BioSTEAM-inherited ones
# =============================================================================

class HXprocess(SanUnit, HXprocess):
    '''
    Similar to :class:`biosteam.units.HXprocess`,
    but can be initilized with :class:`qsdsan.SanStream` and :class:`qsdsan.WasteStream`.

    See Also
    --------
    `biosteam.units.HXprocess <https://biosteam.readthedocs.io/en/latest/units/heat_exchange.html>`_
    '''

    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='Stream', F_BM_default=None,
                 *, U=None, dT=5., T_lim0=None, T_lim1=None,
                 material="Carbon steel/carbon steel",
                 heat_exchanger_type="Floating head",
                 N_shells=2, ft=None,
                 phase0=None,
                 phase1=None,
                 H_lim0=None,
                 H_lim1=None):
        SanUnit.__init__(self, ID, ins, outs, thermo,
                         init_with=init_with, F_BM_default=F_BM_default)

        #: [float] Enforced overall heat transfer coefficent (kW/m^2/K)
        self.U = U

        #: [float] Total heat transfered.
        self.Q = None

        #: Number of shells for LMTD correction factor method.
        self.N_shells = N_shells

        #: User imposed correction factor.
        self.ft = ft

        #: [float] Pinch temperature difference.
        self.dT = dT

        #: [float] Temperature limit of outlet stream at index 0.
        self.T_lim0 = T_lim0

        #: [float] Temperature limit of outlet stream at index 1.
        self.T_lim1 = T_lim1

        #: [float] Temperature limit of outlet stream at index 0.
        self.H_lim0 = H_lim0

        #: [float] Temperature limit of outlet stream at index 1.
        self.H_lim1 = H_lim1

        #: Enforced phase of outlet at index 0
        self.phase0 = phase0

        #: Enforced phase of outlet at index 1
        self.phase1 = phase1

        self.material = material
        self.heat_exchanger_type = heat_exchanger_type


class HXutility(SanUnit, HXutility):
    '''
    Similar to :class:`biosteam.units.HXutility`,
    but can be initilized with :class:`qsdsan.SanStream` and :class:`qsdsan.WasteStream`.

    See Also
    --------
    `biosteam.units.HXutility <https://biosteam.readthedocs.io/en/latest/units/heat_exchange.html>`_
    '''

    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='Stream', F_BM_default=None,
                 *, T=None, V=None, rigorous=False, U=None, H=None,
                 heat_exchanger_type="Floating head",
                 material="Carbon steel/carbon steel",
                 N_shells=2, ft=None, heat_only=None, cool_only=None):
        SanUnit.__init__(self, ID, ins, outs, thermo,
                         init_with=init_with, F_BM_default=F_BM_default)
        self.T = T #: [float] Temperature of outlet stream (K).
        self.V = V #: [float] Vapor fraction of outlet stream.
        self.H = H #: [float] Enthalpy of outlet stream.

        #: [bool] If true, calculate vapor liquid equilibrium
        self.rigorous = rigorous

        #: [float] Enforced overall heat transfer coefficent (kW/m^2/K)
        self.U = U

        #: Number of shells for LMTD correction factor method.
        self.N_shells = N_shells

        #: User imposed correction factor.
        self.ft = ft

        #: [bool] If True, heat exchanger can only heat.
        self.heat_only = heat_only

        #: [bool] If True, heat exchanger can only cool.
        self.cool_only = cool_only

        self.material = material
        self.heat_exchanger_type = heat_exchanger_type


# %%

oilhx_path = ospath.join(data_path,'sanunit_data/_oil_heat_exchanger.tsv')

@price_ratio(default_price_ratio=1)
class OilHeatExchanger(SanUnit):
    '''
    Oil heat exhanger utilizes an organic Rankin cycle. This CHP system is
    used to generate additional electricity that the refinery and/or
    facility can use to decrease the units electrical demand on the
    electrical grid. This type of system is required for ISO 31800
    certification as the treatment unity needs to be energy independent
    when processing faecal sludge.

    The following components should be included in system thermo object for simulation:
    N2O.

    The following impact items should be pre-constructed for life cycle assessment:
    OilHeatExchanger, Pump.

    Parameters
    ----------
    ins : stream obj
        Hot gas.
    outs : stream obj
        Hot gas.

    References
    ----------
    #!!! Reference?

    '''
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with, F_BM_default=1)

        self.construction = (
            Construction('oil_heat_exchanger', linked_unit=self, item='OilHeatExchanger', quantity_unit='ea'),
            Construction('pump', linked_unit=self, item='Pump', quantity_unit='ea'),
            )

        data = load_data(path=oilhx_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)

    _N_ins = 1
    _N_outs = 1

    def _run(self):
        heat_in = self.ins[0]
        heat_out = self.outs[0]
        heat_out.copy_like(heat_in)
        heat_out.phase = 'g'

        # Set temperature
        heat_out.T = self.ohx_temp

        #!!! Below aren't used anywhere
        # Calculate the power that was delivered to the ORC
        self.power_delivery_orc = self.oil_flowrate * ((273.15 + self.oil_temp_out) - (273.15 + self.oil_temp_in)) * self.oil_density * self.oil_specific_heat * (60/1000)  # MJ/hr
        # Calculate losses through pipe
        self.pipe_heat_loss = ((273.15 + self.oil_temp_out) - (273.15 + self.amb_temp)) / (self.oil_r_pipe_k_k_w + self.oil_r_insulation_k_k_w) * 3.6  # MJ/hr


    def _design(self):
        design = self.design_results
        constr = self.construction
        design['OilHeatExchanger'] = constr[0].quantity = 4/200
        design['Pump'] = constr[1].quantity = 2.834/2.27
        self.add_construction()


    def _cost(self):
        self.baseline_purchase_costs['Oil Heat Exchanger'] = self.orc_cost * self.price_ratio
        self.power_utility(self.oil_pump_power-self.oil_electrical_energy_generated) # kWh/hr