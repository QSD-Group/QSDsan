#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems
Copyright (C) 2020, Quantitative Sustainable Design Group

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the UIUC open-source license. Please refer to 
https://github.com/QSD-Group/QSDsan/blob/master/LICENSE.txt
for license details.

Ref:
    [1] Trimmer et al., Navigating Multidimensional Social–Ecological System
        Trade-Offs across Sanitation Alternatives in an Urban Informal Settlement.
        Environ. Sci. Technol. 2020, 54 (19), 12641–12653.
        https://doi.org/10.1021/acs.est.0c03296.

TODO:
    [1] Add sysB and sysC, debugging

#!!! Questions:
    [1] WWTP power consumption very low, only 6.5/0.8 kWh per hour for
        existing/alternative WWTP

'''


# %%

import numpy as np
import biosteam as bst
import qsdsan as qs
from sklearn.linear_model import LinearRegression as LR
from qsdsan import sanunits as su
from qsdsan import WasteStream, SanUnit, ImpactItem, StreamImpactItem, \
    SimpleTEA, LCA
from bwaise import cmps

bst.settings.set_thermo(cmps)
items = ImpactItem._items
GWP = qs.ImpactIndicator._indicators['GWP']
currency = qs.currency = 'USD'
qs.CEPCI = qs.CEPCI_by_year[2018]


# %%

# =============================================================================
# 
# =============================================================================

class WWTPCost(SanUnit):
    '''Lumped CAPEX and electricity cost of the wastewater treatment plant.'''
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 lumped_CAPEX=0., lumped_power=0.):
        SanUnit.__init__(self, ID, ins, outs, thermo)
        self.lumped_CAPEX = lumped_CAPEX
        self.lumped_power = lumped_power
    
    def _run(self):
        self.outs[0].copy_like(self.ins[0])
    
    _BM = {'Wastewater treatment plant': 1}
    
    def _cost(self):
        self.purchase_costs['Wastewater treatment plant'] = self.lumped_CAPEX
        self.power_utility(self.lumped_power)



# %%

# =============================================================================
# Unit Parameters
# =============================================================================

toilet_user = 16 # four people per household, four households per toilet
# Number of people served by the existing plant (sysA and SceC), sewer + sludge
ppl_existing = 4e4 + 416667
# Number of people served by the alternative plant (SceB)
ppl_alternative = 5e4

exchange_rate = 3700 # UGX per USD, triangular of 3600, 3700, 3900
get_exchange_rate = lambda: exchange_rate

# Time take for full degradation, [yr]
tau_deg = 2
# Log reduction at full degradation
log_deg = 3
# Get reduction rate constant k for COD and N, use a function so that k can be
# changed during uncertainty analysis
def get_decay_k(tau_deg=2, log_deg=3):
    k = (-1/tau_deg)*np.log(10**-log_deg)
    return k

max_CH4_emission = 0.25

# Model for tanker truck cost based on capacity (m3)
# price = a*capacity**b -> ln(price) = ln(a) + bln(capacity)
UGX_price_dct = np.array((8e4, 12e4, 20e4, 25e4))
capacities = np.array((3, 4.5, 8, 15))
def get_tanker_truck_cost(capacity):
    # Add 15% additional costs
    price_dct = UGX_price_dct*1.15/get_exchange_rate()
    ln_p = np.log(price_dct)
    ln_cap = np.log(capacities)
    model = LR().fit(ln_cap.reshape(-1,1), ln_p.reshape(-1,1))
    [[predicted]] = model.predict(np.array((np.log(capacity))).reshape(1, -1)).tolist()
    cost = np.exp(predicted)
    return cost

# Nutrient loss during applciation
app_loss = dict.fromkeys(('NH3', 'NonNH3', 'P', 'K', 'Mg', 'Ca'), 0.02)
app_loss['NH3'] = 0.05


# =============================================================================
# price_dct and GWP CFs
# =============================================================================

# Recycled nutrients are sold at a lower price than commercial fertilizers
price_factor = 0.25
get_price_factor = lambda: price_factor

price_dct = {
    'Electricity': 0.17,
    'Concrete': 194,
    'Steel': 2.665,
    'N': 1.507*get_price_factor(),
    'P': 3.983*get_price_factor(),
    'K': 1.333*get_price_factor()
    }

GWP_dct = {
    'Electricity': 0.15,
    'CH4': 28,
    'N2O': 265,
    'N': -5.4,
    'P': -4.9,
    'K': -1.5
    }

bst.PowerUtility.price = price_dct['Electricity']
items['Concrete'].price = price_dct['Concrete']
items['Steel'].price = price_dct['Steel']
# This is universal to all scenarios
e_item = ImpactItem(ID='e_item', functional_unit='kWh', GWP=GWP_dct['Electricity'])


# %%

# =============================================================================
# Scenario A: existing system
# =============================================================================

def batch_creating_streams(prefix):
    stream_dct = {}
    stream_dct['CH4'] = WasteStream(f'{prefix}_CH4', phase='g')
    stream_dct['N2O'] = WasteStream(f'{prefix}_N2O', phase='g')
    stream_dct['liq_N'] = WasteStream(f'{prefix}_liq_N', phase='l', price=price_dct['N'])
    stream_dct['sol_N'] = WasteStream(f'{prefix}_sol_N', phase='l', price=price_dct['N'])
    stream_dct['liq_P'] = WasteStream(f'{prefix}_liq_P', phase='l', price=price_dct['P'])
    stream_dct['sol_P'] = WasteStream(f'{prefix}_sol_P', phase='l', price=price_dct['P'])
    stream_dct['liq_K'] = WasteStream(f'{prefix}_liq_K', phase='l', price=price_dct['K'])
    stream_dct['sol_K'] = WasteStream(f'{prefix}_sol_K', phase='l', price=price_dct['K'])

    item_dct = {}
    item_dct['CH4'] = StreamImpactItem(stream_dct['CH4'], GWP=GWP_dct['CH4'])
    item_dct['N2O'] = StreamImpactItem(stream_dct['N2O'], GWP=GWP_dct['N2O'])
    item_dct['liq_N'] = StreamImpactItem(stream_dct['liq_N'], GWP=GWP_dct['N'])
    item_dct['sol_N'] = StreamImpactItem(stream_dct['sol_N'], GWP=GWP_dct['N'])
    item_dct['liq_P'] = StreamImpactItem(stream_dct['liq_P'], GWP=GWP_dct['P'])
    item_dct['sol_P'] = StreamImpactItem(stream_dct['sol_P'], GWP=GWP_dct['P'])
    item_dct['liq_K'] = StreamImpactItem(stream_dct['liq_K'], GWP=GWP_dct['K'])
    item_dct['sol_K'] = StreamImpactItem(stream_dct['sol_K'], GWP=GWP_dct['K'])
    return stream_dct, item_dct

streamsA, itemsA = batch_creating_streams('A')


#################### Human Inputs ####################
A1 = su.Excretion('A1', outs=('urine', 'feces'))
def refresh_sysA_streams():
    A_streams, A_items = batch_creating_streams('A')
    A1._run()
A1.specification = refresh_sysA_streams

################### User Interface ###################
A2 = su.PitLatrine('A2', ins=(A1-0, A1-1,
                              'toilet_paper', 'flushing_water',
                              'cleansing_water', 'desiccant'),
                      outs=('mixed_waste', 'leachate', '', ''),
                      N_user=toilet_user, N_toilet=ppl_existing/toilet_user,
                      OPEX_over_CAPEX=0.05,
                      decay_k_COD=get_decay_k(tau_deg, log_deg),
                      decay_k_N=get_decay_k(tau_deg, log_deg),
                      max_CH4_emission=max_CH4_emission)

##################### Conveyance #####################
A3 = su.Trucking('A3', ins=A2-0, outs=('transported', 'conveyance_loss'),
                 load_type='mass',
                 distance=5, distance_unit='km',
                 interval=A2.emptying_period, interval_unit='yr',
                 loss_ratio=0.02)
def update_A3_param():
    A3._run()
    truck = A3.single_truck
    truck.load = A3.F_mass_in*A2.emptying_period*365*24/A2.N_toilet
    vol = truck.load/1e3 # Assume the density of water
    A3.fee = get_tanker_truck_cost(vol)
A3.specification = update_A3_param

###################### Treatment ######################
A4 = WWTPCost('A4', ins=A3-0, lumped_CAPEX=18606700, lumped_power=57120/(365*24))

A5 = su.SedimentationTank('A5', ins=A4-0,
                          outs=('liq', 'sol', '', ''),
                          decay_k_COD=get_decay_k(tau_deg, log_deg),
                          decay_k_N=get_decay_k(tau_deg, log_deg),
                          max_CH4_emission=max_CH4_emission)
def A5_cost():
    A5.purchase_costs.clear()
A5._cost = A5_cost

A6 = su.Lagoon('A6', ins=A5-0, outs=('anaerobic_treated', '', ''),
                  design_type='anaerobic',
                  decay_k_N=get_decay_k(tau_deg, log_deg),
                  max_CH4_emission=max_CH4_emission)

A7 = su.Lagoon('A7', ins=A6-0, outs=('facultative_treated', '', ''),
                  design_type='facultative',
                  decay_k_N=get_decay_k(tau_deg, log_deg),
                  max_CH4_emission=max_CH4_emission)

A8 = su.DryingBed('A8', ins=A5-1, outs=('dried_sludge', 'evaporated', '', ''),
                     design_type='unplanted',
                     decay_k_COD=get_decay_k(tau_deg, log_deg),
                     decay_k_N=get_decay_k(tau_deg, log_deg),
                     max_CH4_emission=max_CH4_emission)

treatA = bst.System('treatA', path=(A4, A5, A6, A7, A8))

################## Reuse or Disposal ##################
A9 = su.CropApplication('A8', ins=A7-0, outs=('liquid_fertilizer', 'reuse_loss'),
                        loss_ratio=app_loss)
def adjust_NH3_loss():
    A9._run()
    # Assume the slight higher loss of NH3 does not affect COD,
    # does not matter much since COD not considered in crop application
    A9.outs[0]._COD = A9.outs[1]._COD = A9.ins[0]._COD
A9.specification = adjust_NH3_loss

A10 = su.Mixer('A10', ins=(A2-2, A5-2, A6-1, A7-1, A8-2), outs=streamsA['CH4'])
A10.line = 'CH4 mixer'

A11 = su.Mixer('A11', ins=(A2-3, A5-3, A6-2, A7-2, A8-3), outs=streamsA['N2O'])
A11.line = 'N2O mixer'

A12 = su.ComponentSplitter('A12', ins=A8-0,
                           outs=(streamsA['sol_N'], streamsA['sol_P'], streamsA['sol_K'],
                                 'A_sol_non_fertilizers'),
                           split_keys=(('NH3', 'NonNH3'), 'P', 'K'))

A13 = su.ComponentSplitter('A13', ins=A9-0,
                           outs=(streamsA['liq_N'], streamsA['liq_P'], streamsA['liq_K'],
                                 'A_liq_non_fertilizers'),
                           split_keys=(('NH3', 'NonNH3'), 'P', 'K'))

############### Simulation, TEA, and LCA ###############
sysA = bst.System('sysA', path=(A1, A2, A3, treatA, A9, A10, A11, A12, A13))
sysA.simulate()
# sysA.save_report('results/sysA.xlsx')

teaA = SimpleTEA(system=sysA, discount_rate=0.05, start_year=2018,
                 lifetime=8, uptime_ratio=1, lang_factor=None,
                 annual_maintenance=0, annual_labor=12*3e6*12/get_exchange_rate(),
                 construction_schedule=None)

# 57120 is the annual electricity usage for the whole treatment plant
get_annual_e = lambda: A4.power_utility.rate*teaA._operating_hours

lcaA = LCA(system=sysA, lifetime=8, lifetime_unit='yr', uptime_ratio=1,
           # Assuming all additional WWTP OPEX from electricity
           e_item=get_annual_e()*8)


# %%

# =============================================================================
# Scenario B: anaerobic treatment with existing pit latrines and conveyance
# =============================================================================

# streamsB, itemsB = batch_creating_streams('B')

# ws1 = A2.outs[0].copy('ws1')
# A4 = su.AnaerobicDigestion('A4', ins=ws1,
#                               outs=('treated', 'CH4', 'N2O'),
#                               # tau_previous=A2.emptying_period*365,
#                               decay_k_N=get_decay_k(tau_deg, log_deg),
#                                 max_CH4_emission=max_CH4_emission)
# A4.simulate()
# # A4.show()

# ws2 = A2.outs[0].copy('ws2')
# A5 = su.SludgeSeparator('A5', ins=ws2, outs=('liq', 'sol'))

# A5.simulate()
# # A5.show()




# ws5 = A2.outs[0].copy('ws5')
# A8 = su.DryingBed('A8', ins=ws5, outs=('solid', 'evaporated', 'CH4', 'N2O'),
#                       design_type='unplanted',
#                       decay_k_COD=get_decay_k(tau_deg, log_deg),
#                       decay_k_N=get_decay_k(tau_deg, log_deg),
#                       max_CH4_emission=max_CH4_emission)

# A8.simulate()
# # A8.show()

# ws6 = A2.outs[0].copy('ws6')
# A9 = su.AnaerobicBaffledReactor('A9', ins=ws6, outs=('treated', 'CH4', 'N2O'),
#                                     decay_k_COD=get_decay_k(tau_deg, log_deg),
#                                     max_CH4_emission=max_CH4_emission)

# A9.simulate()
# # A9.show()

# ws7 = A2.outs[0].copy('ws7')
# A10 = su.LiquidTreatmentBed('A10', ins=ws7, outs=('treated', 'CH4', 'N2O'),
#                                 decay_k_COD=get_decay_k(tau_deg, log_deg),
#                                 decay_k_N=get_decay_k(tau_deg, log_deg),
#                                 max_CH4_emission=max_CH4_emission)

# A10.simulate()
# A10.show()

# biogas = WasteStream('biogas', CH4=1)
# AX2 = su.BiogasCombustion('AX2', ins=(biogas, 'air'),
#                               outs=('used', 'lost', 'wasted'),
#                               if_combustion=True,
#                               biogas_loss=0.1, biogas_eff=0.55)



# sysA = bst.System('sysA', path=(A1, A2, A3, AX, AX2))

# sysA.simulate()



# lcaA = LCA(sysA, lifetime=8, lifetime_unit='yr',
#                 e_item=(5000, 'Wh'))









# %%

# =============================================================================
# Summarizing Functions
# =============================================================================

def get_total_inputs(unit):
    if len(unit.ins) == 0: # Excretion units do not have ins
        ins = unit.outs
    else:
        ins = unit.ins
    inputs = {}
    inputs['COD'] = sum(i.COD*i.F_vol/1e3 for i in ins)
    inputs['energy'] = inputs['COD'] * 14e3
    inputs['N'] = sum(i.TN*i.F_vol/1e3 for i in ins)
    inputs['NH3'] = sum(i.imass['NH3'] for i in ins)
    inputs['P'] = sum(i.TP*i.F_vol/1e3 for i in ins)
    inputs['K'] = sum(i.TK*i.F_vol/1e3 for i in ins)
    inputs['Mg'] = sum(i.TMg*i.F_vol/1e3 for i in ins)
    inputs['Ca'] = sum(i.TCa*i.F_vol/1e3 for i in ins)
    for i, j in inputs.items():
        inputs[i] = j*365*24
    return inputs

def get_recovery(ins=None, outs=None, hours=365*24, ppl=1, if_relative=True):
    try: iter(outs)
    except: outs = (outs,)
    non_g = tuple(i for i in outs if i.phase != 'g')
    recovery = {}
    recovery['COD'] = sum(i.COD*i.F_vol/1e3 for i in non_g)
    recovery['energy'] = recovery['COD'] * 14e3
    recovery['N'] = sum(i.TN*i.F_vol/1e3 for i in non_g)
    recovery['NH3'] = sum(i.imass['NH3'] for i in non_g)
    recovery['P'] = sum(i.TP*i.F_vol/1e3 for i in non_g)
    recovery['K'] = sum(i.TK*i.F_vol/1e3 for i in non_g)
    recovery['Mg'] = sum(i.TMg*i.F_vol/1e3 for i in non_g)
    recovery['Ca'] = sum(i.TCa*i.F_vol/1e3 for i in non_g)
    for i, j in recovery.items():
        if if_relative:
            inputs = get_total_inputs(ins)
            recovery[i] /= inputs[i]/hours * ppl
        else:
            recovery[i] /= 1/hours * ppl
    return recovery

def get_direct_emissions(outs=None, hours=365*24, ppl=1):
    try: iter(outs)
    except: outs = (outs,)
    gas = tuple(i for i in outs if i.phase == 'g')
    emission = {}
    emission['direct'] = \
        sum((i.imass['CH4', 'N2O']*(GWP_dct['CH4'], GWP_dct['N2O'])).sum()
            for i in gas)*hours/ppl
    return emission

sys_dct = {
    #!!! UPDATE VALUES AFTER CREATING SYSB AND SYSC
    'ppl': dict(sysA=ppl_existing, sysB=ppl_alternative, sysC=ppl_existing),
    'input_unit': dict(sysA=A1, sysB=A1, sysC=A1),
    'liq_unit': dict(sysA=A13, sysB=A13, sysC=A13),
    'sol_unit': dict(sysA=A12, sysB=A12, sysC=A12),
    'stream_dct': dict(sysA=streamsA, sysB=streamsA, sysC=streamsA),
    'TEA': dict(sysA=teaA, sysB=teaA, sysC=teaA),
    'LCA': dict(sysA=lcaA, sysB=lcaA, sysC=lcaA),
    }

def print_summaries(sys):
    if isinstance(sys, bst.System):
        sys = sys.ID
    ppl = sys_dct['ppl'][sys]
    print(f'\n---------- Summary for {sys} ----------\n')
    tea = sys_dct['TEA'][sys]
    tea.show()
    print('\n')
    lca = sys_dct['LCA'][sys]
    lca.show()
    print('\n')

    print(f'Net cost is {tea.EAC/ppl:.1f} {currency}/cap/yr.')
    print(f'Construction cost is {tea.annualized_CAPEX/ppl:.1f} {currency}/cap/yr.')
    print(f'Operating cost is {tea.AOC/ppl:.1f} {currency}/cap/yr.')
    print('\n')
    
    ind = 'GlobalWarming'
    factor = lca.lifetime * ppl
    total = lca.total_impacts[ind]/factor
    constr = lca.total_construction_impacts[ind]/factor
    trans = lca.total_transportation_impacts[ind]/factor
    ws = lca.get_stream_impacts(stream_items=lca.stream_inventory)[ind]/factor
    other = lca.total_other_impacts[ind]/factor
    print(f'Net emission is {total:.1f} {GWP.unit}/cap/yr.')
    print(f'Construction emission is {constr:.1f} {GWP.unit}/cap/yr.')
    print(f'Transportation emission is {trans:.1f} {GWP.unit}/cap/yr.')
    print(f'Stream emission is {ws:.1f} {GWP.unit}/cap/yr.')
    print(f'Other emission is {other:.1} {GWP.unit}/cap/yr.')
    print('\n')

    input_unit = sys_dct['input_unit'][sys]
    liq_unit = sys_dct['liq_unit'][sys]
    sol_unit = sys_dct['sol_unit'][sys]
    liq = get_recovery(ins=input_unit, outs=liq_unit.outs[:3], ppl=ppl)
    sol = get_recovery(ins=input_unit, outs=sol_unit.outs[:3], ppl=ppl)
    for i in ('N', 'P', 'K'):  
        print(f'Total {i} recovery is {liq[i]+sol[i]:.1%}, '
              f'{liq[i]:.1%} in liquid, '
              f'{sol[i]:.1%} in solid.')
    liq_COD = get_recovery(ins=input_unit, outs=liq_unit.ins, ppl=ppl)['COD']
    sol_COD = get_recovery(ins=input_unit, outs=sol_unit.ins, ppl=ppl)['COD']
    print(f'Total COD recovery is {liq_COD+sol_COD:.1%}, '
          f'{liq_COD:.1%} in liquid, '
          f'{sol_COD:.1%} in solid.')



# %%

# # =============================================================================
# # Scenario C: containaer-based sanitation with existing treatment
# # =============================================================================

# print('\n----------Scenario C----------\n')
# C1 = su.Excretion('C1', outs=('urine', 'feces'), N_user=toilet_user)
# C2 = su.UDDT('C2', ins=(C1-0, C1-1,
#                             'toilet_paper', 'flushing_water',
#                             'cleaning_water', 'desiccant'),
#                 outs=('liquid_waste', 'solid_waste',
#                       'struvite', 'HAP', 'CH4', 'N2O'),
#                 N_user=toilet_user, OPEX_over_CAPEX=0.1,
#                 decay_k_COD=get_decay_k(tau_deg, log_deg),
#                 decay_k_N=get_decay_k(tau_deg, log_deg),
#                 max_CH4_emission=max_CH4_emission)
# C1.simulate()
# C2.simulate()

# truck_cost = {
#     'TankerTruck1': 8e4/get_exchange_rate()*1.15, # 15% additional, per m3
#     'TankerTruck2': 12e4/get_exchange_rate()*1.15, # 15% additional, per m3
#     'TankerTruck3': 2e5/get_exchange_rate()*1.15, # 15% additional, per m3
#     'TankerTruck4': 25e4/get_exchange_rate()*1.15, # 15% additional, per m3
#     'HandCart': 0.01, # per cap/d
#     'CBSTruck': 23e3/get_exchange_rate() # per m3
#     }

# # Assume density is 1 tonne/m3 (as water)
# V = (3, 4.5, 8, 15, 1, 1)
# truck_V = dict.fromkeys(truck_cost.keys())
# for i, j in zip(truck_V.keys(), V):
#     truck_V[i] = j

# truck = 'HandcartAndTruck'
# # Liquid waste
# interval = (C2.collection_period*truck_V[truck])/C2.tank_V
# C3 = su.Trucking('C3', ins=C2-0, outs=('transported_l', 'loss_l'),
#                     load_type='mass', load=truck_V[truck], load_unit='tonne',
#                     distance=5, distance_unit='km',
#                     interval=interval, interval_unit='day',
#                     fee=truck_cost[truck]/truck_V[truck]+0.01*ppl_alternative*interval,
#                     loss_ratio=0.02)

# # Solid waste
# interval = (C2.collection_period*truck_V[truck])/C2.tank_V
# C4 = su.Trucking('C4', ins=C2-1, outs=('transported_s', 'loss_s'),
#                     load_type='mass', load=truck_V[truck], load_unit='tonne',
#                     distance=5, distance_unit='km',
#                     interval=interval, interval_unit='day',
#                     fee=truck_cost[truck]/truck_V[truck]+0.01*ppl_alternative*interval,
#                     loss_ratio=0.02)


# CX = su.CropApplication('CX', ins=WasteStream(), loss_ratio=app_loss)
# def adjust_NH3_loss():
#     CX._run()
#     CX.outs[0]._COD = CX.outs[1]._COD = CX.ins[0]._COD
# CX.specification = adjust_NH3_loss


# SceC = bst.System('SceC', path=(C1, C2, C3, C4, CX))

# SceC.simulate()

