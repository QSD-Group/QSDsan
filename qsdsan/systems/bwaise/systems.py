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

TODOs:
    [1] Recheck COD calculation acorss units once WasteStream is ready

Questions:
    [1] WWTP power consumption very low, only 6.5/0.8 kWh per hour for
        existing/alternative WWTP

'''


# %%

import numpy as np
import biosteam as bst
import qsdsan as qs
from sklearn.linear_model import LinearRegression as LR
from qsdsan import sanunits as su
from qsdsan import WasteStream, ImpactItem, StreamImpactItem, SimpleTEA, LCA
from qsdsan.systems.bwaise._cmps import cmps


# =============================================================================
# Unit parameters
# =============================================================================

bst.settings.set_thermo(cmps)
items = ImpactItem._items
GWP = qs.ImpactIndicator._indicators['GWP']
currency = qs.currency = 'USD'
qs.CEPCI = qs.CEPCI_by_year[2018]
bst.speed_up()

toilet_user = 16 # four people per household, four households per toilet
# Number of people served by the existing plant (sysA and SceC), sewer + sludge
ppl_existing = 4e4 + 416667
# Number of people served by the alternative plant (SceB)
ppl_alternative = 5e4

exchange_rate = 3700 # UGX per USD, triangular of 3600, 3700, 3900
get_exchange_rate = lambda: exchange_rate
discount_rate = 0.05
get_discount_rate = lambda: discount_rate

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
def get_tanker_truck_fee(capacity):
    # Add 15% additional
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
# Prices and GWP CFs
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
    'K': 1.333*get_price_factor(),
    'Biogas': 6500/get_exchange_rate()
    }

GWP_dct = {
    'Electricity': 0.15,
    'CH4': 28,
    'N2O': 265,
    'N': -5.4,
    'P': -4.9,
    'K': -1.5,
    'Biogas': -3
    }

bst.PowerUtility.price = price_dct['Electricity']
items['Concrete'].price = price_dct['Concrete']
items['Steel'].price = price_dct['Steel']
# This is universal to all scenarios
e_item = ImpactItem(ID='e_item', functional_unit='kWh', GWP=GWP_dct['Electricity'])

# =============================================================================
# Universal units and functions
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


# Costs of WWTP units have been considered in the lumped unit
def clear_unit_costs(sys):
    for i in sys.units:
        if isinstance(i, su.LumpedCost):
            continue
        i.purchase_costs.clear()


def adjust_NH3_loss(unit):
    unit._run()
    # Assume the slight higher loss of NH3 does not affect COD,
    # does not matter much since COD not considered in crop application
    unit.outs[0]._COD = unit.outs[1]._COD = unit.ins[0]._COD


# %%

# =============================================================================
# Scenario A (sysA): pit latrine with existing treatment system
# =============================================================================

flowsheetA = bst.Flowsheet('sysA')
bst.main_flowsheet.set_flowsheet(flowsheetA)
streamsA, itemsA = batch_creating_streams('A')

#################### Human Inputs ####################
A1 = su.Excretion('A1', outs=('urine', 'feces'))
def refresh_sysA_streams():
    streamsA, itemsA = batch_creating_streams('A')
    A1._run()
A1.specification = refresh_sysA_streams

################### User Interface ###################
A2 = su.PitLatrine('A2', ins=(A1-0, A1-1,
                              'toilet_paper', 'flushing_water',
                              'cleansing_water', 'desiccant'),
                   outs=('mixed_waste', 'leachate', 'A2_CH4', 'A2_N2O'),
                   N_user=toilet_user, N_toilet=ppl_existing/toilet_user,
                   OPEX_over_CAPEX=0.05,
                   decay_k_COD=get_decay_k(tau_deg, log_deg),
                   decay_k_N=get_decay_k(tau_deg, log_deg),
                   max_CH4_emission=max_CH4_emission)

##################### Conveyance #####################
A3 = su.Trucking('A3', ins=A2-0, outs=('transported', 'conveyance_loss'),
                 load_type='mass', distance=5, distance_unit='km',
                 interval=A2.emptying_period, interval_unit='yr',
                 loss_ratio=0.02)
def update_A3_param():
    A3._run()
    truck = A3.single_truck
    truck.interval = A2.emptying_period*365*24
    truck.load = A3.F_mass_in*truck.interval/A2.N_toilet
    rho = A3.F_mass_in/A3.F_vol_in
    vol = truck.load/rho
    A3.fee = get_tanker_truck_fee(vol)
    A3._design()
A3.specification = update_A3_param

###################### Treatment ######################
A4 = su.LumpedCost('A4', ins=A3-0, cost_item_name='Lumped WWTP',
                   CAPEX=18606700, power=57120/(365*24), lifetime=8)
A4.line = 'Lumped WWTP cost'

A5 = su.SedimentationTank('A5', ins=A4-0,
                          outs=('liq', 'sol', 'A5_CH4', 'A5_N2O'),
                          decay_k_COD=get_decay_k(tau_deg, log_deg),
                          decay_k_N=get_decay_k(tau_deg, log_deg),
                          max_CH4_emission=max_CH4_emission)

A6 = su.Lagoon('A6', ins=A5-0, outs=('anaerobic_treated', 'A6_CH4', 'A6_N2O'),
               design_type='anaerobic',
               decay_k_N=get_decay_k(tau_deg, log_deg),
               max_CH4_emission=max_CH4_emission)

A7 = su.Lagoon('A7', ins=A6-0, outs=('facultative_treated', 'A7_CH4', 'A7_N2O'),
               design_type='facultative',
               decay_k_N=get_decay_k(tau_deg, log_deg),
               max_CH4_emission=max_CH4_emission,
               if_N2O_emission=True)

A8 = su.DryingBed('A8', ins=A5-1, outs=('dried_sludge', 'evaporated',
                                        'A8_CH4', 'A8_N2O'),
                  design_type='unplanted',
                  decay_k_COD=get_decay_k(tau_deg, log_deg),
                  decay_k_N=get_decay_k(tau_deg, log_deg),
                  max_CH4_emission=max_CH4_emission)
treatA = bst.System('treatA', path=(A4, A5, A6, A7, A8))
A8._cost = lambda: clear_unit_costs(treatA)

################## Reuse or Disposal ##################
A9 = su.CropApplication('A9', ins=A7-0, outs=('liquid_fertilizer', 'reuse_loss'),
                        loss_ratio=app_loss)
A9.specification = lambda: adjust_NH3_loss(A9)

A10 = su.Mixer('A10', ins=(A2-2, A5-2, A6-1, A7-1, A8-2), outs=streamsA['CH4'])
A10.line = 'fugative CH4 mixer'

A11 = su.Mixer('A11', ins=(A2-3, A5-3, A6-2, A7-2, A8-3), outs=streamsA['N2O'])
A11.line = 'fugative N2O mixer'

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

teaA = SimpleTEA(system=sysA, discount_rate=get_discount_rate(), start_year=2018,
                 lifetime=8, uptime_ratio=1, lang_factor=None,
                 annual_maintenance=0, annual_labor=12*3e6*12/get_exchange_rate(),
                 construction_schedule=None)

lcaA = LCA(system=sysA, lifetime=8, lifetime_unit='yr', uptime_ratio=1,
           # Assuming all additional WWTP OPEX from electricity
           e_item=lambda: A4.power_utility.rate*(365*24)*8)


# %%

# =============================================================================
# Scenario B (sysB): pit latrine with anaerobic treatment
# =============================================================================

flowsheetB = bst.Flowsheet('sysB')
bst.main_flowsheet.set_flowsheet(flowsheetB)
streamsB, itemsB = batch_creating_streams('B')
streamsB['biogas'] = WasteStream('B_biogas', phase='g', price=price_dct['Biogas'])
items['biogas'] = StreamImpactItem(streamsB['biogas'], GWP=GWP_dct['Biogas'])

#################### Human Inputs ####################
B1 = su.Excretion('B1', outs=('urine', 'feces'))
def refresh_sysB_streams():
    streamsB, itemsB = batch_creating_streams('B')
    streamsB['biogas'] = WasteStream('B_biogas', phase='g', price=price_dct['Biogas'])
    items['biogas'] = StreamImpactItem(streamsB['biogas'], GWP=GWP_dct['Biogas'])
    B1._run()
B1.specification = refresh_sysB_streams

################### User Interface ###################
B2 = su.PitLatrine('B2', ins=(B1-0, B1-1,
                              'toilet_paper', 'flushing_water',
                              'cleansing_water', 'desiccant'),
                   outs=('mixed_waste', 'leachate', 'B2_CH4', 'B2_N2O'),
                   N_user=toilet_user, N_toilet=ppl_alternative/toilet_user,
                   OPEX_over_CAPEX=0.05,
                   decay_k_COD=get_decay_k(tau_deg, log_deg),
                   decay_k_N=get_decay_k(tau_deg, log_deg),
                   max_CH4_emission=max_CH4_emission)

##################### Conveyance #####################
B3 = su.Trucking('B3', ins=B2-0, outs=('transported', 'conveyance_loss'),
                 load_type='mass', distance=5, distance_unit='km',
                 interval=B2.emptying_period, interval_unit='yr',
                 loss_ratio=0.02)
def update_B3_param():
    B3._run()
    truck = B3.single_truck
    truck.interval = B2.emptying_period*365*24
    truck.load = B3.F_mass_in*truck.interval/B2.N_toilet
    rho = B3.F_mass_in/B3.F_vol_in
    vol = truck.load/rho
    B3.fee = get_tanker_truck_fee(vol)
    B3._design()
B3.specification = update_B3_param

###################### Treatment ######################
B4 = su.LumpedCost('B4', ins=B3-0, cost_item_name='Lumped WWTP',
                   CAPEX=337140, power=6854/(365*24), lifetime=10)
B4.line = 'Lumped WWTP cost'

B5 = su.AnaerobicBaffledReactor('B5', ins=B4-0, outs=('treated', 'biogas',
                                                      'B5_CH4', 'B5_N2O'),
                                decay_k_COD=get_decay_k(tau_deg, log_deg),
                                max_CH4_emission=max_CH4_emission)

B6 = su.SludgeSeparator('B6', ins=B5-0, outs=('liq', 'sol'))

B7 = su.LiquidTreatmentBed('B7', ins=B6-0, outs=('treated', 'B7_CH4', 'B7_N2O'),
                           decay_k_COD=get_decay_k(tau_deg, log_deg),
                           decay_k_N=get_decay_k(tau_deg, log_deg),
                           max_CH4_emission=max_CH4_emission)

B8 = su.DryingBed('B8', ins=B6-1, outs=('dried_sludge', 'evaporated',
                                        'B8_CH4', 'B8_N2O'),
                  design_type='planted',
                  decay_k_COD=get_decay_k(tau_deg, log_deg),
                  decay_k_N=get_decay_k(tau_deg, log_deg),
                  max_CH4_emission=max_CH4_emission)

treatB = bst.System('treatB', path=(B4, B5, B6, B7, B8))
B8._cost = lambda: clear_unit_costs(treatB)

################## Reuse or Disposal ##################
B9 = su.CropApplication('B9', ins=B7-0, outs=('liquid_fertilizer', 'reuse_loss'),
                        loss_ratio=app_loss)
B9.specification = lambda: adjust_NH3_loss(B9)

B10 = su.Mixer('B10', ins=(B2-2, B5-2, B7-1, B8-2), outs=streamsB['CH4'])
B10.line = 'fugative CH4 mixer'

B11 = su.Mixer('B11', ins=(B2-3, B5-3, B7-2, B8-3), outs=streamsB['N2O'])
B11.line = 'fugative N2O mixer'

B12 = su.ComponentSplitter('B12', ins=B8-0,
                            outs=(streamsB['sol_N'], streamsB['sol_P'], streamsB['sol_K'],
                                  'B_sol_non_fertilizers'),
                            split_keys=(('NH3', 'NonNH3'), 'P', 'K'))

B13 = su.ComponentSplitter('B13', ins=B9-0,
                            outs=(streamsB['liq_N'], streamsB['liq_P'], streamsB['liq_K'],
                                  'B_liq_non_fertilizers'),
                            split_keys=(('NH3', 'NonNH3'), 'P', 'K'))

B14 = su.BiogasCombustion('B14', ins=(B5-1, 'air'),
                          outs=('used', 'lost', 'wasted'),
                          if_combustion=False,
                          biogas_loss=0.1, biogas_eff=0.55)
B15 = su.Mixer('B15', ins=(B14-0, B14-2), outs=streamsB['biogas'])

############### Simulation, TEA, and LCA ###############
sysB = bst.System('sysB', path=(B1, B2, B3, treatB, B9, B10, B11, B12, B13, B14, B15))

teaB = SimpleTEA(system=sysB, discount_rate=get_discount_rate(), start_year=2018,
                  lifetime=10, uptime_ratio=1, lang_factor=None,
                  annual_maintenance=0,
                  annual_labor=(5*5e6+5*75e4)*12/get_exchange_rate(),
                  construction_schedule=None)

lcaB = LCA(system=sysB, lifetime=10, lifetime_unit='yr', uptime_ratio=1,
           # Assuming all additional WWTP OPEX from electricity
           e_item=lambda: B4.power_utility.rate*(365*24)*10)


# %%

# =============================================================================
# Scenario C (sysC): containaer-based sanitation with existing treatment system
# =============================================================================

flowsheetC = bst.Flowsheet('sysC')
bst.main_flowsheet.set_flowsheet(flowsheetC)
streamsC, itemsC = batch_creating_streams('C')

#################### Human Inputs ####################
C1 = su.Excretion('C1', outs=('urine', 'feces'))
def refresh_sysC_streams():
    streamsC, itemsC = batch_creating_streams('C')
    C1._run()
C1.specification = refresh_sysC_streams

################### User Interface ###################
C2 = su.UDDT('C2', ins=(C1-0, C1-1,
                        'toilet_paper', 'flushing_water',
                        'cleaning_water', 'desiccant'),
             outs=('liq_waste', 'sol_waste',
                   'struvite', 'HAP', 'C2_CH4', 'C2_N2O'),
             N_user=toilet_user, N_toilet=ppl_existing/toilet_user,
             OPEX_over_CAPEX=0.1,
             decay_k_COD=get_decay_k(tau_deg, log_deg),
             decay_k_N=get_decay_k(tau_deg, log_deg),
             max_CH4_emission=max_CH4_emission)

##################### Conveyance #####################
# Liquid waste
get_handtruck_and_truck_fee = \
    lambda vol, ppl: 23e3/get_exchange_rate()*vol+0.01*ppl*C2.collection_period
C3 = su.Trucking('C3', ins=C2-0, outs=('transported_l', 'loss_l'),
                 load_type='mass', distance=5, distance_unit='km',
                 loss_ratio=0.02)

# Solid waste
C4 = su.Trucking('C4', ins=C2-1, outs=('transported_s', 'loss_s'),
                 load_type='mass', load=1, load_unit='tonne',
                 distance=5, distance_unit='km',
                 loss_ratio=0.02)
def update_C3_C4_param():
    C4._run()
    truck3, truck4 = C3.single_truck, C4.single_truck
    hr = truck3.interval = truck4.interval = C2.collection_period*24
    ppl = ppl_existing / C2.N_toilet
    truck3.load = C3.F_mass_in * hr / C2.N_toilet
    truck4.load = C4.F_mass_in * hr / C2.N_toilet
    rho3 = C3.F_mass_in/C3.F_vol_in
    rho4 = C4.F_mass_in/C4.F_vol_in
    C3.fee = get_handtruck_and_truck_fee(truck3.load/rho3, ppl)
    C4.fee = get_handtruck_and_truck_fee(truck4.load/rho4, ppl)
    C3._design()
    C4._design()
C4.specification = update_C3_C4_param

###################### Treatment ######################
C5 = su.LumpedCost('C5', ins=(C3-0, C4-0),
                   cost_item_name='Lumped WWTP',
                   CAPEX=18606700, power=57120/(365*24), lifetime=8)

C6 = su.Lagoon('C6', ins=C5-0, outs=('anaerobic_treated', 'C6_CH4', 'C6_N2O'),
               design_type='anaerobic',
               decay_k_N=get_decay_k(tau_deg, log_deg),
               max_CH4_emission=max_CH4_emission)

C7 = su.Lagoon('C7', ins=C6-0, outs=('facultative_treated', 'C7_CH4', 'C7_N2O'),
               design_type='facultative',
               decay_k_N=get_decay_k(tau_deg, log_deg),
               max_CH4_emission=max_CH4_emission,
               if_N2O_emission=True)

C8 = su.DryingBed('C8', ins=C5-1, outs=('dried_sludge', 'evaporated', 'C8_CH4', 'C8_N2O'),
                 design_type='unplanted',
                 decay_k_COD=get_decay_k(tau_deg, log_deg),
                 decay_k_N=get_decay_k(tau_deg, log_deg),
                 max_CH4_emission=max_CH4_emission)
treatC = bst.System('treatC', path=(C5, C6, C7, C8))
C8._cost = lambda: clear_unit_costs(treatC)

################## Reuse or Disposal ##################
C9 = su.CropApplication('C9', ins=C7-0, outs=('liquid_fertilizer', 'reuse_loss'),
                        loss_ratio=app_loss)
C9.specification = lambda: adjust_NH3_loss(C9)

C10 = su.Mixer('C10', ins=(C2-4, C6-1, C7-1, C8-2), outs=streamsC['CH4'])
C10.line = 'fugative CH4 mixer'

C11 = su.Mixer('C11', ins=(C2-5, C6-2, C7-2, C8-3), outs=streamsC['N2O'])
C11.line = 'fugative N2O mixer'

C12 = su.ComponentSplitter('C12', ins=C8-0,
                           outs=(streamsC['sol_N'], streamsC['sol_P'], streamsC['sol_K'],
                                 'C_sol_non_fertilizers'),
                           split_keys=(('NH3', 'NonNH3'), 'P', 'K'))

C13 = su.ComponentSplitter('C13', ins=C9-0,
                           outs=(streamsC['liq_N'], streamsC['liq_P'], streamsC['liq_K'],
                                 'C_liq_non_fertilizers'),
                           split_keys=(('NH3', 'NonNH3'), 'P', 'K'))

############### Simulation, TEA, and LCA ###############
sysC = bst.System('sysC', path=(C1, C2, C3, C4, treatC, C9, C10, C11, C12, C13))

teaC = SimpleTEA(system=sysC, discount_rate=get_discount_rate(), start_year=2018,
                 lifetime=8, uptime_ratio=1, lang_factor=None,
                 annual_maintenance=0, annual_labor=12*3e6*12/get_exchange_rate(),
                 construction_schedule=None)

lcaC = LCA(system=sysC, lifetime=8, lifetime_unit='yr', uptime_ratio=1,
           # Assuming all additional WWTP OPEX from electricity
           e_item=lambda: C5.power_utility.rate*(365*24)*8)


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
    inputs['N'] = sum(i.TN*i.F_vol/1e3 for i in ins)
    inputs['NH3'] = sum(i.imass['NH3'] for i in ins)
    inputs['P'] = sum(i.TP*i.F_vol/1e3 for i in ins)
    inputs['K'] = sum(i.TK*i.F_vol/1e3 for i in ins)
    hr = 365 * 24
    for i, j in inputs.items():
        inputs[i] = j * hr
    return inputs

def get_recovery(ins=None, outs=None, hr=365*24, ppl=1, if_relative=True):
    try: iter(outs)
    except: outs = (outs,)
    non_g = tuple(i for i in outs if i.phase != 'g')
    recovery = {}
    recovery['COD'] = sum(i.COD*i.F_vol/1e3 for i in non_g)
    recovery['N'] = sum(i.TN*i.F_vol/1e3 for i in non_g)
    recovery['NH3'] = sum(i.imass['NH3'] for i in non_g)
    recovery['P'] = sum(i.TP*i.F_vol/1e3 for i in non_g)
    recovery['K'] = sum(i.TK*i.F_vol/1e3 for i in non_g)
    for i, j in recovery.items():
        if if_relative:
            inputs = get_total_inputs(ins)
            recovery[i] /= inputs[i]/hr * ppl
        else:
            recovery[i] /= 1/hr * ppl
    return recovery

def get_stream_emissions(streams=None, hr=365*24, ppl=1):
    try: iter(streams)
    except: streams = (streams,)
    emission = {}
    factor = hr / ppl
    for i in streams:
        if not i.impact_item: continue
        emission[f'{i.ID}'] = i.F_mass*i.impact_item.CFs['GlobalWarming']*factor
    return emission

sys_dct = {
    'ppl': dict(sysA=ppl_existing, sysB=ppl_alternative, sysC=ppl_existing),
    'input_unit': dict(sysA=A1, sysB=B1, sysC=C1),
    'liq_unit': dict(sysA=A13, sysB=B13, sysC=C13),
    'sol_unit': dict(sysA=A12, sysB=B12, sysC=C12),
    'gas_unit': dict(sysA=None, sysB=B15, sysC=None),
    'stream_dct': dict(sysA=streamsA, sysB=streamsB, sysC=streamsC),
    'TEA': dict(sysA=teaA, sysB=teaB, sysC=teaC),
    'LCA': dict(sysA=lcaA, sysB=lcaB, sysC=lcaC),
    'cache': dict(sysA={}, sysB={}, sysC={}),
    }

def cache_recoveries(sys):
    total_COD = get_total_inputs(sys_dct['input_unit'][sys.ID])['COD']
    ppl = sys_dct['ppl'][sys.ID]
    if sys_dct['gas_unit'][sys.ID]:
        gas_COD = -sys_dct['gas_unit'][sys.ID].outs[0].HHV*365*24/14e3/ppl/total_COD
    else:
        gas_COD = 0
    cache = {
        'liq': get_recovery(ins=sys_dct['input_unit'][sys.ID],
                            outs=sys_dct['liq_unit'][sys.ID].ins,
                            ppl=ppl),
        'sol': get_recovery(ins=sys_dct['input_unit'][sys.ID],
                            outs=sys_dct['sol_unit'][sys.ID].ins,
                            ppl=ppl),
        'gas': dict(COD=gas_COD, N=0, P=0, K=0)
        }
    return cache

def update_cache(sys):
    last_u = sys.path[-1]
    last_u._run()
    sys_dct['cache'][sys.ID] = cache_recoveries(sys)

A13.specification = lambda: update_cache(sysA)
B15.specification = lambda: update_cache(sysB)
C13.specification = lambda: update_cache(sysC)


def get_summarizing_fuctions():
    func_dct = {}
    func_dct['get_annual_cost'] = lambda tea, ppl: tea.EAC/ppl
    func_dct['get_annual_CAPEX'] = lambda tea, ppl: tea.annualized_CAPEX/ppl
    func_dct['get_annual_OPEX'] = lambda tea, ppl: tea.AOC/ppl
    ind = 'GlobalWarming'
    func_dct['get_annual_GWP'] = \
        lambda lca, ppl: lca.total_impacts[ind]/lca.lifetime/ppl
    func_dct['get_constr_GWP'] = \
        lambda lca, ppl: lca.total_construction_impacts[ind]/lca.lifetime/ppl
    func_dct['get_trans_GWP'] = \
        lambda lca, ppl: lca.total_transportation_impacts[ind]/lca.lifetime/ppl  
    func_dct['get_direct_emission_GWP'] = \
        lambda lca, ppl: lca.get_stream_impacts(stream_items=lca.stream_inventory, kind='direct_emission')[ind] \
            /lca.lifetime/ppl
    func_dct['get_offset_GWP'] = \
        lambda lca, ppl: lca.get_stream_impacts(stream_items=lca.stream_inventory, kind='offset')[ind] \
            /lca.lifetime/ppl
    func_dct['get_other_GWP'] = \
        lambda lca, ppl: lca.total_other_impacts[ind]/lca.lifetime/ppl
    for i in ('COD', 'N', 'P', 'K'):
        func_dct[f'get_liq_{i}_recovery'] = \
            lambda sys, i: sys_dct['cache'][sys.ID]['liq'][i]
        func_dct[f'get_sol_{i}_recovery'] = \
            lambda sys, i: sys_dct['cache'][sys.ID]['sol'][i]
        func_dct[f'get_gas_{i}_recovery'] = \
            lambda sys, i: sys_dct['cache'][sys.ID]['gas'][i]
        func_dct[f'get_tot_{i}_recovery'] = \
            lambda sys, i: \
                sys_dct['cache'][sys.ID]['liq'][i] + \
                sys_dct['cache'][sys.ID]['sol'][i] + \
                sys_dct['cache'][sys.ID]['gas'][i]
    return func_dct


def print_summaries(systems):
    try: iter(systems)
    except: systems = (systems, )
    func = get_summarizing_fuctions()
    for sys in systems:
        sys.simulate()
        ppl = sys_dct['ppl'][sys.ID]
        print(f'\n---------- Summary for {sys} ----------\n')
        tea = sys_dct['TEA'][sys.ID]
        tea.show()
        print('\n')
        lca = sys_dct['LCA'][sys.ID]
        lca.show()
        
        unit = f'{currency}/cap/yr'
        print(f'\nNet cost: {func["get_annual_cost"](tea, ppl):.1f} {unit}.')
        print(f'Capital: {func["get_annual_CAPEX"](tea, ppl):.1f} {unit}.')
        print(f'Operating: {func["get_annual_OPEX"](tea, ppl):.1f} {unit}.')
        
        unit = f'{GWP.unit}/cap/yr'
        print(f'\nNet emission: {func["get_annual_GWP"](lca, ppl):.1f} {unit}.')
        print(f'Construction: {func["get_constr_GWP"](lca, ppl):.1f} {unit}.')
        print(f'Transportation: {func["get_trans_GWP"](lca, ppl):.1f} {unit}.')
        print(f'Direct emission: {func["get_direct_emission_GWP"](lca, ppl):.1f} {unit}.')
        print(f'Offset: {func["get_offset_GWP"](lca, ppl):.1f} {unit}.')
        print(f'Other: {func["get_other_GWP"](lca, ppl):.1} {unit}.\n')

        for i in ('COD', 'N', 'P', 'K'):
            print(f'Total {i} recovery is {func[f"get_tot_{i}_recovery"](sys, i):.1%}, '
                  f'{func[f"get_liq_{i}_recovery"](sys, i):.1%} in liquid, '
                  f'{func[f"get_sol_{i}_recovery"](sys, i):.1%} in solid, '
                  f'{func[f"get_gas_{i}_recovery"](sys, i):.1%} in gas.')

def save_all_reports():
    import os
    path = os.path.dirname(os.path.realpath(__file__))
    del os
    for i in (sysA, sysB, sysC, lcaA, lcaB, lcaC):
        if isinstance(i, bst.System):
            i.simulate()
            i.save_report(f'{path}/results/{i.ID}.xlsx')
        else:
            i.save_report(f'{path}/results/{i.system.ID}_lca.xlsx')

__all__ = ('sysA', 'sysB', 'sysC', 'teaA', 'teaB', 'teaC', 'lcaA', 'lcaB', 'lcaC',
           'print_summaries', 'save_all_reports',
           *(i.ID for i in sysA.units),
           *(i.ID for i in sysB.units),
           *(i.ID for i in sysC.units),
           )







