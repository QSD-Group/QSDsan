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
    [1] Recheck unit consistency, power law degradation
    [2] Why existing plant has sewer and sludge population? 

'''


# %%

import numpy as np
import biosteam as bst
import qsdsan as qs
from sklearn.linear_model import LinearRegression as LR
from qsdsan import sanunits as su
from qsdsan import WasteStream, SanUnit, ImpactItem, StreamImpactItem, \
    SimpleTEA, LCA
import bwaise
cmps = bwaise._cmps.cmps

bst.settings.set_thermo(cmps)
e = bst.PowerUtility
e.price = 0.17
currency = qs.currency = 'USD'
qs.CEPCI = qs.CEPCI_by_year[2018]
items = ImpactItem._items
GWP = qs.ImpactIndicator._indicators['GWP']



# %%

class ExistingWWTPCost(SanUnit):
    '''Lumped CAPEX and electricity cost of the wastewater treatment plant.'''
    def _run(self):
        self.outs[0].copy_like(self.ins[0])
    
    _BM = {'Wastewater treatment plant': 1}
    
    def _cost(self):
        self.purchase_costs['Wastewater treatment plant'] = 18606700
        self.power_utility(57120/(365*24)) #!!! Really? Only 6.5 kWh per hour?



# %%

# =============================================================================
# Assumptions
# =============================================================================

toilet_user = 16 # four people per household, four households per toilet
# ppl_existing = 4e4 # number of people served by the existing plant (sysA and SceC)
ppl_existing = 416667
ppl_alternative = 5e4 # number of people served by the alternative plant (SceB)

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
UGX_prices = np.array((8e4, 12e4, 20e4, 25e4))
capacities = np.array((3, 4.5, 8, 15))
def get_tanker_truck_cost(capacity):
    # Add 15% additional costs
    prices = UGX_prices*1.15/get_exchange_rate()
    ln_p = np.log(prices)
    ln_cap = np.log(capacities)
    model = LR().fit(ln_cap.reshape(-1,1), ln_p.reshape(-1,1))
    [[predicted]] = model.predict(np.array((np.log(capacity))).reshape(1, -1)).tolist()
    cost = np.exp(predicted)
    return cost

items['Concrete'].price = 194
items['Steel'].price = 2.665

N_AD_rx = 3

# Nutrient loss during applciation
#!!! Maybe the loss shouldn't be taken into account in cost and emission?
app_loss = dict.fromkeys(('NH3', 'NonNH3', 'P', 'K', 'Mg', 'Ca'), 0.02)
app_loss['NH3'] = 0.05



# %%

# =============================================================================
# Scenario A: existing system
# =============================================================================

fugitive_CH4 = WasteStream('fugitive_CH4', phase='g')
fugitive_N2O = WasteStream('fugitive_N2O', phase='g')
# Recycled nutrients are sold at a lower price than commercial fertilizers
price_factor = 0.25
liq_N = WasteStream('liq_N', phase='l', price=1.507*price_factor)
sol_N = WasteStream('sol_N', phase='l', price=1.507*price_factor)
liq_P = WasteStream('liq_P', phase='l', price=3.983*price_factor)
sol_P = WasteStream('sol_P', phase='l', price=3.983*price_factor)
liq_K = WasteStream('liq_K', phase='l', price=1.333*price_factor)
sol_K = WasteStream('sol_K', phase='l', price=1.333*price_factor)
fertilizers = (liq_N, sol_N, liq_P, sol_P, liq_K, sol_K)

print('\n----------Scenario A----------\n')
A1 = su.Excretion('A1', outs=('urine', 'feces'))

A2 = su.PitLatrine('A2', ins=(A1-0, A1-1,
                              'toilet_paper', 'flushing_water',
                              'cleansing_water', 'desiccant'),
                      outs=('mixed_waste', 'leachate', '', ''),
                      N_user=toilet_user, N_toilet=ppl_existing/toilet_user,
                      OPEX_over_CAPEX=0.05,
                      decay_k_COD=get_decay_k(tau_deg, log_deg),
                      decay_k_N=get_decay_k(tau_deg, log_deg),
                      max_CH4_emission=max_CH4_emission)

A3 = su.Trucking('A3', ins=A2-0, outs=('transported', 'loss'),
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
A3._specification = update_A3_param


A4 = ExistingWWTPCost('A4', ins=A3-0)

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


A9 = su.CropApplication('A8', ins=A7-0, outs=('liquid_fertilizer', 'loss'),
                        loss_ratio=app_loss)
def adjust_NH3_loss():
    A9._run()
    # Assume the slight higher loss of NH3 does not affect COD,
    # does not matter much since COD not considered in crop application
    A9.outs[0]._COD = A9.outs[1]._COD = A9.ins[0]._COD
A9.specification = adjust_NH3_loss

A10 = su.Mixer('A10', ins=(A2-2, A5-2, A6-1, A7-1, A8-2),
                 outs=fugitive_CH4)
A10.line = 'CH4 mixer'

A11 = su.Mixer('A11', ins=(A2-3, A5-3, A6-2, A7-2, A8-3),
                 outs=fugitive_N2O)
A11.line = 'N2O mixer'

A12 = su.ComponentSplitter('A12', ins=A8-0,
                              outs=(sol_N, sol_P, sol_K, 'sol_non_fertilizers'),
                              splits=(('NH3', 'NonNH3'), 'P', 'K'))

A13 = su.ComponentSplitter('A13', ins=A9-0,
                              outs=(liq_N, liq_P, liq_K, 'liq_non_fertilizers'),
                              splits=(('NH3', 'NonNH3'), 'P', 'K'))

sysA = bst.System('sysA',
                  path=(A1, A2, A3, treatA, A9, A10, A11, A12, A13))
sysA.simulate()
# sysA.save_report('results/sysA.xlsx')


# Emissions and product credits
CH4_item = StreamImpactItem(fugitive_CH4, GWP=28)
N2O_item = StreamImpactItem(fugitive_N2O, GWP=265)
liq_N_item = StreamImpactItem(liq_N, GWP=-5.4)
sol_N_item = StreamImpactItem(sol_N, GWP=-5.4)
liq_P_item = StreamImpactItem(liq_P, GWP=-4.9)
sol_P_item = StreamImpactItem(sol_P, GWP=-4.9)
liq_K_item = StreamImpactItem(liq_K, GWP=-1.5)
sol_K_item = StreamImpactItem(sol_K, GWP=-1.5)

teaA = SimpleTEA(system=sysA, discount_rate=0.05, start_year=2018,
                 lifetime=8, uptime_ratio=1, lang_factor=None,
                 annual_maintenance=0, annual_labor=12*3e6*12/get_exchange_rate(),
                 construction_schedule=None)
teaA.show()
print('\n')

e_item = ImpactItem(ID='e_item', functional_unit='kWh', GWP=0.15)
# 57120 is the annual electricity usage for the whole treatment plant
get_annual_e = lambda: A4.power_utility.rate*teaA._operating_hours

lcaA = LCA(system=sysA, lifetime=8, lifetime_unit='yr', uptime_ratio=1,
            # Assuming all additional WWTP OPEX from electricity
            e_item=get_annual_e()*8)
lcaA.show()
print('\n')


# %%

# =============================================================================
# Summaries
# =============================================================================

def get_total_inputs(unit):
    if unit is A1:
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

def get_recovery(unit_in=A1, outs=None, if_relative=True):
    inputs = get_total_inputs(unit_in)
    try: iter(outs)
    except: outs = (outs,)
    liq_sol = tuple(i for i in outs if i.phase != 'g')
    recovery = {}
    recovery['COD'] = sum(i.COD*i.F_vol/1e3 for i in liq_sol)
    recovery['energy'] = recovery['COD'] * 14e3
    recovery['N'] = sum(i.TN*i.F_vol/1e3 for i in liq_sol)
    recovery['NH3'] = sum(i.imass['NH3'] for i in liq_sol)
    recovery['P'] = sum(i.TP*i.F_vol/1e3 for i in liq_sol)
    recovery['K'] = sum(i.TK*i.F_vol/1e3 for i in liq_sol)
    recovery['Mg'] = sum(i.TMg*i.F_vol/1e3 for i in liq_sol)
    recovery['Ca'] = sum(i.TCa*i.F_vol/1e3 for i in liq_sol)
    for i, j in inputs.items():
        if if_relative:
            recovery[i] /= j/(365*24) * ppl_existing
        else:
            recovery[i] /= 1/(365*24) * ppl_existing
    return recovery

def get_emissions(outs):
    try: iter(outs)
    except: outs = (outs,)
    gas = tuple(i for i in outs if i.phase == 'g')
    emission = {}
    emission['direct'] = \
        sum((i.imass['CH4', 'N2O']*(28, 265)).sum() for i in gas)*365*24/ppl_existing
    return emission


def print_summaries():
    get_AOC_cap = lambda: teaA.AOC/ppl_existing
    get_EAC_cap = lambda: teaA.EAC/ppl_existing
    print(f'Without CAPEX, the net cost is {get_AOC_cap():.1f} {currency}/cap/yr.')
    print(f'With CAPEX, the net cost is {get_EAC_cap():.1f} {currency}/cap/yr.')
    
    get_GWP = lambda: lcaA.total_impacts['GlobalWarming']/8/ppl_existing
    print(f'Net emission is {get_GWP():.1f} {GWP.unit}/cap/yr.')

    get_total_N = lambda: \
        (A1.outs[0].imass['NH3', 'NonNH3']+A1.outs[1].imass['NH3', 'NonNH3']).sum()
    get_liq_N_recovery = lambda: liq_N.F_mass/ppl_existing/get_total_N()
    get_sol_N_recovery = lambda: sol_N.F_mass/ppl_existing/get_total_N()
    get_N_recovery = lambda: get_liq_N_recovery()+get_sol_N_recovery()
    print(f'Total N recovery is {get_N_recovery():.1%}, '
          f'{get_liq_N_recovery():.1%} in liquid, '
          f'{get_sol_N_recovery():.1%} in solid.')
    
    get_total_P = lambda: \
        (A1.outs[0].imass['P']+A1.outs[1].imass['P']).sum()
    get_liq_P_recovery = lambda: liq_P.F_mass/ppl_existing/get_total_P()
    get_sol_P_recovery = lambda: sol_P.F_mass/ppl_existing/get_total_P()
    get_P_recovery = lambda: get_liq_P_recovery()+get_sol_P_recovery()
    print(f'Total P recovery is {get_P_recovery():.1%}, '
          f'{get_liq_P_recovery():.1%} in liquid, '
          f'{get_sol_P_recovery():.1%} in solid.')
    
    get_total_K = lambda: \
        (A1.outs[0].imass['K']+A1.outs[1].imass['K']).sum()
    get_liq_K_recovery = lambda: liq_K.F_mass/ppl_existing/get_total_K()
    get_sol_K_recovery = lambda: sol_K.F_mass/ppl_existing/get_total_K()
    get_K_recovery = lambda: get_liq_K_recovery()+get_sol_K_recovery()
    print(f'Total K recovery is {get_K_recovery():.1%}, '
          f'{get_liq_K_recovery():.1%} in liquid, '
          f'{get_sol_K_recovery():.1%} in solid.')
    
    get_COD = lambda stream: stream.COD*stream.F_vol/1e3
    get_total_COD = lambda: get_COD(A1.outs[0])+get_COD(A1.outs[1])
    get_liq_COD_recovery = lambda: get_COD(A13.ins[0])/ppl_existing/get_total_COD()
    get_sol_COD_recovery = lambda: get_COD(A12.ins[0])/ppl_existing/get_total_COD()
    get_COD_recovery = lambda: get_liq_COD_recovery()+get_sol_COD_recovery()
    print(f'Total COD recovery is {get_COD_recovery():.1%}, '
          f'{get_liq_COD_recovery():.1%} in liquid, '
          f'{get_sol_COD_recovery():.1%} in solid.')



# %%

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
#                      design_type='unplanted',
#                      decay_k_COD=get_decay_k(tau_deg, log_deg),
#                      decay_k_N=get_decay_k(tau_deg, log_deg),
#                      max_CH4_emission=max_CH4_emission)

# A8.simulate()
# # A8.show()

# ws6 = A2.outs[0].copy('ws6')
# A9 = su.AnaerobicBaffledReactor('A9', ins=ws6, outs=('treated', 'CH4', 'N2O'),
#                                    decay_k_COD=get_decay_k(tau_deg, log_deg),
#                                    max_CH4_emission=max_CH4_emission)

# A9.simulate()
# # A9.show()

# ws7 = A2.outs[0].copy('ws7')
# A10 = su.LiquidTreatmentBed('A10', ins=ws7, outs=('treated', 'CH4', 'N2O'),
#                                decay_k_COD=get_decay_k(tau_deg, log_deg),
#                                decay_k_N=get_decay_k(tau_deg, log_deg),
#                                max_CH4_emission=max_CH4_emission)

# A10.simulate()
# A10.show()








# biogas = WasteStream('biogas', CH4=1)
# AX2 = su.BiogasCombustion('AX2', ins=(biogas, 'air'),
#                               outs=('used', 'lost', 'wasted'),
#                               if_combustion=True,
#                               biogas_loss=0.1, biogas_eff=0.55)



# sysA = bst.System('sysA', path=(A1, A2, A3, AX, AX2))

# sysA.simulate()


# e_item = ImpactItem(ID='e_item', functional_unit='kWh', GWP=20)

# lcaA = LCA(sysA, lifetime=8, lifetime_unit='yr',
#                e_item=(5000, 'Wh'))



# # %%

# # =============================================================================
# # Scenario B: anaerobic treatment with existing latrines and conveyance
# # =============================================================================


# # %%

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



# %%

# # =============================================================================
# # Try out impact-related classes
# # =============================================================================

# CEDf = ImpactIndicator(ID='CEDf', unit='MJ', method='Cumulative energy demand',
#                        category='Fossil')





