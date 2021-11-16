#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Code is based on Bwaise systems.py developed by:
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/master/LICENSE.txt
for license details.
"""
import numpy as np
import biosteam as bst
import qsdsan as qs
from sklearn.linear_model import LinearRegression as LR
from qsdsan import sanunits as su
from qsdsan import WasteStream, ImpactIndicator, ImpactItem, StreamImpactItem, SimpleTEA, LCA
from _cmps import cmps



#import all units
# from _murt_toilet import MURTToilet
# from _primary_septic_tank import PrimarySepticTank
# from _ultrafiltration import Ultrafiltration
# from _ion_exchange_reclaimer import IonExchangeReclaimer
# from _ECR import ECR
# from _controls import Controls
# from _housing import Housing
# from _system_reclaimer import SystemReclaimer
# from _solar import Solar

# =============================================================================
# Unit parameters
# =============================================================================
bst.settings.set_thermo(cmps)
#items = ImpactItem._items

currency = qs.currency = 'USD'
bst.speed_up()

household_size = 4
get_household_size = lambda: household_size
household_per_toilet = 5
get_household_per_toilet = lambda: household_per_toilet
get_toilet_user = lambda: get_household_size()*get_household_per_toilet()

# Number of people served by the Reclaimer 2.0
ppl = 120

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
get_max_CH4_emission = lambda: max_CH4_emission

# Model for tanker truck cost based on capacity (m3)
# price = a*capacity**b -> ln(price) = ln(a) + bln(capacity)
USD_price_dct = np.array((21.62, 32.43, 54.05, 67.57))
capacities = np.array((3, 4.5, 8, 15))
emptying_fee = 0.15
get_emptying_fee = lambda: emptying_fee
def get_tanker_truck_fee(capacity):
    price_dct = USD_price_dct*(1+get_emptying_fee())
    ln_p = np.log(price_dct)
    ln_cap = np.log(capacities)
    model = LR().fit(ln_cap.reshape(-1,1), ln_p.reshape(-1,1))
    [[predicted]] = model.predict(np.array((np.log(capacity))).reshape(1, -1)).tolist()
    cost = np.exp(predicted)
    return cost

# =============================================================================
# Prices and GWP CFs
# =============================================================================

# Recycled nutrients are sold at a lower price than commercial fertilizers

price_factor = 0.25
get_price_factor = lambda: price_factor

operator_daily_wage = 29
get_operator_daily_wage = lambda: operator_daily_wage

price_dct = {
    'Electricity': 0.06,
    'Concrete': 194,
    'Steel': 2.665,
    'N': 1.507*get_price_factor(),
    'P': 3.983*get_price_factor(),
    'K': 1.333*get_price_factor(),
    'Polymer': 0.75,
    'Resin': 3.335,
    'FilterBag': 4.81,
    'MgOH2': 0,
    'MgCO3': 0.9,
    'H2SO4': 0.3,
    'struvite': 0,
    'salt': 0,
    'HCl_acid': 0,
    'KCl': 0
    }


GWP_dct = {
    'Electricity': 0.69,
    'CH4': 28,
    'N2O': 265,
    'N': -5.4,
    'P': -4.9,
    'K': -1.5,
    'Polymer':2.8,
    'Resin': 1.612,
    'FilterBag': 0.464,
    'MgOH2': 1.176277921,
    'MgCO3': 1.176277921,
    'H2SO4': 0.158899487,
    'struvite': 0,
    'salt': 0.266695553, 
    'HCl_acid': 0.8,
    'KCl': 0.8
    
    }

items = ImpactItem.get_all_items()

if not items.get('Excavation'): # prevent from reloading
    import os
    path = os.path.dirname(os.path.realpath(__file__)) + '/data'
    ImpactIndicator.load_from_file(path+'/impact_indicators.csv')
    item_path = path+'/_impact_item.xlsx'
    ImpactItem.load_from_file(item_path)
    del os

GWP = qs.ImpactIndicator.get_indicator('GWP')

bst.PowerUtility.price = price_dct['Electricity']
items['Concrete'].price = price_dct['Concrete']
items['Steel'].price = price_dct['Steel']

# =============================================================================
# Universal units and functions
# =============================================================================

CH4_item = StreamImpactItem(ID='CH4_item', GWP=GWP_dct['CH4'])
N2O_item = StreamImpactItem(ID='N2O_item', GWP=GWP_dct['N2O'])
N_item = StreamImpactItem(ID='N_item', GWP=GWP_dct['N'])
P_item = StreamImpactItem(ID='P_item', GWP=GWP_dct['P'])
K_item = StreamImpactItem(ID='K_item', GWP=GWP_dct['K'])
e_item = ImpactItem(ID='e_item', functional_unit='kWh', GWP=GWP_dct['Electricity'])
polymer_item = StreamImpactItem(ID='polymer_item', GWP=GWP_dct['Polymer'])
resin_item = StreamImpactItem(ID='resin_item', GWP=GWP_dct['Resin'])
filter_bag_item = StreamImpactItem(ID='filter_bag_item', GWP=GWP_dct['FilterBag'])
MgOH2_item = StreamImpactItem(ID='MgOH2_item', GWP=GWP_dct['MgOH2'])
MgCO3_item = StreamImpactItem(ID='MgCO3_item', GWP=GWP_dct['MgCO3'])
H2SO4_item = StreamImpactItem(ID='H2SO4_item', GWP=GWP_dct['H2SO4'])
struvite_item = StreamImpactItem(ID = 'struvite_item', GWP=GWP_dct['struvite'])
salt_item = StreamImpactItem(ID = 'salt_item', GWP=GWP_dct['salt'])
KCl_item = StreamImpactItem(ID = 'KCl_item', GWP=GWP_dct['KCl'])
#sludge_item = StreamImpactItem(ID = 'sludge_item', GWP=GWP_dct['sludge'])

def batch_create_streams(prefix):
    stream_dct = {}
    stream_dct['CH4'] = WasteStream(f'{prefix}_CH4', phase='g',
                                    stream_impact_item=CH4_item.copy(set_as_source=True))
    stream_dct['N2O'] = WasteStream(f'{prefix}_N2O', phase='g',
                                    stream_impact_item=N2O_item.copy(set_as_source=True))
    stream_dct['liq_N'] = WasteStream(f'{prefix}_liq_N', phase='l', price=price_dct['N'],
                                      stream_impact_item=N_item.copy(set_as_source=True))
    stream_dct['sol_N'] = WasteStream(f'{prefix}_sol_N', phase='l', price=price_dct['N'],
                                      stream_impact_item=N_item.copy(set_as_source=True))
    stream_dct['liq_P'] = WasteStream(f'{prefix}_liq_P', phase='l', price=price_dct['P'],
                                      stream_impact_item=P_item.copy(set_as_source=True))
    stream_dct['sol_P'] = WasteStream(f'{prefix}_sol_P', phase='l', price=price_dct['P'],
                                      stream_impact_item=P_item.copy(set_as_source=True))
    stream_dct['liq_K'] = WasteStream(f'{prefix}_liq_K', phase='l', price=price_dct['K'],
                                      stream_impact_item=K_item.copy(set_as_source=True))
    stream_dct['sol_K'] = WasteStream(f'{prefix}_sol_K', phase='l', price=price_dct['K'],
                                      stream_impact_item=K_item.copy(set_as_source=True))
    stream_dct['polymer'] = WasteStream(f'{prefix}_polymer', phase='s', price=price_dct['Polymer'],
                                      stream_impact_item=resin_item.copy(set_as_source=True))
    stream_dct['resin'] = WasteStream(f'{prefix}_resin', phase='s', price=price_dct['Resin'],
                                      stream_impact_item=resin_item.copy(set_as_source=True))
    stream_dct['filter_bag'] = WasteStream(f'{prefix}_filter_bag', phase='s', price=price_dct['FilterBag'],
                                      stream_impact_item=filter_bag_item.copy(set_as_source=True))
    stream_dct['MgOH2'] = WasteStream(f'{prefix}_MgOH2', phase='s', price=price_dct['MgOH2'],
                                      stream_impact_item=MgOH2_item.copy(set_as_source=True))
    stream_dct['MgCO3'] = WasteStream(f'{prefix}_MgCO3', phase='s', price=price_dct['MgCO3'],
                                      stream_impact_item=MgCO3_item.copy(set_as_source=True))
    stream_dct['H2SO4'] = WasteStream(f'{prefix}_H2SO4', phase='s', price=price_dct['H2SO4'],
                                      stream_impact_item=H2SO4_item.copy(set_as_source=True))
    stream_dct['struvite'] = WasteStream(f'{prefix}_struvite', phase='s', price=price_dct['struvite'],
                                      stream_impact_item=struvite_item.copy(set_as_source=True))
    stream_dct['salt'] = WasteStream(f'{prefix}_salt', phase='s', price=price_dct['salt'],
                                      stream_impact_item=salt_item.copy(set_as_source=True))
    stream_dct['KCl'] = WasteStream(f'{prefix}_KCl', phase='l', price=price_dct['KCl'],
                                      stream_impact_item=KCl_item.copy(set_as_source=True))
   # stream_dct['sludge'] = WasteStream(f'{prefix}_sludge', phase='s', price=price_dct['sludge'],
                                     # impact_item=sludge_item.copy(set_as_source=True))
    
    

    return stream_dct
    

def add_fugitive_items(unit, item):
    unit._run()
    for i in unit.ins:
        i.impact_item = item.copy(set_as_source=True)
        


# %%

# =============================================================================
# Scenario A (sysA): Duke Reclaimer Design 2.0
# =============================================================================

flowsheetA = bst.Flowsheet('sysA')
bst.main_flowsheet.set_flowsheet(flowsheetA)

streamsA = batch_create_streams('A')

#################### Human Inputs ####################
A1 = su.Excretion('A1', outs=('urine','feces'))

################### User Interface ###################
#Reclaimer 2.0 can process ~30L/hr(net), 720L/24 hours of constant operation
#flush volume of 6L per flush determines the number of users would be 120 users
A2 = su.MURTToilet('A2', ins=(A1-0, A1-1,
                              'toilet_paper', 'flushing_water',
                              'cleansing_water', 'desiccant'),
                    outs=('mixed_waste', 'A2_CH4', 'A2_N2O'),
                    decay_k_COD=get_decay_k(tau_deg, log_deg), 
                    decay_k_N=get_decay_k(tau_deg, log_deg),
                    max_CH4_emission=get_max_CH4_emission(),
                    N_user=120/7, N_toilet=7, 
                    if_flushing=True, if_desiccant=False, if_toilet_paper=True,
                    CAPEX = 0,
                    OPEX_over_CAPEX= 0.07)

###################### Treatment ######################
#Septic Tank 
A3 = su.PrimaryReclaimer('A3', ins=(A2-0), 
                    outs=('A3_treated', 'A3_CH4', 'A3_N2O', 'A3_sludge'), 
                    decay_k_COD=get_decay_k(tau_deg, log_deg), 
                    decay_k_N=get_decay_k(tau_deg, log_deg),
                    max_CH4_emission=get_max_CH4_emission())


A4 = su.Ultrafiltration('A4', ins=(A3-0), outs = ('A4_treated', 'retentate'))
                        
A5 = su.IonExchangeReclaimer('A5', ins=(A4-0, streamsA['GAC'], streamsA['KCl']),
                                outs=('A5_treated', 'SpentGAC',streamsA['Conc_NH3']),
                                decay_k_COD=get_decay_k(tau_deg, log_deg), 
                                decay_k_N=get_decay_k(tau_deg, log_deg),
                                max_CH4_emission=get_max_CH4_emission(), if_gridtied=True)

A6 = su.ECR('A6', ins=(A5-0, streamsA['salt'], streamsA['HCl_acid']), 
                    outs = ('A6_treated'),
                    decay_k_COD=get_decay_k(tau_deg, log_deg),)


A7 = su.Mixer('A7', ins=(A3-1, A2-1), outs=streamsA['CH4'])
A7.specification = lambda: add_fugitive_items(A7, CH4_item)
A7.line = 'fugitive CH4 mixer' 
        
A8 = su.Mixer('A8', ins=(A3-2, A2-2), outs=streamsA['N2O'])
A8.specification = lambda: add_fugitive_items(A8, N2O_item)
A8.line = 'fugitive N2O mixer'

################## Other impacts and costs ##################
#A9 = su.Misch('A9', ins=(A6-0), outs = ('A9_out'))



A10 = su.Trucking('A12', ins=A3-3, outs=('transported', 'conveyance_loss'),
                  load_type='mass', distance=5, distance_unit='km',
                  interval=365, interval_unit='d',
                  loss_ratio=0.02)
def update_A10_param():
    A10._run()
    truck = A10.single_truck
    truck.interval = 365*24
    truck.load = A10.F_mass_in*truck.interval
    rho = A10.F_mass_in/A10.F_vol_in
    vol = truck.load/rho
    A10.fee = get_tanker_truck_fee(vol)
    A10._design()
A10.specification = update_A10_param


############### Simulation, TEA, and LCA ###############
sysA = bst.System('sysA', path= (A1, A2, A3, A4, A5, A6, A7, A8, A10))

sysA.simulate()

power = sum([u.power_utility.rate for u in sysA.units])

#!!! update labor to input country specific data and be a distribution
teaA = SimpleTEA(system=sysA, discount_rate=get_discount_rate(), 
                  start_year=2020, lifetime=10, uptime_ratio=1, 
                  lang_factor=None, annual_maintenance=0, 
                  annual_labor=(get_operator_daily_wage() * 1 * 52), construction_schedule=None)

lcaA = LCA(system=sysA, lifetime=10, lifetime_unit='yr', uptime_ratio=1,
            e_item=lambda: power*(365*24)*10)



# %%

# =============================================================================
# Scenario B (sysB): MBR (Anoxic and Facultative Aerobic with MBR)
# =============================================================================


# %%

# =============================================================================
# Summarizing Functions
# =============================================================================

#!!! items need to be updated in sys_dct, system_streams, learning curve assumptions
sys_dct = {
    'ppl':  dict(sysA=300, sysB=300, sysC=300),
    'input_unit': dict(sysA=A1),
    'liq_unit': dict(sysA=None, sysB=None, sysC=None),
    'sol_unit': dict(sysA=None, sysB=None, sysC=None),
    'gas_unit': dict(sysA=None, sysB=None, sysC=None),
    'stream_dct': dict(sysA=streamsA),
    'TEA': dict(sysA=teaA),
    'LCA': dict(sysA=lcaA),
    'cache': dict(sysA={}, sysB={}, sysC={}),
    }

system_streams = {sysA: streamsA}

#learning curve assumptions
percent_CAPEX_to_scale = 0.1 #this number is the decimal of the fraction of scost of pecialty parts/cost of total parts
get_percent_CAPEX_to_scale = lambda: percent_CAPEX_to_scale

number_of_units = 100000

percent_limit = 0.015
get_percent_limit = lambda: percent_limit

learning_curve_percent=0.925
get_learning_curve_percent = lambda: learning_curve_percent

def get_scaled_capital(tea):
    CAPEX_to_scale = tea.annualized_CAPEX * get_percent_CAPEX_to_scale()
    CAPEX_not_scaled = tea.annualized_CAPEX - CAPEX_to_scale
    scaled_limited = CAPEX_to_scale * get_percent_limit()
    b=(np.log(get_learning_curve_percent())/np.log(2))
    scaled_CAPEX_annualized  = (CAPEX_to_scale - scaled_limited)*number_of_units**b+scaled_limited    
    new_CAPEX_annualized = scaled_CAPEX_annualized + CAPEX_not_scaled    
    return new_CAPEX_annualized
    

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

#10 corresponds to year
def get_fugitive_emissions(streams=None, hr=365*24*10):
    try: iter(streams)
    except: streams = (streams,)
    factor = hr
    emission=0
    for i in streams:
        if not i.impact_item: continue
        emission = i.F_mass*i.impact_item.CFs['GlobalWarming']*factor
    return emission

#N2O and CH4 emissions
def get_fugitive_gases(system_num=None):
    fugitive_streams = ('CH4', 'N2O')
    gases = 0
    for i in fugitive_streams:
        gases+=get_fugitive_emissions(system_streams[system_num][i])
    return gases

def cache_recoveries(sys):
    total_COD = get_total_inputs(sys_dct['input_unit'][sys.ID])['COD']
    ppl = sys_dct['ppl'][sys.ID]
    if sys_dct['gas_unit'][sys.ID]:
        gas_mol = sys_dct['gas_unit'][sys.ID].outs[0].imol['CH4']
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

def get_summarizing_fuctions():
    func_dct = {}
    func_dct['get_annual_cost'] = lambda tea, ppl: (tea.EAC-tea.annualized_CAPEX
                                                    +get_scaled_capital(tea))/ppl
    func_dct['get_annual_CAPEX'] = lambda tea, ppl: get_scaled_capital(tea)/ppl
    func_dct['get_annual_energy'] = lambda tea, ppl: sum(i.power_utility.cost for i in tea.units)*tea.operating_hours/ppl
    func_dct['get_annual_OPEX'] = lambda tea, ppl: (tea.AOC/ppl -(sum(i.power_utility.cost for i in tea.units)
                                                    *tea.operating_hours/ppl) - (tea.annual_labor/ppl))
    func_dct['get_annual_labor'] = lambda tea, ppl: tea.annual_labor/ppl
    func_dct['get_annual_sales'] = lambda tea, ppl: (tea.sales/ppl)
    ind = 'GlobalWarming'
    func_dct['get_annual_GWP'] = \
        lambda lca, ppl: lca.total_impacts[ind]/lca.lifetime/ppl
    func_dct['get_constr_GWP'] = \
        lambda lca, ppl: lca.total_construction_impacts[ind]/lca.lifetime/ppl
    func_dct['get_trans_GWP'] = \
        lambda lca, ppl: lca.total_transportation_impacts[ind]/lca.lifetime/ppl  
    func_dct['get_CH4_N2O_GWP'] = \
        lambda sys, lca, ppl: get_fugitive_gases(sys)/lca.lifetime/ppl
    func_dct['get_stream_items_emission_GWP'] = \
        lambda sys, lca, ppl: (lca.get_stream_impacts(stream_items=lca.stream_inventory, kind='direct_emission')[ind]) \
            /lca.lifetime/ppl - (get_fugitive_gases(sys)/lca.lifetime/ppl) #maybe ignore
    func_dct['get_offset_GWP'] = \
        lambda lca, ppl: lca.get_stream_impacts(stream_items=lca.stream_inventory, kind='offset')[ind] \
            /lca.lifetime/ppl
    func_dct['get_other_GWP'] = \
        lambda lca, ppl: lca.total_other_impacts[ind]/lca.lifetime/ppl
    # for i in ('COD', 'N', 'P', 'K'):
    #     func_dct[f'get_liq_{i}_recovery'] = \
    #         lambda sys, i: sys_dct['cache'][sys.ID]['liq'][i]
    #     func_dct[f'get_sol_{i}_recovery'] = \
    #         lambda sys, i: sys_dct['cache'][sys.ID]['sol'][i]
    #     func_dct[f'get_gas_{i}_recovery'] = \
    #         lambda sys, i: sys_dct['cache'][sys.ID]['gas'][i]
    #     func_dct[f'get_tot_{i}_recovery'] = \
    #         lambda sys, i: \
    #             sys_dct['cache'][sys.ID]['liq'][i] + \
    #             sys_dct['cache'][sys.ID]['sol'][i] + \
    #             sys_dct['cache'][sys.ID]['gas'][i]
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
        print(f'\nTotal cost: {func["get_annual_cost"](tea, ppl):.2f} {unit}.')
        print(f'Capital: {func["get_annual_CAPEX"](tea, ppl):.2f} {unit}.')
        print(f'Energy: {func["get_annual_energy"](tea,ppl):.2f} {unit}.')
        print(f'Operating: {func["get_annual_OPEX"](tea, ppl):.2f} {unit}.')
        print(f'Labor: {func["get_annual_labor"](tea, ppl):.2f} {unit}.')
        print(f'Sales: {func["get_annual_sales"](tea, ppl):.2f} {unit}.')
        
        unit = f'{GWP.unit}/cap/yr'
        print(f'\nNet emission: {func["get_annual_GWP"](lca, ppl):.2f} {unit}.')
        print(f'Construction: {func["get_constr_GWP"](lca, ppl):.2f} {unit}.')
        print(f'Transportation: {func["get_trans_GWP"](lca, ppl):.2f} {unit}.')
        print(f'Fugitive gas: {func["get_CH4_N2O_GWP"](sys, lca, ppl):.2f} {unit}.')
        print(f'Stream items emission: {func["get_stream_items_emission_GWP"](sys, lca, ppl):.2f} {unit}.')
        print(f'Offset: {func["get_offset_GWP"](lca, ppl):.2f} {unit}.')
        print(f'Other: {func["get_other_GWP"](lca, ppl):.2} {unit}.\n')

        # for i in ('COD', 'N', 'P', 'K'):
        #     print(f'Total {i} recovery is {func[f"get_tot_{i}_recovery"](sys, i):.1%}, '
        #           f'{func[f"get_liq_{i}_recovery"](sys, i):.1%} in liquid, '
        #           f'{func[f"get_sol_{i}_recovery"](sys, i):.1%} in solid, '
        #           f'{func[f"get_gas_{i}_recovery"](sys, i):.1%} in gas.')

def save_all_reports():
 # import os
    # path = os.path.dirname(os.path.realpath(__file__))
    # path += '/results'
    # if not os.path.isdir(path):
    #     os.path.mkdir(path)
    # del os
    path = '/Users/torimorgan/opt/anaconda3/lib/python3.8/site-packages/exposan/Duke_Reclaimer'
    for i in (sysA,lcaA):
        if isinstance(i, bst.System):
            i.simulate()
            i.save_report(f'{path}/{i.ID}.xlsx')
        else:
            i.save_report(f'{path}/{i.system.ID}_lca.xlsx')


__all__ = ('sysA', 'teaA', 'lcaA',
            'print_summaries', 'save_all_reports',
            *(i.ID for i in sysA.units)
            )