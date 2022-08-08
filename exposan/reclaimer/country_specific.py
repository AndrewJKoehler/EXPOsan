#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Tori Morgan <tvlmorgan@gmail.com>
    Hannah Lohman <hlohman94@gmail.com>
    Yalin Li <mailto.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

from chaospy import distributions as shape
from qsdsan import ImpactItem, PowerUtility
from exposan.utils import general_country_specific_inputs, run_module_country_specific
from exposan import reclaimer as re
from exposan.reclaimer import (
    create_model,
    GWP_dct,
    H_Ecosystems_dct,
    H_Health_dct,
    H_Resources_dct,
    price_dct,
    run_uncertainty,
    results_path,
    )

__all__ = ('create_country_specific_model',)

# Filter out warnings related to uptime ratio
import warnings
warnings.filterwarnings('ignore', message='uptime_ratio')

from copy import copy
old_price_ratio = copy(re.price_ratio)


# %%

# =============================================================================
# Create model for country-specific analysis
# =============================================================================

get_default_uniform = lambda b, ratio: shape.Uniform(lower=b*(1-ratio), upper=b*(1+ratio))
format_key = lambda key: (' '.join(key.split('_'))).capitalize()

def create_country_specific_model(ID, country, country_data=None, model=None):
    model = model.copy() if model is not None else create_model(model_ID=ID, country_specific=True)
    param = model.parameter
    sys = model.system
    sys_stream = sys.flowsheet.stream
    country_data = country_data or general_country_specific_inputs[country]
    param_dct = {p.name: p for p in model.parameters}

    def get_param_name_b_D(key, ratio):
        name = format_key(key)
        p = param_dct.get(name)
        # Throw out this parameter (so that a new one can be added with updated values)
        # if it's already there
        if p: model.parameters = [i for i in model.parameters if i is not p]
        b = country_data[key]
        D = get_default_uniform(b, ratio)
        return name, p, b, D

    ##### Constant parameter settings #####
    # Diet and excretion
    excretion_unit = sys.path[0]
    excretion_unit.p_anim = country_data['p_anim']
    excretion_unit.p_veg = country_data['p_veg']
    excretion_unit.e_cal = country_data['e_cal']
    excretion_unit.waste_ratio = country_data['food_waste_ratio']

    # Price ratio
    i = country_data['price_ratio']
    re.price_ratio = i
    for u in sys.units:
        if hasattr(u, 'price_ratio'):
            u.price_ratio = i

    # Energy price
    PowerUtility.price = country_data['energy_price']

    ##### Uncertain parameter settings #####
    # Wages
    key = 'operator_daily_wage'
    wage_D_ratio = 0.5
    name, p, b, D = get_param_name_b_D(key, wage_D_ratio)
    @param(name=name,
           element=excretion_unit,  # just want to add to the 0th unit of the system
           kind='coupled', units='USD/h',
           baseline=b, distribution=D)
    def set_labor_wages(i):
        for u in sys.units:
            if hasattr(u, '_calc_maintenance_labor_cost'):
                u.wages = i

    # Energy GWP
    key = 'energy_GWP'
    GWP_D_ratio = 0.1
    name, p, b, D = get_param_name_b_D(key, GWP_D_ratio)
    @param(name=name, element='LCA', kind='isolated',
           units='kg CO2-eq/kWh', baseline=b, distribution=D)
    def set_electricity_CF(i):
        GWP_dct['Electricity'] = ImpactItem.get_item('e_item').CFs['GlobalWarming'] = i

    # Energy H_Ecosystems
    key = 'energy_H_Ecosystems'
    name, p, b, D = get_param_name_b_D(key, GWP_D_ratio)
    @param(name=name, element='LCA', kind='isolated',
           units='points/kWh', baseline=b, distribution=D)
    def set_electricity_ecosystems_CF(i):
        H_Ecosystems_dct['Electricity'] = ImpactItem.get_item('e_item').CFs['H_Ecosystems'] = i

    # Energy H_Health
    key = 'energy_H_Health'
    name, p, b, D = get_param_name_b_D(key, GWP_D_ratio)
    @param(name=name, element='LCA', kind='isolated',
           units='points/kWh', baseline=b, distribution=D)
    def set_electricity_health_CF(i):
        H_Health_dct['Electricity'] = ImpactItem.get_item('e_item').CFs['H_Health'] = i

    # Energy H_Resources
    key = 'energy_H_Resources'
    name, p, b, D = get_param_name_b_D(key, GWP_D_ratio)
    @param(name=name, element='LCA', kind='isolated',
           units='points/kWh', baseline=b, distribution=D)
    def set_electricity_resources_CF(i):
        H_Resources_dct['Electricity'] = ImpactItem.get_item('e_item').CFs['H_Resources'] = i

    # LPG
    key = 'LPG_price'
    price_D_ratio = 0.2
    name, p, b, D = get_param_name_b_D(key, price_D_ratio)
    @param(name=name, element='TEA', kind='isolated', units='USD/kg',
           baseline=b, distribution=D)
    def set_LPG_price(i):
        price_dct['LPG'] = sys_stream.LPG.price = i

    # N fertilizer price
    key = 'N_fertilizer_price'
    name, p, b, D = get_param_name_b_D(key, price_D_ratio)
    @param(name=name, element='TEA', kind='isolated', units='USD/kg N',
           baseline=b, distribution=D)
    def set_N_price(i):
        price_dct['N'] = sys_stream.liq_N.price = sys_stream.sol_N.price = i * re.price_factor

    # P fertilizer price
    key = 'P_fertilizer_price'
    name, p, b, D = get_param_name_b_D(key, price_D_ratio)
    @param(name=name, element='TEA', kind='isolated', units='USD/kg P',
           baseline=b, distribution=D)
    def set_P_price(i):
        price_dct['P'] = sys_stream.liq_P.price = sys_stream.sol_P.price = i * re.price_factor

    # K fertilizer price
    key = 'K_fertilizer_price'
    name, p, b, D = get_param_name_b_D(key, price_D_ratio)
    @param(name=name, element='TEA', kind='isolated', units='USD/kg K',
           baseline=b, distribution=D)
    def set_K_price(i):
        price_dct['K'] = sys_stream.liq_K.price = sys_stream.sol_K.price = i * re.price_factor

    # Concentrated NH3 fertilizer price
    key = 'NH3_fertilizer_price'
    name, p, b, D = get_param_name_b_D(key, price_D_ratio)
    @param(name=name, element='TEA', kind='isolated', units='USD/kg N',
           baseline=b, distribution=D)
    def set_con_NH3_price(i):
        price_dct['conc_NH3'] = sys_stream.conc_NH3.price = i * re.price_factor

    return model


def run_country(system_IDs, seed=None, N=1000):
    return run_module_country_specific(
        create_country_specific_model_func=create_country_specific_model,
        run_uncertainty_func=run_uncertainty,
        folder_path=results_path,
        system_IDs=system_IDs,
        seed=seed,
        N=N,
        )


if __name__ == '__main__':
    models = run_country(
        system_IDs=('sysB', 'sysC'),
        # system_IDs=('sysA', 'sysB', 'sysC', 'sysD'),
        seed=5, N=10,
        )