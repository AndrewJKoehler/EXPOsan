# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Ga-Yeong Kim <gayeong1225@gmail.com>
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''
#%% Import everything for now

from qsdsan.utils import ospath, time_printer
from exposan.pm2_batch import (
    results_path,
    create_model,
    sensitive_params,
    )

import numpy as np, pandas as pd
from scipy.optimize import minimize, basinhopping, shgo

# import winsound as sd
# from datetime import datetime

__all__ = ('cali_setup', 'optimizer', 'objective_function')

# start_time = datetime.now()
# print('Started (hh:mm:ss.ms) {}'.format(start_time))

#%%

mdl = create_model(kind='exclude', analysis='cali')   # with uniform distribution of sensitive params (in model.py)
# mdl = create_model(kind='exclude', analysis='cali')   # with uniform distribution of sensitive params (in model.py)

def cali_setup():
    params = []
    baseline = np.array([])
    boundary = []

    for k, v in sensitive_params.items():
        b, units, bounds = v
        params.append(k)
        baseline = np.append(baseline, b)
        boundary.append(bounds)
    boundary = tuple(boundary)

    return params, baseline, boundary

sense_params, opt_params, bnds = cali_setup()         # start from init_guess

'''
sense_params: list of sensitive parameters
opt_params: initial guess of sensitive parameters
bnds: min & max of sensitive parameters
'''

#%%
def optimizer():

    opt = shgo(objective_function, bounds=bnds, iters=5, minimizer_kwargs={'method':'SLSQP', 'ftol':1e-2})

    opt_as_series = pd.Series(opt)

    opt_as_series.to_excel(excel_writer=(ospath.join(results_path, 'calibration_result_exclude_newbase_shgo_kinetic.xlsx')))

#%%
@time_printer
def objective_function(opt_params, *args):

    try:

        # mdl._update_state(opt_params, t_span=(0, 0.25), t_eval = np.arange(0, 0.26, 0.01), method='BDF', state_reset_hook='reset_cache')
        mdl._update_state(opt_params, t_span=(0, 7), t_eval = np.arange(0, 7.01, 0.01), method='BDF', state_reset_hook='reset_cache')

    except:
        return 0.5

    out = [metric() for metric in mdl.metrics]
    obj = np.average(out)

    return obj

optimizer()

# time_elapsed = datetime.now()-start_time
# print('Time elapsed (hh:mm:ss.ms) {}'.format(time_elapsed))

# sd.Beep(2000, 1000)