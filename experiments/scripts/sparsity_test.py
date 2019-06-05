#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
Produce the results of SM Figure 2.

The tested time-series all have the same noise levels, but varying sampling
frequency. We expect more accurate results when the sampling is dense.

Author: Jean-Gabriel Young <jgyou@umich.edu>
"""
import pandas as pd
import pickle
from os.path import join

model_dir = '../../stan_models/models_bin/'
ts_dir = '../time_series/SM/'
res_dir = '../results/SM/'

# Load models
with open(join(model_dir, 'prevalence_model.bin'), 'rb') as f:
    prevalence_model = pickle.load(f)
with open(join(model_dir, 'incidence_model.bin'), 'rb') as f:
    incidence_model = pickle.load(f)

# Loop over time series and sample from the posterior distribution
for ts_file in ['T50_sigma_0p02',
                'T100_sigma_0p02',
                'T250_sigma_0p02']:
    # Read data as CSV
    df = pd.read_csv(join(ts_dir, ts_file + ".txt"), sep='\t')
    # Normalize time
    t = df['t'] / df['t'].iloc[-1]
    # Prepare data dicts
    data_prevalence = {'T': len(t),
                       'ts': t,
                       'Y': df['Y(t)'],
                       # Hyperparameters
                       'scale_sigma': 1,
                       'scale_gamma': 100,
                       'scale_xi_mean': 100,
                       'scale_xi_spread': 1,
                       'N': 8,
                       # Misc.
                       'overshoot': 0.1,
                       'num_steps_beta': 100,
                       'num_steps_y': 100,
                       'max_iter': 25000}
    data_incidence = {'T': len(t),
                      'ts': t,
                      'Z': df['Z(t)'],
                      'population': 100000,
                      'max_Y': 1,
                      # Hyperparameters
                      'scale_sigma': 1,
                      'loc_gamma': 0.1  * df['t'].iloc[-1],
                      'scale_gamma': 0.1,
                      'scale_xi_mean': 100,
                      'scale_xi_spread': 1,
                      'N': 8,
                      # Misc.
                      'overshoot': 0.00,
                      'num_steps_beta': 100,
                      'num_steps_y': 100,
                      'max_iter': 25000}
    # Fit with prevalence model and dump results to disk
    fit = prevalence_model.sampling(data_prevalence,
                                    iter=1000, chains=4,
                                    control={'max_treedepth': 15})
    with open(join(res_dir, 'sprs_prevalence_' + ts_file + '.pck'), 'wb') as f:
        pickle.dump(fit, f)
    # Fit with incidence model and dump results to disk
    fit = incidence_model.sampling(data_incidence,
                                   iter=1000, chains=4,
                                   control={'max_treedepth': 15})
    with open(join(res_dir, 'sprs_incidence_' + ts_file + '.pck'), 'wb') as f:
        pickle.dump(fit, f)
