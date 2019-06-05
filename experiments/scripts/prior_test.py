#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
Produce the results of SM Figure 4-5.

The tested time-series are all identical. We vary the prior on gamma and check
how this affects inference with the incidence model.

Author: Jean-Gabriel Young <jgyou@umich.edu>
"""
import pandas as pd
import pickle
from os.path import join

model_dir = '../../stan_models/models_bin/'
ts_dir = '../time_series/SM/'
res_dir = '../results/SM/'

# Load models
with open(join(model_dir, 'incidence_model.bin'), 'rb') as f:
    incidence_model = pickle.load(f)

ts_file = "T100_sigma_0p01"
# Read data as CSV
df = pd.read_csv(join(ts_dir, ts_file + ".txt"), sep='\t')
# Normalize time
t = df['t'] / df['t'].iloc[-1]
# Loop over priors
for (loc_gamma, scale_gamma) in [(0.20, 10),    # overshoot, vague
                                 (0.05, 10),    # undershoot, vague
                                 (0.10, 10),    # on point, vague
                                 (0.20, 0.1)]:  # overshoot, confident
    # Prepare data dict
    data_incidence = {'T': len(t),
                      'ts': t,
                      'Z': df['Z(t)'],
                      'population': 100000,
                      'max_Y': 1,
                      # Hyperparameters
                      'scale_sigma': 1,
                      'loc_gamma': loc_gamma * df['t'].iloc[-1],
                      'scale_gamma': scale_gamma,
                      'scale_xi_mean': 100,
                      'scale_xi_spread': 1,
                      'N': 8,
                      # Misc.
                      'overshoot': 0.00,
                      'num_steps_beta': 100,
                      'num_steps_y': 100,
                      'max_iter': 15000}
    # Fit with incidence model and dump results to disk
    fit = incidence_model.sampling(data_incidence,
                                   iter=1000, chains=4,
                                   control={'max_treedepth': 15})
    fname = 'gamma_test_loc' + str(loc_gamma).replace('.', 'p') +\
            '_scale' + str(scale_gamma).replace('.', 'p') + '.pck'
    with open(join(res_dir, fname), 'wb') as f:
        pickle.dump(fit, f)
