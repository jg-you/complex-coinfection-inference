#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
Compiles the models with pystan.

Author: Jean-Gabriel Young <jgyou@umich.edu>
"""
from pickle import dump
from pystan import StanModel

with open('stan_models/models_bin/prevalence_model.bin', 'wb') as f:
    model = StanModel('stan_models/models_code/prevalence_model.stan',
                      model_name='prevalence',
                      include_paths='stan_models/models_code')
    dump(model, f)

with open('stan_models/models_bin/incidence_model.bin', 'wb') as f:
    model = StanModel('stan_models/models_code/incidence_model.stan',
                      model_name='incidence',
                      include_paths='stan_models/models_code')
    dump(model, f)
