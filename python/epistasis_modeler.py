
"""
Version of Starr model that recreates epistasis modeling for machine learning cross-validation
Some differences: 1) Starr model uses actual measured single mutation effects when measured for each library
2) Starr  models each library separately and averages log10Ka for final version (but say joint is nearly as good)
3) Some updates in dms_variants and binarymap seem to have altered model output slightly from published results
4) Starr models used raw data for fitting, did not aggregate over same mutations first

df 

Created on 12/23/2022

@author: Thomas Sheffield
"""

# import collections
# import os
# import itertools
# import random
# import tempfile
import time
# import warnings

import pandas as pd

# from plotnine import *

import binarymap
# import dms_variants.codonvarianttable
import dms_variants.globalepistasis
# import dms_variants.plotnine_themes
# import dms_variants.simulate
# from dms_variants.constants import CBPALETTE, CODONS_NOSTOP

# df = pd.read_csv("C:/Users/tysheff/covid_variants/Datasets/SARS-CoV-2-RBD_DMS/results/binding_Kds/binding_Kds.csv")
# df.rename(columns={'log10Ka':'func_score'},inplace=True)
# df['func_score_var'] = df['log10SE']**2

def epistasis_modeler(df, dt):
  
  #ensure no NA or NULL in the R code now
  func_scores = df[pd.notnull(df['func_score'])]
  func_scores.fillna('',inplace=True)
  
  bmap = binarymap.BinaryMap(func_scores, n_pre_col = None, n_post_col = None)
  
  start = time.time()
  model = dms_variants.globalepistasis.MonotonicSplineEpistasisCauchyLikelihood(bmap)
  model.fit()  # do NOT change ftol in normal use, this is just for test
  print(f"fitting took {time.time() - start:.1f} sec.")
  
  #add predicted and latent phenotypes to table by barcode with additional columns used for interpretation, save output file
  dt[['aa_substitutions']] = dt[['aa_substitutions']].fillna(value='')
  
  dt = model.add_phenotypes_to_df(df=dt, substitutions_col='aa_substitutions',latent_phenotype_col='latent_phenotype',observed_phenotype_col='preds',unknown_as_nan=True)
  lat = model.single_mut_effects(phenotype='latent', standardize_range=False)
  obs = model.single_mut_effects(phenotype='observed', standardize_range=False)


  return [dt, lat, obs]


