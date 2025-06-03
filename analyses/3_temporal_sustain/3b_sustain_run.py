
# Following this tutorial: https://github.com/ucl-pond/pySuStaIn/blob/master/notebooks/SuStaInWorkshop.ipynb

#%% Import modules

import os
import shutil
import numpy as np
import matplotlib.pyplot as plt
import pickle
from pathlib import Path
import sklearn.model_selection
import pandas as pd
import pylab
import sys
import pySuStaIn

# Run the SuStaIn disease modeling algorithm for WMH pathophysiology within each of the 3 spatial regions

os.makedirs(f'./results/3b_sustain_run', exist_ok=True)
os.makedirs(f'./visualization/3b_sustain_run', exist_ok=True)

#%% Prepare inputs

# Call script from terminal with the following arguments:

# Arguments:
#   - Cluster #
#   - Total number of fold
#   - fold number for this run
# N_folds = int(10)
# fold_no = int(1)
# clust = int(1)

clust = int(sys.argv[1])
N_folds = int(sys.argv[2])
fold_no = int(sys.argv[3])

print(f'N folds = {N_folds}')
print(f'Fold number = {fold_no}')
print(f"Cluster number {clust}")

df = pd.read_csv("./results/3a_sustain_prep/subj_WMH_patho_clean.tsv", delimiter='\t')

# Keep only micro columns with specific cluster
df = df[df.columns[df.columns.str.contains(f'c{clust}')]]

# Remove WMHvol column
df = df.iloc[:,:-1]

# Double check
for col in df.columns:
    print(f"First 10 values of {col}:")
    print(df[col].iloc[:10].tolist())
    max_value = df[col].max()
    print(f"Maximum value: {max_value}")

micro = ["FA", "MD", "ICVF", "ISOVF", "OD", "T2star", "QSM"]

# Custom z-vals thresholds for each marker and each cluster (see zscore_subject_regions)

if clust==1:
    Z_vals = np.zeros((len(micro), 7))
    Z_vals[[0]] = [0,0,0,0,0,0.5,1] # FA
    Z_vals[[1]] = [0,0,0,0,0.5,1,1.5] # MD
    Z_vals[[2]] = [0,0,0,0,0.5,1,1.5] # ICVF
    Z_vals[[3]] = [0,0,0,0,0,0.5,1] # ISOVF
    Z_vals[[4]] = [0,0,0,0,0,0,0.5] # OD
    Z_vals[[5]] = [0,0,0,0,0,0.5,1] # T2star
    Z_vals[[6]] = [0,0,0,0,0,0,0.5] # QSM
    Z_max  = np.array([2,3,3,2,1,2,1])


if clust==2:
    Z_vals = np.zeros((len(micro), 7))
    Z_vals[[0]] = [0,0,0,0.5,1,1.5,2] # FA
    Z_vals[[1]] = [0,0.5,1,1.5,2,3,4] # MD
    Z_vals[[2]] = [0,0,0.5,1,1.5,2,3] # ICVF
    Z_vals[[3]] = [0,0,0,0.5,1,1.5,2] # ISOVF
    Z_vals[[4]] = [0,0,0,0,0.5,1,1.5] # OD
    Z_vals[[5]] = [0,0,0,0.5,1,1.5,2] # T2star
    Z_vals[[6]] = [0,0,0,0,0,0.5,1] # QSM
    Z_max  = np.array([3,5,4,3,3,3,2])

if clust==3:
    Z_vals = np.zeros((len(micro), 7))
    Z_vals[[0]] = [0,0,0,0.5,1,1.5,2] # FA
    Z_vals[[1]] = [0.5,1,1.5,2,3,4,5] # MD
    Z_vals[[2]] = [0.5,1,1.5,2,3,4,5] # ICVF
    Z_vals[[3]] = [0,0,0,0.5,1,1.5,2] # ISOVF
    Z_vals[[4]] = [0,0,0,0,0.5,1,1.5] # OD
    Z_vals[[5]] = [0,0,0,0.5,1,1.5,2] # T2star
    Z_vals[[6]] = [0,0,0,0,0,0.5,1] # QSM
    Z_max  = np.array([3,7,7,3,3,3,2])

print(f"Zvals = {Z_vals}")
print(f"Zmax = {Z_max}")

N = len(micro)
SuStaInLabels = micro

N_startpoints = 15
N_S_max = 3
N_iterations_MCMC = int(1e4)
output_folder = f'results/3b_sustain_run/c{clust}'
dataset_name = f'WMH_micro_c{clust}'

sustain_input = pySuStaIn.ZscoreSustain(
                              df.values,
                              Z_vals,
                              Z_max,
                              SuStaInLabels,
                              N_startpoints,
                              N_S_max, 
                              N_iterations_MCMC, 
                              output_folder, 
                              dataset_name, 
                              True)

cv = sklearn.model_selection.KFold(N_folds, shuffle=True, random_state=123)
cv_it = cv.split(df.values)

test_ids = []
for train, test in cv_it:
    test_ids.append(test)

#%% Run Sustain with cross-validation

CVIC, loglike_matrix = sustain_input.cross_validate_sustain_model(test_ids, select_fold=fold_no)
