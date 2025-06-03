
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pySuStaIn
import statsmodels.formula.api as smf
from scipy import stats
import sklearn.model_selection
import plotly.express as px
import sys
import pylab

# Bind cross-validated results across folds

os.makedirs(f'./results/3b_sustain_run', exist_ok=True)
os.makedirs(f'./visualization/3b_sustain_run', exist_ok=True)

#%% Prepare inputs

# Arguments:
#   - Cluster #
# clust = int(3)

clust = int(sys.argv[1])

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

# Remove columns with all 0s (added for pos var diagram)
Z_vals = Z_vals[:,(np.any(Z_vals != 0, axis=0))]

print(f"Zvals = {Z_vals}")
print(f"Zmax = {Z_max}")

N = len(micro)
SuStaInLabels = micro

N_startpoints = 15
N_S_max = 3
N_iterations_MCMC = int(1e4)
input_folder = 'results/3b_sustain_run'
dataset_name = f'WMH_micro_c{clust}'

sustain_input = pySuStaIn.ZscoreSustain(
                              df.values,
                              Z_vals,
                              Z_max,
                              SuStaInLabels,
                              N_startpoints,
                              N_S_max, 
                              N_iterations_MCMC, 
                              input_folder, 
                              dataset_name, 
                              True)

#%% Cross-validation

output_folder = 'results/3c_sustain_bind_CV'

# Load results (this won't rerun, but will load the CV results)
N_folds = 10

cv = sklearn.model_selection.KFold(N_folds, shuffle=True, random_state=123)
cv_it = cv.split(df.values)

test_ids = []
for train, test in cv_it:
    test_ids.append(test)

CVIC, loglike_matrix = sustain_input.cross_validate_sustain_model(test_ids)

np.savetxt(f'{output_folder}/c{clust}/CV_CVIC.csv', CVIC, delimiter=',')
np.savetxt(f'{output_folder}/c{clust}/CV_loglike_matrix.csv', loglike_matrix, delimiter=',')

#%% Save results to csv

pk_base = pd.read_pickle(f'./{output_folder}/pickle_files/{dataset_name}_fold0_subtype0.pickle')

samples_sequence_cval = np.zeros([pk_base['samples_sequence'].shape[1], pk_base['samples_sequence'].shape[2]*N_folds])

for s in range(N_S_max):
    for this_s in range(s+1):
        for f in range(N_folds):
            print(f's = {s}; this_s = {this_s}; f = {f}')
            pk = pd.read_pickle(f"./{output_folder}/pickle_files/{dataset_name}_fold{f}_subtype{s}.pickle")
            samples_sequence_cval[:,(f*N_iterations_MCMC):((f+1)*N_iterations_MCMC)] = pk['samples_sequence'][this_s,:,:]
        print(samples_sequence_cval.shape)
        np.savetxt(f'{output_folder}/c{clust}/subtype{s}/samples_sequence_cval_{this_s}.csv', samples_sequence_cval, delimiter=',')

