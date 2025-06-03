
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
# import seaborn as sns
# import joypy
from sklearn.model_selection import train_test_split
from pcntoolkit.normative import estimate, evaluate, predict
from pcntoolkit.util.utils import create_bspline_basis, compute_MSLL
import datatable as dt
import os
from itertools import product
import plotnine as pn
import random
import math
from joblib import Parallel, delayed
import multiprocessing as mp
import pickle
import shutil

random.seed(123)

# Normative modeling with the PCN-toolkit

# PCN-toolkit resources
# Paper: https://www.nature.com/articles/s41596-022-00696-5
# Documentation: https://pcntoolkit.readthedocs.io/en/latest/pages/BLR_normativemodel_protocol.html
# Github: https://github.com/predictive-clinical-neuroscience/PCNtoolkit-demo/tree/main

names = ["dti_FA", "dti_MD", "NODDI_ICVF", "NODDI_ISOVF", "NODDI_OD"]
new_names = ["FA", "MD", "ICVF", "ISOVF", "OD"]

# Run this script for each microstructural marker (number from 0 to 6)
name = names[int(sys.argv[1])]
new_name = new_names[int(sys.argv[1])]
print(new_name)

# Create directories to save files in RAM
os.makedirs('/dev/shm/parent41', exist_ok=True)
os.makedirs(f'/dev/shm/parent41/nm_n{name}', exist_ok=True)
os.makedirs(f'/dev/shm/parent41/nm_n{name}/tmp', exist_ok=True)

# Create output directories
os.makedirs(f'./results/1a_norm_models_BLR', exist_ok=True)
os.makedirs(f'./visualization/1a_norm_models_BLR', exist_ok=True)

# Load UKB data

# ids: final sample of UKB inclusions (32,526 subjects)
ids_ukb = dt.fread("../../QC/inclusions_final.txt").to_pandas()
ids_ukb.rename(columns={'C0': 'ID'}, inplace=True)

# demo: age and sex
demo_ukb = dt.fread("../../../UKB/tabular/df_demo/UKBB_demo_wider.tsv").to_pandas()
demo_ukb = demo_ukb[(demo_ukb['InstanceID'] == 2)]
demo_ukb = demo_ukb[['SubjectID', 'Sex_31_0', 'Age_when_attended_assessment_centre_21003_0']]
demo_ukb.rename(columns={'SubjectID': 'ID', 'Sex_31_0':'Sex', 'Age_when_attended_assessment_centre_21003_0':'Age'}, inplace=True)
demo_ukb['Sex'] = demo_ukb['Sex'].replace({'Female': 0, 'Male': 1}) # 0 = Female, 1 = Male
demo_ukb = pd.merge(ids_ukb, demo_ukb, on='ID', how='inner')

# micro: microstructural values (subject by voxel in UKB space)
micro_ukb = dt.fread(f"../../micro_matrices/ses2_{name}_after_exclusions.tsv").to_pandas()
# label: BISON labels (subject by voxel in UKB space)
label_ukb = dt.fread(f"../../micro_matrices/ses2_Label_after_exclusions.tsv").to_pandas()

print("These datasets should have the same number of rows:")
print(ids_ukb.shape)
print(demo_ukb.shape)
print(micro_ukb.shape)
print(label_ukb.shape)

# Remove all voxels with labels other than NAWM
micro_ukb[label_ukb != 8] = np.nan

# Load ADNI data

# demo: age and sex
demo_adni = dt.fread("../../../ADNI/Analyses/clean_data/results/df_demo.tsv").to_pandas()
demo_adni = demo_adni[['PTID', 'Year', 'age_int', "sex", "dx"]]
demo_adni.rename(columns={'PTID': 'ID', 'sex': 'Sex', 'age_int':'Age'}, inplace=True)
demo_adni['Sex'] = demo_adni['Sex'].replace({'Female': 0, 'Male': 1}) # 0 = Female, 1 = Male

# micro: microstructural values (subject by voxel in UKB space)
micro_adni = dt.fread(f"../../../ADNI/micro_matrices/{new_name}.tsv").to_pandas()
# label: BISON labels (subject by voxel in UKB space)
label_adni = dt.fread(f"../../../ADNI/micro_matrices/Label.tsv").to_pandas()

# Remove all voxels with labels other than NAWM
micro_adni[label_adni != 8] = np.nan

# Select cognitively unimpaired ADNI subjects <81 years old
adni_max_81yo = demo_adni.index[(demo_adni["Age"] <= 81)].tolist()
adni_max_81yo_hc = demo_adni.index[(demo_adni["Age"] <= 81) & (demo_adni["dx"] == "CN")].tolist()

demo_adni_adapt = demo_adni.iloc[adni_max_81yo_hc,:].reset_index(drop=True)
micro_adni_adapt = micro_adni.iloc[adni_max_81yo_hc,:].reset_index(drop=True)
label_adni_adapt = label_adni.iloc[adni_max_81yo_hc,:].reset_index(drop=True)

demo_adni_max81 = demo_adni.iloc[adni_max_81yo,:].reset_index(drop=True)
micro_adni_max81 = micro_adni.iloc[adni_max_81yo,:].reset_index(drop=True)
label_adni_max81 = label_adni.iloc[adni_max_81yo,:].reset_index(drop=True)

demo_adni = demo_adni[['ID', 'Year', 'Age', "Sex"]]
demo_adni_adapt = demo_adni_adapt[['ID', 'Year', 'Age', "Sex"]]
demo_adni_max81 = demo_adni_max81[['ID', 'Year', 'Age', "Sex"]]

print("These datasets should have the same number of rows:")
print(demo_adni_adapt.shape)
print(micro_adni_adapt.shape)
print(label_adni_adapt.shape)

print("These datasets should have the same number of rows:")
print(demo_adni_max81.shape)
print(micro_adni_max81.shape)
print(label_adni_max81.shape)

# Add intercept and B-spline with age to demo (add 5 to min and max, according to PCN tutorial at OHBM)
# https://github.com/predictive-clinical-neuroscience/NM_educational_OHBM24/blob/main/slot1_Fraza/1_fit_normative_models.ipynb)
B = create_bspline_basis(demo_ukb['Age'].min() - 5, demo_ukb['Age'].max() + 5, p=3, nknots = 4)

# UKB
Phi = np.array([B(i) for i in demo_ukb.iloc[:,2]])
intercept_pd = pd.DataFrame(np.ones((demo_ukb.shape[0],1)), columns = ["Intercept"])
Phi_pd = pd.DataFrame(Phi, columns=[f'Bspline_{i+1}' for i in range(Phi.shape[1])])
demo_ukb_pd = pd.concat([demo_ukb.iloc[:, 1:], intercept_pd, Phi_pd], axis=1)

# ADNI
Phi = np.array([B(i) for i in demo_adni.iloc[:,2]])
intercept_pd = pd.DataFrame(np.ones((demo_adni.shape[0],1)), columns = ["Intercept"])
Phi_pd = pd.DataFrame(Phi, columns=[f'Bspline_{i+1}' for i in range(Phi.shape[1])])
demo_adni_pd = pd.concat([demo_adni.iloc[:, 2:], intercept_pd, Phi_pd], axis=1)

demo_adni_pd_adapt = demo_adni_pd.iloc[adni_max_81yo_hc,:].reset_index(drop=True)
demo_adni_pd_max81 = demo_adni_pd.iloc[adni_max_81yo,:].reset_index(drop=True)

# All possibilities of age (integers) and sex, with intercept
possibilities_demo = list(product([0,1], list(range(demo_adni['Age'].min().astype(int), demo_adni['Age'].max().astype(int)+1))))
possibilities_demo = np.hstack((possibilities_demo, np.ones((len(possibilities_demo), 1), dtype=int)))
Phi_demo = np.array([B(i) for i in possibilities_demo[:,1]]) # Add 4th order B-spline to age
possibilities_demo = np.concatenate((possibilities_demo, Phi_demo), axis=1) # Add intercept column to iv

# UKB output file: fit metrics + age- and sex-specific means and SDs
metrics_columns = ['ROI', 'N', 'MSLL', 'EXPV', 'SMSE', 'RMSE', 'Rho', 'R2']
male_means_columns = [f'Male_mean_{i+1}' for i in range(int(demo_ukb['Age'].min())-1, int(demo_ukb['Age'].max()))]
male_sd_columns = [f'Male_sd_{i+1}' for i in range(int(demo_ukb['Age'].min())-1, int(demo_ukb['Age'].max()))]
female_means_columns = [f'Female_mean_{i+1}' for i in range(int(demo_ukb['Age'].min())-1, int(demo_ukb['Age'].max()))]
female_sd_columns = [f'Female_sd_{i+1}' for i in range(int(demo_ukb['Age'].min())-1, int(demo_ukb['Age'].max()))]
all_columns_ukb = metrics_columns + male_means_columns + male_sd_columns + female_means_columns + female_sd_columns

# ADNI output file: age- and sex-specific means and SDs
metrics_columns = ['ROI', 'N']
male_means_columns = [f'Male_mean_{i+1}' for i in range(int(demo_adni_max81['Age'].min())-1, int(demo_adni_max81['Age'].max()))]
male_sd_columns = [f'Male_sd_{i+1}' for i in range(int(demo_adni_max81['Age'].min())-1, int(demo_adni_max81['Age'].max()))]
female_means_columns = [f'Female_mean_{i+1}' for i in range(int(demo_adni_max81['Age'].min())-1, int(demo_adni_max81['Age'].max()))]
female_sd_columns = [f'Female_sd_{i+1}' for i in range(int(demo_adni_max81['Age'].min())-1, int(demo_adni_max81['Age'].max()))]
all_columns_adni = metrics_columns + male_means_columns + male_sd_columns + female_means_columns + female_sd_columns

# Iterate over voxels
def run_nm(i):
    print(i)
    # i = 23045
    os.makedirs(f'/dev/shm/parent41/nm_n{name}/tmp/v{i}', exist_ok=True)
    # Select voxel
    vox_ukb = micro_ukb.iloc[:,i]
    vox_adni_adapt = micro_adni_adapt.iloc[:, i]
    vox_adni_max81 = micro_adni_max81.iloc[:, i]
    # Count number of NAWM labels in voxel
    count_tissue_voxels_ukb = vox_ukb.shape[0] - (vox_ukb.isna().sum().sum())
    count_tissue_voxels_adni_adapt = vox_adni_adapt.shape[0] - (vox_adni_adapt.isna().sum().sum())
    print(f'Number of voxels for specific tissue type in voxel = {count_tissue_voxels_ukb} - {count_tissue_voxels_adni_adapt}')
    res_ukb = pd.DataFrame(np.zeros((1, len(all_columns_ukb))), columns=all_columns_ukb)
    res_adni = pd.DataFrame(np.zeros((1, len(all_columns_adni))), columns=all_columns_adni)
    # If less than 100 NAWM labels, return NaNs
    if count_tissue_voxels_ukb < 100 or count_tissue_voxels_adni_adapt < 5:
        # Results UKB
        res_ukb.iloc[0,:] = np.full(len(all_columns_ukb), np.nan)
        res_ukb['ROI'] = i
        res_ukb['N'] = count_tissue_voxels_ukb
        # Results ADNI
        res_adni.iloc[0,:] = np.full(len(all_columns_adni), np.nan)
        res_adni['ROI'] = i
        res_adni['N'] = count_tissue_voxels_ukb
    # If more or equal to 100 NAWM labels, run normative modeling
    if count_tissue_voxels_ukb >= 100 and count_tissue_voxels_adni_adapt >= 5:
        # Concatenate with demo
        vox_ukb = pd.concat([demo_ukb_pd, vox_ukb], axis=1)
        vox_adni_adapt = pd.concat([demo_adni_pd_adapt, vox_adni_adapt], axis=1)
        vox_adni_max81 = pd.concat([demo_adni_pd_max81, vox_adni_max81], axis=1)
        # Drow rows where not tissue type (NaN)
        vox_ukb = vox_ukb.dropna(subset=[vox_ukb.columns[-1]]).reset_index(drop=True)
        vox_adni_adapt = vox_adni_adapt.dropna(subset=[vox_adni_adapt.columns[-1]]).reset_index(drop=True)
        vox_adni_max81 = vox_adni_max81.dropna(subset=[vox_adni_max81.columns[-1]]).reset_index(drop=True)
        # Save UKB files
        with open(f'/dev/shm/parent41/nm_n{name}/tmp/v{i}/{name}_ukb_cov_bspline_{i}.pkl', 'wb') as file:
            pickle.dump(pd.DataFrame(vox_ukb.iloc[:,:-1].to_numpy(dtype=float)), file)
        with open(f'/dev/shm/parent41/nm_n{name}/tmp/v{i}/{name}_ukb_resp_{i}.pkl', 'wb') as file:
            pickle.dump(pd.DataFrame(vox_ukb.iloc[:,-1].to_numpy(dtype=float)), file)
        # Save "adaptation files" from ADNI
        with open(f'/dev/shm/parent41/nm_n{name}/tmp/v{i}/{name}_adni_adapt_cov_bspline_{i}.pkl', 'wb') as file:
            pickle.dump(pd.DataFrame(vox_adni_adapt.iloc[:,:-1].to_numpy(dtype=float)), file)
        with open(f'/dev/shm/parent41/nm_n{name}/tmp/v{i}/{name}_adni_adapt_resp_{i}.pkl', 'wb') as file:
            pickle.dump(pd.DataFrame(vox_adni_adapt.iloc[:,-1].to_numpy(dtype=float)), file)
        # Save max81 from ADNI
        with open(f'/dev/shm/parent41/nm_n{name}/tmp/v{i}/{name}_adni_max81_cov_bspline_{i}.pkl', 'wb') as file:
            pickle.dump(pd.DataFrame(vox_adni_max81.iloc[:,:-1].to_numpy(dtype=float)), file)
        with open(f'/dev/shm/parent41/nm_n{name}/tmp/v{i}/{name}_adni_max81_resp_{i}.pkl', 'wb') as file:
            pickle.dump(pd.DataFrame(vox_adni_max81.iloc[:,-1].to_numpy(dtype=float)), file)
        # Load files
        ukb_cov_file = os.path.join(f'/dev/shm/parent41/nm_n{name}/tmp/v{i}/{name}_ukb_cov_bspline_{i}.pkl')
        ukb_resp_file = os.path.join(f'/dev/shm/parent41/nm_n{name}/tmp/v{i}/{name}_ukb_resp_{i}.pkl')
        adni_adapt_cov_file = os.path.join(f'/dev/shm/parent41/nm_n{name}/tmp/v{i}/{name}_adni_adapt_cov_bspline_{i}.pkl')
        adni_adapt_resp_file = os.path.join(f'/dev/shm/parent41/nm_n{name}/tmp/v{i}/{name}_adni_adapt_resp_{i}.pkl')
        adni_max81_cov_file = os.path.join(f'/dev/shm/parent41/nm_n{name}/tmp/v{i}/{name}_adni_max81_cov_bspline_{i}.pkl')
        adni_max81_resp_file = os.path.join(f'/dev/shm/parent41/nm_n{name}/tmp/v{i}/{name}_adni_max81_resp_{i}.pkl')
        # Estimate model on UKB data
        yhat, s2, nm, Z, metrics_te = estimate(ukb_cov_file, ukb_resp_file, testresp=ukb_resp_file, testcov=ukb_cov_file, alg = 'blr', optimizer = 'powell', cvfolds = None, savemodel = True, saveoutput = False, standardize = False, outputsuffix=f"n{new_name}v{i}", dir_model=f'/dev/shm/parent41/nm_n{name}/tmp/v{i}')
        # Transfer model to ADNI
        yhat_adni, s2_adni, Z_adni = predict(covfile=adni_max81_cov_file, alg="blr", respfile=adni_max81_resp_file, model_path=f'/dev/shm/parent41/nm_n{name}/tmp/v{i}', inputsuffix=f"n{new_name}v{i}", adaptrespfile = adni_adapt_resp_file, adaptcovfile = adni_adapt_cov_file)
        # remove temporary files
        shutil.rmtree(f'/dev/shm/parent41/nm_n{name}/tmp/v{i}', ignore_errors=True)
        # UKB: Test predict for all possibilities of demographics
        demo_pred_ukb = pd.DataFrame(columns=['Sex', 'Age', 'Mean', 'Var', 'SD'])
        r=0
        for s in range(0,2):
            for a in range(int(demo_ukb['Age'].min()), int(demo_ukb['Age'].max()+1)):
                idx_age_sex = vox_ukb[(vox_ukb['Sex'] == s) & (vox_ukb['Age'] == a)].index.tolist()
                vec = [s, a, np.mean(yhat[idx_age_sex]), np.mean(s2[idx_age_sex]), np.sqrt(np.mean(s2[idx_age_sex]))]
                vec = pd.DataFrame([vec], columns=['Sex', 'Age', 'Mean', 'Var', 'SD'])
                demo_pred_ukb = pd.concat([demo_pred_ukb, vec], ignore_index=True)
                r = r+1
        # ADNI: Test predict for all possibilities of demographics
        demo_pred_adni = pd.DataFrame(columns=['Sex', 'Age', 'Mean', 'Var', 'SD'])
        r=0
        for s in range(0,2):
            for a in range(int(demo_adni_max81['Age'].min()), int(demo_adni_max81['Age'].max()+1)):
                if ((vox_adni_max81['Sex'] == s) & (vox_adni_max81['Age'] == a)).sum() > 0:
                    idx_age_sex = vox_adni_max81[(vox_adni_max81['Sex'] == s) & (vox_adni_max81['Age'] == a)].index.tolist()
                    vec = [s, a, np.mean(yhat_adni[idx_age_sex]), np.mean(s2_adni[idx_age_sex]), np.sqrt(np.mean(s2_adni[idx_age_sex]))]
                    vec = pd.DataFrame([vec], columns=['Sex', 'Age', 'Mean', 'Var', 'SD'])
                    demo_pred_adni = pd.concat([demo_pred_adni, vec], ignore_index=True)
                else:
                    vec = [s, a, np.nan, np.nan, np.nan]
                    vec = pd.DataFrame([vec], columns=['Sex', 'Age', 'Mean', 'Var', 'SD'])
                    demo_pred_adni = pd.concat([demo_pred_adni, vec], ignore_index=True)
                r = r+1
        # Save UKB results to dataframe
        for idx,c in enumerate(res_ukb.columns, start=0):
            # print(c)
            if c == 'ROI':
                res_ukb[c] = i
            if c == 'N':
                res_ukb[c] = count_tissue_voxels_ukb
            if c == 'MSLL':
                res_ukb[c] = metrics_te['MSLL'][0]
            if c == 'EXPV':
                res_ukb[c] = metrics_te['EXPV'][0]
            if c == 'SMSE':
                res_ukb[c] = metrics_te['SMSE'][0]
            if c == 'RMSE':
                res_ukb[c] = metrics_te['RMSE'][0]
            if c == 'Rho':
                res_ukb[c] = metrics_te['Rho'][0]
            if c == 'R2':
                res_ukb[c] = math.pow(metrics_te['Rho'][0],2)            
            if c.startswith("Male_"):
                age = int(res_ukb.columns[idx][-2:])
                sex = 1
                if 'mean' in c:
                    res_ukb[c] = demo_pred_ukb[(demo_pred_ukb['Sex'] == sex) & (demo_pred_ukb['Age'] == age)]['Mean'].values
                if 'sd' in c:
                    res_ukb[c] = demo_pred_ukb[(demo_pred_ukb['Sex'] == sex) & (demo_pred_ukb['Age'] == age)]['SD'].values
            elif c.startswith("Female"):
                age = int(res_ukb.columns[idx][-2:])
                sex = 0
                if 'mean' in c:
                    res_ukb[c] = demo_pred_ukb[(demo_pred_ukb['Sex'] == sex) & (demo_pred_ukb['Age'] == age)]['Mean'].values
                if 'sd' in c:
                    res_ukb[c] = demo_pred_ukb[(demo_pred_ukb['Sex'] == sex) & (demo_pred_ukb['Age'] == age)]['SD'].values
        # Save ADNI results to dataframe
        for idx,c in enumerate(res_adni.columns, start=0):
            # print(c)
            if c == 'ROI':
                res_adni[c] = i
            if c == 'N':
                res_adni[c] = count_tissue_voxels_ukb         
            if c.startswith("Male_"):
                age = int(res_adni.columns[idx][-2:])
                sex = 1
                if 'mean' in c:
                    res_adni[c] = demo_pred_adni[(demo_pred_adni['Sex'] == sex) & (demo_pred_adni['Age'] == age)]['Mean'].values
                if 'sd' in c:
                    res_adni[c] = demo_pred_adni[(demo_pred_adni['Sex'] == sex) & (demo_pred_adni['Age'] == age)]['SD'].values
            elif c.startswith("Female"):
                age = int(res_adni.columns[idx][-2:])
                sex = 0
                if 'mean' in c:
                    res_adni[c] = demo_pred_adni[(demo_pred_adni['Sex'] == sex) & (demo_pred_adni['Age'] == age)]['Mean'].values
                if 'sd' in c:
                    res_adni[c] = demo_pred_adni[(demo_pred_adni['Sex'] == sex) & (demo_pred_adni['Age'] == age)]['SD'].values
    return res_ukb, res_adni

# Run normative modeling for all voxels using parallelized processes
num_cpus = mp.cpu_count()
results_all = Parallel(n_jobs=int(num_cpus-10))(delayed(run_nm)(i) for i in range(micro_ukb.shape[1]))
res_ukb_all, res_adni_all = zip(*results_all)

# Concatenate results
results_ukb_all = pd.concat(res_ukb_all, ignore_index=True)
results_adni_all = pd.concat(res_adni_all, ignore_index=True)

# Save results (only one example included because too large otherwise)
# Performance metrics, predicted means and SDs for every combination of age and sex of both UKB and ADNI datasets
results_ukb_all.to_csv(f"./results/1a_norm_models_BLR/nm_{new_name}_ukb.tsv", sep='\t', index=False, na_rep="NA")
results_adni_all.to_csv(f"./results/1a_norm_models_BLR/nm_{new_name}_adni.tsv", sep='\t', index=False, na_rep="NA")

# Remove temporary directory
shutil.rmtree(f"/dev/shm/parent41/", ignore_errors=True)


