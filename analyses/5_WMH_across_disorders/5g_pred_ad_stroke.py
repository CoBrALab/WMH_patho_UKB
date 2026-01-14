
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import StratifiedKFold, cross_validate, train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import (
    classification_report, confusion_matrix, roc_auc_score, 
    roc_curve, precision_recall_curve, average_precision_score,
    make_scorer, f1_score
)
from sklearn.utils.class_weight import compute_class_weight
import warnings
import datatable as dt
from neuroHarmonize import harmonizationLearn, harmonizationApply
from sklearn.metrics import accuracy_score, precision_score, recall_score, balanced_accuracy_score

warnings.filterwarnings('ignore')

# Set random seed for reproducibility
np.random.seed(42)

#region Load data for spect_clust_k3 and ComBat harmonization with ADNI

# Load UKB stroke data
df_ukb = pd.read_csv(f"../3_temporal_sustain/results/3f_clean_wmh_data/wmh_combined.tsv", delimiter="\t")
df_ukb = df_ukb.drop(columns=['stage_c1', 'stage_c2', 'stage_c3', 'post_ant_subtype', 'post_ant_stage', 'T2star_c1', 'T2star_c2', 'T2star_c3', 'QSM_c1', 'QSM_c2', 'QSM_c3'])

df_dx = dt.fread("../../../UKB/Analyses/clean_firstocc/results/firstocc_categ.tsv").to_pandas()
df_dx = df_dx[df_dx["chapter"].isin([
    "Endocrine, nutritional and metabolic diseases",
    "Mental, Behavioral and Neurodevelopmental disorders",
    "Diseases of the circulatory system",
    "Diseases of the nervous system"
])]
df_dx = df_dx.merge(df_ukb['ID'], how="right")
controls = df_dx[df_dx.isna().any(axis=1)]["ID"].unique()
df_ukb_combat = df_ukb[df_ukb["ID"].isin(controls)].reset_index(drop=True)

stroke_cases = (df_dx.loc[df_dx["icd_code"].isin(["I60", "I61", "I63", "I64"]), ["ID"]].drop_duplicates().reset_index(drop=True))
dementia_cases = (df_dx.loc[df_dx["icd_code"].isin(["F00","F03","G30"]), ["ID"]].drop_duplicates().reset_index(drop=True))

overlap_ids = set(stroke_cases['ID']).intersection(dementia_cases['ID'])
stroke_cases = stroke_cases[~stroke_cases['ID'].isin(overlap_ids)].reset_index(drop=True)
dementia_cases = dementia_cases[~dementia_cases['ID'].isin(overlap_ids)].reset_index(drop=True)

df_stroke = pd.merge(stroke_cases, df_ukb)
df_stroke = df_stroke.drop(columns=['ID'])
df_stroke['Group'] = 'Stroke'
df_stroke['SITE'] = 'UKB'

df_dementia_ukb = pd.merge(dementia_cases, df_ukb)
df_dementia_ukb = df_dementia_ukb.drop(columns=['ID'])
df_dementia_ukb['Group'] = 'cog_imp'
df_dementia_ukb['SITE'] = 'UKB'

# Load ADNI data
df_adni = pd.read_csv(f"../clean_wmh_data/results/spect_clust_k3/wmh_combined_ADNI.tsv", delimiter="\t")
df_adni = df_adni.rename(columns={"Age": "age", "Sex": "sex"})

df_demo = pd.read_csv("../../../ADNI/Analyses/clean_data/results/df_demo.tsv", delimiter="\t")
df_demo = df_demo.rename(columns={"PTID": "ID", "EXAMDATE": "date"})
# Keep last timepoint only
df_demo = (
    df_demo
    .sort_values(['ID'])      # optional unless order matters within ID
    .groupby('ID', as_index=False)
    .tail(1)
    .reset_index(drop=True)
)
# Remove CI due to non AD
df_demo = df_demo[
    (df_demo['dx_mci_cause'].isna() | (df_demo['dx_mci_cause'] != "Due to other")) &
    (df_demo['dx_dementia_cause'].isna() | (df_demo['dx_dementia_cause'] != "Due to other"))
]
# Merge with WMH data
df_demo_merge = df_demo.loc[:,['ID', 'date']]
df_adni = pd.merge(df_demo_merge, df_adni, on=['ID', 'date'])
df_adni = df_adni.dropna().reset_index(drop=True)
df_demo = pd.merge(df_demo, df_adni.loc[:,['ID', 'date']], on=['ID', 'date'])

# Keep only MCI or dementia
df_demo_merge = df_demo[df_demo['dx'].isin(['MCI', 'Dementia'])].reset_index(drop=True)
df_demo_merge = df_demo_merge.loc[:,["ID", "date"]]

df_ci_adni = pd.merge(df_demo_merge, df_adni, how="left", on=['ID', 'date'])
df_ci_adni['Group'] = 'cog_imp'
df_ci_adni['SITE'] = 'ADNI'

df_ci_adni = df_ci_adni.loc[:,df_stroke.columns]

# Select healthy IDs from ADNI
hc_ids = df_demo[(df_demo["dx"] == "CN")]['ID']
df_adni_combat = df_adni[df_adni["ID"].isin(hc_ids)].reset_index(drop=True)
df_adni_combat = df_adni_combat.drop(columns=["date", "Year"])

# ComBat harmonization
df_adni_combat['SITE'] = 'ADNI'
df_ukb_combat['SITE'] = 'UKB'
df_combat = pd.concat([df_adni_combat, df_ukb_combat], ignore_index=True)
df_combat['sex'] = df_combat['sex'].map({'Male': 1, 'Female': 0}).values

exclude_cols = ["ID", "age", "sex", "SITE"]
feature_cols = [c for c in df_combat.columns if c not in exclude_cols]
X_combat = df_combat[feature_cols].to_numpy()
covars = df_combat.loc[:,['SITE', 'age', 'sex']]

# Train
combat_model, df_combat_adj = harmonizationLearn(X_combat, covars, ref_batch="UKB")

# Predict
df_all = pd.concat([df_stroke, df_dementia_ukb, df_ci_adni], axis=0, ignore_index=True)
df_all['sex'] = df_all['sex'].map({'Male': 1, 'Female': 0}).values

exclude_cols = ["age", "sex", 'Group', "SITE"]
feature_cols = [c for c in df_all.columns if c not in exclude_cols]
X_combat = df_all[feature_cols].to_numpy()
covars = df_all.loc[:,['SITE', 'age', 'sex']]

X_clust = harmonizationApply(X_combat, covars, combat_model)
X_clust_vol = X_clust[:,0:3]
y_clust = df_all['Group'].map({'Stroke': 0, 'cog_imp': 1}).values

#endregion

#region Load data for PV/deep and ComBat harmonization with ADNI

# Load UKB stroke data
df_ukb = pd.read_csv(f"../clean_wmh_data/results/pv_deep/wmh_combined.tsv", delimiter="\t")
df_ukb = df_ukb.drop(columns=['pv_deep_T2star_c1', 'pv_deep_T2star_c2', 'pv_deep_QSM_c1', 'pv_deep_QSM_c2'])

df_dx = dt.fread("../../../UKB/Analyses/clean_firstocc/results/firstocc_categ.tsv").to_pandas()
df_dx = df_dx[df_dx["chapter"].isin([
    "Endocrine, nutritional and metabolic diseases",
    "Mental, Behavioral and Neurodevelopmental disorders",
    "Diseases of the circulatory system",
    "Diseases of the nervous system"
])]
df_dx = df_dx.merge(df_ukb['ID'], how="right")
controls = df_dx[df_dx.isna().any(axis=1)]["ID"].unique()
df_ukb_combat = df_ukb[df_ukb["ID"].isin(controls)].reset_index(drop=True)

stroke_cases = (df_dx.loc[df_dx["icd_code"].isin(["I60", "I61", "I63", "I64"]), ["ID"]].drop_duplicates().reset_index(drop=True))
dementia_cases = (df_dx.loc[df_dx["icd_code"].isin(["F00","F03","G30"]), ["ID"]].drop_duplicates().reset_index(drop=True))

overlap_ids = set(stroke_cases['ID']).intersection(dementia_cases['ID'])
stroke_cases = stroke_cases[~stroke_cases['ID'].isin(overlap_ids)].reset_index(drop=True)
dementia_cases = dementia_cases[~dementia_cases['ID'].isin(overlap_ids)].reset_index(drop=True)

df_stroke = pd.merge(stroke_cases, df_ukb)
df_stroke = df_stroke.drop(columns=['ID'])
df_stroke['Group'] = 'Stroke'
df_stroke['SITE'] = 'UKB'

df_dementia_ukb = pd.merge(dementia_cases, df_ukb)
df_dementia_ukb = df_dementia_ukb.drop(columns=['ID'])
df_dementia_ukb['Group'] = 'cog_imp'
df_dementia_ukb['SITE'] = 'UKB'

# Load ADNI data
df_adni = pd.read_csv(f"../clean_wmh_data/results/pv_deep/wmh_combined_ADNI.tsv", delimiter="\t")
df_adni = df_adni.rename(columns={"Age": "age", "Sex": "sex"})

df_demo = pd.read_csv("../../../ADNI/Analyses/clean_data/results/df_demo.tsv", delimiter="\t")
df_demo = df_demo.rename(columns={"PTID": "ID", "EXAMDATE": "date"})
# Keep last timepoint only
df_demo = (
    df_demo
    .sort_values(['ID'])      # optional unless order matters within ID
    .groupby('ID', as_index=False)
    .tail(1)
    .reset_index(drop=True)
)
# Remove CI due to non AD
df_demo = df_demo[
    (df_demo['dx_mci_cause'].isna() | (df_demo['dx_mci_cause'] != "Due to other")) &
    (df_demo['dx_dementia_cause'].isna() | (df_demo['dx_dementia_cause'] != "Due to other"))
]
# Merge with WMH data
df_demo_merge = df_demo.loc[:,['ID', 'date']]
df_adni = pd.merge(df_demo_merge, df_adni, on=['ID', 'date'])
df_adni = df_adni.dropna().reset_index(drop=True)
df_demo = pd.merge(df_demo, df_adni.loc[:,['ID', 'date']], on=['ID', 'date'])

# Keep only MCI or dementia
df_demo_merge = df_demo[df_demo['dx'].isin(['MCI', 'Dementia'])].reset_index(drop=True)
df_demo_merge = df_demo_merge.loc[:,["ID", "date"]]

df_ci_adni = pd.merge(df_demo_merge, df_adni, how="left", on=['ID', 'date'])
df_ci_adni['Group'] = 'cog_imp'
df_ci_adni['SITE'] = 'ADNI'

df_ci_adni = df_ci_adni.loc[:,df_stroke.columns]

# Select healthy IDs from ADNI
hc_ids = df_demo[(df_demo["dx"] == "CN")]['ID']
df_adni_combat = df_adni[df_adni["ID"].isin(hc_ids)].reset_index(drop=True)
df_adni_combat = df_adni_combat.drop(columns=["date", "Year"])

# ComBat harmonization
df_adni_combat['SITE'] = 'ADNI'
df_ukb_combat['SITE'] = 'UKB'
df_combat = pd.concat([df_adni_combat, df_ukb_combat], ignore_index=True)
df_combat['sex'] = df_combat['sex'].map({'Male': 1, 'Female': 0}).values

exclude_cols = ["ID", "age", "sex", "SITE"]
feature_cols = [c for c in df_combat.columns if c not in exclude_cols]
X_combat = df_combat[feature_cols].to_numpy()
covars = df_combat.loc[:,['SITE', 'age', 'sex']]

# Train
combat_model, df_combat_adj = harmonizationLearn(X_combat, covars, ref_batch="UKB")

# Predict
df_all = pd.concat([df_stroke, df_dementia_ukb, df_ci_adni], axis=0, ignore_index=True)
df_all['sex'] = df_all['sex'].map({'Male': 1, 'Female': 0}).values

exclude_cols = ["age", "sex", 'Group', "SITE"]
feature_cols = [c for c in df_all.columns if c not in exclude_cols]
X_combat = df_all[feature_cols].to_numpy()
covars = df_all.loc[:,['SITE', 'age', 'sex']]

X_pvdeep = harmonizationApply(X_combat, covars, combat_model)
X_pvdeep_vol = X_pvdeep[:,0:2]
y_pvdeep = df_all['Group'].map({'Stroke': 0, 'cog_imp': 1}).values

#endregion

#region Load data for no parcellation and ComBat harmonization with ADNI

# Load UKB stroke data
df_ukb = pd.read_csv(f"../clean_wmh_data/results/no_clust/wmh_combined.tsv", delimiter="\t")
df_ukb = df_ukb.drop(columns=['stage_c1', 'T2star_c1', 'QSM_c1'])

df_dx = dt.fread("../../../UKB/Analyses/clean_firstocc/results/firstocc_categ.tsv").to_pandas()
df_dx = df_dx[df_dx["chapter"].isin([
    "Endocrine, nutritional and metabolic diseases",
    "Mental, Behavioral and Neurodevelopmental disorders",
    "Diseases of the circulatory system",
    "Diseases of the nervous system"
])]
df_dx = df_dx.merge(df_ukb['ID'], how="right")
controls = df_dx[df_dx.isna().any(axis=1)]["ID"].unique()
df_ukb_combat = df_ukb[df_ukb["ID"].isin(controls)].reset_index(drop=True)

stroke_cases = (df_dx.loc[df_dx["icd_code"].isin(["I60", "I61", "I63", "I64"]), ["ID"]].drop_duplicates().reset_index(drop=True))
dementia_cases = (df_dx.loc[df_dx["icd_code"].isin(["F00","F03","G30"]), ["ID"]].drop_duplicates().reset_index(drop=True))

overlap_ids = set(stroke_cases['ID']).intersection(dementia_cases['ID'])
stroke_cases = stroke_cases[~stroke_cases['ID'].isin(overlap_ids)].reset_index(drop=True)
dementia_cases = dementia_cases[~dementia_cases['ID'].isin(overlap_ids)].reset_index(drop=True)

df_stroke = pd.merge(stroke_cases, df_ukb)
df_stroke = df_stroke.drop(columns=['ID'])
df_stroke['Group'] = 'Stroke'
df_stroke['SITE'] = 'UKB'

df_dementia_ukb = pd.merge(dementia_cases, df_ukb)
df_dementia_ukb = df_dementia_ukb.drop(columns=['ID'])
df_dementia_ukb['Group'] = 'cog_imp'
df_dementia_ukb['SITE'] = 'UKB'

# Load ADNI data
df_adni = pd.read_csv(f"../clean_wmh_data/results/no_clust/wmh_combined_ADNI.tsv", delimiter="\t")
df_adni = df_adni.rename(columns={"Age": "age", "Sex": "sex", "region_vol[, seq(2, ncol(region_vol))]": "WMHvol_c1"})

df_demo = pd.read_csv("../../../ADNI/Analyses/clean_data/results/df_demo.tsv", delimiter="\t")
df_demo = df_demo.rename(columns={"PTID": "ID", "EXAMDATE": "date"})
# Keep last timepoint only
df_demo = (
    df_demo
    .sort_values(['ID'])      # optional unless order matters within ID
    .groupby('ID', as_index=False)
    .tail(1)
    .reset_index(drop=True)
)
# Remove CI due to non AD
df_demo = df_demo[
    (df_demo['dx_mci_cause'].isna() | (df_demo['dx_mci_cause'] != "Due to other")) &
    (df_demo['dx_dementia_cause'].isna() | (df_demo['dx_dementia_cause'] != "Due to other"))
]
# Merge with WMH data
df_demo_merge = df_demo.loc[:,['ID', 'date']]
df_adni = pd.merge(df_demo_merge, df_adni, on=['ID', 'date'])
df_adni = df_adni.dropna().reset_index(drop=True)
df_demo = pd.merge(df_demo, df_adni.loc[:,['ID', 'date']], on=['ID', 'date'])

# Keep only MCI or dementia
df_demo_merge = df_demo[df_demo['dx'].isin(['MCI', 'Dementia'])].reset_index(drop=True)
df_demo_merge = df_demo_merge.loc[:,["ID", "date"]]

df_ci_adni = pd.merge(df_demo_merge, df_adni, how="left", on=['ID', 'date'])
df_ci_adni['Group'] = 'cog_imp'
df_ci_adni['SITE'] = 'ADNI'

df_ci_adni = df_ci_adni.loc[:,df_stroke.columns]

# Select healthy IDs from ADNI
hc_ids = df_demo[(df_demo["dx"] == "CN")]['ID']
df_adni_combat = df_adni[df_adni["ID"].isin(hc_ids)].reset_index(drop=True)
df_adni_combat = df_adni_combat.drop(columns=["date", "Year"])

# ComBat harmonization
df_adni_combat['SITE'] = 'ADNI'
df_ukb_combat['SITE'] = 'UKB'
df_combat = pd.concat([df_adni_combat, df_ukb_combat], ignore_index=True)
df_combat['sex'] = df_combat['sex'].map({'Male': 1, 'Female': 0}).values

exclude_cols = ["ID", "age", "sex", "SITE"]
feature_cols = [c for c in df_combat.columns if c not in exclude_cols]
X_combat = df_combat[feature_cols].to_numpy()
covars = df_combat.loc[:,['SITE', 'age', 'sex']]

# Train
combat_model, df_combat_adj = harmonizationLearn(X_combat, covars, ref_batch="UKB")

# Predict
df_all = pd.concat([df_stroke, df_dementia_ukb, df_ci_adni], axis=0, ignore_index=True)
df_all['sex'] = df_all['sex'].map({'Male': 1, 'Female': 0}).values

exclude_cols = ["age", "sex", 'Group', "SITE"]
feature_cols = [c for c in df_all.columns if c not in exclude_cols]
X_combat = df_all[feature_cols].to_numpy()
covars = df_all.loc[:,['SITE', 'age', 'sex']]

X_noparc = harmonizationApply(X_combat, covars, combat_model)
X_noparc_vol = X_noparc[:,0]
y_noparc = df_all['Group'].map({'Stroke': 0, 'cog_imp': 1}).values

#endregion

#region Function to do model fit and CV

def fit_repeated_cv(X, y, n_repeats=100, cv=5):
    """
    Perform repeated stratified k-fold cross-validation
    
    Parameters:
    -----------
    X : numpy array
        Feature matrix
    y : numpy array
        Target vector
    n_repeats : int
        Number of times to repeat cross-validation
    cv : int
        Number of folds per repeat
    model_type : str
        'lasso' for Lasso Logistic Regression or 'rf' for Random Forest
    """
    from sklearn.metrics import accuracy_score, precision_score, recall_score, balanced_accuracy_score
    from sklearn.ensemble import RandomForestClassifier
    # Store metrics from all repeats
    all_repeat_metrics = {
        'balanced_accuracy': [],
        'f1': [], 'roc_auc': []
    }
    # Store predictions from first repeat only (for visualization)
    first_repeat_y_true = []
    first_repeat_y_pred = []
    first_repeat_y_pred_proba = []
    print("\n" + "="*60)
    model_name = "LASSO LOGISTIC REGRESSION"
    print(f"REPEATED CROSS-VALIDATION - {model_name}")
    print(f"({n_repeats} x {cv}-Fold Stratified)")
    print("="*60)
    print(f"Total model fits: {n_repeats * cv} = {n_repeats * cv}")
    # Compute class weights
    class_weights = compute_class_weight('balanced', 
                                         classes=np.unique(y), 
                                         y=y)
    class_weight_dict = {0: class_weights[0], 1: class_weights[1]}
    print(f"Class weights for imbalance: {class_weight_dict}")
    for repeat in range(n_repeats):
        # Use different random state for each repeat
        skf = StratifiedKFold(n_splits=cv, shuffle=True, random_state=repeat)
        repeat_metrics = {
            'balanced_accuracy': [],
            'accuracy': [], 'precision': [], 'recall': [], 
            'f1': [], 'roc_auc': []
        }
        for fold, (train_idx, test_idx) in enumerate(skf.split(X, y), 1):
            # Split data
            X_train_fold, X_test_fold = X[train_idx], X[test_idx]
            y_train_fold, y_test_fold = y[train_idx], y[test_idx]
            # Standardize features (fit on train, transform both)
            scaler = StandardScaler()
            X_train_scaled = scaler.fit_transform(X_train_fold)
            X_test_scaled = scaler.transform(X_test_fold)
            # Train model
            model = LogisticRegression(
                penalty='l1',
                C=1.0,
                solver='saga',
                class_weight='balanced',  # HANDLES CLASS IMBALANCE
                max_iter=1000,
                random_state=42,
                n_jobs=-1
            )
            model.fit(X_train_scaled, y_train_fold)
            # Predict
            y_pred = model.predict(X_test_scaled)
            y_pred_proba = model.predict_proba(X_test_scaled)[:, 1]
            # Store predictions from first repeat only
            if repeat == 0:
                first_repeat_y_true.extend(y_test_fold)
                first_repeat_y_pred.extend(y_pred)
                first_repeat_y_pred_proba.extend(y_pred_proba)
            # Calculate metrics for this fold
            repeat_metrics['balanced_accuracy'].append(balanced_accuracy_score(y_test_fold, y_pred))
            repeat_metrics['f1'].append(f1_score(y_test_fold, y_pred, zero_division=0))
            repeat_metrics['roc_auc'].append(roc_auc_score(y_test_fold, y_pred_proba))
        # Store mean metrics for this repeat
        for metric in repeat_metrics:
            all_repeat_metrics[metric].append(np.mean(repeat_metrics[metric]))
        # Progress update every 10 repeats
        if (repeat + 1) % 10 == 0:
            current_auc = np.mean(all_repeat_metrics['roc_auc'])
            print(f"  Completed {repeat + 1}/{n_repeats} repeats... "
                  f"Mean ROC-AUC: {current_auc:.3f}")
    # Convert first repeat predictions to arrays
    first_repeat_y_true = np.array(first_repeat_y_true)
    first_repeat_y_pred = np.array(first_repeat_y_pred)
    first_repeat_y_pred_proba = np.array(first_repeat_y_pred_proba)
    # Print summary statistics across all repeats
    print("\n" + "-"*60)
    print(f"SUMMARY ACROSS {n_repeats} REPEATS")
    print("-"*60)
    for metric in ['balanced_accuracy', 'f1', 'roc_auc']:
        scores = all_repeat_metrics[metric]
        print(f"{metric.upper()}:")
        print(f"  Mean: {np.mean(scores):.4f}")
        print(f"  Std:  {np.std(scores):.4f}")
        print(f"  95% CI: [{np.percentile(scores, 2.5):.4f}, {np.percentile(scores, 97.5):.4f}]")
    return (first_repeat_y_true, first_repeat_y_pred, first_repeat_y_pred_proba, 
            all_repeat_metrics)

#endregion

#region Run

################# Spect_clust_k3

# Volumes + micro
y_true, y_pred, y_pred_proba, fold_metrics = fit_repeated_cv(X_clust, y_clust, cv=5, model_type='lasso')
df_metrics = pd.DataFrame(fold_metrics)
df_metrics.to_csv("./results/5g_pred_ad_stroke/spect_clust_k3_micro_metrics.csv", index=False)

# Volumes only
y_true, y_pred, y_pred_proba, fold_metrics = fit_repeated_cv(X_clust_vol, y_clust, cv=5, model_type='lasso')
df_metrics = pd.DataFrame(fold_metrics)
df_metrics.to_csv("./results/5g_pred_ad_stroke/spect_clust_k3_vol_metrics.csv", index=False)

################# PV/deep

# Volumes + micro
y_true, y_pred, y_pred_proba, fold_metrics = fit_repeated_cv(X_pvdeep, y_pvdeep, cv=5, model_type='lasso')
df_metrics = pd.DataFrame(fold_metrics)
df_metrics.to_csv("./results/5g_pred_ad_stroke/pv_deep_micro_metrics.csv", index=False)

# Volumes only
y_true, y_pred, y_pred_proba, fold_metrics = fit_repeated_cv(X_pvdeep_vol, y_pvdeep, cv=5, model_type='lasso')
df_metrics = pd.DataFrame(fold_metrics)
df_metrics.to_csv("./results/5g_pred_ad_stroke/pv_deep_vol_metrics.csv", index=False)

################# No parcellation

# Volumes + micro
y_true, y_pred, y_pred_proba, fold_metrics = fit_repeated_cv(X_noparc, y_noparc, cv=5, model_type='lasso')
df_metrics = pd.DataFrame(fold_metrics)
df_metrics.to_csv("./results/5g_pred_ad_stroke/noparc_micro_metrics.csv", index=False)

# Volumes only
y_true, y_pred, y_pred_proba, fold_metrics = fit_repeated_cv(X_noparc_vol.reshape(-1,1), y_noparc, cv=5)
df_metrics = pd.DataFrame(fold_metrics)
df_metrics.to_csv("./results/5g_pred_ad_stroke/noparc_vol_metrics.csv", index=False)

#endregion
