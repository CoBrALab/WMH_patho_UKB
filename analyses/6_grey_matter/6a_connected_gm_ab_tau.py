
import neuromaps.datasets as nmd
import neuromaps.images as nmi
import neuromaps.transforms as nmt
import neuromaps.plotting as nmp
import neuromaps.stats as nms
import neuromaps.nulls as nmn
from neuromaps.parcellate import Parcellater
import matplotlib.pyplot as plt
from nilearn.image import smooth_img
import pandas as pd
import numpy as np
import nibabel as nib
import copy
import glob
import os
from joblib import Parallel, delayed
import multiprocessing as mp
import pickle
from pathlib import Path
import datatable as dt

# Compare cortical connectivity profiles of WMH clusters with AB and Tau distributions

os.makedirs(f'./results/6a_connected_gm_ab_tau', exist_ok=True)
os.makedirs(f'./visualization/6a_connected_gm_ab_tau', exist_ok=True)

# Load volumetric parcellation
parc_dkt_sym = nib.load("../../cerebra_parc/CerebrA_cortex.nii")
parc_dkt_asym = nib.load("../../cerebra_parc/CerebrA_cortex_asym.nii")

# Load cortical parcellation for visualization purposes (DKT parcellation is almost identical to CerebrA)
parc_dkt_fsa = nmi.annot_to_gifti((f'../../dkt_atlas/lh.aparc.annot', f'../../dkt_atlas/rh.aparc.annot'))
parcels = np.concatenate([np.unique(parc_dkt_fsa[0].darrays[0].data), np.unique(parc_dkt_fsa[1].darrays[0].data)])
parcels = parcels[parcels != 0]

# Load CerebrA mapping
cerebra_mapping = dt.fread("../../cerebra_parc/CerebrA_LabelDetails_modified.csv").to_pandas()
cerebra_mapping_lr = pd.concat([cerebra_mapping, cerebra_mapping], axis=0)
cerebra_mapping_filter = cerebra_mapping_lr[cerebra_mapping_lr['Included'] == True]

parcels_cerebra = np.unique(np.sort(np.concatenate([cerebra_mapping_filter['RH Label'], cerebra_mapping_filter['LH Labels']])))

# Load DKT34 region names and mapping
dkt_mapping = dt.fread("../../dkt_atlas/dkt_name_mapping.csv").to_pandas()
dkt_mapping_filter = dkt_mapping[dkt_mapping['dkt_surf_34_label'].notna()]

# Function to sample average values of map within volumetric parcels
# Applied to group averages of AB and Tau PET
def proc_ad_map(file, parc_vol, parc_surf, name, name_parc, viz):
    print(name)
    # Load data
    print("Load data")
    out_map = {}
    out_map['name']=name
    out_map['img'] = nib.load(file)
    # Prep surface parcellation
    out_map['surf_parc'] = copy.deepcopy(parc_surf)
    out_map['surf_parc'][0].darrays[0].data = out_map['surf_parc'][0].darrays[0].data.astype(float)
    out_map['surf_parc'][1].darrays[0].data = out_map['surf_parc'][1].darrays[0].data.astype(float)
    # Parcellate data
    print("Parcellate data")
    out_map['list_parc'] = []
    out_map['list_dkt'] = np.full((parcels.shape), np.nan)
    # For every parcel
    for p in range(0, parcels_cerebra.shape[0]):
        # Calculate mean in parcel in volume space
        mean_parcel = np.mean(out_map['img'].get_fdata()[parc_vol.get_fdata() == parcels_cerebra[p]]) if np.any(parc_vol.get_fdata() == parcels_cerebra[p]) else 0
        out_map['list_parc'].append(mean_parcel)
        # Right hemisphere
        if parcels_cerebra[p] < 52:
            # Get corresponding DKT34 region(s)
            dkt_labels_str = cerebra_mapping[cerebra_mapping['RH Label'] == parcels_cerebra[p]]['dkt_cortex_label'].to_list()
            dkt_labels = [int(label.strip()) for label in dkt_labels_str[0].split(',')]
            # Assign to surface
            for p_dkt in dkt_labels:
                mask = parc_dkt_fsa[1].agg_data() == p_dkt + 34
                out_map['surf_parc'][1].darrays[0].data[mask] = mean_parcel
                out_map['list_dkt'][p_dkt-1+34] = mean_parcel
        # Left hemisphere
        if parcels_cerebra[p] >= 52:
            # Get corresponding DKT34 region(s)
            dkt_labels_str = cerebra_mapping[cerebra_mapping['LH Labels'] == parcels_cerebra[p]]['dkt_cortex_label'].to_list()
            dkt_labels = [int(label.strip()) for label in dkt_labels_str[0].split(',')]
            # Assign to surface
            for p_dkt in dkt_labels:
                mask = parc_dkt_fsa[0].agg_data() == p_dkt
                out_map['surf_parc'][0].darrays[0].data[mask] = mean_parcel
                out_map['list_dkt'][p_dkt-1] = mean_parcel
    # Visualize
    if viz == True:
        print("Visualize")
        abs_max = max(abs(out_map['surf_parc'][0].darrays[0].data.max()),
                    abs(out_map['surf_parc'][0].darrays[0].data.min()),
                    abs(out_map['surf_parc'][1].darrays[0].data.max()),
                    abs(out_map['surf_parc'][1].darrays[0].data.min()))
        nmp.plot_surf_template(out_map['surf_parc'], 'fsaverage', '164k')
        plt.savefig(f"./visualization/6a_connected_gm_ab_tau/{name}_{name_parc}.png",bbox_inches = "tight", dpi=300)
        plt.clf()
        plt.close()
        print(f"./visualization/6a_connected_gm_ab_tau/{name}_{name_parc}.png")
    return out_map

# Run for AB and Tau PET averages from PREVENT-AD
other_maps = {}

files_pad = glob.glob('./data/prevent_ad/*_mni152sym09c.nii.gz', recursive=True)
names = [os.path.splitext(os.path.splitext(os.path.basename(f))[0])[0] for f in files_pad]

results = Parallel(n_jobs=-1)(delayed(proc_ad_map)(files_pad[f], parc_dkt_sym, parc_dkt_fsa, f'pad_{names[f]}', 'dkt', viz=False) for f in range(0, len(files_pad)))

for i in range(len(results)):
    other_maps[results[i]['name']] = results[i]

# Visualize cortical connectivity profiles of WMH clusters

# Load connected GM regions
conn_gm = {}
conn_gm['clust1'] = {}
conn_gm['clust2'] = {}
conn_gm['clust3'] = {}

conn_gm['clust1']['sums'] = np.loadtxt(f"./data/conn_gm_cerebra_dkt/sum_clust1.txt").tolist()
conn_gm['clust2']['sums'] = np.loadtxt(f"./data/conn_gm_cerebra_dkt/sum_clust2.txt").tolist()
conn_gm['clust3']['sums'] = np.loadtxt(f"./data/conn_gm_cerebra_dkt/sum_clust3.txt").tolist()

# Normalize
total_conn_gm = np.sum(np.vstack([value['sums'] for value in conn_gm.values()]), axis=0)

conn_gm['clust1']['sums_norm'] = conn_gm['clust1']['sums'] / total_conn_gm
conn_gm['clust2']['sums_norm'] = conn_gm['clust2']['sums'] / total_conn_gm
conn_gm['clust3']['sums_norm'] = conn_gm['clust3']['sums'] / total_conn_gm

# Parcellate
for c in [1,2,3]:
    conn_gm[f'clust{c}']['surf_norm'] = copy.deepcopy(parc_dkt_fsa)
    conn_gm[f'clust{c}']['surf_norm'][0].darrays[0].data = conn_gm[f'clust{c}']['surf_norm'][0].darrays[0].data.astype(float)
    conn_gm[f'clust{c}']['surf_norm'][1].darrays[0].data = conn_gm[f'clust{c}']['surf_norm'][1].darrays[0].data.astype(float)
    conn_gm[f'clust{c}']['list_dkt_norm'] = np.full((parcels.shape), np.nan)
    # Replace labels by values
    for p in range(0, 51):
        if cerebra_mapping.iloc[p, :]['Included'] == True:
            print(p)
            # Get indices for cerebra
            cerebra_right = cerebra_mapping.iloc[p, :]['RH Label'] - 1
            cerebra_left = cerebra_mapping.iloc[p, :]['LH Labels'] - 1
            # Get corresponding DKT34 region(s)
            dkt_labels_str = cerebra_mapping.iloc[p, :]['dkt_cortex_label']
            dkt_labels = [int(label.strip()) for label in dkt_labels_str.split(',')]
            # Get mean in parcel
            mean_parcel_left = conn_gm[f'clust{c}']['sums'][cerebra_left]
            mean_parcel_norm_left = conn_gm[f'clust{c}']['sums_norm'][cerebra_left]
            mean_parcel_right = conn_gm[f'clust{c}']['sums'][cerebra_right]
            mean_parcel_norm_right = conn_gm[f'clust{c}']['sums_norm'][cerebra_right]
            # Assign to surface
            for p_dkt in dkt_labels:
                mask_left = parc_dkt_fsa[0].agg_data() == p_dkt
                mask_right = parc_dkt_fsa[1].agg_data() == p_dkt+34
                conn_gm[f'clust{c}']['surf'][0].darrays[0].data[mask_left] = mean_parcel_left
                conn_gm[f'clust{c}']['surf'][1].darrays[0].data[mask_right] = mean_parcel_right
                conn_gm[f'clust{c}']['surf_norm'][0].darrays[0].data[mask_left] = mean_parcel_norm_left
                conn_gm[f'clust{c}']['surf_norm'][1].darrays[0].data[mask_right] = mean_parcel_norm_right
                conn_gm[f'clust{c}']['list_dkt_norm'][int(p_dkt)-1] = mean_parcel_norm_left
                conn_gm[f'clust{c}']['list_dkt_norm'][int(p_dkt)-1+34] = mean_parcel_norm_right

# Plot
for c in [1,2,3]:
    # Norm
    nmp.plot_surf_template(conn_gm[f'clust{c}']['surf_norm'], 'fsaverage', '164k')
    plt.savefig(f"./visualization/6a_connected_gm_ab_tau/connected_gm_clust{c}_norm.png",bbox_inches = "tight", dpi=300)
    plt.clf()
    plt.close()
    print(f"./visualization/6a_connected_gm_ab_tau/connected_gm_clust{c}_norm.png")

# Spin tests: spatial permutation-based null that handles spatial autocorrelation of cortical GM maps

# Define function
def spin_test(df):
    nulls = nmn.alexander_bloch(data=df, atlas='fsaverage', density='164k',
                n_perm=10000, seed=np.random.choice(10000), parcellation=parc_dkt_fsa)
    return nulls

# Run for other maps
results = Parallel(n_jobs=-1)(delayed(spin_test)(other_maps[list(other_maps.keys())[f]]['list_dkt']) for f in range(0, len(list(other_maps.keys()))))

for f in range(0, len(list(other_maps.keys()))):
    other_maps[list(other_maps.keys())[f]]['nulls_parc'] = results[f]

del results

# Run for connected GM maps
for i, name_i in enumerate(conn_gm.keys()):
    print(f'{i} {name_i}')
    conn_gm[name_i]['nulls_parc'] = spin_test(conn_gm[name_i]['list_dkt_norm'])

# Relate maps together

# Define function
def spin_test_corr(map1, name):
    results_corr = pd.DataFrame(columns=['x', 'y', 'pearson_corr', 'pearson_p'])
    # Run correlations
    for c in range(0,3):
        p_corr, p_pval = nms.compare_images(map1, conn_gm[f'clust{c+1}']['list_dkt_norm'][0:], metric="pearsonr", nulls=other_maps[name]['nulls_parc'])
        results_corr.loc[len(results_corr)] = [name, f'clust{c+1}', p_corr, p_pval]
    return results_corr

# Run for all other maps
other_maps_type = list(other_maps.keys())

results = Parallel(n_jobs=-1)(delayed(spin_test_corr)(other_maps[list(other_maps.keys())[f]]['list_dkt'], other_maps[list(other_maps.keys())[f]]['name']) for f in range(0, len(other_maps_type)))

# Save results
results = pd.concat(results, ignore_index=True)
results.to_csv("./results/6a_connected_gm_ab_tau/corr_between_maps.csv", index=False)


