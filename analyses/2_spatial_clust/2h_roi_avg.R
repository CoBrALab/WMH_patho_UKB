
library(data.table)
library(RMINC)
library(matrixStats)

# Calculate ROI-based statistics per subject
# Average WMH pathophysiology per spatial region

# Individual-level data cannot be shared

dir.create("./results/2h_roi_avg", showWarnings=FALSE)
dir.create("./visualization/2h_roi_avg", showWarnings=FALSE)

names=c("FA", "MD", "ICVF", "ISOVF", "OD", "T2star", "QSM")

mask = mincGetVolume("../../data/WMH_mask.mnc")

# Load parcellation

parc = fread("./results/2d_assign_excluded_voxels/final_parc_k3.tsv")
parc = parc$V1

# Calculate median value of WMH patho in each cluster

# Load label and volumes
label = as.data.frame(fread(paste0("../../micro_matrices/ses2_Label_after_exclusions.tsv"), header=FALSE))
volumes = as.data.frame(fread("../../../UKB/BISON/BISON_volumes.csv"))
inclusions = fread("../../QC/inclusions_final.txt")
inclusions = inclusions$V1

colnames(volumes)[2] = "tp"

volumes = subset(volumes, tp == ses)
volumes = volumes[volumes$ID %in% inclusions, ]

# Initialize results dataframe
results = volumes

# Calculate median in WMHs for each cluster
for (m in 1:length(names)) {
    cat(paste0("\tMicro = ", names[m], "\n"))

    # Load micro for specific ses and remove non-WMH values
    zvals = as.data.frame(fread(paste0("../norm_zscores_subjects/results/WMH_zscores_subj_",names[m],"_NAs.tsv")))
    zvals[!(label == 9)] = NA

    for (i in sort(unique(parc[parc > 0]))) {
        cat(paste0("\t\t\tCluster = ", i, "\n"))

        # Select voxels in cluster
        zvals_clust = subset(zvals, select=which(parc==i))
        # Calculate median for each subject
        med = rowMedians(as.matrix(zvals_clust), na.rm=TRUE)
        results = cbind(results, med)
        colnames(results)[ncol(results)] = paste0("final_parc_k3_", names[m], "_c",i)
        rm(zvals_clust)
        rm(med)

    }

    rm(zvals)
}

# Calculate parcellated WMH volumes in common space
for (i in sort(unique(parc[parc > 0]))) {
    label_clust = subset(label, select=which(parc==i))
    # Count number of WMH voxels in spatial region
    count = rowSums(label_clust == 9, na.rm = TRUE)
    results = cbind(results, count)
    colnames(results)[ncol(results)] = paste0("final_parc_k3_WMHvol_c",i)
}

# Write to tsv
fwrite(results, paste0("./results/2h_roi_avg/subj_WMH_patho.tsv"), col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

