
library(data.table)
library(RMINC)
library(tidyverse)
library(MRIcrotome)
library(viridis)

# Assign voxels that were previously excluded due to low NAWM or WMH prevalence
# to closest spatial cluster using a search area strategy

dir.create("./results/2d_assign_excluded_voxels", showWarnings=FALSE)
dir.create("./visualization/2d_assign_excluded_voxels", showWarnings=FALSE)

# Load WMH and NAWM prevalence
wmh_thresh = 30
nawm_thresh = 5000

mask=mincGetVolume("../../data/WMH_mask.mnc")
prev = as.data.frame(fread("../tissue_prevalence/results/tissue_prevalence_after_exclusions_WMHmask_ses2.tsv"))
nawm_prev = as.numeric(prev[8,])
wmh_prev = as.numeric(prev[9,])

indices_above_thres = intersect(which(nawm_prev > nawm_thresh), which(wmh_prev > wmh_thresh))
logical_above_thres = rep(FALSE, length(nawm_prev))
logical_above_thres = 1:length(logical_above_thres) %in% indices_above_thres

indices_below_thres = union(which(nawm_prev <= nawm_thresh), which(wmh_prev <= wmh_thresh))
logical_below_thres = rep(FALSE, length(nawm_prev))
logical_below_thres = 1:length(logical_below_thres) %in% indices_below_thres

# Visualize left-out voxels
outvol <- mincGetVolume("../../data/UKB_template_2mm.mnc")
outvol[] <- 0
outvol[mask > 0.5] <- logical_below_thres
mincWriteVolume(outvol, paste0("./results/2d_assign_excluded_voxels/voxels_to_assign_clust.mnc"), clobber=TRUE)

anatVol = mincArray(mincGetVolume("../../data/UKB_template_2mm.mnc"))

png(file=paste0("./visualization/2d_assign_excluded_voxels/voxels_to_assign_clust.png"), width=3000, height=3000, pointsize = 80)
sliceSeries(nrow=3, ncol=3, begin=40, end=52, dimension=3) %>%
    #addtitle("On average") %>%
    anatomy(anatVol, low=10, high=140) %>%
    overlay(mincArray(mincGetVolume(paste0("./results/2d_assign_excluded_voxels/voxels_to_assign_clust.mnc"))),
            low=0, high=1, col="red") %>%
    # legend("Cluster") %>%
    draw(layout="row")
dev.off()

# Strategy: in a region of x voxels surrounding the voxel to assign, assign to the cluster that appears most often
# Repeat by iteratively increasing x

library(foreach)
library(doParallel)

dir.create("tmp", showWarnings=FALSE)

# parallelize process
cores = detectCores()
parallelCluster <- makeCluster(cores-5, type = "SOCK", methods = FALSE)
setDefaultCluster(parallelCluster)
registerDoParallel(parallelCluster)

# First: dilate voxels by x voxels (from 1 to 7)
write_voxel_dilated = function(vox_idx) {
    voxel_mask = rep(0, length(logical_below_thres))
    voxel_mask[vox_idx] = 1

    # Write voxel to minc
    outvol <- mincGetVolume("../../data/UKB_template_2mm.mnc")
    outvol[] <- 0
    outvol[mask > 0.5] <- voxel_mask
    RMINC::mincWriteVolume(outvol, paste0("./tmp/vox_to_assign_",vox_idx,".mnc"), clobber=TRUE)

    # Dilate mask by 1-5 voxels
    command=paste0("mincmorph -clobber -3D26 -successive D ./tmp/vox_to_assign_",vox_idx,".mnc ./tmp/vox_to_assign_",vox_idx,"_dil_1.mnc")
    system(command)
    command=paste0("mincmorph -clobber -3D26 -successive DD ./tmp/vox_to_assign_",vox_idx,".mnc ./tmp/vox_to_assign_",vox_idx,"_dil_2.mnc")
    system(command)
    command=paste0("mincmorph -clobber -3D26 -successive DDD ./tmp/vox_to_assign_",vox_idx,".mnc ./tmp/vox_to_assign_",vox_idx,"_dil_3.mnc")
    system(command)
    command=paste0("mincmorph -clobber -3D26 -successive DDDD ./tmp/vox_to_assign_",vox_idx,".mnc ./tmp/vox_to_assign_",vox_idx,"_dil_4.mnc")
    system(command)
    command=paste0("mincmorph -clobber -3D26 -successive DDDDD ./tmp/vox_to_assign_",vox_idx,".mnc ./tmp/vox_to_assign_",vox_idx,"_dil_5.mnc")
    system(command)
    command=paste0("mincmorph -clobber -3D26 -successive DDDDDD ./tmp/vox_to_assign_",vox_idx,".mnc ./tmp/vox_to_assign_",vox_idx,"_dil_6.mnc")
    system(command)
    command=paste0("mincmorph -clobber -3D26 -successive DDDDDDD ./tmp/vox_to_assign_",vox_idx,".mnc ./tmp/vox_to_assign_",vox_idx,"_dil_7.mnc")
    system(command)
}

foreach(i=indices_below_thres, .packages='RMINC') %dopar% {
    write_voxel_dilated(i)
}

# Second: Assign voxels in ventricles to c1 (otherwise under-represented due to lower total amount of voxels)
vent_label = 9
clust_in_vent = 1

bison_ukb = mincGetVolume("../../../UKB/temporary_template/bison_ukbb/RF0_ukbb_bison_Label_2mm.mnc")[mask>0.5]
prev_parc = mincGetVolume("./results/2c_viz_clust/k3_all_clusters.mnc")[mask>0.5]
prev_parc[which(bison_ukb == vent_label)] = clust_in_vent

# Third: get highest prevalence cluster in surrounding area, and repeat until all remaining voxels are assigned to a cluster
highest_prev_clust = function(vox_idx, parc, dil) {
    # Calculate highest prevalence of clusters in surrounding voxels
    mask=mincGetVolume("../../data/WMH_mask.mnc")
    vox_dil = mincGetVolume(paste0("./tmp/vox_to_assign_",vox_idx,"_dil_",dil,".mnc"))[mask>0.5]
    surrounding_clusts = parc[which(vox_dil>0)]
    mode_value = as.numeric(names(table(surrounding_clusts))[which.max(table(surrounding_clusts))])

    if (length(mode_value) == 0) {mode_value=0}

    return(mode_value)
}

iter_assign_clust = function(left_to_assign, prev_parc, d, parc) {
    mask=mincGetVolume("../../data/WMH_mask.mnc")
    assigned_clusts = foreach(i=left_to_assign, .export = "highest_prev_clust", .packages='RMINC', .combine=c) %dopar% {
        highest_prev_clust(i, prev_parc, d)
    }

    # New parcellation: Added voxels that were assigned a cluster
    new_parc = rep(0, length(logical_above_thres))
    new_parc[prev_parc>0] = prev_parc[prev_parc>0]
    new_parc[left_to_assign] = assigned_clusts

    # Write to mnc and visualize
    outvol <- mincGetVolume("../../data/UKB_template_2mm.mnc")
    outvol[] <- 0
    outvol[mask > 0.5] <- new_parc
    mincWriteVolume(outvol, paste0("./tmp/new_parc_d",d,"_",parc,".mnc"), clobber=TRUE)

    color_scale = rainbow(10)
    color_scale = color_scale[1:4]

    # Calculate number of voxels left to assign
    count_left_to_assign = sum(new_parc == 0)

    print(paste0("Dilating ",d," voxels; Left to assign after ",parc," attempts = ", count_left_to_assign))
    prev_parc = mincGetVolume(paste0("./tmp/new_parc_d",d,"_",parc,".mnc"))[mask>0.5]
    left_to_assign = which(prev_parc==0)

    ret = list(count_left_to_assign, left_to_assign, prev_parc)

    return(ret)
}

count_left_to_assign = sum(logical_above_thres==FALSE)
left_to_assign = indices_below_thres
parc=1
d=1

while (length(left_to_assign) > 0) {

    iter_parc = iter_assign_clust(left_to_assign, prev_parc, d, parc)

    count_left_to_assign_old = count_left_to_assign
    count_left_to_assign = iter_parc[[1]]
    left_to_assign = iter_parc[[2]]
    prev_parc = iter_parc[[3]]

    parc = parc + 1

    if (count_left_to_assign == count_left_to_assign_old & d <= 7) { d = d+1 }
    if (count_left_to_assign == count_left_to_assign_old & d > 7) { break }

}

# Export results
final_parc = prev_parc
fread(final_parc, "./results/2d_assign_excluded_voxels/final_parc_k3.tsv", row.names=FALSE, col.names=TRUE, sep="\t")

# Write to mnc
outvol <- mincGetVolume("../../data/UKB_template_2mm.mnc")
outvol[] <- 0
outvol[mask > 0.5] <- final_parc
mincWriteVolume(outvol, paste0("./results/2d_assign_excluded_voxels/final_parc_k3.mnc"), clobber=TRUE)

# Visualize final parcellation

color_scale = viridis(3)

png(file=paste0("./visualization/2d_assign_excluded_voxels/final_parc_k3.png"), width=3000, height=3000, pointsize = 80)
sliceSeries(nrow=3, ncol=3, begin=40, end=52, dimension=3) %>%
    #addtitle("On average") %>%
    anatomy(mincArray(mincGetVolume("../../data/UKB_template_2mm.mnc")), low=10, high=140) %>%
    overlay(mincArray(mincGetVolume(paste0("./results/2d_assign_excluded_voxels/final_parc_k3.mnc"))),
            low=0, high=4, col=color_scale) %>%
    # legend("Cluster") %>%
    draw(layout="row")
dev.off()
