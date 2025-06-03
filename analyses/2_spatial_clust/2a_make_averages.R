
library(RMINC)
library(data.table)

# Make between-subject averages of WMH zscore maps

dir.create("./results/2a_make_averages", showWarnings=FALSE)
dir.create("./visualization/2a_make_averages", showWarnings=FALSE)

names=c( "MD", "ISOVF", "FA", "ICVF", "OD", "T2star", "QSM")

# mask: where WMH prevalence >1
mask=mincGetVolume("../../data/WMH_mask.mnc")
# label: BISON labels (subject by voxel in UKB space)
label=as.data.frame(fread("../../../micro_matrices/ses2_Label_after_exclusions.tsv"))

# Load WMH and NAWM prevalence
prev = as.data.frame(fread("../../../Analyses_nm/tissue_prevalence/results/tissue_prevalence_after_exclusions_WMHmask_ses2.tsv"))
nawm_prev = as.numeric(prev[8,])
wmh_prev = as.numeric(prev[9,])

wmh_thresh = 30
nawm_thresh = 5000

# Iterate for each microstructural marker
for (m in 1:length(names)) {
    print(names[m])

    # Load subject Z-scores
    zvals = as.data.frame(fread(paste0("../1_norm_mod/results/1c_zscore_UKB/WMH_zscores_UKB_",names[m],"_NAs.tsv")))
    zvals[!(label == 9)] = NA

    # Assign 0 to voxels with prevalence NAWM < 5000 and WMH < 30
    indices_above_thres = intersect(which(nawm_prev > nawm_thresh), which(wmh_prev > wmh_thresh))
    logical_above_thres = rep(FALSE, length(nawm_prev))
    logical_above_thres = 1:length(logical_above_thres) %in% indices_above_thres

    # Voxel averages of zscores
    avg = colMeans(zvals, na.rm=TRUE)
    avg[is.na(avg)] = 0
    avg[logical_above_thres==FALSE] = 0

    # Write to tsv
    fwrite(as.data.frame(avg), paste0("./results/2a_make_averages/WMH_zscore_", names[m],"_avg.tsv"), row.names=FALSE, col.names=FALSE, sep='\t', quote=FALSE)

    # Write to mnc
    outvol = mincGetVolume("../../data/UKB_template_2mm.mnc")
    outvol[] = 0
    outvol[mask > 0.5] = avg
    mincWriteVolume(outvol, paste0("./results/2a_make_averages/WMH_zscore_", names[m],"_avg.mnc"), clobber=TRUE)
}


# Visualize

library(grid)
library(gridExtra) 
library(tidyverse)
library(MRIcrotome)
library(magrittr)
library(viridis)
library(RColorBrewer)

color_scale_div_2 = colorRampPalette(brewer.pal(9,"Blues"))(255)
color_scale_div_1 = colorRampPalette(brewer.pal(9,"Reds"))(255)

anatVol = mincArray(mincGetVolume("../../data/UKB_template_2mm.mnc"))

png(file=paste0("./visualization/2a_make_averages/WMH_30_NAWM_5000.png"), width=7500, height=4000, pointsize = 150)
sliceSeries(nrow=4, ncol=1, begin=42, end=52, dimension=3) %>%
    addtitle(names[1]) %>%
    anatomy(anatVol, low=10, high=200) %>%
    overlay(mincArray(mincGetVolume(paste0("./results/2a_make_averages/WMH_zscore_", names[1],"_avg.mnc"))),
        low=0.015, high=6,col=color_scale_div_1, rCol=color_scale_div_2, symmetric=TRUE) %>%
sliceSeries() %>%
    addtitle(names[2]) %>%
    anatomy(anatVol, low=10, high=200) %>%
    overlay(mincArray(mincGetVolume(paste0("./results/2a_make_averages/WMH_zscore_", names[2],"_avg.mnc"))),
        low=0.015, high=6,col=color_scale_div_1, rCol=color_scale_div_2, symmetric=TRUE) %>%
sliceSeries() %>%
    addtitle(names[3]) %>%
    anatomy(anatVol, low=10, high=200) %>%
    overlay(mincArray(mincGetVolume(paste0("./results/2a_make_averages/WMH_zscore_", names[3],"_avg.mnc"))),
        low=0.015, high=6,col=color_scale_div_1, rCol=color_scale_div_2, symmetric=TRUE) %>%
sliceSeries() %>%
    addtitle(names[4]) %>%
    anatomy(anatVol, low=10, high=200) %>%
    overlay(mincArray(mincGetVolume(paste0("./results/2a_make_averages/WMH_zscore_", names[4],"_avg.mnc"))),
        low=0.015, high=6,col=color_scale_div_1, rCol=color_scale_div_2, symmetric=TRUE) %>%
sliceSeries() %>%
    addtitle(names[5]) %>%
    anatomy(anatVol, low=10, high=200) %>%
    overlay(mincArray(mincGetVolume(paste0("./results/2a_make_averages/WMH_zscore_", names[5],"_avg.mnc"))),
        low=0.015, high=6,col=color_scale_div_1, rCol=color_scale_div_2, symmetric=TRUE) %>%
sliceSeries() %>%
    addtitle(names[6]) %>%
    anatomy(anatVol, low=10, high=200) %>%
    overlay(mincArray(mincGetVolume(paste0("./results/2a_make_averages/WMH_zscore_", names[6],"_avg.mnc"))),
        low=0.015, high=6,col=color_scale_div_1, rCol=color_scale_div_2, symmetric=TRUE) %>%
sliceSeries() %>%
    addtitle(names[7]) %>%
    anatomy(anatVol, low=10, high=200) %>%
    overlay(mincArray(mincGetVolume(paste0("./results/2a_make_averages/WMH_zscore_", names[7],"_avg.mnc"))),
        low=0.015, high=6,col=color_scale_div_1, rCol=color_scale_div_2, symmetric=TRUE) %>%
legend("Z-values") %>%
draw()
dev.off()


