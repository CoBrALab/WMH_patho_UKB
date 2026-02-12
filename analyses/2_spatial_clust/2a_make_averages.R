
suppressPackageStartupMessages({
    library(RMINC)
    library(data.table)
    library(matrixStats)
    library(grid)
    library(gridExtra) 
    library(tidyverse)
    library(MRIcrotome)
    library(magrittr)
    library(viridis)
    library(RColorBrewer)
    library(ggplot2)
    library(dplyr)
    library(patchwork)
})

# Make between-subject averages of WMH zscore maps

dir.create("./results/2a_make_averages", showWarnings=FALSE)
dir.create("./visualization/2a_make_averages", showWarnings=FALSE)

names=c( "MD", "ISOVF", "FA", "ICVF", "OD", "T2star", "QSM")

# mask: where WMH prevalence >1
mask=mincGetVolume("../../data/WMH_mask.mnc")
# label: BISON labels (subject by voxel in UKB space)
label=as.data.frame(fread("../../../WMH_micro_spatial/micro_matrices/ses2_Label_after_exclusions.tsv"))

# Load WMH and NAWM prevalence
prev = as.data.frame(fread("../../../WMH_micro_spatial/Analyses_nm/tissue_prevalence/results/tissue_prevalence_after_exclusions_WMHmask_ses2.tsv"))
nawm_prev = as.numeric(prev[8,])
wmh_prev = as.numeric(prev[9,])

# Thresholds for minimum prevalence counts
WMH_MIN_PREVALENCE = 30
NAWM_MIN_PREVALENCE = 5000

wmh_thresh = WMH_MIN_PREVALENCE
nawm_thresh = NAWM_MIN_PREVALENCE

# Iterate for each microstructural marker
for (m in 1:length(names)) {
    print(names[m])

    # Load subject Z-scores
    zvals = as.data.frame(fread(paste0("../../../WMH_micro_spatial/Analyses_nm/norm_zscores_subjects/results/WMH_zscores_subj_",names[m],"_NAs.tsv")))
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
    outvol = mincGetVolume("../../data/UKB_template_T1_2mm.mnc")
    outvol[] = 0
    outvol[mask > 0.5] = avg
    mincWriteVolume(outvol, paste0("./results/2a_make_averages/WMH_zscore_", names[m],"_avg.mnc"), clobber=TRUE)

    # Averages with NAWM > 1000 and WMH > 1

    # Assign 0 to voxels with prevalence NAWM < 1000 and WMH < 1
    indices_above_thres = intersect(which(nawm_prev > 1000), which(wmh_prev > 1))
    logical_above_thres = rep(FALSE, length(nawm_prev))
    logical_above_thres = 1:length(logical_above_thres) %in% indices_above_thres

    # Voxel averages of zscores
    avg = colMeans(zvals, na.rm=TRUE)
    avg[is.na(avg)] = 0
    avg[logical_above_thres==FALSE] = 0

    # Write to tsv
    fwrite(as.data.frame(avg), paste0("./results/2a_make_averages/WMH_zscore_", names[m],"_avg_NAWM1000_WMH1.tsv"), row.names=FALSE, col.names=FALSE, sep='\t', quote=FALSE)

    # Write to mnc
    outvol = mincGetVolume("../../data/UKB_template_T1_2mm.mnc")
    outvol[] = 0
    outvol[mask > 0.5] = avg
    mincWriteVolume(outvol, paste0("./results/2a_make_averages/WMH_zscore_", names[m],"_avg_NAWM1000_WMH1.mnc"), clobber=TRUE)

    # Calculate standard error of the mean
    voxel_sem = colSds(as.matrix(zvals), na.rm = TRUE) / sqrt(colCounts(!is.na(as.matrix(zvals))))
    voxel_sem[is.na(voxel_sem)] = 0
    fwrite(as.data.frame(voxel_sem), paste0("./results/2a_make_averages/WMH_zscore_", names[m],"_sem.tsv"), row.names=FALSE, col.names=FALSE, sep='\t', quote=FALSE)

    outvol = mincGetVolume("../../data/UKB_template_T1_2mm.mnc")
    outvol[] = 0
    outvol[mask > 0.5] = voxel_sem
    mincWriteVolume(outvol, paste0("./results/2a_make_averages/WMH_zscore_", names[m],"_sem.mnc"), clobber=TRUE)

}

# Visualize

color_scale_div_2 = colorRampPalette(brewer.pal(9,"Blues"))(255)
color_scale_div_1 = colorRampPalette(brewer.pal(9,"Reds"))(255)

anatVol = mincArray(mincGetVolume("../../data/UKB_template_T1_2mm.mnc"))

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

# Averages with NAWM > 1000 and WMH > 1

png(file=paste0("./visualization/2a_make_averages/WMH_1_NAWM_1000.png"), width=7500, height=4000, pointsize = 150)
sliceSeries(nrow=4, ncol=1, begin=42, end=52, dimension=3) %>%
    addtitle(names[1]) %>%
    anatomy(anatVol, low=10, high=200) %>%
    overlay(mincArray(mincGetVolume(paste0("./results/2a_make_averages/WMH_zscore_", names[1],"_avg_NAWM1000_WMH1.mnc"))),
        low=0.015, high=6,col=color_scale_div_1, rCol=color_scale_div_2, symmetric=TRUE) %>%
sliceSeries() %>%
    addtitle(names[2]) %>%
    anatomy(anatVol, low=10, high=200) %>%
    overlay(mincArray(mincGetVolume(paste0("./results/2a_make_averages/WMH_zscore_", names[2],"_avg_NAWM1000_WMH1.mnc"))),
        low=0.015, high=6,col=color_scale_div_1, rCol=color_scale_div_2, symmetric=TRUE) %>%
sliceSeries() %>%
    addtitle(names[3]) %>%
    anatomy(anatVol, low=10, high=200) %>%
    overlay(mincArray(mincGetVolume(paste0("./results/2a_make_averages/WMH_zscore_", names[3],"_avg_NAWM1000_WMH1.mnc"))),
        low=0.015, high=6,col=color_scale_div_1, rCol=color_scale_div_2, symmetric=TRUE) %>%
sliceSeries() %>%
    addtitle(names[4]) %>%
    anatomy(anatVol, low=10, high=200) %>%
    overlay(mincArray(mincGetVolume(paste0("./results/2a_make_averages/WMH_zscore_", names[4],"_avg_NAWM1000_WMH1.mnc"))),
        low=0.015, high=6,col=color_scale_div_1, rCol=color_scale_div_2, symmetric=TRUE) %>%
sliceSeries() %>%
    addtitle(names[5]) %>%
    anatomy(anatVol, low=10, high=200) %>%
    overlay(mincArray(mincGetVolume(paste0("./results/2a_make_averages/WMH_zscore_", names[5],"_avg_NAWM1000_WMH1.mnc"))),
        low=0.015, high=6,col=color_scale_div_1, rCol=color_scale_div_2, symmetric=TRUE) %>%
sliceSeries() %>%
    addtitle(names[6]) %>%
    anatomy(anatVol, low=10, high=200) %>%
    overlay(mincArray(mincGetVolume(paste0("./results/2a_make_averages/WMH_zscore_", names[6],"_avg_NAWM1000_WMH1.mnc"))),
        low=0.015, high=6,col=color_scale_div_1, rCol=color_scale_div_2, symmetric=TRUE) %>%
sliceSeries() %>%
    addtitle(names[7]) %>%
    anatomy(anatVol, low=10, high=200) %>%
    overlay(mincArray(mincGetVolume(paste0("./results/2a_make_averages/WMH_zscore_", names[7],"_avg_NAWM1000_WMH1.mnc"))),
        low=0.015, high=6,col=color_scale_div_1, rCol=color_scale_div_2, symmetric=TRUE) %>%
legend("Z-values") %>%
draw()
dev.off()

# Visualize SEM

color_scale = viridis(100)
anatVol = mincArray(mincGetVolume("../../data/UKB_template_T1_2mm.mnc"))

png(file=paste0("./visualization/2a_make_averages/SEM.png"), width=7500, height=4000, pointsize = 150)
sliceSeries(nrow=4, ncol=1, begin=42, end=52, dimension=3) %>%
    addtitle(names[1]) %>%
    anatomy(anatVol, low=10, high=200) %>%
    overlay(mincArray(mincGetVolume(paste0("./results/2a_make_averages/WMH_zscore_", names[1],"_sem.mnc"))),
        low=0.015, high=3,col=color_scale) %>%
sliceSeries() %>%
    addtitle(names[2]) %>%
    anatomy(anatVol, low=10, high=200) %>%
    overlay(mincArray(mincGetVolume(paste0("./results/2a_make_averages/WMH_zscore_", names[2],"_sem.mnc"))),
        low=0.015, high=3,col=color_scale) %>%
sliceSeries() %>%
    addtitle(names[3]) %>%
    anatomy(anatVol, low=10, high=200) %>%
    overlay(mincArray(mincGetVolume(paste0("./results/2a_make_averages/WMH_zscore_", names[3],"_sem.mnc"))),
        low=0.015, high=3,col=color_scale) %>%
sliceSeries() %>%
    addtitle(names[4]) %>%
    anatomy(anatVol, low=10, high=200) %>%
    overlay(mincArray(mincGetVolume(paste0("./results/2a_make_averages/WMH_zscore_", names[4],"_sem.mnc"))),
        low=0.015, high=3,col=color_scale) %>%
sliceSeries() %>%
    addtitle(names[5]) %>%
    anatomy(anatVol, low=10, high=200) %>%
    overlay(mincArray(mincGetVolume(paste0("./results/2a_make_averages/WMH_zscore_", names[5],"_sem.mnc"))),
        low=0.015, high=3,col=color_scale) %>%
sliceSeries() %>%
    addtitle(names[6]) %>%
    anatomy(anatVol, low=10, high=200) %>%
    overlay(mincArray(mincGetVolume(paste0("./results/2a_make_averages/WMH_zscore_", names[6],"_sem.mnc"))),
        low=0.015, high=3,col=color_scale) %>%
sliceSeries() %>%
    addtitle(names[7]) %>%
    anatomy(anatVol, low=10, high=200) %>%
    overlay(mincArray(mincGetVolume(paste0("./results/2a_make_averages/WMH_zscore_", names[7],"_sem.mnc"))),
        low=0.015, high=3,col=color_scale) %>%
legend("Z-values") %>%
draw()
dev.off()

# SEM by voxel prevalence graphs

results_sem = list()
for (m in 1:length(names)) {
    results_sem[[m]] = as.data.frame(fread(paste0("./results/2a_make_averages/WMH_zscore_", names[m],"_sem.tsv")))
    colnames(results_sem[[m]])[1] = "SEM"
    results_sem[[m]]$Marker = names[m]
    results_sem[[m]]$wmh_prev = wmh_prev
    results_sem[[m]]$nawm_prev = nawm_prev
}
results_sem = do.call(rbind, results_sem)
results_sem$Marker = factor(results_sem$Marker, levels=names)

results_sem = results_sem %>% filter(SEM > 0)

color_scale = c(MD = "#04319E", ISOVF = "#8298CF", FA = "#2B520B", ICVF="#448312", OD="#A2C189", T2star="#D36108", QSM="#E9B084")

plt_wmh = ggplot(results_sem, aes(x=wmh_prev, y=SEM, color=Marker, group=Marker)) +
    # geom_point() + 
    geom_smooth() +
    geom_vline(xintercept = 30) +
    scale_x_continuous(name = "WMH prevalence", limits=c(0,100)) +
    scale_y_continuous(name = "SEM") +
    scale_color_manual(name="", values=color_scale) +
    coord_cartesian(xlim = c(0, 100), ylim = c(0, 2)) +
    theme_classic() +
    theme(text=element_text(size=30, color="black"), plot.title=element_text(hjust=0.5, size=40),
      legend.position = "none")

plt_nawm = ggplot(results_sem, aes(x=nawm_prev, y=SEM, color=Marker, group=Marker)) +
    # geom_point() + 
    geom_smooth() +
    geom_vline(xintercept = 5000) +
    scale_x_continuous(name = "NAWM prevalence", limits=c(0,10000)) +
    scale_y_continuous(name = "SEM") +
    scale_color_manual(name="", values=color_scale) +
    coord_cartesian(xlim = c(0, 10000), ylim = c(0, 0.5)) +
    theme_classic() +
    theme(text=element_text(size=30, color="black"), plot.title=element_text(hjust=0.5, size=40),
      legend.position = "none")

wrap_plots(plt_wmh, plt_nawm, nrow=1)

ggsave("./visualization/2a_make_averages/SEM_by_prev.png", width = 14, height=7)
