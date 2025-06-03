
library(RMINC)
library(data.table)
library(ggplot2)

# Visualize the WMH z-scores of 10 subjects with highest WMH volumes in UKB
# Cannot share the data due to confidentiality of individual participants

dir.create("./results/1e_viz_zscores_subjects", showWarnings=FALSE)
dir.create("./visualization/1e_viz_zscores_subjects", showWarnings=FALSE)

# Load WMH volumes
vol = fread("../../../UKB/BISON/BISON_volumes.csv")
vol = subset(as.data.frame(vol), ses==2, select=c("ID", "WMH"))

# Load subject IDs for z-score matrix (+ test if all lists are the same - should be)

names=c("MD", "ISOVF", "FA", "ICVF", "OD", "T2star", "QSM")
names_long = c("dti_MD", "NODDI_ISOVF", "dti_FA", "NODDI_ICVF", "NODDI_OD", "T2star", "QSM")

# Merge with volumes
# inclusions: final sample of UKB inclusions (32,526 subjects)
df = as.data.frame(fread("../../QC/inclusions_final.txt"))
colnames(df) = c("ID")
df = merge(df, vol, by="ID", all=FALSE)

# Select 10 subject IDs with highest WMH volumes
df_order = df[order(df$WMH, decreasing = TRUE), ]
ids_to_viz = df_order$ID[seq(1,10)]

# Make minc files out of WMH z-scores for selected subjects

mask=mincGetVolume("../../data/WMH_mask.mnc")
label=as.data.frame(fread("../../micro_matrices/ses2_Label_after_exclusions.tsv"))

# Make individual minc files for WMH z-scores AND raw micro in WMHs only
# ids: final sample of UKB inclusions (32,526 subjects)
ids = fread(paste0("../../QC/inclusions_final.txt"))
ids = as.numeric(ids$V1)

for(m in 1:length(names)) {
    print(names[m])

    # Load z-score matrix and IDs
    zscores = as.data.frame(fread(paste0("./results/1c_zscore_UKB/WMH_zscores_UKB_",names[m],".tsv")))
    # Remove labels other than WMH
    zscores[!(label == 9)] = 0

    for(i in 1:length(ids_to_viz)) {
        print(ids_to_viz[i])
        
        #### For zscores

        # Select z-score row for ID
        index = which(ids == ids_to_viz[i])
        zscores_id = as.numeric(zscores[index,])

        # Write to mnc
        outvol = mincGetVolume("../../data/UKB_template_2mm.mnc")
        outvol[] = 0
        outvol[mask > 0.5] = zscores_id
        mincWriteVolume(outvol, paste0("./results/1e_viz_zscores_subjects/",ids_to_viz[i],"_",names[m],"_WMH_zscore.mnc"), clobber=TRUE, like=outvol)

        #### For raw micro

        # Add 1000 to QSM values so that all values are positive
        # Mask out non-WMH voxels
        if (names[m] == "QSM") {
            command = paste0("mincmath -clobber -add -const 1000 ../../../UKB/maps_UKB_space/sub-",ids_to_viz[i],"_ses-2_",names_long[m],"_UKB.mnc ./results/1e_viz_zscores_subjects/",ids_to_viz[i],"_",names[m],"_plus1000.mnc")
            system(command)

            command = paste0("minccalc -clobber -expression \"(A[0]==9)?A[1]:0\" ../../../UKB/maps_UKB_space/sub-",ids_to_viz[i],"_ses-2_Label_UKB.mnc ./results/1e_viz_zscores_subjects/",ids_to_viz[i],"_",names[m],"_plus1000.mnc ./results/1e_viz_zscores_subjects/",ids_to_viz[i],"_",names[m],"_WMH_raw.mnc")
            system(command)

        } else {
            command = paste0("minccalc -clobber -expression \"(A[0]==9)?A[1]:0\" ../../../UKB/maps_UKB_space/sub-",ids_to_viz[i],"_ses-2_Label_UKB.mnc ../../../UKB/maps_UKB_space/sub-",ids_to_viz[i],"_ses-2_",names_long[m],"_UKB.mnc ./results/1e_viz_zscores_subjects/",ids_to_viz[i],"_",names[m],"_WMH_raw.mnc")
            system(command)
        }
    }
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
color_scale_raw = viridis(255, option = "magma")

# Visualize subject WMH z-scores overlayed on their FLAIRs

for(i in 1:length(ids_to_viz)) {
    anatVol = mincArray(mincGetVolume(paste0("../../../UKB/maps_UKB_space/sub-",ids_to_viz[i],"_ses-2_flair_UKB.mnc")))
    micro_subj = list()
    micro_subj[[1]] = mincArray(mincGetVolume(paste0("./results/1e_viz_zscores_subjects/",ids_to_viz[i],"_",names[1],"_WMH_zscore.mnc")))
    micro_subj[[2]] = mincArray(mincGetVolume(paste0("./results/1e_viz_zscores_subjects/",ids_to_viz[i],"_",names[2],"_WMH_zscore.mnc")))
    micro_subj[[3]] = mincArray(mincGetVolume(paste0("./results/1e_viz_zscores_subjects/",ids_to_viz[i],"_",names[3],"_WMH_zscore.mnc")))
    micro_subj[[4]] = mincArray(mincGetVolume(paste0("./results/1e_viz_zscores_subjects/",ids_to_viz[i],"_",names[4],"_WMH_zscore.mnc")))
    micro_subj[[5]] = mincArray(mincGetVolume(paste0("./results/1e_viz_zscores_subjects/",ids_to_viz[i],"_",names[5],"_WMH_zscore.mnc")))
    micro_subj[[6]] = mincArray(mincGetVolume(paste0("./results/1e_viz_zscores_subjects/",ids_to_viz[i],"_",names[6],"_WMH_zscore.mnc")))
    micro_subj[[7]] = mincArray(mincGetVolume(paste0("./results/1e_viz_zscores_subjects/",ids_to_viz[i],"_",names[7],"_WMH_zscore.mnc")))

    png(file=paste0("./visualization/1e_viz_zscores_subjects/",ids_to_viz[i],"_WMH_zscore.png"), width=8500, height=4000, pointsize = 150)
    sliceSeries(nrow=4, ncol=1, begin=42, end=52, dimension=3) %>%
        addtitle("FLAIR") %>%
        anatomy(anatVol, low=10, high=200) %>%
    sliceSeries(nrow=4, ncol=1, begin=42, end=52, dimension=3) %>%
        addtitle(names[1]) %>%
        anatomy(anatVol, low=10, high=200) %>%
        overlay(micro_subj[[1]], low=0.015, high=6,col=color_scale_div_1, rCol=color_scale_div_2, symmetric=TRUE) %>%
    sliceSeries(nrow=4, ncol=1, begin=42, end=52, dimension=3) %>%
        addtitle(names[2]) %>%
        anatomy(anatVol, low=10, high=200) %>%
        overlay(micro_subj[[2]], low=0.015, high=6,col=color_scale_div_1, rCol=color_scale_div_2, symmetric=TRUE) %>%
    sliceSeries(nrow=4, ncol=1, begin=42, end=52, dimension=3) %>%
        addtitle(names[3]) %>%
        anatomy(anatVol, low=10, high=200) %>%
        overlay(micro_subj[[3]], low=0.015, high=6,col=color_scale_div_1, rCol=color_scale_div_2, symmetric=TRUE) %>%
    sliceSeries(nrow=4, ncol=1, begin=42, end=52, dimension=3) %>%
        addtitle(names[4]) %>%
        anatomy(anatVol, low=10, high=200) %>%
        overlay(micro_subj[[4]], low=0.015, high=6,col=color_scale_div_1, rCol=color_scale_div_2, symmetric=TRUE) %>%
    sliceSeries(nrow=4, ncol=1, begin=42, end=52, dimension=3) %>%
        addtitle(names[5]) %>%
        anatomy(anatVol, low=10, high=200) %>%
        overlay(micro_subj[[5]], low=0.015, high=6,col=color_scale_div_1, rCol=color_scale_div_2, symmetric=TRUE) %>%
    sliceSeries(nrow=4, ncol=1, begin=42, end=52, dimension=3) %>%
        addtitle(names[6]) %>%
        anatomy(anatVol, low=10, high=200) %>%
        overlay(micro_subj[[6]], low=0.015, high=6,col=color_scale_div_1, rCol=color_scale_div_2, symmetric=TRUE) %>%
    sliceSeries(nrow=4, ncol=1, begin=42, end=52, dimension=3) %>%
        addtitle(names[7]) %>%
        anatomy(anatVol, low=10, high=200) %>%
        overlay(micro_subj[[7]], low=0.015, high=6,col=color_scale_div_1, rCol=color_scale_div_2, symmetric=TRUE) %>%
    legend("Z-values") %>%
    draw()
    dev.off()

    print(paste0("./visualization/1e_viz_zscores_subjects/",ids_to_viz[i],"_WMH_zscore.png"))

}

# Visualize subject raw micro

for(i in 1:length(ids_to_viz)) {
    anatVol = mincArray(mincGetVolume(paste0("../../../UKB/maps_UKB_space/sub-",ids_to_viz[i],"_ses-2_flair_UKB.mnc")))
    micro_subj = list()
    micro_subj[[1]] = mincArray(mincGetVolume(paste0("./results/1e_viz_zscores_subjects/",ids_to_viz[i],"_",names[1],"_WMH_raw.mnc")))
    micro_subj[[2]] = mincArray(mincGetVolume(paste0("./results/1e_viz_zscores_subjects/",ids_to_viz[i],"_",names[2],"_WMH_raw.mnc")))
    micro_subj[[3]] = mincArray(mincGetVolume(paste0("./results/1e_viz_zscores_subjects/",ids_to_viz[i],"_",names[3],"_WMH_raw.mnc")))
    micro_subj[[4]] = mincArray(mincGetVolume(paste0("./results/1e_viz_zscores_subjects/",ids_to_viz[i],"_",names[4],"_WMH_raw.mnc")))
    micro_subj[[5]] = mincArray(mincGetVolume(paste0("./results/1e_viz_zscores_subjects/",ids_to_viz[i],"_",names[5],"_WMH_raw.mnc")))
    micro_subj[[6]] = mincArray(mincGetVolume(paste0("./results/1e_viz_zscores_subjects/",ids_to_viz[i],"_",names[6],"_WMH_raw.mnc")))
    micro_subj[[7]] = mincArray(mincGetVolume(paste0("./results/1e_viz_zscores_subjects/",ids_to_viz[i],"_",names[7],"_WMH_raw.mnc")))

    down_lims = c(0, 0, 0, 0, 0, 0, 900)
    up_lims = c(0.003, 0.5, 0.6, 0.6, 0.4, 120, 1150)

    png(file=paste0("./visualization/1e_viz_zscores_subjects/",ids_to_viz[i],"_WMH_raw.png"), width=8500, height=4000, pointsize = 150)
    sliceSeries(nrow=4, ncol=1, begin=42, end=52, dimension=3) %>%
        addtitle("FLAIR") %>%
        anatomy(anatVol, low=10, high=200) %>%
    sliceSeries(nrow=4, ncol=1, begin=42, end=52, dimension=3) %>%
        addtitle(names[1]) %>%
        anatomy(anatVol, low=10, high=200) %>%
        overlay(micro_subj[[1]], low=down_lims[1], high=up_lims[1],col=color_scale_raw) %>%
    sliceSeries(nrow=4, ncol=1, begin=42, end=52, dimension=3) %>%
        addtitle(names[2]) %>%
        anatomy(anatVol, low=10, high=200) %>%
        overlay(micro_subj[[2]], low=down_lims[2], high=up_lims[2],col=color_scale_raw) %>%
    sliceSeries(nrow=4, ncol=1, begin=42, end=52, dimension=3) %>%
        addtitle(names[3]) %>%
        anatomy(anatVol, low=10, high=200) %>%
        overlay(micro_subj[[3]], low=down_lims[3], high=up_lims[3],col=color_scale_raw) %>%
    sliceSeries(nrow=4, ncol=1, begin=42, end=52, dimension=3) %>%
        addtitle(names[4]) %>%
        anatomy(anatVol, low=10, high=200) %>%
        overlay(micro_subj[[4]], low=down_lims[4], high=up_lims[4],col=color_scale_raw) %>%
    sliceSeries(nrow=4, ncol=1, begin=42, end=52, dimension=3) %>%
        addtitle(names[5]) %>%
        anatomy(anatVol, low=10, high=200) %>%
        overlay(micro_subj[[5]], low=down_lims[5], high=up_lims[5],col=color_scale_raw) %>%
    sliceSeries(nrow=4, ncol=1, begin=42, end=52, dimension=3) %>%
        addtitle(names[6]) %>%
        anatomy(anatVol, low=10, high=200) %>%
        overlay(micro_subj[[6]], low=down_lims[6], high=up_lims[6],col=color_scale_raw) %>%
    sliceSeries(nrow=4, ncol=1, begin=42, end=52, dimension=3) %>%
        addtitle(names[7]) %>%
        anatomy(anatVol, low=10, high=200) %>%
        overlay(micro_subj[[7]], low=down_lims[7], high=up_lims[7],col=color_scale_raw) %>%
    legend("Raw microstructure") %>%
    draw()
    dev.off()

    print(paste0("./visualization/1e_viz_zscores_subjects/",ids_to_viz[i],"_WMH_raw.png"))
}


