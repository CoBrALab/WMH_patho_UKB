
library(data.table)
library(RMINC)
library(grid)
library(gridExtra) 
library(tidyverse)
library(MRIcrotome)
library(magrittr)
library(viridis)
library(RColorBrewer)
library(ggplot2)
library(patchwork)

# Visualize normative models of a few voxels

dir.create("./results/1b_viz_norm_models", showWarnings=FALSE)
dir.create("./visualization/1b_viz_norm_models", showWarnings=FALSE)

names = c("MD", "ISOVF", "FA", "ICVF", "OD", "T2star", "QSM")
names_long = c("dti_MD", "NODDI_ISOVF", "dti_FA", "NODDI_ICVF", "NODDI_OD", "T2star", "QSM")

# mask: where WMH prevalence >1
mask = mincGetVolume("../../data/WMH_mask.mnc")

# inclusions: final sample of UKB inclusions (32,526 subjects)
inclusions = as.data.frame(fread("../../QC/inclusions_final.txt"))
colnames(inclusions) = "ID"

# demo: age and sex
demo = as.data.frame(fread("../../../UKB/tabular/df_demo/UKBB_demo_wider.tsv"))
demo = subset(demo, InstanceID == 2, select=c('SubjectID', 'Sex_31_0', 'Age_when_attended_assessment_centre_21003_0'))
colnames(demo) = c("ID", "Sex", "Age")
demo = merge(inclusions, demo, by="ID", all=FALSE)

# anatVol: UKB template
anatVol = mincArray(mincGetVolume("../../data/UKB_template_2mm.mnc"))

# label: BISON labels (subject by voxel in UKB space)
label = as.data.frame(fread(paste0("../../micro_matrices/ses2_Label_after_exclusions.tsv")))
label = cbind(inclusions$ID, label)
colnames(label)[1] = "ID"

# Selected voxels that had high WMH and NAWm prevalence
vox_to_viz = c(32376, 32734)

micro = list()
nm = list()

# Load micro and NM
for (n in 1:length(names)) {
    print(names[n])
    # Load raw microstructure matrix
    micro[[n]] = as.data.frame(fread(paste0("../../micro_matrices/ses2_",names_long[n],"_after_exclusions.tsv")))
    micro[[n]] = micro[[n]][,vox_to_viz]
    micro[[n]] = cbind(inclusions$ID, micro[[n]])
    colnames(micro[[n]])[1] = "ID"

    # Load normative model outputs
    nm[[n]] = as.data.frame(fread(paste0("./results/1a_norm_models_BLR/nm_",names[n],"_ukb.tsv")))
    nm[[n]] = nm[[n]][vox_to_viz,]
}

# Plotting functions
plot_vox = function(nm_df, df_wmh, df_nawm, name, ymax, ymin) {
    
    male_avg = as.numeric(nm_df[, grep("Male_mean", names(nm_df), value = TRUE)])
    male_sd = as.numeric(nm_df[, grep("Male_sd", names(nm_df), value = TRUE)])
    female_avg = as.numeric(nm_df[, grep("Female_mean", names(nm_df), value = TRUE)])
    female_sd = as.numeric(nm_df[, grep("Female_sd", names(nm_df), value = TRUE)])

    nm_df_2 = data.frame("Age" = seq(45,81), "Male_avg" = male_avg, "Male_SD" = male_sd, "Female_avg" = female_avg, "Female_SD" = female_sd)
    nm_df_2 = nm_df_2 %>%
        pivot_longer(cols = -Age, 
        names_to = c("Sex_nm", ".value"), 
        names_pattern = "(Male|Female)_(avg|SD)")

    nm_df = as.data.frame(nm_df_2)

    plt = ggplot() +
            stat_density_2d(data=df_nawm, aes(x=Age, y=Value, fill = after_stat(density)), geom = "raster", alpha = 0.5, contour = FALSE) + 
            geom_line(data = nm_df, aes(x=Age, y=avg, linetype=Sex_nm), color="black", size=1)
            geom_point(data=df_wmh, aes(x=Age, y=Value, color=Sex), alpha=0.5, size=2) +
            scale_fill_viridis_c(name="              ") +   
            # geom_point(data=df_wmh, aes(x=Age, y=Value, shape=Sex), alpha=0.5, size=2, color="blue") + 
            # scale_shape_discrete(name="Sex (WMH data)") + 
            # geom_line(data = nm_df, aes(x=Age, y=avg, linetype=Sex_nm), color="red", size=1) + 
            # scale_linetype_discrete(name="Sex (NAWM pred)") + 
            ggtitle(names[n]) +
            scale_y_continuous(name="Value", limits=c(ymin, ymax)) +
            scale_x_continuous(name="Age") +
            scale_color_manual(values=c(Female="#f90609", Male="#1006f9")) + 
            theme_classic() + 
            theme(text=element_text(size=30), plot.title=element_text(hjust=0.5, size=30))

    ggsave(paste0(name, "_withWMH.png"), width=10, height=6)
    print(paste0(name, "_withWMH.png"))

    return(plt)
}

# Constant y axis thresholds
raw_up_thresh = c(0.003, 0.85, 0.85, 1, 1, 0.6, 100)
raw_down_thresh = c(0, 0, 0, 0, 0, 0, -100)

# Visualize voxels
for (i in 1:length(vox_to_viz)) {
    print(vox_to_viz[i])
    dir.create(paste0("./visualization/1b_viz_norm_models/vox_",vox_to_viz[i]), showWarnings=FALSE)
    plots = list()

    # Voxel values and nm fit
    for (n in 1:length(names)) {
        print(names[n])

        df = as.data.frame(cbind(demo, micro[[n]][,i+1], label[,vox_to_viz[i]+1]))
        colnames(df)[c(4,5)] = c("Value", "Label")

        # Remove voxels other than NAWM and WMH
        df_wmh = subset(df, Label == 9)
        df_wmh$Label = as.factor(df_wmh$Label)
        df_nawm = subset(df, Label == 8)
        df_nawm$Label = as.factor(df_nawm$Label)

        # Import age and sex trends from NM
        nm_df = nm[[n]][i,]

        # Plots
        plots[[n]]=plot_vox(nm_df, df_wmh, df_nawm, paste0("./visualization/1b_viz_norm_models/vox_",vox_to_viz[i],"/",names[n],"_",vox_to_viz[i]), raw_up_thresh[n], raw_down_thresh[n])
    }
    wrap_plots(plots, ncol=3)
    ggsave(paste0("./visualization/1b_viz_norm_models/vox_",vox_to_viz[i],"/all_",vox_to_viz[i],".png"), width=24, height=15)
}

# Visualize where voxels are located
for (i in 1:length(vox_to_viz)) {
    # Save mask of voxel
    outvol <- mincGetVolume("../../data/UKB_template_2mm.mnc")
    outvol[] <- 0
    vox = c(rep(0, vox_to_viz[i] - 1), 1, rep(0, (length(wmh_prev) - (vox_to_viz[i]))))
    outvol[mask > 0.5] <- vox
    mincWriteVolume(outvol, paste0("./visualization/1b_viz_norm_models/vox_",vox_to_viz[i],"/vox_",vox_to_viz[i],".mnc"), clobber=TRUE)

    # Image visualizing where voxel is located
    vox_mask = mincArray(mincGetVolume(paste0("./visualization/1b_viz_norm_models/vox_",vox_to_viz[i],"/vox_",vox_to_viz[i],".mnc")))
    slices = which(vox_mask[] == 1, arr.ind = TRUE)
    png(file=paste0("./visualization/1b_viz_norm_models/vox_",vox_to_viz[i],"/vox_",vox_to_viz[i],".png"), width=6000, height=2500, pointsize = 80)
    sliceSeries(nrow=1, ncol=1, begin=slices[1], end=slices[1], dimension=1) %>%
        anatomy(anatVol, low=10, high=140) %>%
        overlay(vox_mask, low=0, high=1, col="red") %>%
    sliceSeries(nrow=1, ncol=1, begin=slices[2], end=slices[2], dimension=2) %>%
        anatomy(anatVol, low=10, high=140) %>%
        overlay(vox_mask, low=0, high=1, col="red") %>%
    sliceSeries(nrow=1, ncol=1, begin=slices[3], end=slices[3], dimension=3) %>%
        anatomy(anatVol, low=10, high=140) %>%
        overlay(vox_mask, low=0, high=1, col="red") %>%
    draw()
    dev.off()
    print(paste0("./visualization/1b_viz_norm_models/vox_",vox_to_viz[i],"/vox_",vox_to_viz[i],".png"))

    command = paste0("convert -gravity center ./visualization/1b_viz_norm_models/vox_",vox_to_viz[i],"/vox_",vox_to_viz[i],".png ./visualization/1b_viz_norm_models/vox_",vox_to_viz[i],"/all_",vox_to_viz[i],".png +append ./visualization/1b_viz_norm_models/vox_",vox_to_viz[i],"/vox_all_",vox_to_viz[i],".png")
    system(command)
}

