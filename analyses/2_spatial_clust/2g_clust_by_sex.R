
library(RMINC)
library(data.table)
library(ggplot2)
library(ClusterR)
library(grid)
library(gridExtra) 
library(tidyverse)
library(MRIcrotome)
library(magrittr)
library(viridis)
library(RColorBrewer)
library(kernlab)
library(parallel)
library(Spectrum)
library(scales)

# Calculate spatial clusters by sex

dir.create("./results/2g_clust_by_sex", showWarnings=FALSE)
dir.create("./visualization/2g_clust_by_sex", showWarnings=FALSE)

# Load data
# demo: age and sex
demo = as.data.frame(fread("../../../UKB/tabular/df_demo/UKBB_demo_wider.tsv"))
demo = subset(demo, InstanceID == 2 & is.na(Age_when_attended_assessment_centre_21003_0) == FALSE, select=c("SubjectID", "Sex_31_0", "Age_when_attended_assessment_centre_21003_0"))
colnames(demo) = c("ID", "sex", "age")

# inclusions: final sample of UKB inclusions (32,526 subjects)
inc = as.data.frame(fread("../../QC/inclusions_final.txt"))
colnames(inc) = "ID"
df = merge(inc, demo, by="ID", all.x=TRUE)
df$sex = as.factor(df$sex)

names=c( "MD", "ISOVF", "FA", "ICVF", "OD", "T2star", "QSM")

# mask: where WMH prevalence >1
mask = mincGetVolume("../../data/WMH_mask.mnc")
# label: BISON labels (subject by voxel in UKB space)
label = as.data.frame(fread(paste0("../../micro_matrices/ses2_Label_after_exclusions.tsv")))

# Load WMH and NAWM prevalence
prev = as.data.frame(fread("../tissue_prevalence/results/tissue_prevalence_after_exclusions_WMHmask_ses2.tsv"))
nawm_prev = as.numeric(prev[8,])
wmh_prev = as.numeric(prev[9,])

wmh_thresholds = 15
nawm_thresh = 5000

# Make WMH prevalence files by WMH volume severity (with MD zscore file only sinve biggest amplitude)
zvals = as.data.frame(fread(paste0("../norm_zscores_subjects/results/WMH_zscores_subj_",names[2],"_NAs.tsv")))
zvals[!(label == 9)] = NA

for (i in c("Female", "Male")) {
    print(paste0("sex = ",i))
    # Make z-score matrix
    zvals_sex = zvals[which(df$sex == i),]

    # Get WMH prevalence for WMH severity subset
    wmh_prev_sex = colSums(zvals_sex != 0, na.rm=TRUE)
    write.table(wmh_prev_sex, paste0("./results/2g_clust_by_sex/WMH_prev_",i,".tsv"), row.names=FALSE, col.names=FALSE, quote=FALSE)
}

# Calculate averages
for (m in 1:length(names)) {
    print(names[m])

    # Load subject Z-scores
    zvals = as.data.frame(fread(paste0("../norm_zscores_subjects/results/WMH_zscores_subj_",names[m],"_NAs.tsv")))
    zvals[!(label == 9)] = NA

    # For every sex
    for (i in c("Female", "Male")) {
        print(paste0("Sex = ",i))

        # Make z-score matrix
        zvals_sex = zvals[which(df$sex == i),]

        # Get WMH prevalence for sex
        wmh_prev_sex = fread(paste0("./results/2g_clust_by_sex/WMH_prev_",i,".tsv"))
        wmh_prev_sex = as.numeric(wmh_prev_sex$V1)

        # Assign 0 to voxels with prevalence NAWM < x and WMH < y
        wmh_thresh = wmh_thresholds

        indices_above_thres = intersect(which(nawm_prev > nawm_thresh), which(wmh_prev_sex >= wmh_thresh))
        logical_above_thres = rep(FALSE, length(nawm_prev))
        logical_above_thres = 1:length(logical_above_thres) %in% indices_above_thres

        # Voxel averages of zscores
        avg = colMeans(zvals_sex, na.rm=TRUE)
        avg[is.na(avg)] = 0
        avg[logical_above_thres==FALSE] = 0

        # Write to tsv
        write.table(avg, paste0("./results/2g_clust_by_sex/sex_",i,"_WMH_zscore_", names[m],"_avg.tsv"), row.names=FALSE, col.names=FALSE, sep='\t', quote=FALSE)

        # Write to mnc
        outvol = mincGetVolume("../../../UKB/temporary_template/avg.020_2mm.mnc")
        outvol[] = 0
        outvol[mask > 0.5] = avg
        mincWriteVolume(outvol, paste0("./results/2g_clust_by_sex/sex_",i,"_WMH_zscore_", names[m],"_avg.mnc"), clobber=TRUE)
    
        rm (zvals_sex)
    }
    rm(zvals)
}

# Visualize averages

color_scale_div_2 = colorRampPalette(brewer.pal(9,"Blues"))(255)
color_scale_div_1 = colorRampPalette(brewer.pal(9,"Reds"))(255)

anatVol = mincArray(mincGetVolume("../../data/UKB_template_2mm.mnc"))

names=c("MD", "ISOVF", "FA", "ICVF", "OD", "T2star", "QSM")

for (i in c("Female", "Male")) {
    png(file=paste0("./visualization/2g_clust_by_sex/sex_",i,"_WMH_15_NAWM_5000.png"), width=7500, height=4000, pointsize = 150)
    sliceSeries(nrow=4, ncol=1, begin=42, end=52, dimension=3) %>%
        addtitle(names[1]) %>%
        anatomy(anatVol, low=10, high=200) %>%
        overlay(mincArray(mincGetVolume(paste0("./results/2g_clust_by_sex/sex_",i,"_WMH_zscore_", names[1],"_avg.mnc"))),
            low=0.015, high=6,col=color_scale_div_1, rCol=color_scale_div_2, symmetric=TRUE) %>%
    sliceSeries() %>%
        addtitle(names[2]) %>%
        anatomy(anatVol, low=10, high=200) %>%
        overlay(mincArray(mincGetVolume(paste0("./results/2g_clust_by_sex/sex_",i,"_WMH_zscore_", names[2],"_avg.mnc"))),
            low=0.015, high=6,col=color_scale_div_1, rCol=color_scale_div_2, symmetric=TRUE) %>%
    sliceSeries() %>%
        addtitle(names[3]) %>%
        anatomy(anatVol, low=10, high=200) %>%
        overlay(mincArray(mincGetVolume(paste0("./results/2g_clust_by_sex/sex_",i,"_WMH_zscore_", names[3],"_avg.mnc"))),
            low=0.015, high=6,col=color_scale_div_1, rCol=color_scale_div_2, symmetric=TRUE) %>%
    sliceSeries() %>%
        addtitle(names[4]) %>%
        anatomy(anatVol, low=10, high=200) %>%
        overlay(mincArray(mincGetVolume(paste0("./results/2g_clust_by_sex/sex_",i,"_WMH_zscore_", names[4],"_avg.mnc"))),
            low=0.015, high=6,col=color_scale_div_1, rCol=color_scale_div_2, symmetric=TRUE) %>%
    sliceSeries() %>%
        addtitle(names[5]) %>%
        anatomy(anatVol, low=10, high=200) %>%
        overlay(mincArray(mincGetVolume(paste0("./results/2g_clust_by_sex/sex_",i,"_WMH_zscore_", names[5],"_avg.mnc"))),
            low=0.015, high=6,col=color_scale_div_1, rCol=color_scale_div_2, symmetric=TRUE) %>%
    sliceSeries() %>%
        addtitle(names[6]) %>%
        anatomy(anatVol, low=10, high=200) %>%
        overlay(mincArray(mincGetVolume(paste0("./results/2g_clust_by_sex/sex_",i,"_WMH_zscore_", names[6],"_avg.mnc"))),
            low=0.015, high=6,col=color_scale_div_1, rCol=color_scale_div_2, symmetric=TRUE) %>%
    sliceSeries() %>%
        addtitle(names[7]) %>%
        anatomy(anatVol, low=10, high=200) %>%
        overlay(mincArray(mincGetVolume(paste0("./results/2g_clust_by_sex/sex_",i,"_WMH_zscore_", names[7],"_avg.mnc"))),
            low=0.015, high=6,col=color_scale_div_1, rCol=color_scale_div_2, symmetric=TRUE) %>%
    legend("Z-values") %>%
    draw()
    dev.off()

    print(paste0("./visualization/sex_",i,"_WMH_15_NAWM_5000.png"))
}

# Cluster average maps WMH z-scores by sex

# For each sex
for (sex_clust in c("Female", "Male")) {

    # Import data
    maps=c("FA", "MD", "ICVF", "ISOVF", "OD", "T2star", "QSM")

    wmh_thresh = 15
    nawm_thresh = 5000

    df = fread(paste0("../results/2g_clust_by_sex/sex_",sex_clust,"_WMH_zscore_",maps[1],"_avg.tsv"))
    for (m in 2:length(maps)) {
        df = cbind(df, fread(paste0("../results/2g_clust_by_sex/sex_",sex_clust,"_WMH_zscore_",maps[m],"_avg.tsv")))
    }
    colnames(df) = maps

    # Mask out voxels lower than prevalence threshold
    mask=mincGetVolume("../../data/WMH_mask.mnc")

    indices_above_thres = which(rowSums(df) != 0)
    df_mask = df[indices_above_thres,]

    set.seed(123)
    
    # Transpose input matrix (to get voxels * subject to cluster voxels)
    df_clust = as.data.frame(t(df_mask))

    # Run spectral clustering with the Spectrum package
    # Set k=3 to match main clustering solution on whole sample
    k=3
    spect = Spectrum(df_clust, method=3, runrange=FALSE, fixk = k, showpca=TRUE, fontsize=8, dotsize=2)

    # Save clustering results
    clust_final = cbind(df_mask, spect[[k]]$assignments)
    colnames(clust_final)[ncol(clust_final)] = "cluster"
    write.csv(clust_final, paste0("./results/2g_clust_by_sex/",sex_clust,"_clust_final_k",k,".csv"), row.names = FALSE)

}


# Visualize clusters

# Different choices of k
for (sex_clust in c("Female", "Male")) {

    # Mask out voxels lower than prevalence threshold
    mask=mincGetVolume("../../data/WMH_mask.mnc")
    df = fread(paste0("./results/2g_clust_by_sex/sev_",i,"_WMH_zscore_", names[m],"_avg.tsv"))
    indices_above_thres = which(rowSums(df) != 0)

    logical_above_thres = rep(FALSE, dim(df)[1])
    logical_above_thres = 1:length(logical_above_thres) %in% indices_above_thres
    
    cluster_number = 3
    k=3
    clust_number = cluster_number
    clust_final = read.csv(paste0("./results/2g_clust_by_sex/",sex_clust,"_clust_final_k",k,".csv"))

    # Re-order clusters manually
    if (s == "Female") {clust_final$cluster = ifelse(clust_final$cluster==1, 1, ifelse(clust_final$cluster==2, 3, ifelse(clust_final$cluster==3, 2,0)))}
    if (s == "Male") {clust_final$cluster = ifelse(clust_final$cluster==1, 3, ifelse(clust_final$cluster==2, 1, ifelse(clust_final$cluster==3, 2,0)))}

    # Add back masked out voxels
    clust_no_mask = matrix(data=0,nrow=length(logical_above_thres), ncol=ncol(clust_final))

    valid_idx = 1
    for (i in 1:length(logical_above_thres)) {
    if(logical_above_thres[i]==TRUE) {
        clust_no_mask[i,] = as.numeric(clust_final[valid_idx,])
        valid_idx = valid_idx+1
        }
    }

    clust_no_mask = as.data.frame(clust_no_mask)
    colnames(clust_no_mask) = colnames(clust_final)

    # Combine cluster mnc files
    outvol <- mincGetVolume("../../data/UKB_template_2mm.mnc")
    outvol[] <- 0
    outvol[mask > 0.5] <- clust_no_mask$cluster
    mincWriteVolume(outvol, paste0("./results/2g_clust_by_sex/",sex_clust,"_all_clusters.mnc"), clobber=TRUE)

    # Spatial clustering images  
    anatVol = mincArray(mincGetVolume("../../data/UKB_template_2mm.mnc"))

    color_scale = viridis(clust_number)

    png(file=paste0("./visualization/2g_clust_by_sex/",sex_clust,"_clusters.png"), width=3000, height=3000, pointsize = 80)
    sliceSeries(nrow=3, ncol=3, begin=40, end=52, dimension=3) %>%
        anatomy(anatVol, low=10, high=140) %>%
        overlay(mincArray(mincGetVolume(paste0("./results/2g_clust_by_sex/",sex_clust,"_all_clusters.mnc"))),
                low=0, high=clust_number+1, col=color_scale) %>%
        legend("Cluster") %>%
        draw(layout="row")
    dev.off()

    # Microstructural distributions by WMH cluster

    # By cluster
    micro_means_cluster = as.data.frame(matrix(nrow=0, ncol=3))
    colnames(micro_means_cluster) = c("micro", "means", "cluster")

    for (c in 1:clust_number){
        print(c)

        # Select micro data in cluster
        micro_maps = subset(clust_final, cluster==c, select=c(maps))
        micro_2 = as.data.frame(matrix(nrow=0, ncol=2))
        colnames(micro_2) = c("val", "micro")

        for (m in 1:length(maps)){
            micro_1 = as.data.frame(micro_maps[,m])
            colnames(micro_1) = "val"
            micro_1$micro = maps[m]
            micro_2 = rbind(micro_2, micro_1)
        }

        micro_2$micro = factor(micro_2$micro, levels = maps)

        # Calculate median values of distributions
        micro_means = micro_2 %>%
            group_by(micro) %>%
            summarize(mean=median(val))
        micro_means$cluster=c

        micro_means_cluster = rbind(micro_means_cluster, micro_means)

        # Violin plots
        color_scale = c(MD = "#04319E", ISOVF = "#8298CF", FA = "#2B520B", ICVF="#448312", OD="#A2C189", T2star="#D36108", QSM="#E9B084")

        ggplot(micro_2, aes(x=micro, y=val, fill=micro)) + 
            geom_hline(yintercept=0) +
            geom_violin(trim=TRUE, scale="width") + 
            geom_boxplot(width=0.2) + 
            geom_label(data=micro_means, aes(y=11, label=round(mean, 2), color=micro), size=10, fill="white", label.size=0) + 
            theme_classic() + 
            labs(x="", y="Z-scores") + 
            ggtitle(paste0("Cluster ",c)) + 
            scale_y_continuous(breaks=seq(-6,10,by=2), limits=c(-6,11)) + 
            scale_fill_manual(name="", values=color_scale) +
            scale_color_manual(name="", values=color_scale) +
            theme(text=element_text(size=30, color="black"), plot.title=element_text(hjust=0.5, size=40),
            legend.position = "none")
        ggsave(paste0("./visualization/2g_clust_by_sex/",sex_clust,"_allMicro_cluster",c,".png"),
                width = length(maps)+3, height = 10, dpi=300, units="in")

    }  
}

# Calculate dice overlaps with original clusters on whole-sample

# Load data
clust_parc = mincGetVolume("./results/2c_viz_clust/k3_all_clusters.mnc")

other_parc = list()
other_parc[['Male']] = mincGetVolume("./results/2g_clust_by_sex/Male_all_clusters.mnc")
other_parc[['Female']] = mincGetVolume("./results/2g_clust_by_sex/Female_all_clusters.mnc")

clust_parc_labels = c("Periventricular" = 1, "Posterior" = 2, "Anterior" = 3)

# Calculate dice overlap between data-driven clusters and other parcellations

dice_coef <- function(A, B) {
  intersection <- sum(A & B)
  union <- sum(A) + sum(B)
  return(2 * intersection / union)
}

results = data.frame(
  Other_parc = character(),
  spect_clust_num = numeric(),  
  other_parc_num = numeric(),
  spect_clust_label = character(),  
  other_parc_label = character(),
  dice = numeric(),
  stringsAsFactors = FALSE
)

# For every sex
for (o in 1:length(other_parc)) {
    print(names(other_parc)[o])

    # Select voxels that are included in both parcellations
    clust_masked = clust_parc[which(clust_parc[]>0 & other_parc[[o]][] > 0)]
    parc_masked = other_parc[[o]][which(clust_parc[]>0 & other_parc[[o]][] > 0)]

    # For every parcel in original parcellation
    for (c in unique(clust_masked)) {
        clust_label = names(clust_parc_labels[which(clust_parc_labels == c)])
        print(clust_label)

        # For every parcel in parcellation by sex
        for (p in unique(parc_masked)) {

            # Isolate voxels
            clust_masked_indiv = ifelse(clust_masked == c, 1, 0)
            parc_masked_indiv = ifelse(parc_masked == p, 1, 0)

            parc_label = names(clust_parc_labels[which(clust_parc_labels == p)])

            print(parc_label)

            # Calculate dice coefficient
            dice_pair = dice_coef(clust_masked_indiv, parc_masked_indiv)

            # save results
            new_row = data.frame(Other_parc = names(other_parc)[o],
                                spect_clust_num = c,
                                other_parc_num = p,
                                spect_clust_label = clust_label,
                                other_parc_label = parc_label,
                                dice = dice_pair,
                                stringsAsFactors = FALSE)

            # Add the new row to the existing data frame
            results = rbind(results, new_row)  
        }
    }
}

fwrite(results, "./results/2g_clust_by_sex/overlap_clust_whole_cohort.tsv", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

# Visualize in table

results = as.data.frame(fread("./results/2g_clust_by_sex/overlap_clust_whole_cohort.tsv"))

results %>% mutate(Other_parc = factor(Other_parc),
                    other_parc_label = factor(other_parc_label, levels=c("Periventricular","Posterior","Anterior")),
                    spect_clust_label = factor(spect_clust_label, levels=c("Anterior", "Posterior", "Periventricular"))) %>%
ggplot(aes(x=other_parc_label, y=spect_clust_label, fill=dice)) +
    geom_tile(color = "white", lwd = 0.1, linetype = 1) +
    geom_text(aes(label = round(dice, 2)), vjust = 1, size = 5, color = "white") +
    scale_fill_gradientn(colors = viridis_pal(direction=1)(9), limits=c(0, 1), breaks=c(0, 1), name="Dice") +
    facet_grid(. ~ Other_parc, scales="free", space="free") +
    theme_light() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 15),
        legend.text = element_text(size=15), legend.title = element_text(size=15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        strip.text = element_text(size=15, color="black"), strip.background.x = element_rect(fill="#ffffff"))
ggsave("./visualization/2g_clust_by_sex/overlap_clust_whole_cohort.png", height=5, width=13)

