
# Load libraries

library(RMINC)
library(kernlab)
library(parallel)
library(Spectrum)
library(ggplot2)
library(scales)
library(data.table)
library(dplyr)
library(grid)
library(gridExtra)
library(tidyverse)
library(MRIcrotome)
library(magrittr)
library(viridis)

# Visualization tools for spatial WMH clusters

dir.create("./results/2c_viz_clust", showWarnings=FALSE)
dir.create("./visualization/2c_viz_clust", showWarnings=FALSE)

maps = c("MD", "ISOVF", "FA", "ICVF", "OD", "T2star", "QSM")

# Mask out voxels lower than prevalence threshold
mask=mincGetVolume("../../data/WMH_mask.mnc")
prev = as.data.frame(fread("../tissue_prevalence/results/tissue_prevalence_after_exclusions_WMHmask_ses2.tsv"))
nawm_prev = as.numeric(prev[8,])
wmh_prev = as.numeric(prev[9,])

wmh_thresh = 30
nawm_thresh = 5000

indices_above_thres = intersect(which(nawm_prev > nawm_thresh), which(wmh_prev > wmh_thresh))
logical_above_thres = rep(FALSE, length(nawm_prev))
logical_above_thres = 1:length(logical_above_thres) %in% indices_above_thres

# Different choices of k
clust_list = c(2,3,4,5)

anatVol = mincArray(mincGetVolume("../../data/UKB_template_2mm.mnc"))

for (cluster_number in clust_list) {
  
  dir.create(paste0("./visualization/2c_viz_clust/k",cluster_number), showWarnings=FALSE)

  clust_number = cluster_number
  clust_final = read.csv(paste0("./results/2b_spectral_clust/clust_final_k",clust_number,".csv"))
  clust_final$cluster = as.factor(clust_final$cluster)

  # Re-order clusters by ascending order of MD mean
  md_means = clust_final %>%
      group_by(cluster) %>%
      summarize(md_avg = mean(MD))

  clust_order = order(md_means$md_avg)
  level_mapping <- setNames(levels(clust_final$cluster)[clust_order], levels(clust_final$cluster))

  clust_final$cluster = as.numeric(level_mapping[clust_final$cluster])

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
  mincWriteVolume(outvol, paste0("./results/2b_spectral_clust/k",clust_number,"_all_clusters.mnc"), clobber=TRUE)

  # Spatial clustering images

  color_scale = viridis(clust_number)
  
  # Clustered voxel locations
  png(file=paste0("./visualization/2c_viz_clust/k",clust_number,"_all_clusters.png"), width=3000, height=3000, pointsize = 80)
  sliceSeries(nrow=3, ncol=3, begin=40, end=52, dimension=3) %>%
    anatomy(anatVol, low=10, high=140) %>%
    overlay(mincArray(mincGetVolume(paste0("./results/k",clust_number,"_all_clusters.mnc"))),
            low=0, high=clust_number+1, col=color_scale) %>%
    draw(layout="row")
  dev.off()

  # Microstructural distributions by WMH cluster
  
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
    ggsave(paste0("./visualization/2c_viz_clust/k",clust_number,"_allMicro_cluster",c,".png"),
           width = length(maps)+3, height = 10, dpi=300, units="in")
  }  
}


