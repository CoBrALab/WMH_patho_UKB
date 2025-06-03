# Load libraries

library(RMINC)
library(kernlab)
library(parallel)
library(Spectrum)
library(ggplot2)
library(scales)
library(data.table)

# Use spectral clustering on averaged WMH pathophysiology maps

dir.create("./results/2b_spectral_clust", showWarnings=FALSE)
dir.create("./visualization/2b_spectral_clust", showWarnings=FALSE)

# Import data

maps=c("FA", "MD", "ICVF", "ISOVF", "OD", "T2star", "QSM")

wmh_thresh = 30
nawm_thresh = 5000

df = fread(paste0("./results/2a_make_averages/WMH_zscore_",maps[1],"_avg.tsv"))

for (m in 2:length(maps)) {
  df = cbind(df, fread(paste0("./results/2a_make_averages/WMH_zscore_",maps[m],"_avg.tsv")))
}

colnames(df) = maps

# Mask out voxels lower than prevalence threshold
mask=mincGetVolume("../../data/WMH_mask.mnc")
prev = as.data.frame(fread("../tissue_prevalence/results/tissue_prevalence_after_exclusions_WMHmask_ses2.tsv"))
nawm_prev = as.numeric(prev[8,])
wmh_prev = as.numeric(prev[9,])

indices_above_thres = intersect(which(nawm_prev > nawm_thresh), which(wmh_prev > wmh_thresh))
logical_above_thres = rep(FALSE, length(nawm_prev))
logical_above_thres = 1:length(logical_above_thres) %in% indices_above_thres

df_mask = df[indices_above_thres,]

set.seed(123)

# Transpose input matrix (to get voxels * subject to cluster voxels)
df_clust = as.data.frame(t(df_mask))

# Run spectral clustering with the Spectrum package
range_of_ks = 10

spect = Spectrum(df_clust, method=1, runrange=TRUE, krangemax = range_of_ks, showpca=TRUE, fontsize=8, dotsize=2)

# Make elbow plot
elbow=as.data.frame(spect[[2]]$eigenvector_analysis)

write.csv(elbow, "./results/2b_spectral_clust/eigengap.csv", row.names = FALSE)

ggplot(elbow, aes(x = K, y = evals)) +
  geom_line(color="red") + 
  geom_point() +
  scale_x_continuous(breaks = 1:max(elbow$K), name="Eigenvector") + 
  scale_y_continuous(name="Eigenvalue") + 
  theme(text = element_text(size=30)) 
ggsave("./visualization/2b_spectral_clust/nclust_elbow.png", width = 12, height=6)

# Save clustering results
for (k in 1:(range_of_ks-1)){
  print(k+1)
  clust_final = cbind(df_mask, spect[[k+1]]$assignments)
  colnames(clust_final)[ncol(clust_final)] = "cluster"
  
  write.csv(clust_final, paste0("./results/2b_spectral_clust/clust_final_k",k+1,".csv"), row.names = FALSE)
}

