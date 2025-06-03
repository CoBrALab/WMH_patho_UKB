
library(data.table)
library(RMINC)
library(tidyverse)
library(ggplot2)
library(patchwork)

# Clean and combine WMH data
#   - log-transform WMH volumes
#   - add SuStaIn stages

dir.create("./results/3f_clean_wmh_data", showWarnings=FALSE)
dir.create("./visualization/3f_clean_wmh_data", showWarnings=FALSE)

# Load data
inclusions = as.data.frame(fread("../../QC/inclusions_final.txt"))
colnames(inclusions) = "ID"

demo = as.data.frame(fread("../../../UKB/tabular/df_demo/UKBB_demo_wider.tsv"))
demo = subset(demo, InstanceID == 2 & is.na(Age_when_attended_assessment_centre_21003_0) == FALSE, select=c("SubjectID", "Sex_31_0", "Age_when_attended_assessment_centre_21003_0"))
colnames(demo) = c("ID", "sex", "age")

demo = merge(inclusions, demo, by="ID")

# Append SuStaIn stages of different spatial clusters
df_1 = as.data.frame(fread(paste0("../sustain_spect_clust_k3_parc/results/c1/subtype0/comb_subj_cval_ses2.csv")))
df_2 = as.data.frame(fread(paste0("../sustain_spect_clust_k3_parc/results/c2/subtype0/comb_subj_cval_ses2.csv")))
df_3 = as.data.frame(fread(paste0("../sustain_spect_clust_k3_parc/results/c3/subtype0/comb_subj_cval_ses2.csv")))

sustain = as.data.frame(cbind(df_1$ID, df_1$ml_stage_cval, df_2$ml_stage_cval, df_3$ml_stage_cval))
colnames(sustain) = c("ID", "stage_c1", "stage_c2", "stage_c3")

# Median WMH pathophysiology in clusters

df = as.data.frame(fread("../2_spatial_clust/results/2h_roi_avg/subj_WMH_patho.tsv"))

micro = df[, c(1, seq(12, 32))]
region_vol = df[, c(1, seq(33,35))]

colnames(micro) = gsub("spect_clust_k3_", "", colnames(micro))

# WMH volumes (calculated in common space) in clusters
colnames(region_vol) = gsub("spect_clust_k3_", "", colnames(region_vol))

# log-transform
region_vol[,c(2,3,4)] = region_vol[,c(2,3,4)]+1 # Add 1 to not do log(0)
region_vol[,c(2,3,4)] = apply(region_vol[,c(2,3,4)], 2, log)

# Combine

df = merge(demo, region_vol, by="ID")
df = merge(df, sustain, by="ID")
df = merge(df, sustain_post_ant, by="ID")
df = merge(df, micro, by="ID")

fwrite(df, "./results/3f_clean_wmh_data/wmh_combined.tsv", row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
