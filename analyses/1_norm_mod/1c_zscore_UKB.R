
library(RMINC)
library(data.table)
library(matrixStats)
library(dplyr)
# library(tidyverse)

# Calculate subject-wise z-score maps
# Cannot share the data due to confidentiality of individual participants

print("Load data")

dir.create("./results/1c_zscore_UKB", showWarnings=FALSE)
dir.create("./visualization/1c_zscore_UKB", showWarnings=FALSE)

# Arg1: Name of micro
# Arg2: Mask
args = commandArgs(trailingOnly=TRUE)

# args=c(1, "../../data/WMH_mask.mnc")

names = c("FA", "MD", "ICVF", "ISOVF", "OD", "T2star", "QSM")
names_long = c("dti_FA", "dti_MD", "NODDI_ICVF", "NODDI_ISOVF", "NODDI_OD", "T2star", "QSM")

tissue=c(8,9) # Only WMH and NAWM
tissue_names=c('Ventricules', 'CSF', 'Cerebellum_GM', 'Cerebellum_WM', 'Brainstem', 'Subcortical_GM', 'Cortical_GM', 'Cerebral_NAWM', 'WMH')

# Run this script for each microstructural marker (number from 0 to 6)
print(names[as.numeric(args[1])])

# inclusions: final sample of UKB inclusions (32,526 subjects)
inclusions = as.data.frame(fread("../../QC/inclusions_final.txt"))
colnames(inclusions) = c("ID")

# demo: age and sex
demo = as.data.frame(fread("../../../UKB/tabular/df_demo/UKBB_demo_wider.tsv"))
demo = subset(demo, InstanceID == 2, select=c('SubjectID', 'Sex_31_0', 'Age_when_attended_assessment_centre_21003_0'))
colnames(demo) = c("ID", "Sex", "Age")
demo = merge(inclusions, demo, by="ID", all.x=TRUE)

# Load MRI data
# micro: microstructural values (subject by voxel in UKB space)
micro = as.data.frame(fread(paste0("../../micro_matrices/ses2_",names_long[as.numeric(args[1])],"_after_exclusions.tsv")))
# label: BISON labels (subject by voxel in UKB space)
label = as.data.frame(fread(paste0("../../micro_matrices/ses2_Label_after_exclusions.tsv")))
# mask: voxels in UKB space with at least 1 WMH label
mask=mincGetVolume(args[2])

# Remove microstructural values of tissue types other than NAWM and WMH
micro[!(label == 8 | label == 9)] = NA

micro = cbind(demo, micro)

# Load age- and sex-specific NAWM means and SDs
nm = as.data.frame(fread(paste0("./results/1a_norm_models_BLR/nm_",names[as.numeric(args[1])],"_ukb.tsv")))

print("Zscore")

# Function to z-score with age- and sex-specific NAWM averages
nm_zscore = function(x) {
    mean_col = as.numeric(which(colnames(nm)==paste0(x[2], "_mean_", x[3])))
    sd_col = as.numeric(which(colnames(nm)==paste0(x[2], "_sd_", x[3])))
    mean_id = nm[,mean_col]
    sd_id = nm[,sd_col]

    wmh_z = ifelse(is.na(x[seq(4,length(x))]) | is.na(mean_id) | is.na(sd_id), NA, (as.numeric(x[seq(4,length(x))]) - as.numeric(mean_id)) / as.numeric(sd_id))

    return(wmh_z)
}

# Do separately on two halves, otherwise running into memory issues

half_n <- floor(nrow(micro) / 2)

micro_1 <- micro %>% slice(1:half_n)
micro_2 <- micro %>% slice((half_n + 1):nrow(micro))
rm(micro)

micro_1_z = t(apply(micro_1, 1, nm_zscore))
micro_2_z = t(apply(micro_2, 1, nm_zscore))

micro_1_z[is.infinite(micro_1_z)] = NA
micro_2_z[is.infinite(micro_2_z)] = NA

micro_z = rbind(micro_1_z, micro_2_z)
rm(micro_1_z)
rm(micro_2_z)

print("Write to tsv")

# Write with NAs
fwrite(as.data.frame(micro_z), paste0("./results/1c_zscore_UKB/WMH_zscores_UKB_",names[as.numeric(args[1])],"_NAs.tsv"), sep="\t")

# Write with 0s
micro_z[is.na(micro_z)] = 0

fwrite(as.data.frame(micro_z), paste0("./results/1c_zscore_UKB/WMH_zscores_UKB_",names[as.numeric(args[1])],".tsv"))

