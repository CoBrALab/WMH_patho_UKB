
# Load libraries

library(RMINC)
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

# Quantify dice overlap between our data-driven clusters with other holistic parcellations

dir.create("./results/2e_overlap_other_parc", showWarnings=FALSE)
dir.create("./visualization/2e_overlap_other_parc", showWarnings=FALSE)

# Load data
clust_parc = mincGetVolume("./results_2_spectral_clust/k3/all_clusters.mnc")

other_parc = list()
# PV/deep: Based on threshold of 8mm from the ventricle
other_parc[['pv_deep']] = mincGetVolume("../../data/WMH_parc_pv_deep.mnc")
# By lobe
other_parc[['lobar']] = mincGetVolume("../../data/WMH_parc_lobar.mnc")
# Cerebral arterial territories: From https://doi.org/10.1038/s41597-022-01923-0
other_parc[['vascular']] = mincGetVolume("../../data/WMH_parc_vascular.mnc")

other_parc[['pv_deep']][] = round(other_parc[['pv_deep']][], 1)
other_parc[['lobar']][] = round(other_parc[['lobar']][], 1)
other_parc[['vascular']][] = round(other_parc[['vascular']][], 1)

clust_parc_labels = c("Periventricular" = 1, "Posterior" = 2, "Anterior" = 3)

other_parc_labels = list()
other_parc_labels[['pv_deep']] = c("Periventricular" = 2, "Deep" = 1, "OOB" = 0)
other_parc_labels[['lobar']] = c("Frontal" = 1, "Temporal" = 3, "Parietal" = 5, "Occipital" = 9, "Cerebellum" = 7, "Brainstem" = 11, "OOB" = 0)
other_parc_labels[['vascular']] = c("ACA" = 1, "MCA" = 3, "PCA" = 5, "VB" = 7, "Ventricules" = 9, "OOB" = 0)

# Calculate dice overlap between data-driven clusters and other parcellations

# Dice coefficient function
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

clust_masked = clust_parc[which(clust_parc[]>0)]

# For each other parcellation
for (o in 1:length(other_parc)) {
    print(names(other_parc)[o])
    parc_masked = other_parc[[o]][which(clust_parc[]>0)]

    # For each cluster in data-driven clusters
    for (c in unique(clust_masked)) {
        clust_label = names(clust_parc_labels[which(clust_parc_labels == c)])
        print(clust_label)

        # For each parcel in other parcellation
        for (p in unique(parc_masked)) {

            # Isolate voxels
            clust_masked_indiv = ifelse(clust_masked == c, 1, 0)
            parc_masked_indiv = ifelse(parc_masked == p, 1, 0)

            parc_label = names(other_parc_labels[[o]][which(other_parc_labels[[o]] == p)])

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

# Save results
fwrite(results, "./results/2e_overlap_other_parc/k3_overlap_other_parc.tsv", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

# Visualize in table

results = as.data.frame(fread("./results/2e_overlap_other_parc/k3_overlap_other_parc.tsv"))

# Remove overlap with non-brain
results %>% 
    filter(!other_parc_label %in% c("OOB", "Ventricules")) %>%
    filter(Other_parc!="pv_deep_10mm") %>%
    mutate(Other_parc = factor(Other_parc, levels=c("pv_deep_8mm", "lobar", "vascular"), labels=c("PV/deep", "Lobar", "Vascular"))) %>%
    mutate(spect_clust_label = factor(spect_clust_label, levels=c("Anterior", "Posterior", "Periventricular"))) %>%
# Plot
  ggplot(aes(x=other_parc_label, y=spect_clust_label, fill=dice)) +
    geom_tile(color = "white", lwd = 0.1, linetype = 1) +
    geom_text(aes(label = round(dice, 2)), vjust = 1, size = 5, color = "white") +
    scale_fill_gradientn(colors = viridis_pal(direction=1)(9), limits=c(0, 1), breaks=c(0, 1), name="Dice") +
    # facet_wrap(~ Other_parc, nrow=1) +
    facet_grid(. ~ Other_parc, scales="free", space="free") +
    theme_light() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 15),
      legend.text = element_text(size=15), legend.title = element_text(size=15),
      axis.text.y = element_text(size = 15),
      axis.title.x = element_blank(), axis.title.y = element_blank(),
      strip.text = element_text(size=15, color="black"), strip.background.x = element_rect(fill="#ffffff"))
  ggsave("./visualization/2e_overlap_other_parc/k3_overlap_other_parc.png", height=5, width=13)

