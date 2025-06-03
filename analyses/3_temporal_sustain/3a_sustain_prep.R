

library(data.table)
library(RMINC)
library(tidyverse)
library(ggplot2)
library(ggridges)

# Prep data to use as sustain inputs
#   - Invert values for biomarkers which have negative values = more abnormality
#   - 0 where <5 voxels in region
#   - Calculate empirical abnormality thresholds

dir.create("./results/3a_sustain_prep", showWarnings=FALSE)
dir.create("./visualization/3a_sustain_prep", showWarnings=FALSE)

names=c("FA", "MD", "ICVF", "ISOVF", "OD", "T2star", "QSM")

mask=mincGetVolume("../../data/WMH_mask.mnc")

# Define parcellations

parc = fread("../clust_avg_spect/results/k3/final_parc_k3.tsv")
parc = parc$V1

parc_names = c("PV", "Posterior", "Anterior")

# Load median z-scores in regions
df = as.data.frame(fread(paste0("../2_spatial_clust/results/2h_roi_avg/subj_WMH_patho.tsv")))

# Invert values for FA, OD, ICVF (so that more = bad for all markers)
df[, grep("FA", colnames(df))] = df[, grep("FA", colnames(df))] * -1
df[, grep("ICVF", colnames(df))] = df[, grep("ICVF", colnames(df))] * -1
df[, grep("OD", colnames(df))] = df[, grep("OD", colnames(df))] * -1

for (i in sort(unique(parc[parc > 0]))) {
    cat(paste0("\t\tCluster = ", i, "\n"))

    # Where WMH volume < 5 voxels, impute WMH patho as 0
    mask = which(df[,grep(paste0("final_parc_k3_WMHvol_c",i), colnames(df))] < 5)
    print(length(mask))
    df[mask,head(grep(paste0("final_parc_k3_.*_c",i), colnames(df)),-1)] = 0
}

# Write (not available due to individual-level data)
fwrite(df, paste0("./results/3a_sustain_prep/subj_WMH_patho_clean.tsv"), col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

# Calculate cumulative probabilities
thresholds <- c(0.5, 1, 1.5, 2, 3, 4, 5, 6, 7)

calc_cumulative_prob <- function(x, threshold) {
    sum(!is.na(x) & x >= threshold) / sum(!is.na(x))
}

cum_prob = df %>% select(matches(paste0("final_parc_k3_.*_c.*"))) %>%
            select(-matches("WMHvol")) %>%
            pivot_longer(cols = everything()) %>%
            mutate(micro = factor(str_extract(name, "FA|MD|ICVF|ISOVF|OD|T2star|QSM"),
                                        levels = c("MD", "ISOVF", "FA", "ICVF", "OD", "T2star", "QSM"))) %>%
            mutate(cluster = str_extract(name, "(c[0-9])$")) %>%
            mutate(cluster = factor(cluster, levels = sort(unique(cluster), decreasing = FALSE), labels = parc_names[[p]])) %>%
            select(-1) %>%
            group_by(micro, cluster) %>%
            summarize_at(vars(value), list(
                cum_prob_z05 = ~calc_cumulative_prob(., thresholds[1]),
                cum_prob_z1 = ~calc_cumulative_prob(., thresholds[2]),
                cum_prob_z15 = ~calc_cumulative_prob(., thresholds[3]),
                cum_prob_z2 = ~calc_cumulative_prob(., thresholds[4]),
                cum_prob_z3 = ~calc_cumulative_prob(., thresholds[5]),
                cum_prob_z4 = ~calc_cumulative_prob(., thresholds[6]),
                cum_prob_z5 = ~calc_cumulative_prob(., thresholds[7]),
                cum_prob_z6 = ~calc_cumulative_prob(., thresholds[8]),
                cum_prob_z7 = ~calc_cumulative_prob(., thresholds[9])
            )) %>%
            ungroup()

fwrite(cum_prob, paste0("./results/3a_sustain_prep/cum_prob.tsv"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

# Visualize cumulative probabilities
cum_prob %>%
    gather(key = "variable", value = "value", -micro, -cluster) %>%
    ggplot(aes(x = variable, y = value, fill = variable)) +
        geom_bar(stat = "identity") +
        geom_hline(yintercept=0.01) +
        scale_fill_discrete(name="Threshold", labels=c(0.5, 1, 1.5, 2, 3, 4, 5, 6, 7)) + 
        scale_x_discrete(labels = c("0.5", "1", "1.5", "2", "3", "4", "5", "6", "7")) +
        scale_y_continuous(name = "Inverse cumulative probability", breaks = c(0.01, 0.1)) +
        coord_cartesian(ylim = c(0,0.1)) +
        facet_grid(cluster~micro, scale="free") +
        labs(x = " ", y = "Cum Prob", fill = "Thresholds") + 
        theme_light() +
        theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),
            strip.text = element_text(size =25), text = element_text(size=20))
    ggsave(paste0("./visualization/3a_sustain_prep/zscores_med_cum_prob.png"), width = 15, height = 10)
    print(paste0("./visualization/3a_sustain_prep/zscores_med_cum_prob.png"))
