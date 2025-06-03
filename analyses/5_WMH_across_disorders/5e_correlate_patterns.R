
library(data.table)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(viridis)
library(reshape2)
library(dendextend)
library(RColorBrewer)
library(cowplot)

# Correlate effect size patterns across analysis schemes (UKB diagnoses, UKB PRS, ADNI diagnoses)

dir.create("./results/5e_correlate_patterns", showWarnings=FALSE)
dir.create("./visualization/5e_correlate_patterns", showWarnings=FALSE)

# Load data

df_ukb_dx = as.data.frame(fread("./results/5b_UKB_disorders/WMH_diff_UKB_icd10.tsv"))

df_ukb_prs = as.data.frame(fread("./results/5c_UKB_PRS/WMH_diff_UKB_PRS.tsv"))

df_adni_dx = as.data.frame(fread("./results/5d_ADNI_dx/WMH_diff_ADNI_dx.tsv"))

# Clean data

df_ukb_dx = df_ukb_dx %>%
    select(c(wmh_var, clust, categ, beta_groupHC)) %>%
    rename(group="categ", effect="beta_groupHC") %>%
    mutate(group = factor(group, levels=c(
        "ischemic_heart_disease", "stroke", "dementia_no_vasc"
        ), labels=c(
            "UKB - IHD", "UKB - Stroke", "UKB - dementia"
        ))) %>%
    glimpse()

df_ukb_prs = df_ukb_prs %>%
    select(WMH_measure, cluster, PRS, lm_beta) %>%
    rename(wmh_var = "WMH_measure", clust="cluster", group="PRS", effect="lm_beta") %>%
    mutate(group = factor(group, levels=c(
        "cardiovascular_disease", "ischaemic_stroke", "alzheimer_s_disease"
    ), labels=c(
        "UKB - CVD PRS", "UKB - Stroke PRS", "UKB - AD PRS"
    ))) %>%
    mutate(clust = factor(clust, levels=c("c1", "c2", "c3"), labels=c("Periventricular", "Posterior", "Anterior"))) %>%
    glimpse()

df_adni_dx = df_adni_dx %>%
    select(WMH_var, clust, coef_mci, coef_ad) %>%
    pivot_longer(cols=c(coef_mci, coef_ad), names_to="group", values_to="effect") %>%
    rename(wmh_var = "WMH_var") %>%
    mutate(group = factor(group, levels=c(
        "coef_mci", "coef_ad"
    ), labels=c(
        "ADNI - MCI", "ADNI - AD"
    ))) %>%
    mutate(clust = factor(clust, levels=c("c1", "c2", "c3"), labels=c("Periventricular", "Posterior", "Anterior"))) %>%
    glimpse()

# Merge data
df_merge = rbind(df_ukb_dx, df_ukb_prs, df_adni_dx)

# Calculate correlations
results_corr = list()
i = 1
for (xvar in 1:length(unique(df_merge$group))) {
    df_x = df_merge %>% filter(group == unique(df_merge$group)[xvar])
    df_y = df_merge %>% filter(!group == unique(df_merge$group)[xvar])

    for (yvar in 1:length(unique(df_y$group))) {
        df_xy = merge(df_x, df_y %>% filter(group == unique(df_y$group)[yvar]), by=c("wmh_var", "clust"))

        corr_pearson = cor.test(df_xy$effect.x, df_xy$effect.y, method="pearson")
        corr_spear = cor.test(df_xy$effect.x, df_xy$effect.y, method="spearman")

        results_corr[[i]] = data.frame(
            x_var = unique(df_merge$group)[xvar],
            y_var = unique(df_y$group)[yvar],
            coef_pearson = corr_pearson$estimate,
            pval_pearson = corr_pearson$p.value
        )
        i = i + 1
    }
}

results_corr = do.call(rbind, results_corr)

# FDR correction
results_corr = results_corr %>%
    mutate(
        pval_pearson_fdr = p.adjust(pval_pearson, method="fdr")
    ) %>%
    glimpse()

fwrite(results_corr, "./results/5e_correlate_patterns/effects_pattern_correlations.tsv", row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

# Plot clustered correlation matrix

set.seed(123)
results_corr_clust = acast(results_corr, x_var ~ y_var, value.var = "coef_pearson")
results_corr_clust[is.na(results_corr_clust)] = 0
results_corr_clust = (results_corr_clust + t(results_corr_clust)) / 2

# Perform hierarchical clustering
effect_clusters = hclust(dist(results_corr_clust))
ordered_vars = rownames(results_corr_clust)[effect_clusters$order]

results_corr_clust = results_corr[, c("x_var", "y_var", "coef_pearson", "pval_pearson")]
results_corr_clust$x_var <- factor(results_corr_clust$x_var, levels = ordered_vars)
results_corr_clust$y_var <- factor(results_corr_clust$y_var, levels = ordered_vars)

# Order correlation matrix
corr_mat_clust = ggplot(results_corr_clust, aes(x = x_var, y = y_var, fill = coef_pearson)) +
    geom_tile(color = "white") +
    scale_fill_gradientn(name="Pearson r", colors=color_scale, values = scales::rescale(c(-1, 0, 1)), limits=c(-1,1)) +
    geom_text(aes(label = sprintf("%.2f", coef_pearson)), color = "black", size = 5) +
    theme_light() +
    theme(
        text = element_text(size=20),
        axis.text.x = element_text(angle = 30, hjust = 1),
        axis.title = element_blank(),
        legend.position = "none"
    )

# Dendogram
dend <- as.dendrogram(effect_clusters)
dendrogram <- ggplot(dend) +
    theme_classic() +
    theme(
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.margin = margin(t = 0, r = 0, b = -16, l = 0, unit = "pt")
    )

# Combine plots
combined_plot <- plot_grid(
  dendrogram + theme(plot.margin = margin(t = 0, r = 11.5, b = -16, l = 150, unit = "pt")), 
  corr_mat_clust + theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")), 
  ncol = 1, 
  rel_heights = c(1, 4)
)

ggsave("./visualization/5e_correlate_patterns/corr_matrix_clust_dend.png", width=7, height=9)

