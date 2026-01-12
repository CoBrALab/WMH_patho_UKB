
library(ggplot2)
library(data.table)
library(dplyr)
library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(ggpubr)
library(patchwork)
library(multcomp)
library(rstatix)

# Load data

df_clust_micro = as.data.frame(fread("results/5g_pred_ad_stroke/spect_clust_k3_micro_metrics.csv"))
df_clust_micro$run = "clust_vol_micro"
df_clust_micro = df_clust_micro %>% mutate(Fold = row_number())

df_clust_vol = as.data.frame(fread("results/5g_pred_ad_stroke/spect_clust_k3_vol_metrics.csv"))
df_clust_vol$run = "clust_vol"
df_clust_vol = df_clust_vol %>% mutate(Fold = row_number())

df_pvdeep_micro = as.data.frame(fread("results/5g_pred_ad_stroke/pv_deep_micro_metrics.csv"))
df_pvdeep_micro$run = "pvdeep_vol_micro"
df_pvdeep_micro = df_pvdeep_micro %>% mutate(Fold = row_number())

df_pvdeep_vol = as.data.frame(fread("results/5g_pred_ad_stroke/pv_deep_vol_metrics.csv"))
df_pvdeep_vol$run = "pvdeep_vol"
df_pvdeep_vol = df_pvdeep_vol %>% mutate(Fold = row_number())

df_noparc_micro = as.data.frame(fread("results/5g_pred_ad_stroke/noparc_micro_metrics.csv"))
df_noparc_micro$run = "noparc_vol_micro"
df_noparc_micro = df_noparc_micro %>% mutate(Fold = row_number())

df_noparc_vol = as.data.frame(fread("results/5g_pred_ad_stroke/noparc_vol_metrics.csv"))
df_noparc_vol$run = "noparc_vol"
df_noparc_vol = df_noparc_vol %>% mutate(Fold = row_number())

df_prs = as.data.frame(fread("results/5g_pred_ad_stroke/all_prs_metrics.csv"))
df_prs$run = factor(df_prs$run, levels=c("noparc_vol", "pvdeep_vol", "clust_vol", "noparc_vol_micro", "pvdeep_vol_micro", "clust_vol_micro"))

# Bind and format

df = rbind(df_clust_micro, df_clust_vol, df_pvdeep_micro, df_pvdeep_vol, df_noparc_micro, df_noparc_vol)

df = df %>%
    mutate(
        run = factor(run, levels=c("noparc_vol", "pvdeep_vol", "clust_vol", "noparc_vol_micro", "pvdeep_vol_micro", "clust_vol_micro"))
    )

# Plot performance metrics

metrics = c("roc_auc", "balanced_accuracy", "f1")
metrics_clean = c("AUROC", "Balanced accuracy", "F1")

plts = list()
plts_prs = list()

for (m in 1:length(metrics)) {

    df_metric = df %>%
        rename(Metric = all_of(metrics[m])) %>%
        select(Metric, run, Fold)

    df_prs_metric = df_prs %>%
        rename(Metric = all_of(metrics[m])) %>%
        select(Metric, run)

    df_summary = df_metric %>%
        group_by(run) %>%
        summarise(
            mean_metric = mean(Metric),
            sd_metric = sd(Metric),
            .groups = "drop"
        )

    fwrite(df_summary, paste0("./results/5g_pred_ad_stroke/all_",metrics[m],"_mean_sd.tsv"), row.names = FALSE, col.names=TRUE, quote=FALSE, sep="\t")

    # Calculate significant differences with Wilcoxon test
    stat_test = df_metric %>%
        wilcox_test(
            Metric ~ run,
            paired = TRUE
        ) %>%
        adjust_pvalue(method = "bonferroni") %>%
        add_significance("p.adj")

    fwrite(stat_test, paste0("./results/5g_pred_ad_stroke/all_",metrics[m],"_stats.tsv"), row.names = FALSE, col.names=TRUE, quote=FALSE, sep="\t")

    # Plot
    plts[[m]] = ggplot(df_metric, aes(
            x=run,
            y=Metric,
            color = run
        )) +
        geom_violin(trim=TRUE, scale="width") + 
        geom_boxplot(width=0.2) +
        scale_x_discrete(name="", labels=c("No parc (vol)", "PV/deep (vol)", "Cluster (vol)", "No parc (vol + micro)", "PV/deep (vol + micro)", "Clusters (vol + micro)")) + 
        scale_y_continuous(name=metrics_clean[m]) +
        scale_color_manual(values = c("#808080", "#d99a9a", "#92a8df", "#000000","#b23535", "#2451bf"), guide="none") +
        ggtitle(metrics_clean[m]) +
        theme_light() +
        theme(
            axis.text.x = element_text(angle = 30, hjust = 1), text = element_text(size=20),
            plot.title = element_text(hjust = 0.5)
            )
    ggsave(paste0("./visualization/5g_pred_ad_stroke/",metrics[m],".png"), width=5, height=5)
    print(paste0("./visualization/5g_pred_ad_stroke/",metrics[m],".png"))

    plts_prs[[m]] = ggplot(df_metric, aes(
            x=run,
            y=Metric,
            color = run
        )) +
        geom_point(data=df_prs_metric, size=5) +
        scale_x_discrete(name="", labels=c("No parc (vol)", "PV/deep (vol)", "Cluster (vol)", "No parc (vol + micro)", "PV/deep (vol + micro)", "Clusters (vol + micro)")) + 
        scale_y_continuous(name=metrics_clean[m]) +
        scale_color_manual(values = c("#808080", "#d99a9a", "#92a8df", "#000000","#b23535", "#2451bf"), guide="none") +
        ggtitle(metrics_clean[m]) +
        theme_light() +
        theme(
            axis.text.x = element_text(angle = 30, hjust = 1), text = element_text(size=20),
            plot.title = element_text(hjust = 0.5)
            )
    ggsave(paste0("./visualization/5g_pred_ad_stroke/",metrics[m],"_prs.png"), width=5, height=5)
    print(paste0("./visualization/5g_pred_ad_stroke/",metrics[m],"_prs.png"))
}

wrap_plots(c(plts, plts_prs), nrow = 2)
ggsave(paste0("./visualization/5g_pred_ad_stroke/all_metrics.png"), width=6*length(plts), height=14)

