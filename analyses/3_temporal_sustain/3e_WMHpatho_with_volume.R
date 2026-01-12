
library(RMINC)
library(data.table)
library(matrixStats)
library(tidyverse)
library(ggplot2)
library(splines)
library(effects)
library(viridis)
library(scales)
library(patchwork)
library(grid)

# Associations between WMH pathophysiology and WMH volume 
# (as an empirical temporal estimate of WMH progression)

dir.create("./results/3e_WMHpatho_with_volume", showWarnings=FALSE)
dir.create("./visualization/3e_WMHpatho_with_volume", showWarnings=FALSE)

names=c("FA", "MD", "ICVF", "ISOVF", "OD", "T2star", "QSM")

# Define functions

# All micro z-score in one plot (all positive abnormality)
all_micro_abs = function(merged_df, spline_df, color_scale, png_name) {

    merged_df$FA = merged_df$FA*-1
    merged_df$ICVF = merged_df$ICVF*-1
    merged_df$OD = merged_df$OD*-1

    plot = ggplot(data=merged_df, aes(x=log(WMH))) + 
        stat_smooth(aes(y=merged_df[[names[1]]], colour=names[1]), method = lm, formula = y ~ bs(x,spline_df), se=TRUE) + 
        stat_smooth(aes(y=merged_df[[names[2]]], colour=names[2]), method = lm, formula = y ~ bs(x,spline_df), se=TRUE) + 
        stat_smooth(aes(y=merged_df[[names[3]]], colour=names[3]), method = lm, formula = y ~ bs(x,spline_df), se=TRUE) + 
        stat_smooth(aes(y=merged_df[[names[4]]], colour=names[4]), method = lm, formula = y ~ bs(x,spline_df), se=TRUE) + 
        stat_smooth(aes(y=merged_df[[names[5]]], colour=names[5]), method = lm, formula = y ~ bs(x,spline_df), se=TRUE) + 
        stat_smooth(aes(y=merged_df[[names[6]]], colour=names[6]), method = lm, formula = y ~ bs(x,spline_df), se=TRUE) + 
        stat_smooth(aes(y=merged_df[[names[7]]], colour=names[7]), method = lm, formula = y ~ bs(x,spline_df), se=TRUE) + 
        geom_hline(yintercept = 0) + 
        scale_colour_manual(name="Micro", values=color_scale) + 
        scale_y_continuous(breaks=seq(-10,10,by=2), limits=c(-5,9), name="Z-score") + 
        scale_x_continuous(breaks=seq(3,12,by=2), limits=c(3,12), name="log(WMH volume)") + 
        labs(colour = "Micro") + 
        guides(color=guide_legend(override.aes=list(fill=NA))) + 
        theme_classic() + 
        ggtitle("All markers\n(Positive abnormality)") + 
        theme(text=element_text(size=30, color="black"), plot.title=element_text(hjust=0.5, size=40))
    ggsave(png_name, width = 10, height = 8, dpi=300, units="in")

    return(plot)
}


# Micro by biological sensitivity (fluid, fiber, myelin/iron)
by_sensitivity = function(merged_df, spline_df, color_scale, png_name) {
        
    fluid = ggplot(data=merged_df, aes(x=log(WMH))) + 
        geom_hline(yintercept = 0) + 
        geom_point(aes(y=merged_df[[names[2]]], colour=names[2]), alpha=0.1, size=0.1) + 
        geom_point(aes(y=merged_df[[names[4]]], colour=names[4]), alpha=0.1, size=0.1) + 
        stat_smooth(aes(y=merged_df[[names[2]]], colour=names[2]), method = lm, formula = y ~ bs(x,spline_df), se=TRUE, size=1) + 
        stat_smooth(aes(y=merged_df[[names[4]]], colour=names[4]), method = lm, formula = y ~ bs(x,spline_df), se=TRUE, size=1) + 
        labs(y="Z-score", x="log(WMH volume)") + 
        scale_colour_manual(name="", values=color_scale[c(1,2)]) + 
        scale_y_continuous(breaks=seq(-10,10,by=2), limits=c(-5,9)) + 
        scale_x_continuous(breaks=seq(3,12,by=1), limits=c(3,12), name="log(WMH volume)") + 
        labs(colour = "") + 
        theme_classic() + 
        ggtitle("Fluid-sensitive") + 
        theme(text=element_text(size=30, color="black"), plot.title=element_text(hjust=0.5, size=40),
                legend.position = c(0.15, 0.85)) + 
        guides(color = guide_legend(override.aes = list(size = 10, fill=NA)))
    ggsave(paste0(png_name, "_fluid.png"), width = 10, height = 8, dpi=300, units="in")

    fiber = ggplot(data=merged_df, aes(x=log(WMH))) + 
        geom_hline(yintercept = 0) + 
        geom_point(aes(y=merged_df[[names[1]]], colour=names[1]), alpha=0.1, size=0.1) + 
        geom_point(aes(y=merged_df[[names[3]]], colour=names[3]), alpha=0.1, size=0.1) + 
        geom_point(aes(y=merged_df[[names[5]]], colour=names[5]), alpha=0.1, size=0.1) + 
        stat_smooth(aes(y=merged_df[[names[1]]], colour=names[1]), method = lm, formula = y ~ bs(x,spline_df), se=TRUE, size=1) + 
        stat_smooth(aes(y=merged_df[[names[3]]], colour=names[3]), method = lm, formula = y ~ bs(x,spline_df), se=TRUE, size=1) + 
        stat_smooth(aes(y=merged_df[[names[5]]], colour=names[5]), method = lm, formula = y ~ bs(x,spline_df), se=TRUE, size=1) + 
        labs(y="", x="log(WMH volume)") + 
        scale_colour_manual(name="", values=color_scale[c(3,4,5)], limits=c(names[1], names[3], names[5])) + 
        scale_y_continuous(breaks=seq(-10,10,by=2), limits=c(-5,9)) + 
        scale_x_continuous(breaks=seq(3,12,by=2), limits=c(3,12), name="log(WMH volume)") + 
        labs(colour = "") + 
        theme_classic() + 
        ggtitle("Fiber-sensitive") + 
        theme(text=element_text(size=30, color="black"), plot.title=element_text(hjust=0.5, size=40),
                legend.position = c(0.15, 0.85)) + 
        guides(color = guide_legend(override.aes = list(size = 10, fill=NA)))
    ggsave(paste0(png_name, "_fiber.png"), width = 10, height = 8, dpi=300, units="in")

    iron_myelin = ggplot(data=merged_df, aes(x=log(WMH))) + 
        geom_hline(yintercept = 0) + 
        geom_point(aes(y=merged_df[[names[6]]], colour=names[6]), alpha=0.1, size=0.1) + 
        geom_point(aes(y=merged_df[[names[7]]], colour=names[7]), alpha=0.1, size=0.1) + 
        stat_smooth(aes(y=merged_df[[names[6]]], colour=names[6]), method = lm, formula = y ~ bs(x,spline_df), se=TRUE, size=1) + 
        stat_smooth(aes(y=merged_df[[names[7]]], colour=names[7]), method = lm, formula = y ~ bs(x,spline_df), se=TRUE, size=1) + 
        geom_hline(yintercept = 0) + 
        labs(y="", x="log(WMH volume)") + 
        scale_colour_manual(name="", values=color_scale[c(6,7)]) + 
        scale_y_continuous(breaks=seq(-10,10,by=2), limits=c(-5,9)) + 
        scale_x_continuous(breaks=seq(3,12,by=2), limits=c(3,12), name="log(WMH volume)") + 
        labs(colour = "") + 
        theme_classic() + 
        ggtitle("Myelin and iron-sensitive") + 
        theme(text=element_text(size=30, color="black"), plot.title=element_text(hjust=0.5, size=40),
                legend.position = c(0.15, 0.85)) + 
        guides(color = guide_legend(override.aes = list(size = 10, fill=NA)))
    ggsave(paste0(png_name, "_myelin_iron.png"), width = 10, height = 8, dpi=300, units="in")

    plots = list(fluid, fiber, iron_myelin)
    return(plots)

}

# Run

# Load median WMH z-scores per subject
# df = as.data.frame(fread("../2_spatial_clust/results/2h_roi_avg/subj_WMH_patho.tsv"))
df = as.data.frame(fread("../../../WMH_micro_spatial/Analyses_nm/zscore_subject_regions/results/spect_clust_k3/subj_WMH_patho_ses2.tsv"))

total_clusters = 3

for (c in 1:total_clusters) {
    df_clust = df[, c(1, grep(paste0("c",c), names(df)))]
    df_clust = df_clust[, seq(1, ncol(df_clust)-1)]
    colnames(df_clust) = c("ID", names)
    
    df_clust = cbind(df_clust, df$WMH)
    colnames(df_clust)[length(colnames(df_clust))] = "WMH"

    color_scale = c(MD = "#04319E", ISOVF = "#8298CF", FA = "#2B520B", ICVF="#448312", OD="#A2C189", T2star="#D36108", QSM="#E9B084")
    spline_df = 4

    # All micro z-score in one plot
    all_micro_plot_abs = all_micro_abs(df_clust, spline_df, color_scale, paste0("./visualization/3e_WMHpatho_with_volume/c",c,"_all_abs.png"))

    # By type of tissue sensitivity
    by_sensitivity_plot = by_sensitivity(df_clust, spline_df, color_scale, paste0("./visualization/3e_WMHpatho_with_volume/c",c,"_by_sens"))

    wrap_plots(by_sensitivity_plot)
    ggsave(paste0("./visualization/3e_WMHpatho_with_volume/c",c,"_by_sens_all.png"), width = 30, height = 8, dpi=300, units="in")

}

# Plots with sustain

spline_df = 4

df_sustain_c1 = as.data.frame(fread("../../../WMH_micro_spatial/Analyses_nm/sustain_spect_clust_k3_parc/results/c1/subtype0/comb_subj_cval_ses2.csv"))
df_sustain_c2 = as.data.frame(fread("../../../WMH_micro_spatial/Analyses_nm/sustain_spect_clust_k3_parc/results/c2/subtype0/comb_subj_cval_ses2.csv"))
df_sustain_c3 = as.data.frame(fread("../../../WMH_micro_spatial/Analyses_nm/sustain_spect_clust_k3_parc/results/c3/subtype0/comb_subj_cval_ses2.csv"))

df_volumes = df_sustain_c1 %>%
    select(ID, spect_clust_k3_WMHvol_c1, spect_clust_k3_WMHvol_c2, spect_clust_k3_WMHvol_c3) %>%
    rename(WMHvol_c1 = "spect_clust_k3_WMHvol_c1", WMHvol_c2 = "spect_clust_k3_WMHvol_c2", WMHvol_c3 = "spect_clust_k3_WMHvol_c3") %>%
    pivot_longer(cols=-ID, names_to = "Cluster", values_to = "WMHvol") %>%
    mutate(Cluster = factor(Cluster, levels=c("WMHvol_c1", "WMHvol_c2", "WMHvol_c3"), labels=c("Periventricular", "Posterior", "Anterior"))) %>%
    glimpse()

df_sustain = df_sustain_c1 %>%
    mutate(Cluster = "Periventricular") %>%
    bind_rows(df_sustain_c2 %>% mutate(Cluster = "Posterior")) %>%
    bind_rows(df_sustain_c3 %>% mutate(Cluster = "Anterior")) %>%    
    select(ID, ml_stage_cval, Cluster) %>%
    mutate(Cluster = factor(Cluster, levels=c("Periventricular", "Posterior", "Anterior"))) %>%
    left_join(df_volumes, by=c("ID", "Cluster")) %>%
    mutate(WMHvol = log(WMHvol + 1)) %>%
    glimpse()

cor_1 = cor.test(
    df_sustain %>% filter(Cluster == "Periventricular") %>% pull(ml_stage_cval),
    df_sustain %>% filter(Cluster == "Periventricular") %>% pull(WMHvol),
    method = "spearman"
    )

plot1 = ggplot(data=df_sustain %>% filter(Cluster == "Periventricular"), aes(x=WMHvol, y=ml_stage_cval)) + 
    stat_smooth(method = lm, formula = y ~ bs(x,spline_df), se=TRUE, color="black") + 
    geom_point(alpha=0.1, size=0.1) + 
    scale_y_continuous(breaks=seq(0,30,by=5), limits=c(0,31), name="Stage") + 
    scale_x_continuous(breaks=seq(0,12,by=2), limits=c(0,12), name="log(WMH volume)") + 
    theme_classic() + 
    ggtitle("Periventricular") + 
    theme(text=element_text(size=30, color="black"), plot.title=element_text(hjust=0.5, size=40)) +
    annotate(
        "text",
        x = 8, y = 5,
        label = paste0("rho = ", round(cor_1$estimate, 2)),
        hjust = 0, vjust = 1,
        size=10,
        color = "black"
    )

cor_2 = cor.test(
    df_sustain %>% filter(Cluster == "Posterior") %>% pull(ml_stage_cval),
    df_sustain %>% filter(Cluster == "Posterior") %>% pull(WMHvol),
    method = "spearman"
    )

plot2 = ggplot(data=df_sustain %>% filter(Cluster == "Posterior"), aes(x=WMHvol, y=ml_stage_cval)) + 
    stat_smooth(method = lm, formula = y ~ bs(x,spline_df), se=TRUE, color="black") + 
    geom_point(alpha=0.1, size=0.1) + 
    scale_y_continuous(breaks=seq(0,30,by=5), limits=c(0,31), name="Stage") + 
    scale_x_continuous(breaks=seq(0,12,by=2), limits=c(0,12), name="log(WMH volume)") + 
    theme_classic() + 
    ggtitle("Posterior") + 
    theme(text=element_text(size=30, color="black"), plot.title=element_text(hjust=0.5, size=40)) +
    annotate(
        "text",
        x = 8, y = 5,
        label = paste0("rho = ", round(cor_2$estimate, 2)),
        hjust = 0, vjust = 1,
        size=10,
        color = "black"
    )

cor_3 = cor.test(
    df_sustain %>% filter(Cluster == "Anterior") %>% pull(ml_stage_cval),
    df_sustain %>% filter(Cluster == "Anterior") %>% pull(WMHvol),
    method = "spearman"
    )

plot3 = ggplot(data=df_sustain %>% filter(Cluster == "Anterior"), aes(x=WMHvol, y=ml_stage_cval)) + 
    stat_smooth(method = lm, formula = y ~ bs(x,spline_df), se=TRUE, color="black") + 
    geom_point(alpha=0.1, size=0.1) + 
    scale_y_continuous(breaks=seq(0,30,by=5), limits=c(0,31), name="Stage") + 
    scale_x_continuous(breaks=seq(0,12,by=2), limits=c(0,12), name="log(WMH volume)") + 
    theme_classic() + 
    ggtitle("Anterior") + 
    theme(text=element_text(size=30, color="black"), plot.title=element_text(hjust=0.5, size=40)) +
    annotate(
        "text",
        x = 8, y = 5,
        label = paste0("rho = ", round(cor_3$estimate, 2)),
        hjust = 0, vjust = 1,
        size=10,
        color = "black"
    )

wrap_plots(list(plot1, plot2, plot3))
ggsave("./visualization/3e_WMHpatho_with_volume/sustain_with_volume.png", width = 30, height = 8, dpi=300, units="in")
