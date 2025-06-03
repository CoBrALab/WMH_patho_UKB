
library(data.table)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(patchwork)

# Plot spatial correlations in a graph

# Load data

df = as.data.frame(fread("./results/6a_connected_gm_ab_tau/corr_between_maps.csv"))

df = df %>%
    mutate(
        x = factor(x),
        y = factor(y, levels=c("clust1", "clust2", "clust3"), labels=c("Periventricular", "Posterior", "Anterior"))
    ) %>%
    glimpse()

# Clean df and bonferroni correction

df_pad = df %>%
    mutate(
        x = factor(x, levels=c(
            "preventad_AB_ABneg_ses02_dkt",
            "preventad_AB_ABpos_ses02_dkt",
            "preventad_tau_ABneg_ses02_dkt",
            "preventad_tau_ABpos_ses02_dkt"
        ), labels=c(
            "Amyloid PET - AB-",
            "Amyloid PET - AB+",
            "Tau PET - AB-",
            "Tau PET - AB+"
        ))
    ) %>%
    mutate(
        pearson_p_bonf = p.adjust(pearson_p, "bonferroni")
    ) %>%
    mutate(
        pearson_psig = if_else(pearson_p < 0.05, 1, 0),
        pearson_p_bonfsig = if_else(pearson_p_bonf < 0.05, 1, 0),
        pearson_sig = ifelse(pearson_p_bonf < 0.05, 2, ifelse(pearson_p<0.05, 1, 0))
    ) %>%
    glimpse()

fwrite(df_pad, "./results/6a_connected_gm_ab_tau/corr_between_maps_bonf.csv")

# Bar graph, color by AD map, group by WMH cluster
ggplot(df_pad, aes(x=y, y=pearson_corr, group=x, fill=x)) +
    geom_col(position = position_dodge(0.8), width=0.8) +
    geom_vline(xintercept=seq(0,length(levels(df_pad$y)))+0.5 ,color="black", alpha=0.2) +
    geom_hline(yintercept=0, color="black", size=1.5) +
    scale_y_continuous(name="\n\nPearson r", breaks=c(-0.6,-0.4,-0.2,0,0.2,0.4,0.6)) +
    scale_x_discrete(name="") +
    scale_fill_manual(name="", values=c( "#a4e3ce", "#3e7d68", "#8c89ae", "#403a78")) +
    theme_classic() + 
    theme(
        axis.text.x = element_text(angle=30, hjust=1),
        text = element_text(size=30), plot.title = element_text(hjust = 0.5, size=20, face="bold"),
        legend.title = element_blank(), legend.position="right", legend.text=element_text(size=15),
        strip.text = element_text(size = 20, color="black", face="bold"), strip.background = element_rect(fill="#ffffff")
    ) 
ggsave("./visualization/6a_connected_gm_ab_tau/pearson_bar.png", width=10, height=10)

