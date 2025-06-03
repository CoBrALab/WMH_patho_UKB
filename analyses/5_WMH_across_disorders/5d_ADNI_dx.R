
library(data.table)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(RMINC)
library(MRIcrotome)
library(viridis)
library(RColorBrewer)

# Compare the WMHs of cases (ADNI participants with MCI or AD diagnosis) and controls (cognitively unimpaired) with linear models

dir.create("./results/5d_ADNI_dx", showWarnings=FALSE)
dir.create("./visualization/5d_ADNI_dx", showWarnings=FALSE)

# Load data

df_demo = as.data.frame(fread("../../../ADNI/Analyses/clean_data/results/df_demo.tsv"))

df_demo = df_demo %>%
    rename(ID="PTID", Sex = "sex", Age = "age_int", date="EXAMDATE") %>%
    glimpse()

df_wmh = as.data.frame(fread("../clean_wmh_data/results/spect_clust_k3/wmh_combined_ADNI.tsv"))

df = merge(df_demo, df_wmh %>% select(-c(Sex, Age, Year)), by=c("ID", "date"))

df = df %>%
    mutate(across(everything(), ~ na_if(., ""))) %>%
    select(where(~ !all(is.na(.)))) %>%
    mutate(
        date=as.Date(date),
        Sex=factor(Sex, levels=c("Male", "Female")),
        date_birth=as.Date(date_birth),
        dx=factor(dx, levels=c("CN", "MCI", "Dementia")),
        dx_mci_cause = factor(dx_mci_cause, levels=c("Due to AD", "Due to other")),
        dx_dementia_cause = factor(dx_dementia_cause, levels=c("Due to AD", "Due to other"))
    ) %>%
    # Remove non-AD cases
    filter(is.na(dx_mci_cause) | dx_mci_cause != "Due to other",
           is.na(dx_dementia_cause) | dx_dementia_cause != "Due to other") %>%
    glimpse()

# Keep last tp to maximize cognitively impaired participants
df_last_tp = df %>%
    group_by(ID) %>%
    slice_tail(n = 1) %>%
    ungroup

# Associations with dx

# Run linear models

# Columns of WMH measures
y_vars = seq(29, 46)

results = list()
# For each WMH measure
for (yvar in 1:length(y_vars)) {
    print(colnames(df_last_tp)[y_vars[yvar]])
    var_clust = strsplit(colnames(df_last_tp)[y_vars[yvar]], "_")

    df_tmp = df_last_tp %>%
        select(colnames(df_last_tp)[y_vars[yvar]], dx, Age, Sex) %>%
        filter(complete.cases(.))

    df_tmp_count = df_tmp %>%
        group_by(dx) %>%
        summarize(n=n())

    formula = as.formula(paste0("scale(",colnames(df_last_tp)[y_vars[yvar]], ") ~ dx + Age + Sex"))

    lm_dx = summary(lm(formula, data=df_tmp))

    results[[yvar]] = data.frame(
        WMH_var = var_clust[[1]][1],
        clust = var_clust[[1]][2],
        n_hc = df_tmp_count %>% filter(dx == "CN") %>% pull(n),
        n_mci = df_tmp_count %>% filter(dx == "MCI") %>% pull(n),
        n_ad = df_tmp_count %>% filter(dx == "Dementia") %>% pull(n),
        coef_mci = lm_dx$coefficients["dxMCI","Estimate"],
        coef_ad = lm_dx$coefficients["dxDementia","Estimate"],
        pval_mci = lm_dx$coefficients["dxMCI","Pr(>|t|)"],
        pval_ad = lm_dx$coefficients["dxDementia","Pr(>|t|)"]
    )
}

results = do.call(rbind, results)

# FDR correction

results = results %>%
    mutate(
        pval_mci_fdr = p.adjust(pval_mci, method="fdr"),
        pval_ad_fdr = p.adjust(pval_ad, method="fdr")
    )

# Write to tsv
fwrite(results, "./results/5d_ADNI_dx/WMH_diff_ADNI_dx.tsv", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

# Visualize on brain

vars_to_viz = c("coef_mci", "coef_ad")

names=c("WMHvol", "MD", "ISOVF", "FA", "ICVF", "OD")

anatVol = mincArray(mincGetVolume("../../data/UKB_template_2mm.mnc"))
viridis_scale = viridis(n=255, option="B")
color_scale_div_2 = colorRampPalette(brewer.pal(9,"Blues"))(255)
color_scale_div_1 = colorRampPalette(brewer.pal(9,"Reds"))(255)

for (var in vars_to_viz) {

    # Results to mnc
    effects_mnc = list()
    mnc_parc = mincGetVolume("../../2_spatial_clust/results/2d_assign_excluded_voxels/final_parc_k3.mnc")

    for (n in 1:length(names)) {
        to_mnc = results %>%
            filter(WMH_var == names[n]) %>%
            pull(var)

        mnc_tmp = mnc_parc
        mnc_tmp[][which(mnc_tmp[]==1)] = to_mnc[1]
        mnc_tmp[][which(mnc_tmp[]==2)] = to_mnc[2]
        mnc_tmp[][which(mnc_tmp[]==3)] = to_mnc[3]

        effects_mnc[[n]]=mincArray(mnc_tmp)
    }

    # Visualize results on brain map
    name_file = paste0("./visualization/5d_ADNI_dx/dx_",var, ".png")
    png(file=name_file, width=5100, height=1000, pointsize = 150)
    sliceSeries(nrow=1, ncol=1, begin=52, end=52, dimension=3) %>%
        addtitle(names[1]) %>%
        anatomy(anatVol, low=10, high=200) %>%
        overlay(effects_mnc[[1]],
            low=0.0001, high=0.7, col=color_scale_div_1, rCol=color_scale_div_2, symmetric=TRUE) %>%
    sliceSeries(nrow=1, ncol=1, begin=52, end=52, dimension=3) %>%
        addtitle(names[2]) %>%
        anatomy(anatVol, low=10, high=200) %>%
        overlay(effects_mnc[[2]],
            low=0.0001, high=0.7, col=color_scale_div_1, rCol=color_scale_div_2, symmetric=TRUE) %>%
    sliceSeries(nrow=1, ncol=1, begin=52, end=52, dimension=3) %>%
        addtitle(names[3]) %>%
        anatomy(anatVol, low=10, high=200) %>%
        overlay(effects_mnc[[3]],
            low=0.0001, high=0.7, col=color_scale_div_1, rCol=color_scale_div_2, symmetric=TRUE) %>%
    sliceSeries(nrow=1, ncol=1, begin=52, end=52, dimension=3) %>%
        addtitle(names[4]) %>%
        anatomy(anatVol, low=10, high=200) %>%
        overlay(effects_mnc[[4]],
            low=0.0001, high=0.7, col=color_scale_div_1, rCol=color_scale_div_2, symmetric=TRUE) %>%
    sliceSeries(nrow=1, ncol=1, begin=52, end=52, dimension=3) %>%
        addtitle(names[5]) %>%
        anatomy(anatVol, low=10, high=200) %>%
        overlay(effects_mnc[[5]],
            low=0.0001, high=0.7, col=color_scale_div_1, rCol=color_scale_div_2, symmetric=TRUE) %>%
    sliceSeries(nrow=1, ncol=1, begin=52, end=52, dimension=3) %>%
        addtitle(names[6]) %>%
        anatomy(anatVol, low=10, high=200) %>%
        overlay(effects_mnc[[6]],
            low=0.0001, high=0.7, col=color_scale_div_1, rCol=color_scale_div_2, symmetric=TRUE) %>%
    legend("Effect size") %>%
    draw()
    dev.off()
    print(name_file)
}

