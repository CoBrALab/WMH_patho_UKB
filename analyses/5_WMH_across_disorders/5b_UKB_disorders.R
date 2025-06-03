
library(data.table)
library(effsize)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(scales)
library(forcats)
library(cowplot)
library(patchwork)
library(MRIcrotome)
library(magrittr)
library(viridis)
library(RColorBrewer)
library(RMINC)

# Compare the WMHs of cases (UKB ICD-10 diagnoses) and controls with linear models

dir.create("./results/5b_UKB_disorders", showWarnings=FALSE)
dir.create("./visualization/5b_UKB_disorders", showWarnings=FALSE)

# Load data

firstocc = as.data.frame(fread("../../../UKB/Analyses/clean_firstocc/results/firstocc_categ.tsv"))

inclusions = as.data.frame(fread("../../QC/inclusions_final.txt"))
colnames(inclusions) = "ID"
firstocc = merge(inclusions, firstocc, by="ID", all.x=TRUE)
firstocc[,c(3,4,7,8)] = as.data.frame(lapply(firstocc[,c(3,4,7,8)], as.factor))
firstocc[,c(2,5)] = as.data.frame(lapply(firstocc[,c(2,5)], as.Date))

wmh = as.data.frame(fread("../3_temporal_sustain/results/3f_clean_wmh_data/wmh_combined.tsv"))
wmh$sex = as.factor(wmh$sex)

# Manually-defined categories of ICD-10 codes
icd_codes_list <- readRDS("../../../UKB/Analyses/clean_firstocc/results/icd_codes_list.rds")

# Make control group: no lifetime dx on endocrine, mental, nervous, circulatory systems
# Exceptions
controls_exc = firstocc %>% filter(chapter %in% c("Endocrine, nutritional and metabolic diseases", "Mental, Behavioral and Neurodevelopmental disorders", "Diseases of the nervous system", "Diseases of the circulatory system"))

controls_exc_dx = table(controls_exc$dx_name)
controls_exc_dx = controls_exc_dx[order(controls_exc_dx)]
fwrite(as.data.frame(controls_exc_dx), "dx_prevalence.csv")

# dx to keep in controls (high prev and low impact on brain):
# depressive episode (F32; n=4163), migraine (G43; n=2632), other anxiety disorders (F41; n=2632), haemorrhoids (I84; n=1859)
controls_exc = controls_exc %>% filter(!icd_code %in% c("F32", "G43", "F41", "I84"))

controls = inclusions$ID[which(!inclusions$ID %in% controls_exc$ID)]

# Differences betwen cases and controls controling for age and sex

hc_df = as.data.frame(wmh[which(wmh$ID %in% controls),])
hc_df$Group = "HC"

results = list()
i=1

for (dx in 1:length(icd_codes_list)) {
    print(names(icd_codes_list)[dx])

    # Case dataframe: ICD-10 code in disease category and participant NOT in controls
    dx_ids = firstocc %>% filter(icd_code %in% icd_codes_list[[dx]]) %>% pull(ID)
    dx_df = as.data.frame(wmh[which(!wmh$ID %in% controls),])
    dx_df = dx_df[which(dx_df$ID %in% dx_ids),]
    n_dx_tmp = nrow(dx_df)

    dx_df$Group="DX"
    ncol_dx_df = ncol(dx_df)

    firstocc_df = firstocc %>% filter(icd_code %in% icd_codes_list[[dx]])
    dx_df = merge(dx_df, firstocc_df, by="ID")

    # dx_df can contains duplicates because of commorbidities in the same category
    # Still want to visualize the breakdown of all dx, but get the total n from unique values
    dx_df_unique = dx_df %>% distinct(ID, .keep_all = TRUE)

    hc_dx_df = rbind(hc_df, dx_df[,seq(1,ncol_dx_df)])
    hc_dx_df$Group = as.factor(hc_dx_df$Group)

    hc_dx_df_unique = rbind(hc_df, dx_df_unique[,seq(1,ncol_dx_df)])
    hc_dx_df_unique$Group = as.factor(hc_dx_df_unique$Group)
    
    # For every WMH measure
    for (w in c(seq(4,9), seq(12, 32))) {

        wmh_var_clust = strsplit(colnames(hc_dx_df_unique)[w], "_")
        wmh_var = wmh_var_clust[[1]][1]
        clust = wmh_var_clust[[1]][2]

        # Run linear model
        formula = as.formula(paste0("scale(", colnames(hc_dx_df_unique)[w], ") ~ Group + age + sex"))
        linmod = summary(lm(formula, data = hc_dx_df_unique))
        # Save results
        results[[i]] = data.frame(wmh_var = wmh_var, clust=clust,
                                categ=names(icd_codes_list)[dx],
                                n_dx = n_dx_tmp, n_hc = nrow(hc_df),
                                beta_groupHC = linmod$coefficients[2,1],
                                pval_groupHC = linmod$coefficients[2,4])
    
        i=i+1
    }
}

# Significance and FDR-correction
results = do.call(rbind, results)

results = results %>%
    mutate(categ = fct_reorder(categ, n_dx, .desc = TRUE, .fun = max)) %>%
    mutate(pval_fdr = p.adjust(pval_groupHC, method="fdr")) %>%
    mutate(pval_nsig_05 = ifelse(pval_groupHC<=0.05, TRUE, FALSE)) %>%
    mutate(pval_fdr_nsig_05 = ifelse(pval_fdr<=0.05, TRUE, FALSE))

fwrite(results, "./results/5b_UKB_disorders/WMH_diff_UKB_icd10.tsv", row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

# Calculate average betas for each region across markers

avg_beta_region = results %>% group_by(categ, clust) %>%
    summarize(avg = mean(beta_groupHC))

fwrite(as.data.frame(avg_beta_region), "./results/5b_UKB_disorders/avg_beta_region.tsv", row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

# Load results
metrics = c("WMHvol", "stage", "MD", "ISOVF", "FA", "ICVF", "OD", "T2star", "QSM")

results = as.data.frame(fread("./results/5b_UKB_disorders/WMH_diff_UKB_icd10.tsv"))
results = results %>% mutate(wmh_var = factor(wmh_var, levels=metrics), clust = factor(clust, levels=c("c1", "c2", "c3"), labels = c("Periventricular", "Posterior", "Anterior")),
                            categ = factor(categ))

# Visualize effect sizes on brain

names=c("WMHvol", "stage", "MD", "ISOVF", "FA", "ICVF", "OD", "T2star", "QSM")

anatVol = mincArray(mincGetVolume("../../data/UKB_template_2mm.mnc"))
viridis_scale = viridis(n=255, option="B")
color_scale_div_2 = colorRampPalette(brewer.pal(9,"Blues"))(255)
color_scale_div_1 = colorRampPalette(brewer.pal(9,"Reds"))(255)

for (dx_i in 1:length(levels(results$categ))) {

    # Results to mnc
    effects_mnc = list()
    mnc_parc = mincGetVolume("../../2_spatial_clust/results/2d_assign_excluded_voxels/final_parc_k3.mnc")

    for (n in 1:length(names)) {

        to_mnc = results %>%
            filter(wmh_var == names[n] & categ == levels(results$categ)[dx_i]) %>%
            pull(beta_groupHC)

        mnc_tmp = mnc_parc
        mnc_tmp[][which(mnc_tmp[]==1)] = to_mnc[1]
        mnc_tmp[][which(mnc_tmp[]==2)] = to_mnc[2]
        mnc_tmp[][which(mnc_tmp[]==3)] = to_mnc[3]

        effects_mnc[[n]]=mincArray(mnc_tmp)
    }

    # Visualize results on brain map
    name_file = paste0("./visualization/5b_UKB_disorders/",levels(results$categ)[dx_i], ".png")
    png(file=name_file, width=7500, height=1000, pointsize = 150)
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
    sliceSeries(nrow=1, ncol=1, begin=52, end=52, dimension=3) %>%
        addtitle(names[7]) %>%
        anatomy(anatVol, low=10, high=200) %>%
        overlay(effects_mnc[[7]],
            low=0.0001, high=0.7, col=color_scale_div_1, rCol=color_scale_div_2, symmetric=TRUE) %>%
    sliceSeries(nrow=1, ncol=1, begin=52, end=52, dimension=3) %>%
        addtitle(names[8]) %>%
        anatomy(anatVol, low=10, high=200) %>%
        overlay(effects_mnc[[8]],
            low=0.0001, high=0.7, col=color_scale_div_1, rCol=color_scale_div_2, symmetric=TRUE) %>%
    sliceSeries(nrow=1, ncol=1, begin=52, end=52, dimension=3) %>%
        addtitle(names[9]) %>%
        anatomy(anatVol, low=10, high=200) %>%
        overlay(effects_mnc[[9]],
            low=0.0001, high=0.7, col=color_scale_div_1, rCol=color_scale_div_2, symmetric=TRUE) %>%
    legend("Effect size") %>%
    draw()
    dev.off()
    print(name_file)

}
