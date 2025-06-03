
library(data.table)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(scales)
library(forcats)
library(cowplot)
library(MatchIt)
library(patchwork)
library(MRIcrotome)
library(RMINC)
library(viridis)
library(RColorBrewer)
library(effsize)

# Compare the WMHs of cases (top 1% of PRS for certain disease) and controls (bottom 50% of PRS) with linear models

dir.create("./results/5c_UKB_PRS", showWarnings=FALSE)
dir.create("./visualization/5c_UKB_PRS", showWarnings=FALSE)

# Load and clean PRS data

prs = as.data.frame(fread("../../../UKB/tabular/PRS/ukb_category_301.csv"))

colnames(prs) <- gsub("-0.0", "", colnames(prs), fixed = TRUE)

field_ids = as.data.frame(fread("../../../UKB/tabular/field.txt"))
field_ids = field_ids %>% select(c("field_id", "title"))
field_ids$title <- gsub("[ ,'-]", "_", field_ids$title)
field_ids$title <- gsub("[()]", "", field_ids$title)

for (i in 2:length(colnames(prs))) {
    print(colnames(prs)[i])
    field_ids_tmp = field_ids %>% filter(field_id == colnames(prs)[i])
    colnames(prs)[i] = field_ids_tmp[,2]
}

colnames(prs) = gsub("Standard_PRS_for_", "", colnames(prs))
colnames(prs) = gsub("_[^_]+$", "", colnames(prs))

# Merge with demo and keep only IDs included

demo = as.data.frame(fread("../clean_wmh_data/results/demo.tsv"))

wmh = as.data.frame(fread("../3_temporal_sustain/results/3f_clean_wmh_data/wmh_combined.tsv"))

ids = as.data.frame(wmh$ID)
colnames(ids) = "ID"
colnames(prs)[1] = "ID"

prs = merge(merge(ids, demo, by="ID"), prs, by="ID")

# Remove NAs
prs = prs[,!(colSums(is.na(prs)) == nrow(prs))]

prs = prs[complete.cases(prs),]

# save 
fwrite(prs, "./results/5c_UKB_PRS/prs.tsv", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

# Load data

prs = as.data.frame(fread("./results/5c_UKB_PRS/prs.tsv"))
prs$sex = as.factor(prs$sex)

wmh = as.data.frame(fread("../3_temporal_sustain/results/3f_clean_wmh_data/wmh_combined.tsv"))
wmh = wmh %>% mutate(sex = factor(sex))

# Compare extreme groups

results_prs = as.data.frame(matrix(ncol=8, nrow=0))
colnames(results_prs) = c("PRS", "WMH_measure", "cluster", "lm_beta", "lm_pval")

i=1
# For each PRS of different diseases...
for (p in 1:ncol(prs)) {
    print(colnames(prs)[p])
    plots = list()

    g=1
    # For each WMH measure
    for (w in 4:ncol(wmh)) {

    df_tmp = merge(prs, wmh, by="ID")

    df_tmp = df_tmp %>% 
        select(c(colnames(wmh)[w], colnames(prs)[p])) %>%
        mutate(group = factor(ifelse(pull(., 2) < median(pull(., 2)), "low",
                        ifelse(pull(., 2) > as.numeric(quantile(pull(., 2), probs=0.99)), "high",
                        "mid")))) %>%
        filter(group != "mid") %>%
        mutate(group = factor(group, levels=c("high", "low")), WMH = scale(pull(., 1)))

        linmod = lm(formula, df_tmp %>% mutate(group = factor(group, levels=c("low", "high"))))
        linmod = summary(linmod)

        str = str_split(colnames(wmh)[w], "_", simplify = TRUE)

        results_prs[i,] = c(colnames(prs)[p], str[1], str[2], linmod$coefficients[2,1], linmod$coefficients[2,4])

        i = i + 1
    }
}

# Significance and FDR-correction
results_prs = results_prs %>%
        mutate(pval_fdr = p.adjust(lm_pval, method="fdr")) %>%
        mutate(pval_nsig_05 = ifelse(lm_pval<=0.05, TRUE, FALSE)) %>%
        mutate(pval_fdr_nsig_05 = ifelse(pval_fdr<=0.05, TRUE, FALSE))

fwrite(results_prs, "./results/5c_UKB_PRS/WMH_diff_UKB_PRS.tsv", row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

# Calculate average betas for each region across markers for dx in paper

avg_beta_region = results_prs %>% group_by(PRS, cluster) %>%
    summarize(avg = mean(lm_beta))

fwrite(as.data.frame(avg_beta_region), "./results/5c_UKB_PRS/avg_beta_region.tsv", row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

# Load results
results_prs = as.data.frame(fread("./results/5c_UKB_PRS/WMH_diff_UKB_PRS.tsv"))
results_prs[,c(1,2,3)] = as.data.frame(lapply(results_prs[,c(1,2,3)], as.factor))

# Visualize on brain

names=c("WMHvol", "stage", "MD", "ISOVF", "FA", "ICVF", "OD", "T2star", "QSM")

anatVol = mincArray(mincGetVolume("../../data/UKB_template_2mm.mnc"))
viridis_scale = viridis(n=255, option="B")
color_scale_div_2 = colorRampPalette(brewer.pal(9,"Blues"))(255)
color_scale_div_1 = colorRampPalette(brewer.pal(9,"Reds"))(255)

# Prep function
viz_on_brain = function(df, out_name) {

    # For each PRS of different diseases
    for (i in 1:length(levels(df$PRS))) {

        effects_mnc = list()
        mnc_parc = mincGetVolume("../../2_spatial_clust/results/2d_assign_excluded_voxels/final_parc_k3.mnc")

        # For each WMH measure
        for (n in 1:length(names)) {

            to_mnc = df %>%
                filter(WMH_measure == names[n] & PRS == levels(df$PRS)[i]) %>%
                select("lm_beta")
            to_mnc = as.numeric(to_mnc[,1])

            mnc_tmp = mnc_parc
            mnc_tmp[][which(mnc_tmp[]==1)] = to_mnc[1]
            mnc_tmp[][which(mnc_tmp[]==2)] = to_mnc[2]
            mnc_tmp[][which(mnc_tmp[]==3)] = to_mnc[3]

            effects_mnc[[n]]=mincArray(mnc_tmp)
        }

        # Visualize results on brain map
        png(file=paste0("./visualization/5b_UKB_disorders/",out_name,"_",levels(df$PRS)[i],".png"), width=7500, height=1000, pointsize = 150)
        sliceSeries(nrow=1, ncol=1, begin=52, end=52, dimension=3) %>%
            addtitle(names[1]) %>%
            anatomy(anatVol, low=10, high=200) %>%
            overlay(effects_mnc[[1]],
                low=0.0001, high=0.4, col=color_scale_div_1, rCol=color_scale_div_2, symmetric=TRUE) %>%
        sliceSeries(nrow=1, ncol=1, begin=52, end=52, dimension=3) %>%
            addtitle(names[2]) %>%
            anatomy(anatVol, low=10, high=200) %>%
            overlay(effects_mnc[[2]],
                low=0.0001, high=0.4,col=color_scale_div_1, rCol=color_scale_div_2, symmetric=TRUE) %>%
        sliceSeries(nrow=1, ncol=1, begin=52, end=52, dimension=3) %>%
            addtitle(names[3]) %>%
            anatomy(anatVol, low=10, high=200) %>%
            overlay(effects_mnc[[3]],
                low=0.0001, high=0.4,col=color_scale_div_1, rCol=color_scale_div_2, symmetric=TRUE) %>%
        sliceSeries(nrow=1, ncol=1, begin=52, end=52, dimension=3) %>%
            addtitle(names[4]) %>%
            anatomy(anatVol, low=10, high=200) %>%
            overlay(effects_mnc[[4]],
                low=0.0001, high=0.4,col=color_scale_div_1, rCol=color_scale_div_2, symmetric=TRUE) %>%
        sliceSeries(nrow=1, ncol=1, begin=52, end=52, dimension=3) %>%
            addtitle(names[5]) %>%
            anatomy(anatVol, low=10, high=200) %>%
            overlay(effects_mnc[[5]],
                low=0.0001, high=0.4,col=color_scale_div_1, rCol=color_scale_div_2, symmetric=TRUE) %>%
        sliceSeries(nrow=1, ncol=1, begin=52, end=52, dimension=3) %>%
            addtitle(names[6]) %>%
            anatomy(anatVol, low=10, high=200) %>%
            overlay(effects_mnc[[6]],
                low=0.0001, high=0.4,col=color_scale_div_1, rCol=color_scale_div_2, symmetric=TRUE) %>%
        sliceSeries(nrow=1, ncol=1, begin=52, end=52, dimension=3) %>%
            addtitle(names[7]) %>%
            anatomy(anatVol, low=10, high=200) %>%
            overlay(effects_mnc[[7]],
                low=0.0001, high=0.4,col=color_scale_div_1, rCol=color_scale_div_2, symmetric=TRUE) %>%
        sliceSeries(nrow=1, ncol=1, begin=52, end=52, dimension=3) %>%
            addtitle(names[8]) %>%
            anatomy(anatVol, low=10, high=200) %>%
            overlay(effects_mnc[[8]],
                low=0.0001, high=0.4,col=color_scale_div_1, rCol=color_scale_div_2, symmetric=TRUE) %>%
        sliceSeries(nrow=1, ncol=1, begin=52, end=52, dimension=3) %>%
            addtitle(names[9]) %>%
            anatomy(anatVol, low=10, high=200) %>%
            overlay(effects_mnc[[9]],
                low=0.0001, high=0.4,col=color_scale_div_1, rCol=color_scale_div_2, symmetric=TRUE) %>%
        legend("Effect size") %>%
        draw()
        dev.off()
    }
}

# Run
viz_on_brain(results_prs, "prs_")

