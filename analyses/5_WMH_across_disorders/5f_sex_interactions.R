
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

dir.create("./results/5f_sex_interactions", showWarnings=FALSE)
dir.create("./visualization/5f_sex_interactions", showWarnings=FALSE)

#region Diagnosis

# Load data

firstocc = as.data.frame(fread("../../../UKB/Analyses/clean_firstocc/results/firstocc_categ.tsv"))

inclusions = as.data.frame(fread("../../../WMH_micro_spatial/QC/inclusions_final.txt"))
colnames(inclusions) = "ID"
firstocc = merge(inclusions, firstocc, by="ID", all.x=TRUE)
firstocc[,c(3,4,7,8)] = as.data.frame(lapply(firstocc[,c(3,4,7,8)], as.factor))
firstocc[,c(2,5)] = as.data.frame(lapply(firstocc[,c(2,5)], as.Date))

# wmh = as.data.frame(fread("../3_temporal_sustain/results/3f_clean_wmh_data/wmh_combined.tsv"))
wmh = as.data.frame(fread("../../../WMH_micro_spatial/Analyses_nm/clean_wmh_data/results/spect_clust_k3/wmh_combined.tsv"))
wmh$sex = as.factor(wmh$sex)

# Manually-defined categories of ICD-10 codes
icd_codes_list <- readRDS("../../../UKB/Analyses/clean_firstocc/results/icd_codes_list.rds")
icd_codes_list = icd_codes_list[c(
    "tian_ischemic_heart_disease", "tian_stroke", "cust_dementia_no_vasc"
)]
names(icd_codes_list) = c("ischemic_heart_disease", "stroke", "dementia_no_vasc")

# Make control group: no lifetime dx on endocrine, mental, nervous, circulatory systems
# Exceptions
controls_exc = firstocc %>% filter(chapter %in% c("Endocrine, nutritional and metabolic diseases", "Mental, Behavioral and Neurodevelopmental disorders", "Diseases of the nervous system", "Diseases of the circulatory system"))

# dx to keep in controls (high prev and low impact on brain):
# depressive episode (F32; n=4163), migraine (G43; n=2632), other anxiety disorders (F41; n=2632), haemorrhoids (I84; n=1859)
controls_exc = controls_exc %>% filter(!icd_code %in% c("F32", "G43", "F41", "I84"))

controls = inclusions$ID[which(!inclusions$ID %in% controls_exc$ID)]

# Differences betwen cases and controls with sex interactions, controlling for age

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
    hc_dx_df_unique$Group = factor(hc_dx_df_unique$Group, levels=c("HC", "DX"))
    hc_dx_df_unique$sex = factor(hc_dx_df_unique$sex, levels=c("Male", "Female"))

    # For every WMH measure
    for (w in c(seq(4,9), seq(12, 32))) {

        wmh_var_clust = strsplit(colnames(hc_dx_df_unique)[w], "_")
        wmh_var = wmh_var_clust[[1]][1]
        clust = wmh_var_clust[[1]][2]

        # Run linear model
        formula = as.formula(paste0("scale(", colnames(hc_dx_df_unique)[w], ") ~ Group * sex + age"))
        linmod = summary(lm(formula, data = hc_dx_df_unique))
        # Save results
        results[[i]] = data.frame(wmh_var = wmh_var, clust=clust,
                                categ=names(icd_codes_list)[dx],
                                n_dx = n_dx_tmp, n_hc = nrow(hc_df),
                                beta = linmod$coefficients["GroupDX:sexFemale","Estimate"],
                                pval = linmod$coefficients["GroupDX:sexFemale","Pr(>|t|)"])
    
        i=i+1
    }
}

# Significance and FDR-correction
results = do.call(rbind, results)

results = results %>%
    mutate(categ = fct_reorder(categ, n_dx, .desc = TRUE, .fun = max)) %>%
    mutate(pval_fdr = p.adjust(pval, method="fdr")) %>%
    mutate(pval_nsig_05 = ifelse(pval<=0.05, TRUE, FALSE)) %>%
    mutate(pval_fdr_nsig_05 = ifelse(pval_fdr<=0.05, TRUE, FALSE))

fwrite(results, "./results/5f_sex_interactions/sex_interactions_dx.tsv", row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

# Load results
metrics = c("WMHvol", "stage", "MD", "ISOVF", "FA", "ICVF", "OD", "T2star", "QSM")

results = as.data.frame(fread("./results/5f_sex_interactions/sex_interactions_dx.tsv"))
results = results %>% mutate(wmh_var = factor(wmh_var, levels=metrics), clust = factor(clust, levels=c("c1", "c2", "c3"), labels = c("Periventricular", "Posterior", "Anterior")),
                            categ = factor(categ))

# Visualize effect sizes on brain

names=c("WMHvol", "stage", "MD", "ISOVF", "FA", "ICVF", "OD", "T2star", "QSM")

anatVol = mincArray(mincGetVolume("../../data/UKB_template_T1_2mm.mnc"))
viridis_scale = viridis(n=255, option="B")
color_scale_div_2 = colorRampPalette(brewer.pal(9,"Blues"))(255)
color_scale_div_1 = colorRampPalette(brewer.pal(9,"Reds"))(255)

for (dx_i in 1:length(levels(results$categ))) {

    # Results to mnc
    effects_mnc = list()
    mnc_parc = mincGetVolume("../2_spatial_clust/results/2d_assign_excluded_voxels/final_parc_k3.mnc")

    for (n in 1:length(names)) {

        to_mnc = results %>%
            filter(wmh_var == names[n] & categ == levels(results$categ)[dx_i]) %>%
            pull(beta)

        mnc_tmp = mnc_parc
        mnc_tmp[][which(mnc_tmp[]==1)] = to_mnc[1]
        mnc_tmp[][which(mnc_tmp[]==2)] = to_mnc[2]
        mnc_tmp[][which(mnc_tmp[]==3)] = to_mnc[3]

        effects_mnc[[n]]=mincArray(mnc_tmp)
    }

    # Visualize results on brain map
    name_file = paste0("./visualization/5f_sex_interactions/dx_",levels(results$categ)[dx_i], ".png")
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

#endregion

#region PRS

# Load data

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

prs = prs %>% 
    filter(eid %in% inclusions$ID) %>%
    select(eid, alzheimer_s_disease, cardiovascular_disease, ischaemic_stroke) %>%
    rename(ID = "eid")

# Compare extreme groups

results_prs = as.data.frame(matrix(ncol=5, nrow=0))
colnames(results_prs) = c("PRS", "WMH_measure", "cluster", "lm_beta", "lm_pval")

i=1
# For each PRS of different diseases...
for (p in 2:ncol(prs)) {
    print(colnames(prs)[p])

    # For each WMH measure
    for (w in c(seq(4,9), seq(12, 32))) {

    df_tmp = merge(prs, wmh, by="ID")

    df_tmp = df_tmp %>% 
        select(ID, age, sex, all_of(colnames(wmh)[w]), all_of(colnames(prs)[p])) %>%
        rename(WMH = colnames(wmh)[w], PRS = colnames(prs)[p]) %>%
        filter(!is.na(PRS)) %>%
        mutate(group = factor(ifelse(PRS < median(PRS), "low",
                        ifelse(PRS > as.numeric(quantile(PRS, probs=0.99)), "high",
                        "mid")))) %>%
        filter(group != "mid") %>%
        mutate(
            group = factor(group, levels=c("low", "high")),
            sex = factor(sex, levels=c("Male", "Female")),
            WMH = scale(WMH)
        )

        formula = as.formula(paste("scale(WMH) ~ group * sex"))

        linmod = lm(formula, df_tmp %>% mutate(group = factor(group, levels=c("low", "high"))))
        linmod = summary(linmod)

        str = str_split(colnames(wmh)[w], "_", simplify = TRUE)

        results_prs[i,] = c(colnames(prs)[p], str[1], str[2], linmod$coefficients["grouphigh:sexFemale","Estimate"], linmod$coefficients["grouphigh:sexFemale","Pr(>|t|)"])

        i = i + 1
    }
}

# Significance and FDR-correction
results_prs = results_prs %>%
        mutate(pval_fdr = p.adjust(lm_pval, method="fdr")) %>%
        mutate(pval_nsig_05 = ifelse(lm_pval<=0.05, TRUE, FALSE)) %>%
        mutate(pval_fdr_nsig_05 = ifelse(pval_fdr<=0.05, TRUE, FALSE))

fwrite(results_prs, "./results/5f_sex_interactions/sex_interactions_prs.tsv", row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

# Load results
results_prs = as.data.frame(fread("./results/5f_sex_interactions/sex_interactions_prs.tsv"))
results_prs[,c(1,2,3)] = as.data.frame(lapply(results_prs[,c(1,2,3)], as.factor))

# Visualize on brain

names=c("WMHvol", "stage", "MD", "ISOVF", "FA", "ICVF", "OD", "T2star", "QSM")

anatVol = mincArray(mincGetVolume("../../data/UKB_template_T1_2mm.mnc"))
viridis_scale = viridis(n=255, option="B")
color_scale_div_2 = colorRampPalette(brewer.pal(9,"Blues"))(255)
color_scale_div_1 = colorRampPalette(brewer.pal(9,"Reds"))(255)

# For each PRS of different diseases

df = results_prs
out_name = "prs_"
for (i in 1:length(levels(df$PRS))) {

    effects_mnc = list()
    mnc_parc = mincGetVolume("../2_spatial_clust/results/2d_assign_excluded_voxels/final_parc_k3.mnc")

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
    png(file=paste0("./visualization/5f_sex_interactions/",out_name,levels(df$PRS)[i],".png"), width=7500, height=1000, pointsize = 150)
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

#endregion

#region ADNI

df_demo = as.data.frame(fread("../../../ADNI/Analyses/clean_data/results/df_demo.tsv"))

df_demo = df_demo %>%
    rename(ID="PTID", Sex = "sex", Age = "age_int", date="EXAMDATE") %>%
    glimpse()

# df_wmh = as.data.frame(fread("../clean_wmh_data/results/spect_clust_k3/wmh_combined_ADNI.tsv"))
df_wmh = as.data.frame(fread("../../../WMH_micro_spatial/Analyses_nm/clean_wmh_data/results/spect_clust_k3/wmh_combined_ADNI.tsv"))

df = merge(df_demo, df_wmh %>% select(-c(Sex, Age, Year)), by=c("ID", "date"))

df = df %>%
    mutate(across(where(is.character), ~ na_if(., ""))) %>%
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

    formula = as.formula(paste0("scale(",colnames(df_last_tp)[y_vars[yvar]], ") ~ dx * Sex + Age"))

    lm_dx = summary(lm(formula, data=df_tmp))

    results[[yvar]] = data.frame(
        WMH_var = var_clust[[1]][1],
        clust = var_clust[[1]][2],
        n_hc = df_tmp_count %>% filter(dx == "CN") %>% pull(n),
        n_mci = df_tmp_count %>% filter(dx == "MCI") %>% pull(n),
        n_ad = df_tmp_count %>% filter(dx == "Dementia") %>% pull(n),
        coef_mci = lm_dx$coefficients["dxMCI:SexFemale","Estimate"],
        coef_ad = lm_dx$coefficients["dxDementia:SexFemale","Estimate"],
        pval_mci = lm_dx$coefficients["dxMCI:SexFemale","Pr(>|t|)"],
        pval_ad = lm_dx$coefficients["dxDementia:SexFemale","Pr(>|t|)"]
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
fwrite(results, "./results/5f_sex_interactions/sex_interactions_ADNI.tsv", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

# Visualize on brain

vars_to_viz = c("coef_mci", "coef_ad")

names=c("WMHvol", "MD", "ISOVF", "FA", "ICVF", "OD")

anatVol = mincArray(mincGetVolume("../../data/UKB_template_T1_2mm.mnc"))
viridis_scale = viridis(n=255, option="B")
color_scale_div_2 = colorRampPalette(brewer.pal(9,"Blues"))(255)
color_scale_div_1 = colorRampPalette(brewer.pal(9,"Reds"))(255)

for (var in vars_to_viz) {

    # Results to mnc
    effects_mnc = list()
    mnc_parc = mincGetVolume("../2_spatial_clust/results/2d_assign_excluded_voxels/final_parc_k3.mnc")

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
    name_file = paste0("./visualization/5f_sex_interactions/adni_",var, ".png")
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



#endregion
