
library(data.table)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(scales)
library(MRIcrotome)
library(magrittr)
library(viridis)
library(RColorBrewer)
library(RMINC)
library(ggsignif)
library(patchwork)
library(effects)

# Sex differences in WMH pathophysiology

dir.create("./results/4a_sex_differences", showWarnings=FALSE)
dir.create("./visualization/4a_sex_differences", showWarnings=FALSE)

df = as.data.frame(fread("../3_temporal_sustain/results/3f_clean_wmh_data/wmh_combined.tsv"))
df = df %>% 
    mutate(sex = factor(sex, levels=c("Male", "Female")))

y_vars = colnames(df)[-c(1,2,3)]

# Model 1: y ~ sex + age
# Model 2: y ~ sex + age + WMHvol (of corresponding cluster)

#to get cluster info later
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

#loop through models

models = list()

sex_data <- data.frame(
  Model = numeric(),
  wmh_var = character(),
  cluster = character(),
  pval = numeric(),
  beta_coef = numeric(),
  stringsAsFactors = FALSE
)

for (model_num in 1:2){
    print(paste0("Model: ",model_num))

    models[[model_num]] = list()

    for (y_var in 1:length(y_vars)){

        #extract cluster info so that relevant volme information can be extracted
        cluster = substrRight(y_vars[y_var], 2)
        if (model_num == 3 | model_num == 5){
            vol = grep(cluster,c("WMHvol_c1","WMHvol_c2","WMHvol_c3"),value=TRUE)
        }

        # Run linear model
        if (model_num == 1) {
            formula = as.formula(paste0("scale(",y_vars[y_var],") ~ sex + age"))
        } else if (model_num == 2) {
            formula = as.formula(paste0("scale(",y_vars[y_var],") ~ sex + age + ",vol))
        }

        model <- lm(formula, data = df)

        summary_lm <- summary(model)

        models[[model_num]][[y_vars[y_var]]] = model

        # For models 1; save main sex effect
        if (model_num == 1) {
            print("model 1")
            p_value_interaction <- coef(summary_lm)[2,'Pr(>|t|)']
            beta_coefficient <- coef(summary_lm)[2,'Estimate']
        }

        # For model 2; save main sex effect if predictor not WMH volume
        if (model_num == 3 & grepl("WMHvol", y_vars[y_var]) > 0) {
            print("model 3; Volume")
            p_value_interaction <- NA
            beta_coefficient <- NA
        } else if (model_num == 3 & grepl("WMHvol", y_vars[y_var]) == 0) {
            print("model 3; Not volume")
            p_value_interaction <- coef(summary_lm)[2,'Pr(>|t|)']
            beta_coefficient <- coef(summary_lm)[2,'Estimate']
        }

        #remove _cluster from the variable name
        variable_cleaned <- gsub("_c\\d+", "", paste0(y_vars[y_var]))

        #get cluster of relevant variable as an int
        new_row <- data.frame(Model = model_num, wmh_var= variable_cleaned, cluster = cluster,pval = p_value_interaction, beta_coef = beta_coefficient ,stringsAsFactors = FALSE)

        # Add the new row to the existing data frame
        sex_data <- rbind(sex_data, new_row)
    }
}

# Correct p-values
sex_data$pval_fdr = p.adjust(sex_data$pval, method="fdr")

fwrite(sex_data, "./results/sex_diff_effects.tsv", col.names = TRUE, quote = FALSE, row.names=FALSE, sep="\t")

# Visualize coefficients on brain

sex_data = as.data.frame(fread("./results/sex_diff_effects.tsv"))
sex_data[,c(1,2,3)] = lapply(sex_data[,c(1,2,3)], as.factor)

metrics = c("pval_fdr", "beta_coef")

names=c("WMHvol", "stage", "MD", "ISOVF", "FA", "ICVF", "OD", "T2star", "QSM")

anatVol = mincArray(mincGetVolume("../../data/UKB_template_2mm.mnc"))
viridis_scale = viridis(n=255, option="B")
color_scale_div_2 = colorRampPalette(brewer.pal(9,"Blues"))(255)
color_scale_div_1 = colorRampPalette(brewer.pal(9,"Reds"))(255)

for (i in 1:length(metrics)) {
    print(metrics[i])

    for (model_num in 1:5) {
    
        # Results to mnc
        test_df = sex_data %>% filter(Model == model_num)
        if (nrow(test_df) > 0) {

            mnc_tmp = list()
            for (n in 1:length(names)) {
                print(names[n])

                to_mnc = test_df %>%
                    filter(wmh_var == names[n]) %>%
                    select(!!sym(metrics[i]))
                to_mnc = as.numeric(to_mnc[,1])

                mnc_tmp[[n]] = mincGetVolume("../../2_spatial_clust/results/2d_assign_excluded_voxels/final_parc_k3.mnc")
                mnc_tmp[[n]][][which(mnc_tmp[[n]][]==1)] = to_mnc[1]
                mnc_tmp[[n]][][which(mnc_tmp[[n]][]==2)] = to_mnc[2]
                mnc_tmp[[n]][][which(mnc_tmp[[n]][]==3)] = to_mnc[3]
                mnc_tmp[[n]] = mincArray(mnc_tmp[[n]])
            }

            high_thresh = ifelse(metrics[i] == "pval_fdr", 0.05, ifelse(metrics[i] == "beta_coef", 0.2))

            # Visualize results on brain map
            png(file=paste0("./visualization/4a_sex_differences/",metrics[i], "_model",model_num, ".png"), width=7500, height=1000, pointsize = 150)
            sliceSeries(nrow=1, ncol=1, begin=52, end=52, dimension=3) %>%
                addtitle(names[1]) %>%
                anatomy(anatVol, low=10, high=200) %>%
                overlay(mnc_tmp[[1]],
                    low=0, high=high_thresh, col=color_scale_div_1, rCol=color_scale_div_2, symmetric=TRUE) %>%
            sliceSeries(nrow=1, ncol=1, begin=52, end=52, dimension=3) %>%
                addtitle(names[2]) %>%
                anatomy(anatVol, low=10, high=200) %>%
                overlay(mnc_tmp[[2]],
                        low=0, high=high_thresh, col=color_scale_div_1, rCol=color_scale_div_2, symmetric=TRUE) %>%
            sliceSeries(nrow=1, ncol=1, begin=52, end=52, dimension=3) %>%
                addtitle(names[3]) %>%
                anatomy(anatVol, low=10, high=200) %>%
                overlay(mnc_tmp[[3]],
                        low=0, high=high_thresh, col=color_scale_div_1, rCol=color_scale_div_2, symmetric=TRUE) %>%
            sliceSeries(nrow=1, ncol=1, begin=52, end=52, dimension=3) %>%
                addtitle(names[4]) %>%
                anatomy(anatVol, low=10, high=200) %>%
                overlay(mnc_tmp[[4]],
                        low=0, high=high_thresh, col=color_scale_div_1, rCol=color_scale_div_2, symmetric=TRUE) %>%
            sliceSeries(nrow=1, ncol=1, begin=52, end=52, dimension=3) %>%
                addtitle(names[5]) %>%
                anatomy(anatVol, low=10, high=200) %>%
                overlay(mnc_tmp[[5]],
                        low=0, high=high_thresh, col=color_scale_div_1, rCol=color_scale_div_2, symmetric=TRUE) %>%
            sliceSeries(nrow=1, ncol=1, begin=52, end=52, dimension=3) %>%
                addtitle(names[6]) %>%
                anatomy(anatVol, low=10, high=200) %>%
                overlay(mnc_tmp[[6]],
                        low=0, high=high_thresh, col=color_scale_div_1, rCol=color_scale_div_2, symmetric=TRUE) %>%
            sliceSeries(nrow=1, ncol=1, begin=52, end=52, dimension=3) %>%
                addtitle(names[7]) %>%
                anatomy(anatVol, low=10, high=200) %>%
                overlay(mnc_tmp[[7]],
                        low=0, high=high_thresh, col=color_scale_div_1, rCol=color_scale_div_2, symmetric=TRUE) %>%
            sliceSeries(nrow=1, ncol=1, begin=52, end=52, dimension=3) %>%
                addtitle(names[8]) %>%
                anatomy(anatVol, low=10, high=200) %>%
                overlay(mnc_tmp[[8]],
                        low=0, high=high_thresh, col=color_scale_div_1, rCol=color_scale_div_2, symmetric=TRUE) %>%
            sliceSeries(nrow=1, ncol=1, begin=52, end=52, dimension=3) %>%
                addtitle(names[9]) %>%
                anatomy(anatVol, low=10, high=200) %>%
                overlay(mnc_tmp[[9]],
                        low=0, high=high_thresh, col=color_scale_div_1, rCol=color_scale_div_2, symmetric=TRUE) %>%
            legend("Effect size") %>%
            draw()
            dev.off()
        }
    }
}



