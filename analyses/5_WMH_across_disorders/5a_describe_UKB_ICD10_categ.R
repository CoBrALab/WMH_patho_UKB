
library(data.table)
library(effsize)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(scales)
library(forcats)
library(cowplot)
library(patchwork)
library(magrittr)
library(viridis)
library(RColorBrewer)
library(RMINC)
library(ggrepel)
library(cowplot)
library(ggpubr)
library(stringr)

# Describe and plot UKB disease categories based on ICD-10 codes (category 1712; https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=1712)

dir.create("./results/5a_describe_UKB_ICD10_categ", showWarnings=FALSE)
dir.create("./visualization/5a_describe_UKB_ICD10_categ", showWarnings=FALSE)

# Load data

# cleand first-occurences dataframe
firstocc = as.data.frame(fread("../../../UKB/Analyses/clean_firstocc/results/firstocc_categ.tsv"))

inclusions = as.data.frame(fread("../../QC/inclusions_final.txt"))
colnames(inclusions) = "ID"

demo = as.data.frame(fread("../../../UKB/tabular/df_demo/UKBB_demo_wider.tsv"))
demo = subset(demo, InstanceID == 2, select=c('SubjectID', 'Sex_31_0', 'Age_when_attended_assessment_centre_21003_0'))
colnames(demo) = c("ID", "Sex", "Age")
demo = merge(inclusions, demo, by="ID", all.x=TRUE)
demo$Sex = factor(demo$Sex, levels=c("Male", "Female"))

wmh = as.data.frame(fread("../3_temporal_sustain/results/3f_clean_wmh_data/wmh_combined.tsv"))
wmh$sex = as.factor(wmh$sex)

firstocc = merge(inclusions, firstocc, by="ID", all.x=TRUE)
firstocc = merge(firstocc, demo, by="ID", all.x=TRUE)
firstocc = firstocc %>% mutate(
        MRI_date = as.Date(MRI_date),
        icd_code = as.factor(icd_code),
        dx_name = as.factor(dx_name),
        DX_date = as.Date(DX_date),
        chapter = as.factor(chapter),
        categ = as.factor(categ),
        Sex = as.factor(Sex)
    )

# Make control group: no lifetime dx on endocrine, mental, nervous, circulatory systems
# Exceptions
controls_exc = firstocc %>% filter(chapter %in% c("Endocrine, nutritional and metabolic diseases", "Mental, Behavioral and Neurodevelopmental disorders", "Diseases of the nervous system", "Diseases of the circulatory system"))

controls_exc_dx = table(controls_exc$dx_name)
controls_exc_dx = controls_exc_dx[order(controls_exc_dx)]

# dx to keep in controls (high prevalence and low impact on brain):
# depressive episode (F32; n=4163), migraine (G43; n=2632), other anxiety disorders (F41; n=2632), haemorrhoids (I84; n=1859)
controls_exc = controls_exc %>% filter(!icd_code %in% c("F32", "G43", "F41", "I84"))

controls = inclusions$ID[which(!inclusions$ID %in% controls_exc$ID)]

# Manually-defined categories of ICD-10 codes
icd_codes_list <- readRDS("../../../UKB/Analyses/clean_firstocc/results/icd_codes_list.rds")

hc_df = firstocc %>% filter(ID %in% controls) %>% distinct(ID, .keep_all=TRUE)
hc_df$Group = "HC"

# Make plots for each category

dx_hc_colors = c(Controls="")

icd_dx_name = firstocc %>% distinct(icd_code, .keep_all=TRUE) %>%
    select(icd_code, dx_name, chapter, categ)

for (dx in 1:length(icd_codes_list)) {
    print(names(icd_codes_list)[dx])

    plots = list()
    i = 1

    # Case dataframe: ICD-10 code in disease category and participant NOT in controls
    dx_df = firstocc %>% filter(icd_code %in% icd_codes_list[[dx]] & !ID %in% controls)
    
    dx_df$Group="DX"
    ncol_dx_df = ncol(dx_df)

    # dx_df can contains duplicates because of commorbidities in the same category
    # Still want to visualize the breakdown of all dx, but get the total n from unique values
    dx_df_unique = dx_df %>% distinct(ID, .keep_all = TRUE)

    hc_dx_df = rbind(hc_df, dx_df)
    hc_dx_df = hc_dx_df %>% mutate(Group = factor(Group, levels=c("HC", "DX"), labels=c("Controls", "Diagnosis")))

    # Age distributions between cases and controls
    age_dist = ggplot(hc_dx_df, aes(x=Age, color=Group)) + 
        geom_density(size=1.5) +
        scale_color_viridis_d(labels = c("HC", "DX"), option="D", begin = 0.2, end=0.8, direction = -1, name="Group") +
        scale_y_continuous(name="Density") +
        theme_classic() +
        theme(text = element_text(size=15), plot.title = element_text(hjust=0.5, size=10))
    
    # Sex distributions between cases and controls
    sex_perc <- hc_dx_df %>%
        group_by(Group, Sex) %>%
        summarise(count = n()) %>%
        group_by(Group) %>%
        mutate(proportion = count / sum(count))

    sex_dist = ggplot(sex_perc, aes(x=Group, y=proportion, fill=factor(Sex, levels=c("Female", "Male")))) +
        geom_col() + 
        scale_x_discrete(labels=c("HC", "DX"), name=" ") +
        scale_y_continuous(name="Proportion") +
        scale_fill_discrete(name = "Sex") +
        theme_classic() +
        theme(text = element_text(size=15), plot.title = element_text(hjust=0.5))

    # Visualize breakdown of diagnoses included in category
    labels_dx_pie = dx_df %>% count(dx_name) %>%
        mutate(csum = rev(cumsum(rev(n))), 
                pos = n/2 + lead(csum, 1),
                pos = if_else(is.na(pos), n/2, pos)) %>%
        left_join(icd_dx_name)

    labels_legend = stringr::str_wrap(gsub("_", " ",
        paste0(labels_dx_pie$icd_code, " (n=", labels_dx_pie$n, ")\n", labels_dx_pie$dx_name, "\n")),
        width = 30  # Adjust width for better wrapping
        )

    dx_pie = ggplot(labels_dx_pie, aes(x="", y=n, fill=dx_name)) +
            geom_col(width = 1, color = 1) +
            scale_fill_manual(values = hue_pal()(length(labels_legend)), labels = labels_legend) +
            coord_polar(theta = "y") +
            geom_label_repel(aes(y = pos, label = n), size = 4.5, nudge_x = 1, show.legend = FALSE) +
            guides(fill = guide_legend(title = "")) +
            theme(axis.ticks = element_blank(), axis.title = element_blank(), axis.text = element_blank(), panel.background = element_rect(fill = "white"),
                    legend.position = "right", legend.direction = "vertical", legend.text=element_text(size=10),
                    aspect.ratio = 1)

    # Distribution of days between mri and dx
    time_mri_dx = ggplot(hc_dx_df %>% filter(Group == "Diagnosis"), aes(x=days_mri_dx/365.25)) +
        geom_histogram(color="black", fill="#0075ff") + 
        geom_vline(xintercept=0, size=1.5, color="black") +
        scale_y_continuous(name="Count") +
        scale_x_continuous(name = "Years (DX - MRI)") +
        theme_classic() +
        theme(text = element_text(size=15), plot.title = element_text(hjust=0.5))

    # Merge plots and write

    layout <- c(
            area(1, 1, 1, 1),  # First plot
            area(1, 2, 1, 2),  # Second plot
            area(1, 3, 1, 3),  # Third plot
            area(1, 4, 1, 4)   # Fourth plot (twice as wide)
        )

    wrap_plots(age_dist, sex_dist, time_mri_dx, dx_pie, design = layout)
    ggsave(paste0("./visualization/5a_describe_UKB_ICD10_categ/",names(icd_codes_list)[dx],".png"), width=17, height=3*length(plots))

}

