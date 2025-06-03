
library(data.table)
library(ggplot2)
library(tidyverse)
library(viridis)
library(patchwork)
library(grid)
library(Matrix)
library(RMINC)
library(gridExtra) 
library(MRIcrotome)
library(magrittr)
library(RColorBrewer)

# Visualize SuStaIn outputs
# Customized plotting with ggplot

dir.create("./results/3d_sustain_viz", showWarnings=FALSE)
dir.create("./visualization/3d_sustain_viz", showWarnings=FALSE)

# Prepare functions

# Function that takes in the samples_sequence sustain output and outputs a df with the zscore sequence and probabilities
get_zscore_sequence = function(df, maps, Z_vals) {
    df = df + 1 # Add 1 because stages start at 0

    N_stages = ncol(df)
        
    # Calculate max number of z-score events
    count_non_zero = function(row) { sum(row != 0) }
    non_zero_counts = apply(Z_vals, 1, count_non_zero)
    max_zscores = max(non_zero_counts)
    print(paste0("Max z-score = ", max_zscores))

    # Calculate z-score events per marker for each sample sequence
    micro_zevents = list()
    for (m in 1:length(maps)) {
        micro_zevents[[m]] = matrix(0, nrow=nrow(df), ncol=ncol(df))
    }

    # Order of Z-score events
    order = matrix(0, nrow(Z_vals), ncol(Z_vals))
    i=0
    for(c in 1:ncol(Z_vals)) {
        for(r in 1:nrow(Z_vals)) {
            if (Z_vals[r,c] != 0) {
                i=i+1
                order[r,c] = i
            } else (order[r,c] = 0)
        }
    }

    for (i in 1:nrow(df)) {
        # print(i)
        for (m in 1:length(maps)) {
            # pos = seq(from = m, by = length(maps), length.out=4)
            pos = order[m,][order[m,] != 0]
            vals = which(df[i,] %in% pos)
            
            for (v in 1:length(vals)) {
                micro_zevents[[m]][i,][vals[v]] = v
            }
        }
    }

    # Number of times each z-score event appear at each stage position, for each marker
    micro_zevents_prev = list()

    for (m in 1:length(maps)) {
        micro_zevents_prev[[m]] = matrix(0, nrow=max_zscores, ncol=ncol(df))
    }

    for (m in 1:length(maps)) {
        for (i in 1:nrow(micro_zevents_prev[[m]])) {
            occurrences = apply(micro_zevents[[m]], 2, function(column) sum(column == i))
            micro_zevents_prev[[m]][i,] = occurrences
        }
    }

    # Put in long format for ggplot and calculate uncertainty
    df_to_plot = as.data.frame(matrix(nrow=0, ncol=4))
    colnames(df_to_plot) = c("micro", "stage", "zscores", "perc")

    for (m in 1:length(maps)) {
        zevents = as.data.frame(which(micro_zevents_prev[[m]] != 0, arr.ind = TRUE))
        new_zevents = zevents

        # Remove zevents for one zevents if more than one zevent happen at the same stage for this marker (keep highest prob)
        stages = as.data.frame(table(zevents[, 2]))
        for (i in 1:nrow(zevents)) {
            # Check if stage appears more than once
            try(
            if (stages$Freq[which(stages[,1] == zevents[, 2][i])] > 1) {
                print(zevents[, 2][i])
                idx = which(zevents[, 2] == zevents[, 2][i]) # Indices in zevents of events at same stage
                first = micro_zevents_prev[[m]][zevents$row[idx[1]],zevents$col[idx[1]]]
                second = micro_zevents_prev[[m]][zevents$row[idx[2]],zevents$col[idx[2]]]

                if(first > second) {new_zevents[idx[2],] = NA} 
                if(first < second) {new_zevents[idx[1],] = NA}
                else if(first == second) {new_zevents[idx[2],] = NA}

            }
            )
        }
        new_zevents = na.omit(new_zevents)

        new_zevents_vec = rep(0, N_stages)
        new_zevents_vec[new_zevents[, 2]] = new_zevents[, 1]

        mat = as.data.frame(matrix(0, ncol=4, nrow=N_stages))
        colnames(mat) = c("micro", "stage", "zscores", "perc")
        mat$micro = rep(maps[m], N_stages)
        mat$stage = seq(1,N_stages)
        mat$zscores = new_zevents_vec

        perc_vec = rep(0, N_stages)

        for (i in 1:nrow(new_zevents)) {
            perc_vec[new_zevents[, 2]][i] = micro_zevents_prev[[m]][new_zevents$row[i],new_zevents$col[i]] / nrow(df)
        }
        mat$perc = perc_vec

        df_to_plot = rbind(df_to_plot, mat)
    }

    df_to_plot$micro = as.factor(df_to_plot$micro)
    df_to_plot$stage = as.factor(df_to_plot$stage)
    df_to_plot$zscores = as.factor(df_to_plot$zscores)
    df_to_plot$perc = as.numeric(df_to_plot$perc)
    df_to_plot[df_to_plot[, "zscores"] == 0, "zscores"] = NA

    return(df_to_plot)
}

# Plot positional variance diagram with probabilities
plot_pos_var = function(df, fig_name, title, maps, maps_level, N_stages) {

    z_scores_all = c(0.5,1,1.5,2,3,4,5)
    z_scores = z_scores_all[seq(1,max(as.numeric(levels(df$zscores))))]

    color_scale = viridis_pal(alpha = 1, begin = 0, end = 0.8, direction = -1, option = "magma")(length(z_scores_all))
    color_scale = color_scale[seq(1,length(z_scores))]

    plt = ggplot(df, aes(x=stage, y=factor(micro, level = maps_level), fill=zscores)) + 
        geom_tile(aes(alpha=perc)) + 
        geom_segment(data=data.frame(x=rep(0, length(maps)), y=seq(1, length(maps)), xend = rep(N_stages+0.5, length(maps)), yend = seq(1, length(maps))),
                        aes(x=x, y=y, xend=xend, yend=yend), size=0.3, inherit.aes=FALSE) + 
        scale_fill_manual(labels = z_scores, na.translate = FALSE, values=color_scale) +
        scale_alpha_continuous(range = c(0, 1)) + 
        theme_classic() + 
        labs(y=" ", x="SuStaIn Stage") + 
        guides(fill=guide_legend(title="Z-scores"), alpha=guide_legend(title="Prob")) + 
        ggtitle(title) + 
        theme(axis.text.x = element_text(size=13), axis.text.y = element_text(size=20),
            axis.title=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15),
            panel.border = element_rect(colour = "black", fill=NA, size=1), axis.line = element_blank(),
            legend.box = "horizontal", plot.title=element_text(hjust=0.5, size=20))
    ggsave(fig_name, width=12, height=length(maps)/2, dpi=300)

    print(fig_name)

    return(plt)
}

# Jack curves
plot_jack_curves = function(df, fig_name, title, maps) {
    all_seq = as.data.frame(matrix(0, nrow=0, ncol=ncol(df)))
    colnames(all_seq) = colnames(df)

    for (n in 1:length(maps)) {

        # Sequence for one biomarker
        micro_seq = subset(df, micro == maps[n])
        # Add stage 0
        micro_seq = rbind(data.frame("micro" = maps[n], "stage" = 0, "zscores" = NA, "perc" = 0), micro_seq)

        # Winner-take-all
        uncertain_events = subset(micro_seq, perc<1 & is.na(zscores)==FALSE)
        uncertain_z = unique(uncertain_events$zscores)

        for (i in 1:length(uncertain_z)) {
            uncertain_events_z = subset(uncertain_events, zscores==uncertain_z[i])
            row_winner = rownames(uncertain_events_z)[which.max(uncertain_events_z$perc)]
            row_losers = rownames(uncertain_events_z)[-which.max(uncertain_events_z$perc)]

            for(r in row_losers) {
                micro_seq[r,"zscores"] = NA
            }
        }

        # Cumulative z-score
        cum_zscore = 0
        for (i in 1:nrow(micro_seq)) {
            cum_zscore = ifelse(is.na(micro_seq$zscores[i]) == FALSE, cum_zscore+1, cum_zscore)
            micro_seq$zscores[i] = cum_zscore
        }

        all_seq = rbind(all_seq, micro_seq)
    }

    # Exact z-score thresholds (not just order)
    all_seq$zscores =   ifelse(all_seq$zscores == 1, 0.5,
                        ifelse(all_seq$zscores == 2, 1,
                        ifelse(all_seq$zscores == 3, 1.5,
                        ifelse(all_seq$zscores == 4, 2,
                        ifelse(all_seq$zscores == 5, 3,
                        ifelse(all_seq$zscores == 6, 4,
                        ifelse(all_seq$zscores == 7, 5, 0)))))))

    row.names(all_seq) = seq(1, nrow(all_seq))

    all_seq$stage = factor(all_seq$stage, levels = seq(0, max(as.numeric(all_seq$stage))))
    all_seq$micro = factor(all_seq$micro)

    color_scale = c(MD = "#04319E", ISOVF = "#8298CF", FA = "#2B520B", ICVF="#448312", OD="#A2C189", T2star="#D36108", QSM="#E9B084")

    # Plot
    plot = ggplot(all_seq, aes(x=stage, y=zscores, group = micro, color=micro)) + 
        geom_step(size=2) +
        scale_colour_manual(name="", values=color_scale) +
        scale_x_discrete(name = "Stage", breaks=unique(as.numeric(all_seq$stage)[as.numeric(all_seq$stage) %% 5 == 0])) +
        scale_y_continuous(name="Z-scores", breaks=c(0.5, 1, 1.5, 2, 3, 4, 5), limits=c(0,5)) +
        ggtitle(title) + 
        guides(color = guide_legend(override.aes = list(size = 5))) +
        theme_classic() + 
        theme(text = element_text(size=30))
    ggsave(paste0(fig_name, ".png"), width=10, height=8, dpi=300)
    
    return(plot)
}

# Run

clusts = seq(1,3)

N_S_max = 3

out_dir = './results/3c_sustain_bind_CV'

names = c("FA", "MD", "ICVF", "ISOVF", "OD", "T2star", "QSM")

for (clust in clusts) {

    cat(paste0("\n--------------------\nCluster = ",clust,"\n--------------------\n"))

    viz_dir = paste0('./visualization/3d_sustain_viz/c',clust)
    dir.create(viz_dir, showWarnings=FALSE)

    # Z_vals inputs
    
    if (clust==1) {

        Z_vals = matrix(0, nrow=length(names), ncol=7)
        Z_vals[1,] = c(0,0,0,0,0,0.5,1)
        Z_vals[2,] = c(0,0,0,0,0.5,1,1.5)
        Z_vals[3,] = c(0,0,0,0,0.5,1,1.5)
        Z_vals[4,] = c(0,0,0,0,0,0.5,1)
        Z_vals[5,] = c(0,0,0,0,0,0,0.5)
        Z_vals[6,] = c(0,0,0,0,0,0.5,1)
        Z_vals[7,] = c(0,0,0,0,0,0,0.5)

        Z_max = c(2,3,3,2,1,2,1)
    }

    if (clust==2) {
        Z_vals = matrix(0, nrow=length(names), ncol=7)
        Z_vals[1,] = c(0,0,0,0.5,1,1.5,2)
        Z_vals[2,] = c(0,0.5,1,1.5,2,3,4)
        Z_vals[3,] = c(0,0,0.5,1,1.5,2,3)
        Z_vals[4,] = c(0,0,0,0.5,1,1.5,2)
        Z_vals[5,] = c(0,0,0,0,0.5,1,1.5)
        Z_vals[6,] = c(0,0,0,0.5,1,1.5,2)
        Z_vals[7,] = c(0,0,0,0,0,0.5,1)

        Z_max  = c(3,5,4,3,3,3,2)
    }

    if (clust==3) {
        Z_vals = matrix(0, nrow=length(names), ncol=7)
        Z_vals[1,] = c(0,0,0,0.5,1,1.5,2)
        Z_vals[2,] = c(0.5,1,1.5,2,3,4,5)
        Z_vals[3,] = c(0.5,1,1.5,2,3,4,5)
        Z_vals[4,] = c(0,0,0,0.5,1,1.5,2)
        Z_vals[5,] = c(0,0,0,0,0.5,1,1.5)
        Z_vals[6,] = c(0,0,0.5,1,1.5,2,2)
        Z_vals[7,] = c(0,0,0,0,0,0.5,1)

        Z_max = c(3,7,7,3,3,3,2)
    }


    # Load sustain outputs and calculate zscore sequences
    names = c("FA", "MD", "ICVF", "ISOVF", "OD", "T2star", "QSM")

    ss = list()
    df_ss = list()

    for (s in 0:(N_S_max-1)) {
        ss[[s+1]] = list()
        df_ss[[s+1]] = list()
        for (k in 0:s) {
            print(paste0("Subtypes = ",s,"; k = ",k))
            ss[[s+1]][[k+1]] = t(fread(paste0(out_dir, "/c",clust,"/subtype",s,"/samples_sequence_cval_",k,".csv")))
            df_ss[[s+1]][[k+1]] = get_zscore_sequence(ss[[s+1]][[k+1]], names, Z_vals)
        }
    }

    # Plot positional variance diagram with probabilities
    plots = list()
    for (s in 0:(N_S_max-1)) {
        plots[[s+1]] = list()
        for (k in 0:s) {
            print(paste0("Subtypes = ",s,"; k = ",k))
            plots[[s+1]][[k+1]] = plot_pos_var(df_ss[[s+1]][[k+1]], paste0(viz_dir, "/new_pos_var_CV_s",s,"_k",k,"_cval.png"), " ", names, ncol(ss[[s+1]][[k+1]]))
        }
        wrap_plots(plots[[s+1]], ncol=1)
        ggsave(paste0(viz_dir, "/new_pos_var_CV_s",s,"_cval_all.png"), height = 3*(s+1), width = 12)
    }


    # Jack curves
    plots = list()
    for (s in 0:(N_S_max-1)) {
        plots[[s+1]] = list()
        for (k in 0:s) {
            print(paste0("Subtypes = ",s,"; k = ",k))
            plots[[s+1]][[k+1]] = plot_jack_curves(df_ss[[s+1]][[k+1]], paste0(viz_dir, "/jack_curves_s",s,"_k",k), " ", names)
        }
        wrap_plots(plots[[s+1]], ncol=1)
        ggsave(paste0(viz_dir, "/jack_curves_s",s,"_all.png"), height = 8*(s+1), width = 10) # Not cross-validated
    }

}



