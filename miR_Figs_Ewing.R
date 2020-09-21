###This file produces all the figures and tables for the paper.

#######top######
library(reshape2)
library(readr)
library(dplyr)
library(purrr)
library(gtools)
library(tibble)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(forcats)
library(xtable)
library(gridExtra)
library(ggrepel)
library(cowplot)
source('./miR_Analysis_Functions.R')

#load mRNA-protein correlations from the cell paper. 
corr_df <- read_csv("./data_files/Protein_RNA_correlation.csv")

#Load the cleaned CCLE file
load('./data_files/CCLE_df.Rda')

#Import output from main Network Potential Pipeline
load('./data_files/final_df.Rda')

#Load processed miR_targeting data file
load("./data_files/targeting_mat.Rda")

#Load evaluated miR cocktails
load("./data_files/MiRCocktail_output.Rda")
cocktail_df <- output[[2]]
cocktail_df$numtargets <- 10
targets_df <- output[[1]]
rm(output)

load("./data_files/MiRCocktail5_output.Rda")
cocktail5_df <- output[[2]]
cocktail5_df$numtargets <- 5
targets5_df <- output[[1]]
rm(output)

load("./data_files/MiRCocktail15_output.Rda")
cocktail15_df <- output[[2]]
cocktail15_df$numtargets <- 15
targets15_df <- output[[1]]
rm(output)

cocktail_df <- bind_rows(cocktail_df, cocktail5_df,cocktail15_df)

#Load miR candidates matrix
load("./data_files/MiRcandidates_matrix.Rda")

###################Expression and Network Potential Overview Fig################

#This code creates a dataframe with total gibbs calculated for each cell line -experiment pairing
#Commented out code aggregates even farther to just have average gibbs for each ES cell line and 
#Then even farther to just have the global average gibbs/standard deviation
final_df_fig3 <- final_df %>%
  group_by(cell_line, experiment_num) %>% 
  summarise(gibbs_sum = sum(gibbs), 
            expression_sum = sum(expression)) %>% ungroup() %>%
  mutate(cell_line = fct_reorder(cell_line, gibbs_sum, .desc = TRUE)) 

final_df_fig3$experiment_num <- c("one","two","three","one","two","three",
                                  "one","two","three","one","two","three",
                                  "one","two","three","one","two","three")
colnames(final_df_fig3)[colnames(final_df_fig3) == "experiment_num"] <- "replicate"

#Also need to prep data to visualize description of G_i
#Goin to try this filtering out the observations with zero expression
#zero network potential - that should be one eto one. 
hist_df_allgenes <- final_df %>% 
  filter(expression > 0, gibbs < 0) %>% 
  group_by(gene_name) %>% summarise(gibbs = mean(gibbs), 
                                    expression = mean(expression))


#Make a histogram of the distribution of G_i and also mRNA expression
png(filename = "./figures/histogram_allgenes_boxplot.png",
    width = 900, height = 700)
g1 <- ggplot(data = hist_df_allgenes, aes(expression)) + theme_classic() +
  geom_histogram(bins = 70, color = "red3", fill = "darkred") +
  xlab("mRNA expression (normalized count)") + ylab("") + 
  labs(tag = "A") + 
  theme(text = element_text(size = 18),
        plot.tag = element_text(face = "bold"))

g2 <- ggplot(data = hist_df_allgenes, aes(gibbs)) + theme_classic() +
  geom_histogram(bins = 70, color = "blue", fill = "steelblue") + 
  xlab("Network potential") + ylab("") + labs(tag = "B") + 
  theme(text = element_text(size = 18),
        plot.tag = element_text(face = "bold"))  

#code for boxplot of total expression for each cell line
g3 <- ggplot(data = final_df_fig3, aes(cell_line, expression_sum))  + 
  geom_boxplot() + theme_classic() +
  theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust = 0),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 18),
        plot.tag = element_text(face = "bold")) +
  scale_y_continuous(limits = c(106800, 108600), 
                     breaks = seq(106800, 108600, by = 400)) +
  ylab("Total mRNA expression") + xlab("") +
  labs(tag = "C")
# Code for boxplot of total network potential for each cell line
g4 <- ggplot(data = final_df_fig3, aes(cell_line, gibbs_sum))  + 
  geom_boxplot() + theme_classic() +
  theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust = 0),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 18),
        plot.tag = element_text(face = "bold")) +
  scale_y_continuous(limits = c(-343000, -337000), 
                     breaks = seq(-343000 -337000, by = 1000)) +
  ylab("Network potential") + xlab("") +
  labs(tag = "D")

grid.arrange(g1, g2, g3, g4, nrow = 2)
dev.off()




######################mRNA overview Figure###########################

#Calculate average deltaGibbs and average of several other data objects for each 
#cell_line-Gene_name pair #also can specify what we want the deltaGibbs cutoff 
#to be and whether we want to limit to just genes implicated causally in cancer
final_df_mean <- final_df %>% 
  mean_cellline_calc(topx = TRUE, numgenes = 50, justcancer = FALSE, 
                     keep_cell_line = FALSE) 
final_df_mean_corr <- final_df %>% 
  mean_cellline_calc(topx = FALSE, numgenes = 10, justcancer = FALSE, 
                     keep_cell_line = FALSE) 

#Okay lets prep data to be graphed
corr_df1 <- filter(corr_df, `Gene Symbol` %in% final_df_mean$gene_name)
colnames(corr_df1) <- c("gene_name", "pearson_corr", "spearman_corr")
final_df_mean <- final_df_mean %>% left_join(corr_df1)


#create a vector that tells GGPlot to color code the labels genes associated with cancer
final_df_mean$text_format <- "black"
final_df_mean$text_format[final_df_mean$housekeeping == "yes"] <- "blue"
final_df_mean$text_format[final_df_mean$cancer_associated == "yes"] <- "red"
color_vec <- final_df_mean$text_format
final_df_mean <- final_df_mean %>% 
  filter(!is.na(se_deltaGibbs)) %>% 
  mutate(gene_name = fct_reorder(gene_name, mean_deltaGibbs, .desc = TRUE))

#Get a corr_df for the plot that is overlaid on the main figure (top 1000 genes by network potential). 
#CCLE_df2 <- filter(CCLE_df, gene_name %in% final_df_mean_corr$gene_name)

corr_df_hist <- filter(corr_df, `Gene Symbol` %in% final_df_mean_corr$gene_name)
colnames(corr_df_hist) <- c("gene_name", "pearson_corr", "spearman_corr")
final_df_mean_hist <- left_join(final_df_mean_corr, corr_df_hist)
#Prep cancer only data to be graphed
final_df_mean_cancer <- final_df %>% 
  mean_cellline_calc(topx = TRUE, numgenes = 50, justcancer = TRUE, 
                     keep_cell_line = FALSE)

corr_df2 <- filter(corr_df, `Gene Symbol` %in% final_df_mean_cancer$gene_name)
colnames(corr_df2) <- c("gene_name", "pearson_corr", "spearman_corr")
final_df_mean_cancer <- final_df_mean_cancer %>% left_join(corr_df2)

final_df_mean_cancer <- final_df_mean_cancer %>% 
  filter(!is.na(se_deltaGibbs)) %>% 
  mutate(gene_name = fct_reorder(gene_name, mean_deltaGibbs, .desc = TRUE))

#Make the plots
g1 <- ggplot(data = final_df_mean, 
             aes(gene_name, mean_deltaGibbs)) + geom_point(size = 0.9) + 
  geom_errorbar(aes(x = gene_name, 
                    ymin = mean_deltaGibbs - 2*se_deltaGibbs,
                    ymax = mean_deltaGibbs + 2*se_deltaGibbs)) +
  geom_rect(aes(xmin = as.numeric(gene_name) - 0.5, 
                xmax = as.numeric(gene_name) + 0.5,
                ymin = -100, ymax = 0, fill = pearson_corr),
            show.legend = TRUE) + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = -60, vjust = 0.5, hjust = 0, size = 12),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(size = 20),
        plot.tag = element_text(size = 20, face = "bold"),
        plot.margin = unit(c(0.2,1,0.2,0.2), "cm")) +
  ylab(expression(paste(Delta, "Network Potential"))) + 
  theme(axis.text.x = element_text(color = color_vec)) + 
  scale_color_manual(values = colors) + 
  scale_y_continuous(expand = c(0,0)) + 
  scale_fill_viridis_c() + 
  labs(tag = "A", fill = "Pearson Coefficient")

g2 <- ggplot(final_df_mean_hist, aes(x = pearson_corr)) + 
  geom_histogram(fill = 'steelblue2', color = 'navy',bins = 35) + 
  geom_label_repel(aes(label = ifelse(mean_deltaGibbs > 450, gene_name, ""),
                       y = 100),
                   box.padding = 0.4) +
  labs(x = "Pearson Correlation Coefficient", tag = "B") +
  theme_classic() + 
  theme(text = element_text(size = 18),
        plot.tag = element_text(size = 20, face = "bold"))

# g2 <- ggscatter(CCLE_df2, x = 'log_expression', y = 'protein_expression',
#                 size = 1.75, color = "steelblue", add = "reg.line",) +
#   stat_cor(label.y = 9, 
#            aes(label = ..rr.label..), size = 5.5, 
#            method = "spearman") +
#   facet_wrap(vars(gene_name)) + 
#   xlab("mRNA Expression") + 
#   ylab("Protein Expression") + 
#   theme(text = element_text(size = 18),
#         plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + 
#   labs(tag = "B") 
#   
  

#Draw the plots.
png(filename = "./figures/mRNA_Targetoverviewplot.png",
    width = 900, height = 600)
g1 + annotation_custom(
  ggplotGrob(g2), xmin = 20, xmax = 50, ymin = 500, ymax = 2040
)
dev.off()

##Make an mRNA target overview limited to genes associated in cancer
png(filename = "./figures/mRNA_Targetoverviewplot_cancer.png",
    width = 900, height = 600)
g <- ggplot(data = final_df_mean_cancer, 
            aes(gene_name, mean_deltaGibbs)) 
g + geom_point(size = 0.9) + 
  geom_errorbar(aes(x = gene_name, 
                    ymin = mean_deltaGibbs - 2*se_deltaGibbs,
                    ymax = mean_deltaGibbs + 2*se_deltaGibbs)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = -60, vjust = 0.5, hjust = 0),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 18),
        plot.margin = unit(c(0.2,1,0.2,0.2), "cm")) +
  geom_rect(aes(xmin = as.numeric(gene_name) - 0.5, 
                xmax = as.numeric(gene_name) + 0.5,
                ymin = -50, ymax = 0, fill = pearson_corr),
            show.legend = TRUE) + 
  ylab(expression(paste(Delta, "Network Potential"))) + 
  xlab("Gene Name")  + 
  scale_y_continuous(expand = c(0,0)) + 
  scale_fill_viridis_c()
dev.off()

##########################miR Cocktail Fig######################################

#Just going to present the example top cocktail from one cell line + the worst cocktail from the same cell line
#Prepping the dataframe for plotting will use some of the code from the make cocktail function
cell_line_plot <- "A673"
main_targetnum <- 10

targets <- unique(unlist(targets15_df))
#cell_line_targets <- unlist(targets_df[,cell_line_plot])
#pull out the non-target vector
non_targets_df <- filter(final_df, housekeeping == "yes")

#For plotting purposes need to pair down to 30 housekeeping genes
#Bit of a challenge because there are actually too many 
non_targets <- unique(non_targets_df$gene_name)
targeting <- filter(overall_targeting_mat.df, gene_name %in% c(non_targets, targets))
targeting$targets_value <- 0

#Set up the targeting matrix for the loss function
#This section allows the user to fiddle with both the priority of hitting targets vs. non-targets and the assumed 
#repression achieved from a given miR. Defaults are in the function call.
targeting$targets_value[
  targeting$housekeeping == "yes" & targeting$targets_logical ==1] <- 1
targeting$targets_value[
  targeting$housekeeping == "no" & targeting$targets_logical ==1] <- -1
#Create a new value for weighted delta gibbs
targeting <- targeting %>% 
  mutate(weighted_target = targets_value*mean_deltaGibbs)

#fix the housekeeping variable here..
targeting$housekeeping[targeting$housekeeping == "yes"] <- 1
targeting$housekeeping[is.na(targeting$housekeeping)] <- 0
targeting$housekeeping <- as.numeric(targeting$housekeeping)

#Actually I just want to rank every cocktail by cell line/condition
#split-apply-combine
cell_line_vec <- unique(cocktail_df$cell_line)
condition_vec <- unique(cocktail_df$numtargets)
cocktail_df2 <- cocktail_df %>% filter(numtargets ==0)
for(i in 1:length(cell_line_vec)) {
  for(j in 1:length(condition_vec)) {
    cocktail_ij <- cocktail_df %>% 
      filter(cell_line == cell_line_vec[i], numtargets == condition_vec[j]) %>%
      arrange(-desc(loss))
    cocktail_ij$rank <- 1:nrow(cocktail_ij)
    cocktail_df2 <- bind_rows(cocktail_df2, cocktail_ij)
  }
}


top_cocktail <- cocktail_df2 %>% 
  filter(rank == 1, cell_line == cell_line_plot,
         numtargets == main_targetnum) %>% 
  select(1:3) %>% unlist(.)

#need to convert housekeeping to logical to rank gene name by that condition
top_cocktail_df <- targeting %>% filter(miR %in% top_cocktail, 
                                        cell_line == cell_line_plot) %>%
  mutate(gene_name = fct_reorder(gene_name, housekeeping, .desc = FALSE))

#removing extra white space in the fig by getting rid of all non-targeted genes.
nottargeted <- top_cocktail_df %>% group_by(gene_name) %>%
  summarise(sum_target = sum(targets_logical)) %>% 
  filter(sum_target == 0)
top_cocktail_df <- filter(top_cocktail_df, !(gene_name %in% nottargeted$gene_name))

top_cocktail_dfhist <- top_cocktail_df %>% group_by(gene_name) %>% 
  filter(housekeeping == TRUE) %>% 
  summarise(num_hits = sum(targets_logical))
top_cocktail_dfhist$status <- "best cocktail"

#Prep the bottom cocktail dfs for plotting
#Need max_rank to be in a different statement because different numbers of cocktails were 
#evaluated for the different conditions. 
bottom_cocktail <- cocktail_df2 %>% 
  filter(cell_line == cell_line_plot,
         numtargets == main_targetnum) %>% 
  filter(rank == max(rank)) %>%
  select(1:3) %>% unlist(.)

bottom_cocktail_df <- targeting %>% 
  filter(miR %in% bottom_cocktail, cell_line == cell_line_plot) %>% 
  mutate(gene_name = fct_reorder(gene_name, housekeeping, .desc = FALSE))

#same story as above. 
nottargeted <- bottom_cocktail_df %>% group_by(gene_name) %>%
  summarise(sum_target = sum(targets_logical)) %>% 
  filter(sum_target == 0)
bottom_cocktail_df <- filter(bottom_cocktail_df, !(gene_name %in% nottargeted$gene_name))

bottom_cocktail_dfhist <- bottom_cocktail_df %>% group_by(gene_name) %>% 
  filter(housekeeping == 1) %>% 
  summarise(num_hits = sum(targets_logical))
bottom_cocktail_dfhist$status <- "worst cocktail"
cocktail_dfhist <- bind_rows(top_cocktail_dfhist, bottom_cocktail_dfhist)

#Lets plot the number of times a given miR shows up in the bottom 10 or top 10 cocktails (averaged across cell lines)
mir_rank_df <- cocktail_df2 %>% 
  filter(numtargets == 10) %>% 
  group_by(V1,V2,V3) %>% 
  summarise(mean_rank = mean(rank)) %>% 
  pivot_longer(cols = c(V1,V2,V3), values_to = "miR", 
               names_to = "cocktail_position") %>% 
  arrange(-desc(mean_rank))

mir_rank_dftop <- mir_rank_df %>% slice_head(n = 30)
mir_rank_dftop$cocktail_cat <- "top 10"

mir_rank_dfbottom <- mir_rank_df %>% slice_tail(n = 30)
mir_rank_dfbottom$cocktail_cat <- "bottom 10"

mir_rank_df <- bind_rows(mir_rank_dftop, mir_rank_dfbottom)

mir_rank_df <- mir_rank_df %>% 
  group_by(miR, cocktail_cat) %>% 
  summarise(n = n()) %>% ungroup() %>%
  mutate(miR = fct_reorder2(.f = miR, .x = cocktail_cat, .y = n, .desc = TRUE))
#Start of plotting code for main figure
g1 <- ggplot(data = top_cocktail_df, 
             aes(x = gene_name, y = miR, fill = targets_value)) + 
  geom_tile() + 
  scale_fill_gradient2(low = "red", mid = "white", high = "blue") + 
  theme(axis.text.x=element_text(angle = -60, vjust = 0.5, hjust = 0, size = 10), 
        axis.title.x = element_blank(),
        legend.position = "none", 
        text = element_text(size = 20),
        plot.margin = unit(c(0.2,1.5,0.2,0.2), "cm")) + 
  ylab("miRNA") 
g2 <- ggplot(data = top_cocktail_dfhist, aes(num_hits)) + 
  geom_histogram(bins =4, fill = "steelblue") + 
  scale_x_continuous(breaks = c(0,1,2,3), limits = c(1,3)) + 
  theme_classic() + 
  theme(text = element_text(size = 16)) + 
  xlab("Number of hits")
g3 <- ggplot(data = bottom_cocktail_df, 
             aes(x = gene_name, y = miR, fill = targets_value)) + 
  geom_tile() + 
  scale_fill_gradient2(low = "red", mid = "white", high = "blue") + 
  theme(axis.text.x=element_text(angle = -60, vjust = 0.5, hjust = 0, size = 7.5),
        axis.title.x = element_blank(),
        legend.position = "none", 
        text = element_text(size = 20),
        plot.margin = unit(c(0.2,1.5,0.2,0.2), "cm")) + 
  ylab("miRNA") 
g4 <- ggplot(data = bottom_cocktail_dfhist, aes(num_hits)) + 
  geom_histogram(bins =4, fill = "steelblue") + 
  scale_x_continuous(breaks = c(0,1,2,3)) + 
  theme_classic() + 
  theme(text = element_text(size = 16)) + 
  xlab("Number of hits") 
g5 <- ggplot(data = mir_rank_df, aes(miR, n, fill =cocktail_cat)) + 
  geom_col() + 
  theme_classic() + 
  theme(axis.text.x=element_text(angle = -60, vjust = 0.5, hjust = 0, size = 10),
        text = element_text(size = 16),
        plot.margin = unit(c(0.2,1.5,0.2,0.2), "cm"),
        legend.title = element_blank()) + 
  labs(y = "Number of appearances", x = "miRNA")
##Start plot - main text
png(filename = "./figures/miRCocktailPlot.png", height = 800, width = 1000) 
#g3 is s histogram showing the distribution of hits.
plot_grid(g1,g2,g3, g4, g5, nrow = 3, ncol = 2, rel_widths = c(1, 0.5),
          align = "hv", axis = "btl", labels = "AUTO", label_size = 20, 
          label_x = c(0,0.1,0,0.1,0))
dev.off()

#Lets make a new figure with the best fig by condition in a different panel. -this is also the same
#for the 15 target condition, lets plot every cell line
# g_list <- list()
# for(i in 1:length(cell_line_vec)) {
#   top_cocktail <- cocktail_df2 %>% 
#     filter(rank == 1, cell_line == cell_line_vec[i],
#            numtargets == 15) %>% 
#     select(1:3) %>% unlist(.)
#   
#   #need to convert housekeeping to logical to rank gene name by that condition
#   top_cocktail_df <- targeting %>% 
#     filter(miR %in% top_cocktail, cell_line == cell_line_vec[i]) %>%
#     mutate(gene_name = fct_reorder(gene_name, housekeeping, .desc = FALSE))
#   
#   #Create the histogram dataset
#   top_cocktail_dfhist <- top_cocktail_df %>% group_by(gene_name) %>% 
#     filter(housekeeping == 1) %>% 
#     summarise(num_hits = sum(targets_logical))
#   top_cocktail_dfhist$status <- "best cocktail"
#   
#   #filter out non-targeted genes to avoid whitespace on the heatmap
#   nottargeted <- top_cocktail_df %>% group_by(gene_name) %>%
#     summarise(sum_target = sum(targets_logical)) %>% 
#     filter(sum_target == 0)
#   top_cocktail_df <- filter(top_cocktail_df, !(gene_name %in% nottargeted$gene_name))
#   
#   #same story for bottom cocktail
#   bottom_cocktail <- cocktail_df2 %>% 
#     filter(cell_line == cell_line_plot,
#            numtargets == 15) %>% 
#     filter(rank == max(rank)) %>%
#     select(1:3) %>% unlist(.)
#   
#   bottom_cocktail_df <- targeting %>% 
#     filter(miR %in% bottom_cocktail, cell_line == cell_line_vec[i]) %>% 
#     mutate(gene_name = fct_reorder(gene_name, housekeeping, .desc = FALSE))
#   #make the bottom_cocktail hist and put them together
#   bottom_cocktail_dfhist <- bottom_cocktail_df %>% group_by(gene_name) %>% 
#     filter(housekeeping == 1) %>% 
#     summarise(num_hits = sum(targets_logical))
#   bottom_cocktail_dfhist$status <- "worst cocktail"
#   cocktail_dfhist <- bind_rows(top_cocktail_dfhist, bottom_cocktail_dfhist)
#   #get rid of not targeted rows for bottom cocktail
#   nottargeted <- bottom_cocktail_df %>% group_by(gene_name) %>%
#     summarise(sum_target = sum(targets_logical)) %>% 
#     filter(sum_target == 0)
#   bottom_cocktail_df <- filter(bottom_cocktail_df, !(gene_name %in% nottargeted$gene_name))
#   
#   ###Actual plotting code###
#   g1 <- ggplot(data = top_cocktail_df, 
#                aes(x = gene_name, y = miR, fill = targets_value)) + 
#     geom_tile() + 
#     scale_fill_gradient2(low = "red", mid = "white", high = "blue") + 
#     theme(axis.text.x=element_text(angle = -60, vjust = 0.5, hjust = 0, size = 10), 
#           legend.position = "none", 
#           text = element_text(size = 16),
#           plot.margin = unit(c(0.2,1.5,0.2,0.2), "cm")) + 
#     ylab("miRNA") +
#     xlab("gene") + labs(tag  = cell_line_vec[i], title = "best predicted cocktail") 
#   
#   g2 <- ggplot(data = bottom_cocktail_df, 
#                aes(x = gene_name, y = miR, fill = targets_value)) + 
#     geom_tile() + 
#     scale_fill_gradient2(low = "red", mid = "white", high = "blue") + 
#     theme(axis.text.x=element_text(angle = -60, vjust = 0.5, hjust = 0, size = 10), 
#           legend.position = "none", 
#           text = element_text(size = 16),
#           plot.margin = unit(c(0.2,1.5,0.2,0.2), "cm")) + 
#     ylab("miRNA") + 
#     xlab("gene") + labs(tag = "B", title = "worst predicted cocktail")
#   
#   g3 <- ggplot(data = cocktail_dfhist, aes(num_hits)) + 
#     geom_histogram(bins =4, fill = "steelblue") + 
#     scale_x_continuous(breaks = c(0,1,2,3)) + 
#     theme_classic() + 
#     theme(text = element_text(size = 14)) + 
#     xlab("number of microRNA that target a given housekeeping gene") +
#     labs(tag = "C") + facet_wrap(vars(status)) 
#   
#   g_list[[i]] <- arrangeGrob(g1,g2,g3, nrow = 1, ncol = 3, widths = c(1,1,0.3))
# }
# 
# png(filename = "./figures/miRCocktailPlot_numtargets.png", height = 1000, width = 1500) 
# plot_grid(g_list[[1]],g_list[[2]], g_list[[3]], 
#           g_list[[4]], g_list[[5]],g_list[[6]], 
#           nrow = 6, ncol = 1)
# dev.off()
#######################miR Overview Figure#########################################

#Set up miR_table to be graphed
#need to aggregate down to gene-miR level (average away cell lines.)
#The data is only processed to include cell lines because it was easier to do the miR
#cocktail selection with that data structure. 
#This takes some time.

summary_tbl <- overall_targeting_mat.df %>% 
  group_by(miR, gene_name, cell_line) %>% 
  summarise(deltaGibbs_mean = mean(mean_deltaGibbs), 
            targets_logical = unique(targets_logical), 
            weighted_deltaGibbs = mean(weighted_deltaGibbs), 
            sd_deltaGibbs = sd(mean_deltaGibbs),
            housekeeping = unique(housekeeping)) 

#We are going to replace NAs in standard deviation with 0. The NAs arise because our
#process of selecting the top 5% of genes for each cell lines results in a slightly 
#different set of genes for each cell line. In some cases, a given gene will not be in the top 
#5 % for each cell line. Our treatment of SD implicitly imputes the gibb_diff to be the mean
#of the deltaGibbs for the cell_lines that were represented for a given gene. For genes that 
#only showed up in the top 5% for one cell line, the SD is 0 because we are implicitly 
#imputing the deltaGibbs for the other 5 cell lines to be the value for the one cell line that we have.
summary_tbl$sd_deltaGibbs[is.na(summary_tbl$sd_deltaGibbs)] <- 0
summary_tbl <- mutate(summary_tbl, 
                      weighted_sd_deltaGibbs = sd_deltaGibbs * targets_logical)
#Now we can calculate our top miRs
summary_tbl2 <- summary_tbl %>% group_by(miR, cell_line) %>%
  summarise(num_target = sum(targets_logical), 
            gibbs_target = sum(weighted_deltaGibbs),
            sd_target = sum(weighted_sd_deltaGibbs),
            num_housekeeping = sum(housekeeping == "yes" & targets_logical ==TRUE)) %>% 
  filter(num_target > 0) %>% arrange(desc(gibbs_target))

miR_table2 <- summary_tbl2 %>% 
  filter(gibbs_target > 10000) %>% ungroup() %>%
  mutate(miR = fct_reorder(miR, gibbs_target, .desc = TRUE))
summary_tbl2 <- summary_tbl2 %>% group_by(miR) %>% 
  summarise(num_target = unique(num_target),
            gibbs_target = mean(gibbs_target))
mir_rank_df <- filter(mir_rank_df, n>1)

png(filename = "./figures/miROverviewPlot.png", width = 750, height = 600)

g1 <- ggplot(data = miR_table2, aes(miR, (gibbs_target))) +
  geom_boxplot(color = "steelblue", fill = "steelblue") + 
  geom_rect(aes(xmin = as.numeric(miR) - 0.5, 
                xmax = as.numeric(miR) + 0.5,
                ymin = 9000, ymax = 9750, fill = num_housekeeping),
            show.legend = TRUE) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = -60, vjust = 0.5, hjust = 0),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20),
        plot.margin = unit(c(0.2,1.5,0.2,0.2), "cm"),
        plot.tag = element_text(face = "bold")) +
  scale_y_continuous(limits = c(9000, 19000),
                     expand = c(0,0), 
                     breaks = seq(10000,19000, by = 1000)) +
  scale_fill_viridis_c() +
  ylab("Projected disruption in network potential") + 
  xlab("miRNA") + 
  labs(tag = "A",
       fill = "Number of\nhousekeeping\ngenes\ntargeted")
g2 <- ggscatter(summary_tbl2, x = 'num_target', y = 'gibbs_target', add = "reg.line",
                size = 1.75, color = "steelblue") +
  stat_cor(label.y = 17500, 
           aes(label = ..rr.label..), size = 5.5) +
  geom_label_repel(aes(
    label = ifelse(miR %in% mir_rank_df$miR, miR,""), 
    color = ifelse(miR %in% mir_rank_dftop$miR, "blue", "red")
  ), 
  box.padding = 1, size = 4.5, 
  nudge_x = 15) + 
  scale_x_continuous(breaks = seq(0,150,by = 25)) +
  xlab("Number of putative mRNA targets") + 
  ylab("") + 
  theme(text = element_text(size = 18),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"),
        plot.tag = element_text(face = "bold"),
        legend.position = "none") + 
  labs(tag = "B")
g1 + annotation_custom(
  ggplotGrob(g2), xmin = 4.5, xmax = 27.5, ymin = 12400, ymax = 19650)
#grid.arrange(g1,g2, nrow =2)
dev.off()


########################Tables#####################

#The xtable code generates Latex tables
xtable(targets_matrix)
xtable(miR_cocktail_matrix)

#Need to make a table of the 144 distinct transcripts repressed
#need to limit to just the "target set" of the top 5% of mRNA by network potential
top_mir <- overall_targeting_mat.df %>% 
  filter(miR == "miR-3613-3p", targets_logical == 1)
top_mir_targets <- unique(as.character(top_mir$gene_name))
top_mir_targets_df <- enframe(top_mir_targets, value = "miR-3613-3p targets") %>% 
  select(-name)
xtable(top_mir_targets_df)


