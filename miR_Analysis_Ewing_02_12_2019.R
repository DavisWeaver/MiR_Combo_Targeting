
##Generate a list of the best miR candidates based on the total predicted Gibbs
library(reshape2)
library(org.Hs.eg.db)
library(readr)
library(dplyr)
library(purrr)
library(tibble)
library(tidyr)
library(stringr)
library(ggplot2)
library(forcats)
library(xtable)
library(circlize)
library(gridExtra)
source('Mir_Analysis_Functions_09_17_2019.R')

##Load file with miRNAs and all the identified mRNAs they match
if(!file.exists("C:/Users/dtw43/Documents/MIR_Combo_Targeting/data files/targeting_mat.Rda")){
  overall_targeting_mat <- calc_target_matrix()
  save(overall_targeting_mat, file = "C:/Users/dtw43/Documents/MIR_Combo_Targeting/data files/targeting_mat.Rda")
} 
load("C:/Users/dtw43/Documents/MIR_Combo_Targeting/data files/targeting_mat.Rda")

#Say if we want to do the analysis on cancer-related genes only or not. 
cancer_logical = FALSE

##Load and process miR expression data
load("C:/Users/dtw43/Documents/MIR_Combo_Targeting/data files/averaged_miR_values.Rda")
averaged_miR_values.df <- as_tibble(averaged_miR_values, rownames = "miR_names")

#Import output from python pipeline
final_df <- read_csv('C:/Users/dtw43/Documents/MIR_Combo_Targeting/data files/EwingsDatasetFinal_11_12_2019.csv')
ewing_full_df <- read_csv('C:/Users/dtw43/Documents/MIR_Combo_Targeting/data files/EwingsDatasetFullClean_11_12_2019.csv')

#Add logical variables from the cosmic database regarding the relationship of the
#identified genes to cancer, see functions script for more documentation
final_df <- add_cosmic_vars(df = final_df)
ewing_full_df <- add_cosmic_vars(df = ewing_full_df)

#Add entrez_id variable to final_df, using bimap interface to link Gene_Name
#variable to entrez_id variable, see functions script for more documentation
final_df <- add_entrez_id(final_df)
ewing_full_df <- add_entrez_id(ewing_full_df)

#Calculate average gibbs_diff and average of several other data objects for each 
#cell_line-Gene_name pair #also can specify what we want the gibbs_diff cutoff 
#to be and whether we want to limit to just genes implicated causally in cancer
final_df_mean <- final_df %>% 
  mean_cellline_calc(topx = TRUE, numgenes = 200, justcancer = cancer_logical, 
                     keep_cell_line = FALSE) 

#This is going to be to determine how many cell lines the top 10 genes fall into.
final_df_mean_barchart <- final_df %>% 
  mean_cellline_calc(topx = TRUE, numgenes = 20, justcancer = TRUE, 
                     keep_cell_line = TRUE) 
num_cell_lines <- 6 #Just defining the number of cell lines so I can calculate a proportion

#aggregate out cell line and calculate fraction of cell lines per gene
final_df_mean_barchart <- final_df_mean_barchart %>% group_by(Gene_Name) %>% 
  summarise(percent_cell_lines = n()/6)

#Okay lets prep data to be graphed
final_df_mean <- final_df_mean %>% 
  filter(!is.na(se_gibbs_diff)) %>% 
  mutate(Gene_Name = fct_reorder(Gene_Name, mean_gibbs_diff, .desc = TRUE))

#Prep cancer only data to be graphed
final_df_mean_cancer <- final_df %>% 
  mean_cellline_calc(topx = TRUE, numgenes = 50, justcancer = TRUE, 
                     keep_cell_line = FALSE) 
final_df_mean_cancer <- final_df_mean_cancer %>% 
  filter(!is.na(se_gibbs_diff)) %>% 
  mutate(Gene_Name = fct_reorder(Gene_Name, mean_gibbs_diff, .desc = TRUE))


#need to prep full data to be graphed later as well
#ewing_full_df_final <- ewing_full_df %>% select(-X1) %>% 
# group_by(cell_line, experiment, condition) %>%
#summarise(gibbs_sum = sum(gibbs)) #%>% #group_by(cell_line, condition) %>%
#summarise(gibbs_mean = mean(gibbs_sum), gibbs_sd =sd(gibbs_sum))
#%>% group_by(cell_line, condition) %>% 
#summarise(gibbs_mean = mean(gibbs_sum))

#This code creates a dataframe with total gibbs calculated for each cell line -experiment pairing
#Commented out code aggregates even farther to just have average gibbs for each ES cell line and 
#THen even farther to just have the global average gibbs/standard deviation
ewing_full_df_finalcontrol <- ewing_full_df %>% select(-X1) %>% 
  group_by(cell_line, experiment) %>% 
  summarise(gibbs_sum = sum(gibbs), 
            expression_sum = sum(expression)) 
ewing_full_df_finalcontrol <- ewing_full_df_finalcontrol %>% group_by() %>%
  mutate(cell_line = fct_reorder(cell_line, gibbs_sum, .desc = TRUE)) 

#group_by(cell_line) %>% summarise(gibbs_mean = mean(gibbs_sum), 
#gibbs_sd = sd(gibbs_sum))
#ewing_full_df_average <- ewing_full_df_finalcontrol %>%
#summarise(gibbs_mean2 = mean(gibbs_mean),
#   gibbs_sd = sd(gibbs_mean))
#changing experiment number to replicate 1,2,3
ewing_full_df_finalcontrol$experiment <- c("one","two","three","one","two","three",
                                           "one","two","three","one","two","three",
                                           "one","two","three","one","two","three")
colnames(ewing_full_df_finalcontrol)[colnames(ewing_full_df_finalcontrol) == "experiment"] <- "replicate"

#Also need to prep data to visualize description of G_i
#Goin to try this filtering out the observations with zero expression
#/zero network potential - that should be one to one. 
hist_df_allgenes <- ewing_full_df %>% 
  filter(expression > 0, gibbs < 0) %>% 
  group_by(Gene_Name) %>% summarise(gibbs = mean(gibbs), 
                                    expression = mean(expression))
hist_df_targetgenes <- final_df %>% 
  group_by(Gene_Name) %>% summarise(gibbs = mean(gibbs), 
                                    expression = mean(expression))
##########################okay lets dive into the miRs##########################
##turn the targeting matrix into a maneuverable dataframe
overall_targeting_mat.df <- clean_target_matrix(cancer_only = cancer_logical)

miR_candidates_matrix <- Generate_MiR_Candidates(num_miRs = 10, 
                                                 cancer_only = cancer_logical)
#Use xtable to generate latex code for a nice table
#xtable(miR_candidates_matrix)

#Generate a list of predicted miR cocktails for each cell line
output <- Generate_MiR_Cocktail(num_miRs = 5, num_targets = 10,
                                fraction_target = 0.2, 
                                cancer_only = FALSE, 
                                keep_cell_line = TRUE)
targets_matrix <- output[[1]]

#The xtable code generates Latex tables
#xtable(targets_matrix)
miR_cocktail_matrix <- output[[2]]
#xtable(miR_cocktail_matrix)
#see functions script for documentation

#Now generate a series of dataframes to prep for chord diagram - first with 5 targets
chord_df_5 <-  PrepForChordDiagram(targeting_df = overall_targeting_mat.df, 
                                   targets = targets_matrix, 
                                   cocktail = miR_cocktail_matrix)

##Now with 10 targets
output <- Generate_MiR_Cocktail(num_miRs = 5, num_targets = 10,
                                fraction_target = 0.1, 
                                cancer_only = FALSE, 
                                keep_cell_line = FALSE)
targets_matrix <- output[[1]]
miR_cocktail_matrix <- output[[2]]
chord_df_10 <-  PrepForChordDiagram(targeting_df = overall_targeting_mat.df, 
                                    targets = targets_matrix, 
                                    cocktail = miR_cocktail_matrix)

##Now with 15 targets
output <- Generate_MiR_Cocktail(num_miRs = 5, num_targets = 15,
                                fraction_target = 0.1, 
                                cancer_only = FALSE, 
                                keep_cell_line = FALSE)
targets_matrix <- output[[1]]
miR_cocktail_matrix <- output[[2]]
chord_df_15 <-  PrepForChordDiagram(targeting_df = overall_targeting_mat.df, 
                                    targets = targets_matrix, 
                                    cocktail = miR_cocktail_matrix)

#also want to make a cord diagram with just one miR
#I think I can use my existing functions to do this
output <- Generate_MiR_Cocktail(num_miRs = 1, num_targets = 100, 
                                fraction_target = 0.1, cancer_only = FALSE, 
                                keep_cell_line = FALSE)
targets_matrix <- output[[1]]
miR_cocktail_matrix <- output[[2]]
chord_df_singlemiR <- PrepForChordDiagram(targeting_df = overall_targeting_mat.df, 
                                          targets = targets_matrix, 
                                          cocktail = miR_cocktail_matrix)
#Want to add a variable for whether or not the MiR is upregulated or downregulated
#Also want to add a variable for the total gibbs_diff of the targeted genes -
#can do this easy with summarise. 
miR_expression_df <- averaged_miR_values.df %>% 
  gather(key = "cell_line", value = "expression", -miR_names) %>% 
  group_by(miR_names) %>% summarise(mean_expression = mean(expression)) 
colnames(miR_expression_df)[colnames(miR_expression_df) == "miR_names"] <- "miR"
miR_expression_df <- mutate(miR_expression_df, miR = gsub("hsa-", "", miR))


#########################Code for making figures and tables ###################
#Chord Diagram 
png("C:/Users/dtw43/Documents/MIR_Combo_Targeting/figures/miR_mapping_chord_Diagram.png", 
    height = 1000, width = 1500, res = 250)
labels <- c("A", "B")
par(mfrow=c(1,2), cex = 0.55, las = 2)
chordDiagram(chord_df_5, scale = TRUE)
put.fig.letter(labels[1], x = 0.05, y = 0.8, size = 1.5, hue = "steel blue")
chordDiagram(chord_df_10, scale = TRUE, link.overlap = TRUE)
put.fig.letter(labels[2], x = 0.05, y = 0.8, size = 1.5, hue = "steel blue")
#chordDiagram(chord_df_15)
dev.off()

#Make single miR chord diagram for intro
png("C:/Users/dtw43/Documents/MIR_Combo_Targeting/figures/single_miR_chord_Diagram.png", 
    height = 1000, width = 1500, res = 250)
par(cex = 0.9)
chordDiagram(chord_df_singlemiR, scale = TRUE)
dev.off()


#Gotta make a table of the frequency/ proportion of oncogenes and tumor supressor genes
gene_type_tbl <- table(final_df_mean$oncogene, final_df_mean$TSG,
                       dnn = c("oncogene", "TSG"))


##Make the mRNA overview figure 
#create a vector that tells GGPlot to color code the labels genes associated with cancer
final_df_mean$text_format <- "black"
final_df_mean$text_format[final_df_mean$cancer_associated == "yes"] <- "red"
color_vec <- final_df_mean$text_format
colors <- c("yes" = "red", "no" = "black")
png(filename = "C:/Users/dtw43/Documents/MIR_Combo_Targeting/figures/mRNA_Targetoverviewplot.png",
    width = 900, height = 600)
g <- ggplot(data = final_df_mean, 
            aes(Gene_Name, mean_gibbs_diff, color = cancer_associated)) 
g + geom_point(size = 0.9) + 
  geom_errorbar(aes(x = Gene_Name, 
                    ymin = mean_gibbs_diff - 2*se_gibbs_diff,
                    ymax = mean_gibbs_diff + 2*se_gibbs_diff)) +
  theme(axis.text.x = element_text(angle = -60, vjust = 0.5, hjust = 0),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 18),
        plot.margin = unit(c(0.2,1,0.2,0.2), "cm")) +
  ylab(expression(paste(Delta, "Network Potential"))) + 
  xlab("Gene Name") + 
  theme(axis.text.x = element_text(color = color_vec),
        legend.position = "none") + 
  scale_color_manual(values = colors)
dev.off()

##Make an mRNA target overview limited to genes associated in cancer
png(filename = "C:/Users/dtw43/Documents/MIR_Combo_Targeting/figures/mRNA_Targetoverviewplot_cancer.png",
    width = 900, height = 600)
g <- ggplot(data = final_df_mean_cancer, 
            aes(Gene_Name, mean_gibbs_diff)) 
g + geom_point(size = 0.9) + 
  geom_errorbar(aes(x = Gene_Name, 
                    ymin = mean_gibbs_diff - 2*se_gibbs_diff,
                    ymax = mean_gibbs_diff + 2*se_gibbs_diff)) +
  theme(axis.text.x = element_text(angle = -60, vjust = 0.5, hjust = 0),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 18),
        plot.margin = unit(c(0.2,1,0.2,0.2), "cm")) +
  ylab(expression(paste(Delta, "Network Potential"))) + 
  xlab("Gene Name") 
dev.off()
#Make a histogram of the distribution of G_i and also mRNA expression
png(filename = "C:/Users/dtw43/Documents/MIR_Combo_Targeting/figures/histogram_allgenes_boxplot.png",
    width = 900, height = 700)
g1 <- ggplot(data = hist_df_allgenes, aes(expression)) + theme_bw() +
  geom_histogram(bins = 70, color = "red3", fill = "darkred") +
  xlab("mRNA expression (normalized count)") + ylab("") + 
  labs(tag = "A") + theme(text = element_text(size = 18))
  
g2 <- ggplot(data = hist_df_allgenes, aes(gibbs)) + theme_bw() +
  geom_histogram(bins = 70, color = "blue", fill = "steelblue") + 
  xlab("Network potential") + ylab("") + labs(tag = "B") + 
  theme(text = element_text(size = 18))  

#code for boxplot of total expression for each cell line
g3 <- ggplot(data = ewing_full_df_finalcontrol, aes(cell_line, expression_sum))  + 
  geom_boxplot() + theme_bw() +
  theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust = 0),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 18)) +
  scale_y_continuous(limits = c(108300, 110000), 
                     breaks = seq(108300, 110000, by = 400)) +
  ylab("Total mRNA expression") + xlab("") +
  labs(tag = "C")
# Code for boxplot of total network potential for each cell line
g4 <- ggplot(data = ewing_full_df_finalcontrol, aes(cell_line, gibbs_sum))  + 
  geom_boxplot() + theme_bw() +
  theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust = 0),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 18)) +
  scale_y_continuous(limits = c(-347500, -341000), 
                     breaks = seq(-347500, -341000, by = 1000)) +
  ylab("Network potential") + xlab("") +
  labs(tag = "D")

grid.arrange(g1, g2, g3, g4, nrow = 2)
dev.off()

#Set up miR_table to be graphed
#need to aggregate down to gene-miR level (average away cell lines.)
#The data is only processed to include cell lines because it was easier to do the miR
#cocktail selection with that data structure. 

summary_tbl <- overall_targeting_mat.df %>% 
  group_by(miR, Gene_Name) %>% 
  summarise(gibbs_diff_mean = mean(mean_gibbs_diff), 
            targets_logical = unique(targets_logical), 
            weighted_gibbs_diff = mean(weighted_gibbs_diff), 
            sd_gibbs_diff = sd(mean_gibbs_diff), 
            num_cell_lines = n()) 

#We are going to replace NAs in standard deviation with 0. The NAs arise because our
#process of selecting the top 5% of genes for each cell lines results in a slightly 
#different set of genes for each cell line. In some cases, a given gene will not be in the top 
#5 % for each cell line. Our treatment of SD implicitly imputes the gibb_diff to be the mean
#of the gibbs_diff for the cell_lines that were represented for a given gene. For genes that 
#only showed up in the top 5% for one cell line, the SD is 0 because we are implicitly 
#imputing the gibbs_diff for the other 5 cell lines to be the value for the one cell line that we have.
summary_tbl$sd_gibbs_diff[is.na(summary_tbl$sd_gibbs_diff)] <- 0
summary_tbl <- mutate(summary_tbl, 
                      weighted_sd_gibbs_diff = sd_gibbs_diff * targets_logical)
#Now we can calculate our top miRs
summary_tbl2 <- summary_tbl %>% group_by(miR) %>%
  summarise(num_target = sum(targets_logical), 
            gibbs_target = sum(weighted_gibbs_diff),
            sd_target = sum(weighted_sd_gibbs_diff)) %>% 
  filter(num_target > 0) %>% arrange(desc(gibbs_target))

miR_table <- summary_tbl2 %>% ungroup() %>%
  left_join(miR_expression_df) %>% 
  slice(0:20) %>%
  mutate(miR = fct_reorder(miR, gibbs_target, .desc = TRUE))

miR_table2 <- summary_tbl2 %>% left_join(miR_expression_df) %>% 
  filter(gibbs_target > 10000) %>%
  mutate(miR = fct_reorder(miR, gibbs_target, .desc = TRUE))
png(filename = "C:/Users/dtw43/Documents/MIR_Combo_Targeting/figures/miROverviewPlot.png", 
    width = 900, height = 600)
g <- ggplot(data = miR_table, aes(miR, gibbs_target)) 
g + geom_col(color = "steelblue", fill = "steelblue", width = 0.8) + 
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = 0),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20)) +
  geom_errorbar(aes(ymin = gibbs_target - 1.96*sd_target, 
                ymax = gibbs_target + 1.96*sd_target)) +
  scale_y_continuous(limits = c(0, 19000), breaks = seq(0, 19000, by = 3000)) +
  ylab("Projected disruption in Network Potential") + 
  xlab("miR")
dev.off()

#Need to make a table of the 144 distinct transcripts repressed
#need to limit to just the "target set" of the top 5% of mRNA by network potential
top_mir <- overall_targeting_mat.df %>% 
  filter(miR == "miR-3613-3p", targets_logical == 1)
top_mir_targets <- unique(as.character(top_mir$Gene_Name))
top_mir_targets_df <- enframe(top_mir_targets, value = "miR-3613-3p targets") %>% 
  select(-name)
xtable(top_mir_targets_df)
