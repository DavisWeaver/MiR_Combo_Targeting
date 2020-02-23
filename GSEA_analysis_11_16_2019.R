#Use fast_GSEA
#Download gene sets using msigdbr
library(fgsea)
library(msigdbr)
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
source('Mir_Analysis_Functions_09_17_2019.R')

#Round results to make tables more readable
options(scipen = 0)
options(digits = 2)

#Import Data from python pipeline
#final_df <- read_csv('C:/Users/dtw43/Documents/MIR_Combo_Targeting/data files/EwingsDatasetFinal_11_12_2019.csv')
ewing_full_df <- read_csv('C:/Users/dtw43/Documents/MIR_Combo_Targeting/data files/EwingsDatasetFullClean_11_12_2019.csv')
#clean data
final_df <- add_cosmic_vars(df = final_df)
final_df <- add_entrez_id(final_df)
ewing_full_df <- ewing_full_df %>% add_cosmic_vars() %>% add_entrez_id
#calculate top x genes to conduct the GSEA on.
final_df_mean <- final_df %>%
  mean_cellline_calc(topx = TRUE, numgenes = 200, justcancer = FALSE,
                     keep_cell_line = FALSE)
ewing_full_df_mean <- ewing_full_df %>% group_by(Gene_Name) %>%
  summarise(gibbs = mean(gibbs), expression = mean(expression))
  
#Pull out list of gene names for GSEA
gene_set_test <- ewing_full_df_mean$gibbs
gene_set_test <- abs(gene_set_test) #take absolute value of gibbs so the algorithm knows how to rank
names(gene_set_test) <- ewing_full_df_mean$Gene_Name
#Download background gene sets from msigdb
#First using msigdb category H
pathways_hallmarks <- gmtPathways("C:/Users/dtw43/Documents/MIR_Combo_Targeting/data files/h.all.v7.0.symbols.GMT")
#second using msigdb category C6 - cancer
pathways_cancer <- gmtPathways("C:/Users/dtw43/Documents/MIR_Combo_Targeting/data files/c6.all.v7.0.symbols.GMT")

#perform analysis and store results 
results_cancer <- fgsea(pathways_cancer, gene_set_test, nperm = 100)
results_hallmarks <- fgsea(pathways_hallmarks, gene_set_test, nperm = 100)

#Get ready to make a table with hallmarks results
results_hallmarks <- results_hallmarks %>% 
  mutate(pathway = gsub("HALLMARK_", "", pathway)) %>% 
  select(-nMoreExtreme, -leadingEdge) %>% filter(padj < 0.05) %>%
  arrange(-desc(padj))

xtable(results_hallmarks)
