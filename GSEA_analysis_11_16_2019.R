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
source('~/MIR_Combo_Targeting/Code/miRNA_Targeting/Mir_Analysis_Functions_09_17_2019.R')

#Round results to make tables more readable
options(scipen = 0)
options(digits = 2)

#Import Data from python pipeline
#final_df <- read_csv('C:/Users/dtw43/Documents/MIR_Combo_Targeting/data files/EwingsDatasetFinal_11_12_2019.csv')
load('C:/Users/dtw43/Documents/MIR_Combo_Targeting/data_files/gibbsdf_final_EWS.Rda')
final_df <- gibbs_df_final %>% add_cosmic_vars() %>% add_entrez_id()

final_df_mean <- final_df %>% group_by(gene_name) %>%
  summarise(gibbs = mean(gibbs), 
            expression = mean(expression),
            deltaGibbs = mean(deltaGibbs))
  
#Pull out list of gene names for GSEA
gene_set_test <- final_df_mean$deltaGibbs
gene_set_test <- abs(gene_set_test) #take absolute value of gibbs so the algorithm knows how to rank
names(gene_set_test) <- final_df_mean$gene_name
#Download background gene sets from msigdb
#First using msigdb category H
pathways_hallmarks <- gmtPathways("C:/Users/dtw43/Documents/MIR_Combo_Targeting/data_files/h.all.v7.0.symbols.GMT")
#second using msigdb category C6 - cancer
pathways_cancer <- gmtPathways("C:/Users/dtw43/Documents/MIR_Combo_Targeting/data_files/c6.all.v7.0.symbols.GMT")

pathways_mirbio <- read_tsv("C:/Users/dtw43/Documents/MIR_Combo_Targeting/data_files/miRNA_biogenesis_genes_list.txt",
                            col_names = c("gene_name", "entrez_id"))
pathway_mirbio <- list(pathways_mirbio$gene_name)
names(pathway_mirbio) <- "mir_biogenesis"
pathways_hallmarks <- c(pathways_hallmarks, pathway_mirbio)

#Get a pipeline 
#perform analysis and store results 
results_cancer <- fgsea(pathways_cancer, gene_set_test, nperm = 500)
results_hallmarks <- fgsea(pathways_hallmarks, gene_set_test, nperm = 500)
#Get ready to make a table with hallmarks results
results_hallmarks <- results_hallmarks %>% 
  mutate(pathway = gsub("HALLMARK_", "", pathway)) %>% 
  filter(padj < 0.05) %>%
  arrange(-desc(padj)) %>% 
  select(-leadingEdge)

xtable(results_hallmarks)
