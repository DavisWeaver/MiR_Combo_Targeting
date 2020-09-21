#This file does the following:
#1. Adds some new variables to the full gibbs data frame from the Gibbs diff calculation. 
#2. Generate a list of the best miR candidates based on the total predicted disrupted Gibbs
#3. Generate miR cocktails that maximize on-target effects and minimize effects on housekeeping genes.
#4. Saves the post-processed final_df, mir_candidates_matrix, targets_df (important protein targets), 
#and cocktails_df for use by the plotting script. 
##
library(reshape2)
library(CePa)
library(biomaRt)
library(org.Hs.eg.db)
library(EnsDb.Hsapiens.v79)
library(readr)
library(dplyr)
library(purrr)
library(gtools)
library(tibble)
library(tidyr)
library(stringr)
library(forcats)
source('./miR_Analysis_Functions.R')

##Load file with miRNAs and all the identified mRNAs they match + clean for further analysis
#see functions script for more documentation about what is being done to this dataset.
if(!file.exists("./data_files/targeting_mat.Rda")){
  overall_targeting_mat <- calc_target_matrix()
  overall_targeting_mat.df <- clean_target_matrix(cancer_only = FALSE)
  save(overall_targeting_mat.df, file = "./data_files/targeting_mat.Rda")
} 
load("./data_files/targeting_mat.Rda")


#Import output from main Network Potential Pipeline
load('./data_files/gibbsdf_final_EWS.Rda')

#Add logical variables from the cosmic database regarding the relationship of the
#identified genes to cancer, see functions script for more documentation
final_df <- add_cosmic_vars(df = gibbs_df_final)
#ewing_full_df <- add_cosmic_vars(df = ewing_full_df)

#Add entrez_id variable to final_df, using bimap interface to link gene_name
#variable to entrez_id variable, see functions script for more documentation
final_df <- add_entrez_id(final_df)
save(final_df, file = './data_files/final_df.Rda')

#Download gene expression data from CCLE
#Download cell line information from CCLE 
CCLE_meta <- read_tsv(file = "./data_files/Cell_lines_annotations_20181226.txt")
CCLE_rna_df <- read.gct(file = "./data_files/CCLE_RNAseq_genes_rpkm_20180929.gct")
CCLE_rna_df <- CCLE_rna_df %>% as_tibble(rownames = "ensemble_id")

#Download normalized protein data
CCLE_protein_df <- read_csv(file = "./data_files/protein_quant_current_normalized.csv")
CCLE_df <- clean_CCLE(CCLE_rna_df, CCLE_protein_df, CCLE_meta)
save(CCLE_df, file = './data_files/CCLE_df.Rda')

##########################okay lets dive into the miRs##########################
miR_candidates_matrix <- Generate_MiR_Candidates(targeting_df = overall_targeting_mat.df,
                                                 df = final_df,
                                                 num_miRs = 10,
                                                 cancer_only = FALSE)
save(miR_candidates_matrix, file = "./data_files/MiRcandidates_matrix.Rda")
#Use xtable to generate latex code for a nice table
#xtable(miR_candidates_matrix)

#Generate a list of predicted miR cocktails for each cell line - this is a time and memory-intensive
#function that uses multi-threading. Make sure to set ncores < the total number of cores 
#on your machine or I am pretty sure your computer will catch on fire. the max used number of cores is 
#6 (because we are executing all 6 cell lines simultaneously.)

#see functions script for documentation
output <- Generate_MiR_Cocktail(targeting_df = overall_targeting_mat.df,
                                df = final_df,
                                num_miRs = 3, num_targets = 10,
                                cancer_only = FALSE, targetweight = 1, 
                                offtargetweight = 10, ncores = 6)

#save
save(output, file = "./data_files/MiRCocktail_output.Rda")


#Repeat for a few different conditions
output <- Generate_MiR_Cocktail(targeting_df = overall_targeting_mat.df,
                                df = final_df,
                                num_miRs = 3, num_targets = 5,
                                cancer_only = FALSE, targetweight = 1, 
                                offtargetweight = 10, ncores = 6)

#save
save(output, file = "./data_files/MiRCocktail5_output.Rda")

#Repeat for a few different conditions
output <- Generate_MiR_Cocktail(targeting_df = overall_targeting_mat.df,
                                df = final_df,
                                num_miRs = 3, num_targets = 15,
                                cancer_only = FALSE, targetweight = 1, 
                                offtargetweight = 10, ncores = 6)
#save
save(output, file = "./data_files/MiRCocktail15_output.Rda")
