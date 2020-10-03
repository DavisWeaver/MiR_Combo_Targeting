library(tidyr)
library(dplyr)
library(stringr)
#The filepaths are absolute in this because this file is written to be run on 
#a computing cluster.
source('./NP_Pipeline_MiR.R')
source('./NP_Funcs_MiR.R')
df <- read.csv('./data_files/rld_Counts.csv') %>% as_tibble()

#Do some quick and dirty data cleaning
colnames(df)[1] <- "gene_name"
df <- df %>% gather(experiment_id, expression, -gene_name) #convert from wide to long format
experiment_breakout <- str_split(df$experiment_id, pattern = "_", 
                                 simplify = TRUE) %>% as_tibble() %>%
  select(1:3)
colnames(experiment_breakout) <- c("experiment_num", "cell_line", "condition") #breakout the experiment barcode to constituent pieces
df_ready <- cbind(df, experiment_breakout) %>% 
  filter(condition == "DMSO") %>% mutate(gene_name = as.character(gene_name)) #Only working on control cell lines

#run pipeline to save new data file to cluster.
gibbs_pipeline(df = df_ready, filter_proportion = 1, ncores = 3)  #default is 1 less than the machine - can do more on cluster


