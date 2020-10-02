clean_RNAseq <- function(df) {
  library(dplyr)
  library(tidyr)
  #now gotta get rid of gene name so I can log2 scale the entire numeric_matrix
  df <- df %>% log2() %>% 
    as_tibble(rownames = "gene_name") %>%
    gather(df, barcode, expression,-gene_name)
  return(df)
}

gibbs_pipeline <- function(df, filter_proportion = 0.05, 
                           common_subgraph = FALSE, ncores = 7) {
  require(igraph)
  require(arules)
  require(doParallel)
  require(foreach)
  source('./NP_Funcs_MiR.R')
  
  # cancer_type <- paste0(cancer_type, collapse = "")
  # load(
  #   paste0(
  #     "C:/Users/dtw43/Documents/MIR_Combo_Targeting/data files/rnaseq_",
  #     cancer_type, 
  #     ".Rda"
  #   )
  # )
  #Gotta get rid of negative/0 expression values
  df$expression[df$expression <0] <- 0
  #For debugging - make this dataset way way smaller. 
  #Biogrid Data
  if(!dir.exists("./biogrid")) {
    download.file("https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-3.5.171/BIOGRID-ORGANISM-3.5.171.tab2.zip",
                  destfile = "./biogrid.zip")
    unzip('./biogrid.zip', exdir = 'biogrid_unzip')
    file.remove("./biogrid.zip")
  }
  biogrid <-
    read.delim(
     "./biogrid_unzip/BIOGRID-ORGANISM-Homo_sapiens-3.5.171.tab2.txt",
      header = TRUE
    )
  
  biogrid <-
    biogrid[, 8:9] # isolate Official Symbol Interactor A & B columns
  
  #clean the bioGrid
  output = cleanBiogrid(df, biogrid)
  network <- output[[1]]
  networkDF <- output[[2]]
  genesInBoth <- output[[3]]
  gridDF <- output[[4]]
  #rm(list = c("output","rnaseqnorm_df")) #keep the memory free.
  
  #have to do something about the negative expression values/the -inf expression values
  #Add expression values to graph
  output <- addExpression(df = networkDF, grid = gridDF)
  expressionNetwork <- output[[1]]
  cellLines <- output[[2]] 
  cellLines <- cellLines[2:length(cellLines)] #get rid of the gene_name column
  #wideDF <- output[[3]]
  #wideExpressionDF[wideExpressionDF != 0] <- 0 # make blank dataframe for calculations
  # make blank matrix to store calculations
  #rm(output) #keep memory free
  gibbs_list <- list()
  # Create List of Neighbors for Each Node
  
  neighborList <-
    lapply(1:length(genesInBoth),
           getNeighbors, network = expressionNetwork)
  names(neighborList) <- genesInBoth
  
  # Calculate Gibbs Free Energy Value for Each Node
  #Start at 2 because cellLines[1] is gene_name
  for (exp in 1:length(cellLines)) {
    exp_df <- filter(networkDF, experiment_id == cellLines[exp])
    #tic("gene")
    expExpression <-
      get.vertex.attribute(expressionNetwork, cellLines[exp]) # expression for all vertices (genes) in one experiment
    names(expExpression) <- genesInBoth
    expExpression <- expExpression[names(expExpression) %in% exp_df$gene_name] #subset on just the genes for which we have data
    gibbs_vec <- vector(mode= "double")
    for (i in 1:length(exp_df$gene_name)) {
      gibbs_vec[i] <-
        calcGibbs(exp_df$gene_name[i], expressionVector = expExpression, 
                  neighborlist = neighborList)
    }
    exp_df$gibbs <- gibbs_vec
    gibbs_list[[exp]] <- exp_df
  }
  gibbs_df <- bind_rows(gibbs_list)
  
  #rm(gridDF, network, output, exp_df, rnaseqnorm_df, networkDF)
  #Clear up some RAM
  #Putting the rest of the gibbs pipeline back in####
  # Create List of Genes with Lowest Gibbs for Each Cell Line #
  #Going to filter the whole network based on gibbs so don't need this top percentage
  # Going to need to create a new subgraph for each sample
  #gibbs_df_long <- gather(gibbs_df, key = "experiment_id", value = "gibbs", -gene_name)
  #gibbs_df_long <- gibbs_df_long %>% mutate(experiment_id = as.character(experiment_id)) %>%
  #  mutate(gibbs = as.numeric(gibbs)) %>% filter(gibbs < 0) #%>% #Don't need t
  
  
  
  # For Each Experiment, Calculate New Network Gibbs with Each Target Gene Removed #
  #
  #get a experiment_id_vector to iterate through
  experiment_id_vec <- unique(gibbs_df$experiment_id)
  #Lets try to multi-core this step?
  c1 <- makeCluster(ncores)
  registerDoParallel(c1)
  gibbs_df_final <- 
    foreach(exp = 1:length(experiment_id_vec), .combine = rbind, 
            .packages = c('dplyr', 'igraph')) %dopar% {
              source('./NP_Funcs_MiR.R')
              CalcNewGibbs(gibbs_df, experiment_id_vec[exp], 
                           network = expressionNetwork)
            }
  
  #for (exp in 1:length(experiment_id_vec)) {
  # gibbs_df_exp <- CalcNewGibbs(gibbs_df_disc, experiment_id_ind = experiment_id_vec[exp], 
  # network = expressionNetwork)
  #gibbs_df_final[[exp]] <- gibbs_df_exp
  #}
  #gibbs_df_final <- bind_rows(gibbs_df_final)
  save(gibbs_df_final, file = paste0("./data_files/gibbsdf_final_EWS.Rda")) 
}

#Do a log transform on this stuff (logbase2 and then put it in)
#write_csv(rnaseq.df,
#         path = "C:/Users/dtw43/Documents/MIR_Combo_Targeting/data files/clean_TCGA_SARC.csv")