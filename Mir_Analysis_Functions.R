
calc_target_matrix <- function(){
  #convert these to entrez
  library(org.Hs.eg.db)
  library(reshape2)
  setwd("G:/My Drive/MIR_Combo_Targeting/data_files")
  
  
  ## Bimap interface:
  x <- org.Hs.egSYMBOL
  # Get the gene symbol that are mapped to an entrez gene identifiers
  mapped_genes <- mappedkeys(x)
  # Convert to a list
  converted_symbols <- as.list(x[mapped_genes])
  converted_symbols_map <- melt(converted_symbols)
  tmp <- converted_symbols_map$L1
  converted_symbols_map <- converted_symbols_map[,1]
  names(converted_symbols_map) <- tmp
  
  final_df <- read_csv('G:/My Drive/MIR_Combo_Targeting/data_files/EwingsDatasetFinal.csv')  
  gene_list = unique(final_df$gene_name)
  
  #convert these to entrez
  
  gene_names_entrez <- c()
  gene_names <- c()
  
  gene_names_entrez <- names(converted_symbols_map)[match(gene_list, converted_symbols_map)]
  gene_names <- c(gene_names,gene_list[!is.na(gene_names_entrez)])
  gene_names_entrez <- gene_names_entrez[!is.na(gene_names_entrez)]
  
  
  names(gene_names_entrez) <- gene_names#gene_names
  
  load('all_miRNA_targets.rda')
  
  overall_targeting_mat <- matrix(0, nrow=length(gene_names_entrez),ncol=length(names(targets_list)))
  
  rownames(overall_targeting_mat) <- gene_names_entrez
  colnames(overall_targeting_mat) <- names(targets_list)
  
  for(miR_name in names(targets_list)){
    for(gene_name in gene_names_entrez){
      overall_targeting_mat[gene_name, miR_name] <- (gene_name %in% targets_list[[miR_name]])		
    }
  }
  
  head(colnames(overall_targeting_mat)[order(colSums(-overall_targeting_mat))])
  head(-sort(-colSums(overall_targeting_mat)))
  
  ind = 1
  ind.max <- ordercolSums(-overall_targeting_mat)[ind]
  rownames(overall_targeting_mat)[which(overall_targeting_mat[,ind.max]==1)]
  return(overall_targeting_mat)
}

clean_CCLE <- function(df_rna, df_protein, df_meta) {
  #Pull out a vector of the cell lines we care about
  df_meta <- df_meta %>% 
    dplyr::filter(type %in% c("Ewings_Sarcoma", "Ewings_sarcoma",
                              "ewings_sarcoma", "ewings_Sarcoma"))
  ewing_cell_lines <- df_meta$CCLE_ID
  
  #Processing the RNA seq data.
  # 1. Convert from ensembl.gene to gene.symbol
  df_rna <- df_rna %>% 
    mutate(ensemble_id = str_replace(ensemble_id, '\\...$', ""))
  ensemble_vec <- df_rna$ensemble_id
  geneIDs <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ensemble_vec,
                               keytype = "GENEID",
                               columns = c("SYMBOL","GENEID"))
  colnames(geneIDs) <- c("gene_name", "ensemble_id")
  
  
  df_rna <- df_rna %>% 
    left_join(geneIDs) %>%
    pivot_longer(-c(ensemble_id, gene_name), names_to = "cell_line", 
                 values_to = "expression") %>% 
    mutate(log_expression = log2(as.numeric(expression)),
           expression = as.numeric(expression)) %>% 
    filter(!is.infinite(log_expression))
    
  #clean the protein dataset
  df_protein <- df_protein %>% dplyr::select(-contains("Peptides")) %>% 
    pivot_longer(contains('TenPx'), names_to = c("cell_line", "tenplex_number"), 
                 names_sep = "_TenPx", values_to = "protein_expression") %>% 
    dplyr::filter(!is.na(protein_expression)) %>% 
    group_by(Gene_Symbol, cell_line) %>% 
    summarise(protein_expression = mean(protein_expression))
  colnames(df_protein)[1]<- "gene_name"

  #Put them together
  CCLE_df <- left_join(df_protein, df_rna, by = c("cell_line", "gene_name")) %>% 
    dplyr::filter(!is.na(log_expression))
  return(CCLE_df)
}
add_cosmic_vars <- function(df){#Lets download COSMIC's list of genes associated in cancer
  cosmic_df <- read_csv('G:/My Drive/MIR_Combo_Targeting/data_files/Census_allMon_2019.csv')
  oncogene_df <- cosmic_df %>% filter(grepl('oncogene', get('Role in Cancer')))
  oncogene_list <- oncogene_df$`Gene Symbol`
  tumorsuppressor_df <- cosmic_df %>% filter(grepl('TSG', get('Role in Cancer')))
  TSG_list <- tumorsuppressor_df$'Gene Symbol'
  
  #create a variable for whether or not the identified genes in our dataset are oncogenes
  yesoncogene_df <- filter(df, gene_name %in% oncogene_list)
  nooncogene_df <- filter(df, !(gene_name %in% oncogene_list))
  yesoncogene_df$oncogene <- "yes"
  nooncogene_df$oncogene <- "no"
  df <- rbind(yesoncogene_df, nooncogene_df)
  
  #Create variable for whether or not identified genes in our dataset are tumor suppressor genes
  yesTSG_df <- filter(df, gene_name %in% TSG_list)
  noTSG_df <- filter(df, !(gene_name %in% TSG_list))
  yesTSG_df$TSG <- "yes"
  noTSG_df$TSG <- "no"
  df <- rbind(yesTSG_df, noTSG_df)
  
  #create variable for whether or not identified genes are associated with cancer at all
  yescancer_df <- filter(df, gene_name %in% cosmic_df$'Gene Symbol')
  nocancer_df <- filter(df, !(gene_name %in% cosmic_df$'Gene Symbol'))
  yescancer_df$cancer_associated <- "yes"
  nocancer_df$cancer_associated <- "no"
  df <- rbind(yescancer_df, nocancer_df)
  
  #Also need to create a variable describing whether a given gene is an "essential housekeeping gene"
  housekeeping_df <- read_csv("G:/My Drive/MIR_Combo_Targeting/data_files/Housekeeping_GenesHuman.csv")
  yeshousekeeping_df <- filter(df, gene_name %in% housekeeping_df$Gene.name)
  nohousekeeping_df <- filter(df, !(gene_name %in% housekeeping_df$Gene.name))
  yeshousekeeping_df$housekeeping <- "yes"
  nohousekeeping_df$housekeeping <- "no"
  df <- rbind(yeshousekeeping_df, nohousekeeping_df)
  return(df)
}

#Add entrez_id to final_df
#Add column for entrez id
## Bimap interface:
add_entrez_id <- function(df) {
  require(org.Hs.eg.db)
  require(dplyr)
  gene_list <- unique(df$gene_name)
  x <- org.Hs.egSYMBOL
  # Get the gene symbol that are mapped to an entrez gene identifiers
  mapped_genes <- mappedkeys(x)
  # Convert to a list
  converted_symbols <- as.list(x[mapped_genes])
  converted_symbols_map <- melt(converted_symbols)
  rm(list = c("converted_symbols", "x", "mapped_genes")) #Keep memory free
  colnames(converted_symbols_map) <- c("gene_name", "entrez_id")
  #convert factor to numeric
  converted_symbols_map <- converted_symbols_map %>% 
    mutate(entrez_id = as.numeric(entrez_id), 
           gene_name = as.character(gene_name)) %>% 
    filter(gene_name %in% gene_list)
  
  #add entrez_id column
  df <- df %>% left_join(converted_symbols_map) %>% filter(!is.na(entrez_id))
  return(df)
}

#This function aggregates our parent dataframe by cell_line, returning an 
#average value for any continuous variable and whatever the value was for 
#any variables that don't differ by experiment. It also pairs down the dataset 
#to just the genes we want to deal with, if we want just cancer associated genes 
#for example
mean_cellline_calc <- function(df, topx = TRUE, numgenes = 50,
                               justcancer = FALSE, keep_cell_line = FALSE) {
  
  if(keep_cell_line == FALSE) {
    df <- df %>% group_by(gene_name) %>% 
      summarise(
        mean_deltaGibbs = mean(deltaGibbs), 
        se_deltaGibbs = sd(deltaGibbs) / sqrt(n()),
        mean_expression = mean(expression),
        mean_gibbs = mean(gibbs),
        oncogene = unique(oncogene),
        TSG = unique(TSG),
        cancer_associated = unique(cancer_associated),
        housekeeping = unique(housekeeping)
      ) %>% arrange(desc(mean_deltaGibbs))
  } else {
    df <- df %>% group_by(cell_line, gene_name) %>% 
      summarise(
        mean_deltaGibbs = mean(deltaGibbs), 
        se_deltaGibbs = sd(deltaGibbs) / sqrt(n()),
        mean_expression = mean(expression),
        mean_gibbs = mean(gibbs),
        oncogene = unique(oncogene),
        TSG = unique(TSG),
        cancer_associated = unique(cancer_associated),
        housekeeping = unique(housekeeping)
      ) %>% arrange(desc(mean_deltaGibbs))
  }
  if(justcancer == TRUE) {
    df <- df %>% filter(cancer_associated == "yes")
  }
  if(topx == TRUE){
    genes_target <- unique(df$gene_name)[1:numgenes]
    df <- df %>% filter(gene_name %in% genes_target)
  } 
  return(df)
}

clean_target_matrix <- function(targeting_matrix = overall_targeting_mat, 
                                df = final_df, cancer_only = FALSE){
  
  ## Bimap interface:
  x <- org.Hs.egSYMBOL
  # Get the gene symbol that are mapped to an entrez gene identifiers
  mapped_genes <- mappedkeys(x)
  # Convert to a list
  converted_symbols <- as.list(x[mapped_genes])
  converted_symbols_map <- melt(converted_symbols)
  colnames(converted_symbols_map) <- c("gene_name", "entrez_id")
  converted_symbols_map <- mutate(converted_symbols_map, 
                                  entrez_id = as.numeric(entrez_id))
  #Convert data from wide to long so that each row is an observation and each column a variable
  targeting_df <- targeting_matrix %>% 
    as_tibble(rownames = "entrez_id") %>% 
    gather(key = "miR",value = "targets_logical", -entrez_id) %>% 
    filter(entrez_id %in% df$entrez_id) %>% 
    mutate(entrez_id = as.numeric(entrez_id)) %>%
    left_join(converted_symbols_map) %>% filter(!is.na(entrez_id))
  
  #Now need to add average gibbs diff for each gene to this as well. 
  df_mean <- df %>% group_by(entrez_id, cell_line) %>%
    summarise(mean_deltaGibbs = mean(deltaGibbs), 
              sd_deltaGibbs = sd(deltaGibbs),  
              cancer_associated = unique(cancer_associated),
              housekeeping = unique(housekeeping)) 
  
  #add weighted gibbs_df as well
  targeting_df <- targeting_df %>% left_join(df_mean) %>% 
    mutate(weighted_deltaGibbs = targets_logical * mean_deltaGibbs)
  if(cancer_only == TRUE) {
    targeting_df <- filter(targeting_df, cancer_associated == "yes")
  }
  
  #Remove hsa- prefix from the miR names
  targeting_df <- targeting_df %>% mutate(miR = gsub("hsa-", "", miR))
  
  return(targeting_df)
  
}
#Pulling this out of the analysis code - generate a list of top targets for 
#each cell line 
Generate_MiR_Candidates <- function(targeting_df = overall_targeting_mat.df,
                                    df = final_df,
                                    num_miRs = 5,
                                    cancer_only = FALSE) {
  cell_line_list <- unique(df$cell_line)
  miRs_df <- matrix(nrow = num_miRs, ncol = length(cell_line_list))
  for(i in 1:length(cell_line_list)) {
    targeting_df_cellline <- filter(targeting_df, cell_line == cell_line_list[i])
    if(cancer_only == TRUE) {
      summary_tbl <- targeting_df_cellline %>% 
        group_by(miR) %>% filter(cancer_associated == "yes") %>%
        summarise(num_target = sum(targets_logical), 
                  gibbs_target = sum(weighted_deltaGibbs)) %>% 
        filter(num_target > 0) %>% arrange(desc(gibbs_target)) %>% 
        slice_head(n = num_miRs)
    } else {
      summary_tbl <- targeting_df_cellline %>% 
        group_by(miR) %>% summarise(num_target = sum(targets_logical), 
                                    gibbs_target = sum(weighted_deltaGibbs)) %>% 
        filter(num_target > 0) %>% arrange(desc(gibbs_target)) %>% 
        slice_head(n = num_miRs)
      miRs_df[,i] <- summary_tbl$miR
    }
  }
  colnames(miRs_df) <- cell_line_list
  return(miRs_df)
}
#The purpose of this function is to identify a cocktail of miRs that all 
#preferentially downregulate several genes of interest while exhibiting limited 
#shared targets outside of the genes of interest -- We are going to limit the number of miRs to 3 so that a brute force approach is possible
#going to need to toss this to the cluster...
Generate_MiR_Cocktail <- function(targeting_df = overall_targeting_mat.df,
                                  df = final_df,
                                  num_miRs = 3, num_targets = 10,
                                  cancer_only = FALSE,
                                  targetweight = 1,
                                  offtargetweight = 1,
                                  fraction_target = 0.2,
                                  ncores = 6){ 
  require(arules)
  require(doParallel)
  require(foreach)
  if(num_miRs >3) {
    return("num_miRs must be <=3")
  }
  
  #Generate a matrix of the top mRNA by network_potential for each cell line
  cell_line_list <- unique(df$cell_line)
  targets_df <- matrix(nrow = num_targets, ncol = length(cell_line_list))
  for(i in 1:length(cell_line_list)) {
    cell_line_df <- filter(df, cell_line == cell_line_list[i])
    if(cancer_only == TRUE) {
      cell_line_df <- cell_line_df %>% group_by(gene_name) %>%
        filter(cancer_associated == "yes") %>%
        summarise(mean_deltaGibbs = mean(deltaGibbs)) %>% 
        arrange(desc(mean_deltaGibbs)) %>% slice_head(n = num_targets)
    } else{
      cell_line_df <- cell_line_df %>% group_by(gene_name) %>%
        filter(housekeeping == "no") %>% 
        summarise(mean_deltaGibbs = mean(deltaGibbs)) %>% 
        arrange(desc(mean_deltaGibbs)) %>% slice_head(n = num_targets)
    }
    targets_df[,i] <- cell_line_df$gene_name
  }
  colnames(targets_df) <- cell_line_list
  targets_df <- as_tibble(targets_df)
  
 

  c1 <- makeCluster(ncores)
  registerDoParallel(c1)
  cocktail_df <- 
    foreach(i = 1:length(cell_line_list), .combine = rbind,
            .packages = c('tidyverse', 'gtools')) %dopar% {
              
              #Initialize loss function
              #loss function is calculating predicted network disruption where protein targets are considered negative and housekeeping genes are 
              #considered positive. Sum up all the hits and theres your one number. As an added wrinkle, Any mRNA that is only targeted by one miR in the 
              #cocktail is considered to be unaffected (radiation therapy idea). 2 hits = 40% repression, 3 hits = 60% repression
              check_combination <- function(cocktail){ 
                #Need to generate a targeting matrix for the cocktail
                #iterate through the cocktail to generate a new targeting matrix for just 
                #cocktail components - need to do it this way to create duplicates if necessary
                cocktail_df <- list()
                for(i in 1:length(cocktail)) {
                  cocktail_df[[i]] <- filter(targeting_matrix, miR == cocktail[i])
                }
                cocktail_df <- bind_rows(cocktail_df)
                #need to iterate through the columns now to tally up the points
                point_vec <- vector(mode = "numeric", length = ncol(cocktail_df))
                for(i in 2:ncol(cocktail_df)) {
                  if(sum(cocktail_df[,i]) == max(cocktail_df[,i])){ #This step should ensure that single hits are considered as nothing
                    point_vec[i] <- 0
                  } else {
                    point_vec[i] <- sum(cocktail_df[,i])
                  }
                }
                return(sum(point_vec)) #return the sum of the point_vec: approximates N_targets_repressed - N_housekeeping_repressed
              }
              
              cell_line_targets <- unlist(targets_df[,cell_line_list[i]])
              #pull out the non-target vector
              non_targets_df <- filter(df, housekeeping == "yes")
              non_targets <- unique(non_targets_df$gene_name)
              targeting <- filter(targeting_df, cell_line == cell_line_list[i],
                                  gene_name %in% c(non_targets, cell_line_targets))
              targeting$targets_value <- 0
              
              #Set up the targeting matrix for the loss function
              #This section allows the user to fiddle with both the priority of hitting targets vs. non-targets and the assumed 
              #repression achieved from a given miR. Defaults are in the function call.
              targeting$targets_value[
                targeting$housekeeping == "yes" & 
                  targeting$targets_logical ==1] <- fraction_target*offtargetweight
              targeting$targets_value[
                targeting$housekeeping == "no" & 
                  targeting$targets_logical ==1] <- -fraction_target*targetweight
              #Create a new value for weighted delta gibbs
              targeting <- targeting %>% 
                mutate(weighted_target = targets_value*mean_deltaGibbs)
              
              #put together a targeting matrix
              targeting_matrix <- pivot_wider(targeting, id_cols = miR,
                                              names_from = gene_name, 
                                              values_from = weighted_target)
              #need to pair down the miRs_vec --> too long to brute force right now sadly
              #Start by limiting the mir_vec to miRs that target at least one of our target genes.
              targeting_best <- targeting %>% 
                filter(gene_name %in% cell_line_targets, targets_logical == 1) %>% 
                group_by(miR) %>% summarise(n = n()) %>% filter(n>1)
              best_miR_vec <- unique(targeting_best$miR)
              #My god this thing SCALES. Be very careful how many miR combinations you 
              #try to search.
              #top_miRs <- targeting %>% group_by(miR) %>% 
              #  summarise(gibbs_impact = sum(weighted_target)) %>% 
              #  arrange(-desc(gibbs_impact)) %>% slice_head(n = 50)
              #top_miRs <- top_miRs$miR
              #do the combinatorics - this takes a bit
              cocktails <- combinations(length(best_miR_vec), num_miRs, 
                                        best_miR_vec)
              loss_vec <- vector(mode = "numeric", length = nrow(cocktails))
              #ptm <- proc.time()
              for (j in 1:nrow(cocktails)) {
                test_cocktail <- cocktails[j,1:num_miRs]
                loss_vec[j] <- check_combination(test_cocktail)
              }
              #proc.time() - ptm
              cocktails_df <- as_data_frame(cocktails)
              cocktails_df$loss <- loss_vec
              cocktails_df$cell_line <- cell_line_list[i]
              cocktails_df
            }
  #Isolate miRs that target a certain fraction of those mRNA
  #calculate the raw number of mRNAs that must be targeted
  # cocktails_df <- matrix(nrow = num_miRs, ncol = length(cell_line_list))
  # for(i in 1:length(cell_line_list)) {
  #   cell_line_targets <- unlist(targets_df[,cell_line_list[i]])
  #   #
  #   #filter targeting_df for just hits on our target genes
  #   miR_candidates_df <- targeting_df %>% 
  #     filter(gene_name %in% cell_line_targets, targets_logical == 1, 
  #            cell_line == cell_line_list[i])
  #   
  #   #identify miRs with the requisite number of hits and then cut down to 
  #   #generate an optimal cocktail - should I include a tie-breaking rule? 
  #   #Ideally this section would select miRs from different families, 
  #   #when available, to minimize toxicities
  #   miR_candidates <- miR_candidates_df %>% group_by(miR) %>% 
  #     summarise(num_hits = n()) %>% filter(num_hits >= number_targeted) %>% 
  #     arrange(desc(num_hits)) %>% slice(1:num_miRs)
  #   cocktails_df[,i] <- miR_candidates$miR
  # }
  # colnames(cocktails_df) <- cell_line_list
  output <- list(targets_df, cocktail_df)
  return(output)
}

PrepForChordDiagram <- function(targeting_df = overall_targeting_mat.df, 
                                targets = targets_matrix, 
                                cocktail = miR_cocktail_matrix) {
  targeting_df <-  targeting_df %>%
    filter(gene_name %in% targets$gene_name, targets_logical == 1, 
           miR %in% cocktail$miR) %>% 
    group_by(miR, gene_name) %>%
    summarise(targets_logical = max(targets_logical)) %>% 
    pivot_wider(names_from = gene_name, values_from = targets_logical) %>% 
    as.data.frame()
  targeting_df[is.na(targeting_df)] <- 0
  rownames(targeting_df) <- targeting_df$miR
  targeting_df <- targeting_df %>% ungroup() %>% select(-miR) %>%
    as.matrix()
  return(targeting_df)
}

#Found this at https://waterprogramming.wordpress.com/2015/12/02/easy-labels-for-multi-panel-plots-in-r/
put.fig.letter <- function(label, location="topleft", x=NULL, y=NULL, 
                           offset=c(0, 0), size = 1, hue = "black", ...) {
  if(length(label) > 1) {
    warning("length(label) > 1, using label[1]")
  }
  if(is.null(x) | is.null(y)) {
    coords <- switch(location,
                     topleft = c(0.015,0.98),
                     topcenter = c(0.5525,0.98),
                     topright = c(0.985, 0.98),
                     bottomleft = c(0.015, 0.02), 
                     bottomcenter = c(0.5525, 0.02), 
                     bottomright = c(0.985, 0.02),
                     c(0.015, 0.98) )
  } else {
    coords <- c(x,y)
  }
  this.x <- grconvertX(coords[1] + offset[1], from="nfc", to="user")
  this.y <- grconvertY(coords[2] + offset[2], from="nfc", to="user")
  text(labels=label[1], x=this.x, y=this.y, xpd=T, cex = size, col = hue,...)
}
