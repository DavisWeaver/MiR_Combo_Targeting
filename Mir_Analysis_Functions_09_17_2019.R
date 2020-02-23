
calc_target_matrix <- function(){
  #convert these to entrez
  library(org.Hs.eg.db)
  library(reshape2)
  setwd("C:/Users/dtw43/Documents/MIR_Combo_Targeting/data files")
  
  
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
  
  final_df <- read_csv('C:/Users/dtw43/Documents/MIR_Combo_Targeting/data files/EwingsDatasetFinal.csv')  
  gene_list = unique(final_df$Gene_Name)
  
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

add_cosmic_vars <- function(df){#Lets download COSMIC's list of genes associated in cancer
  cosmic_df <- read_csv('C:/Users/dtw43/Documents/MIR_Combo_Targeting/data files/Census_allMon_2019.csv')
  oncogene_df <- cosmic_df %>% filter(grepl('oncogene', get('Role in Cancer')))
  oncogene_list <- oncogene_df$`Gene Symbol`
  tumorsuppressor_df <- cosmic_df %>% filter(grepl('TSG', get('Role in Cancer')))
  TSG_list <- tumorsuppressor_df$'Gene Symbol'
  
  #create a variable for whether or not the identified genes in our dataset are oncogenes
  yesoncogene_df <- filter(df, Gene_Name %in% oncogene_list)
  nooncogene_df <- filter(df, !(Gene_Name %in% oncogene_list))
  yesoncogene_df$oncogene <- "yes"
  nooncogene_df$oncogene <- "no"
  df <- rbind(yesoncogene_df, nooncogene_df)
  
  #Create variable for whether or not identified genes in our dataset are tumor suppressor genes
  yesTSG_df <- filter(df, Gene_Name %in% TSG_list)
  noTSG_df <- filter(df, !(Gene_Name %in% TSG_list))
  yesTSG_df$TSG <- "yes"
  noTSG_df$TSG <- "no"
  df <- rbind(yesTSG_df, noTSG_df)
  
  #create variable for whether or not identified genes are associated with cancer at all
  yescancer_df <- filter(df, Gene_Name %in% cosmic_df$'Gene Symbol')
  nocancer_df <- filter(df, !(Gene_Name %in% cosmic_df$'Gene Symbol'))
  yescancer_df$cancer_associated <- "yes"
  nocancer_df$cancer_associated <- "no"
  df <- rbind(yescancer_df, nocancer_df)
  return(df)
}

#Add entrez_id to final_df
#Add column for entrez id
## Bimap interface:
add_entrez_id <- function(df) {
  require(org.Hs.eg.db)
  require(dplyr)
  gene_list <- unique(df$Gene_Name)
  x <- org.Hs.egSYMBOL
  # Get the gene symbol that are mapped to an entrez gene identifiers
  mapped_genes <- mappedkeys(x)
  # Convert to a list
  converted_symbols <- as.list(x[mapped_genes])
  converted_symbols_map <- melt(converted_symbols)
  rm(list = c("converted_symbols", "x", "mapped_genes")) #Keep memory free
  colnames(converted_symbols_map) <- c("Gene_Name", "entrez_id")
  #convert factor to numeric
  converted_symbols_map <- converted_symbols_map %>% 
    mutate(entrez_id = as.numeric(entrez_id), 
           Gene_Name = as.character(Gene_Name)) %>% 
    filter(Gene_Name %in% gene_list)
  
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
    df <- df %>% group_by(Gene_Name) %>% 
      summarise(
        mean_gibbs_diff = mean(gibbs_diff), 
        se_gibbs_diff = sd(gibbs_diff) / sqrt(n()),
        mean_expression = mean(expression),
        mean_gibbs = mean(gibbs),
        oncogene = unique(oncogene),
        TSG = unique(TSG),
        cancer_associated = unique(cancer_associated)
      ) %>% arrange(desc(mean_gibbs_diff))
  } else {
    df <- df %>% group_by(cell_line, Gene_Name) %>% 
      summarise(
        mean_gibbs_diff = mean(gibbs_diff), 
        se_gibbs_diff = sd(gibbs_diff) / sqrt(n()),
        mean_expression = mean(expression),
        mean_gibbs = mean(gibbs),
        oncogene = unique(oncogene),
        TSG = unique(TSG),
        cancer_associated = unique(cancer_associated)
      ) %>% arrange(desc(mean_gibbs_diff))
  }
  if(justcancer == TRUE) {
    df <- df %>% filter(cancer_associated == "yes")
  }
  if(topx == TRUE){
    genes_target <- unique(df$Gene_Name)[1:numgenes]
    df <- df %>% filter(Gene_Name %in% genes_target)
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
  colnames(converted_symbols_map) <- c("Gene_Name", "entrez_id")
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
    summarise(mean_gibbs_diff = mean(gibbs_diff), 
              sd_gibbs_diff = sd(gibbs_diff),  
              cancer_associated = unique(cancer_associated)) 
  
  #add weighted gibbs_df as well
  targeting_df <- targeting_df %>% left_join(df_mean) %>% 
    mutate(weighted_gibbs_diff = targets_logical * mean_gibbs_diff)
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
                  gibbs_target = sum(weighted_gibbs_diff)) %>% 
        filter(num_target > 0) %>% arrange(desc(gibbs_target)) %>% slice(0:num_miRs)
    } else {
      summary_tbl <- targeting_df_cellline %>% 
        group_by(miR) %>% summarise(num_target = sum(targets_logical), 
                                    gibbs_target = sum(weighted_gibbs_diff)) %>% 
        filter(num_target > 0) %>% arrange(desc(gibbs_target)) %>% slice(0:num_miRs)
    }
    miRs_df[,i] <- summary_tbl$miR
  }
  colnames(miRs_df) <- cell_line_list
  return(miRs_df)
}
#The purpose of this function is to identify a cocktail of miRs that all 
#preferentially downregulate several genes of interest while exhibiting limited 
#shared targets outside of the genes of interest -- this may not be possible
Generate_MiR_Cocktail <- function(targeting_df = overall_targeting_mat.df,
                                  df = final_df,
                                  num_miRs = 5, num_targets = 5, 
                                  fraction_target = 0.6, cancer_only = FALSE,
                                  keep_cell_line = TRUE){ 
  
  #calculate number targeted at the top 
  number_targeted <- round(num_targets*fraction_target)
  #Generate a matrix of the top mRNA by gibb's free energy for each cell line
  if(keep_cell_line == TRUE){
    cell_line_list <- unique(df$cell_line)
    targets_df <- matrix(nrow = num_targets, ncol = length(cell_line_list))
    for(i in 1:length(cell_line_list)) {
      cell_line_df <- filter(df, cell_line == cell_line_list[i])
      if(cancer_only == TRUE) {
        cell_line_df <- cell_line_df %>% group_by(Gene_Name) %>%
          filter(cancer_associated == "yes") %>%
          summarise(mean_gibbs_diff = mean(gibbs_diff)) %>% 
          arrange(desc(mean_gibbs_diff)) %>% slice(1:num_targets)
      } else{
        cell_line_df <- cell_line_df %>% group_by(Gene_Name) %>%
          summarise(mean_gibbs_diff = mean(gibbs_diff)) %>% 
          arrange(desc(mean_gibbs_diff)) %>% slice(1:num_targets)
      }
      targets_df[,i] <- cell_line_df$Gene_Name
    }
    colnames(targets_df) <- cell_line_list
    targets_df <- as_tibble(targets_df)
    
    #Isolate miRs that target a certain fraction of those mRNA
    #calculate the raw number of mRNAs that must be targeted
    cocktails_df <- matrix(nrow = num_miRs, ncol = length(cell_line_list))
    for(i in 1:length(cell_line_list)) {
      cell_line_targets <- unlist(targets_df[,cell_line_list[i]])
      
      #filter targeting_df for just hits on our target genes
      miR_candidates_df <- targeting_df %>% 
        filter(Gene_Name %in% cell_line_targets, targets_logical == 1, 
               cell_line == cell_line_list[i])
      
      #identify miRs with the requisite number of hits and then cut down to 
      #generate an optimal cocktail - should I include a tie-breaking rule? 
      #Ideally this section would select miRs from different families, 
      #when available, to minimize toxicities
      miR_candidates <- miR_candidates_df %>% group_by(miR) %>% 
        summarise(num_hits = n()) %>% filter(num_hits >= number_targeted) %>% 
        arrange(desc(num_hits)) %>% slice(1:num_miRs)
      cocktails_df[,i] <- miR_candidates$miR
    }
    colnames(cocktails_df) <- cell_line_list
  } else {
    if(cancer_only == TRUE) {
      df <- df %>% group_by(Gene_Name) %>%
        filter(cancer_associated == "yes") %>%
        summarise(mean_gibbs_diff = mean(gibbs_diff)) %>% 
        arrange(desc(mean_gibbs_diff)) %>% slice(1:num_targets)
    } else {
      df <- df %>% group_by(Gene_Name) %>%
        summarise(mean_gibbs_diff = mean(gibbs_diff)) %>% 
        arrange(desc(mean_gibbs_diff)) %>% slice(1:num_targets)
    }
    targets <- df$Gene_Name
    targets_df <- df
    
    
    #filter targeting_df for just hits on our target genes
    miR_candidates_df <- targeting_df %>%
      filter(Gene_Name %in% targets, targets_logical == 1) %>% 
      group_by(miR, Gene_Name) %>%
      summarise(targets_logical = max(targets_logical),
                mean_gibbs_diff = mean(mean_gibbs_diff))
    
     
     
    
    
    #need to get rid of 
    #identify miRs with the requisite number of hits and then cut down to 
    #generate an optimal cocktail - should I include a tie-breaking rule? 
    #Ideally this section would select miRs from different families, 
    #when available, to minimize toxicities
    miR_candidates <- miR_candidates_df %>% group_by(miR) %>% 
      summarise(num_hits = n()) %>% filter(num_hits >= number_targeted) %>% 
      arrange(desc(num_hits)) %>% slice(1:num_miRs)
    cocktails_df <- miR_candidates
  }
  output <- list(targets_df, cocktails_df)
  return(output)
}

PrepForChordDiagram <- function(targeting_df = overall_targeting_mat.df, 
                                targets = targets_matrix, 
                                cocktail = miR_cocktail_matrix) {
  targeting_df <-  targeting_df %>%
    filter(Gene_Name %in% targets$Gene_Name, targets_logical == 1, 
           miR %in% cocktail$miR) %>% 
    group_by(miR, Gene_Name) %>%
    summarise(targets_logical = max(targets_logical))
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
