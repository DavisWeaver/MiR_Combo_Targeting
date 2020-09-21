cleanBiogrid <- function(expressionDF, gridDF) {
  # Find Unique Shared Genes in Biogrid & Sample
  geneA <-
    as.character(unique(gridDF$Official.Symbol.Interactor.A))
  geneB <-
    as.character(unique(gridDF$Official.Symbol.Interactor.B))
  
  geneList <- unique(c(geneA, geneB))
  
  
  genesWithExpression <-
    as.character(unique(expressionDF$gene_name))
  
  genesInBoth <- sort(intersect(geneList, genesWithExpression))
  
  # Drop Unnecessary Edges
  x <- gridDF$Official.Symbol.Interactor.A %in% genesInBoth
  y <- gridDF$Official.Symbol.Interactor.B %in% genesInBoth
  
  grid <- cbind(x, y)
  grid[grid == FALSE] <- NA
  gridDF <- gridDF[complete.cases(grid), ]
  
  # Construct the Network
  g <-
    graph_from_data_frame(gridDF, genesInBoth, directed = FALSE)
  g <-
    igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE) ##Seems like there is some conflict in the packages we loaded so just gonna specify igraph as the source for simplify
  
  networkGenes <- V(g)$name # same number as python script
  
  expressionDF <-
    expressionDF[expressionDF$gene_name %in% networkGenes, ] # this results in ~2000 more genes than the python script (python = 601920, R = 603684)
  
  return(list(g, expressionDF, genesInBoth, gridDF))
  
}

#Function to add mRNA expression values to pipeline.
addExpression <- function(df, grid = gridDF) {
  #
  wideDF <- df %>% select(-experiment_num, -cell_line, -condition) %>% 
    spread(key = experiment_id, value = expression)
  
  # Output a Graph with Attached Expression Data
  h <-
    graph_from_data_frame(grid, vertices = wideDF, directed = FALSE)
  
  cellLines <- colnames(wideDF)
  
  return(list(h, cellLines, wideDF))
}

getNeighbors <- function(gene, network) {
  neighborGenes <- as.vector(neighbors(network, gene))
  return(unique(neighborGenes))
  
}


calcGibbs <-
  function(gene, expressionVector, neighborlist = neighborList) {
    #shouldn't have named variables in the function that need to be global variables for it to run
    neighborhoodConc <- sum(expressionVector[neighborlist[[gene]]])
    if(!is.na(expressionVector[gene])) {
      if (expressionVector[gene] > 0) {  
        gibbsval <-
          expressionVector[gene] * log(expressionVector[gene] / (neighborhoodConc + expressionVector[gene]))
      } else {
        gibbsval <- 0
      }
    } else {
      gibbsval <- NA
    }
    return(gibbsval)
    
  }

CalcNewGibbs <- function(df, experiment_id_ind, network, common_subgraph = FALSE){
  gibbs_df_exp <- filter(df, experiment_id == experiment_id_ind)
  if (common_subgraph == FALSE) {
    #need to create subgraph
    #subdataframe for each experiment_id
    #Get a list of node names to calculate the new subgraph
    nodes_exp_names <- gibbs_df_exp$gene_name
    nodes_exp <- as.numeric(V(network)[nodes_exp_names]) #have to do this step because the subgraph function requires numeric index
    subgraph_exp <- induced_subgraph(network, nodes_exp)
    
    
    #Need new neighbors for the subgraph 
    neighbors_subgraph <- lapply(1:length(nodes_exp_names),
                                 FUN = getNeighbors, network = subgraph_exp
    )
    names(neighbors_subgraph) <- nodes_exp_names
    #calculate betweenness centrality in here so we don't have to re-calc it a 1000 times w/ just one subgrpah
    bet_cent <- betweenness(subgraph_exp)
  }
  #get expression values for the subgraph -
  expExpression <- get.vertex.attribute(subgraph_exp, experiment_id_ind) # expression for all vertices (genes) in one experiment
  names(expExpression) <- nodes_exp_names
  #need to calculate the "old gibbs" of the subgraph
  oldgibbs_exp <- sum(sapply(nodes_exp_names, calcGibbs,
                             expressionVector = expExpression,
                             neighborlist = neighbors_subgraph))
  
  
  newgibbs_vector <- vector(mode = 'numeric') #initialize gibbs vector
  
  for (i in nodes_exp_names) {
    copy_expExpression <- expExpression
    copy_expExpression[i] <- 0 #gotta make a copy so we don't gradually set the whole network to zero
    newGibbsVal <- sum(sapply(nodes_exp_names, calcGibbs,
                              expressionVector = copy_expExpression,
                              neighborlist = neighbors_subgraph))
    
    newgibbs_vector[i] <- newGibbsVal
  }
  deltaGibbs <- newgibbs_vector - oldgibbs_exp
  #deltaGibbs <- deltaGibbs[deltaGibbs != 0] i THINK HIS MADE IT BREAK ACTUALLY
  
  gibbs_df_exp$deltaGibbs <- deltaGibbs
  return(gibbs_df_exp)
  
}
