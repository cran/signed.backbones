utils::globalVariables(c("i", "j"))
extract = function(edgelist, directed = TRUE, significance_threshold = '15pc', vigor_threshold = 0.1, 
                   return_weights = FALSE, return_significance = FALSE, 
                   max_iteration = 100, precision = 10e-8){
  
  edgelist[,c(1, 2)] <- lapply(edgelist[,c(1, 2)], as.character)
  
  if (directed) {
    edgelist <- edgelist
  } else {
    edgelist <- rbind(edgelist, setNames(edgelist[,c(2,1,3)], names(edgelist)))
  }
  
  
  edgedf <- calculateStatistics(edgelist, max_iteration, precision)
  
  if (length(significance_threshold) == 2){
    if (is.character(significance_threshold[1]) & is.character(significance_threshold[2])){
      if (all(substr(significance_threshold, nchar(significance_threshold) - 2 + 1, nchar(significance_threshold)) == 'pc')){
        lpc <- as.numeric(substr(significance_threshold[1], 1, nchar(significance_threshold[1])-2)[[1]])
        upc <- as.numeric(substr(significance_threshold[2], 1, nchar(significance_threshold[2])-2)[[1]])
        if (lpc < 0 | lpc > 100 | upc < 0 | upc > 100){
          stop("If significance_threshold is a vector and if filtering is based on ranking, both elements must be within ['0pc','100pc'] range.")         
        } else {
          lb <- quantile(edgedf$StdDist,lpc/100)[[1]]
          ub <- quantile(edgedf$StdDist,1-upc/100)[[1]]
          significance_filter_index <- (edgedf$StdDist <= lb | edgedf$StdDist >= ub)
        }
        
      } else {
        stop("If significance_threshold is a vector and if filtering is based on ranking, both elements must end with 'pc', for instance, ('15pc', '20pc').")
      }
    } else{
      if (significance_threshold[1] > 0) {
        stop('If significance_threshold is a vector and if filtering is not based on ranking, the first element should be non-positive.')
      }
      if (significance_threshold[2] < 0) {
        stop('If significance_threshold is a vector and if filtering is not based on ranking, the second element should be non-negative.')
      }
      significance_filter_index <- (edgedf$StdDist <= significance_threshold[1] | edgedf$StdDist >= significance_threshold[2])
    }
  } else if (length(significance_threshold) == 1) {
    if (is.character(significance_threshold)) {
      if (substr(significance_threshold, nchar(significance_threshold) - 2 + 1, nchar(significance_threshold)) == 'pc') {
        pc <- as.numeric(substr(significance_threshold[1], 1, nchar(significance_threshold[1])-2)[[1]])
        if (pc < 0 | pc > 100){
          stop("If filtering is based on ranking, the value must be within ['0pc','100pc'] range.")
        } else {
          b <- quantile(abs(edgedf$StdDist), 1-pc/100)[[1]]
          significance_filter_index <- abs(edgedf$StdDist) >= b
        }
        
      } else {
        stop("If filtering is based on ranking, the value must end with 'pc', for instance, '15pc'.")
      }
    } else {
      significance_filter_index <- abs(edgedf$StdDist) >= significance_threshold
    }
  } else {
    stop('significance_threshold should be a scalar or a vector of length 2.')
  }
  
  
  if (length(vigor_threshold) == 2){
    if (vigor_threshold[1] > 0) {
      stop('If vigor_threshold is a vector, the first element should be non-positive.')
    }
    if (vigor_threshold[2] < 0) {
      stop('If vigor_threshold is a vector, the second element should be non-negative.')	
    }
    
    vigor_filter_index = (edgedf$LiftScore <= vigor_threshold[1] | edgedf$LiftScore >= vigor_threshold[2])
  } else if (length(vigor_threshold) == 1) {
    vigor_filter_index = (abs(edgedf$LiftScore) >= vigor_threshold)
  } else {
    stop('vigor_threshold should be a scalar or a vector of length 2.')
  }
  
  
  filter_index <- (vigor_filter_index & significance_filter_index)
  edges <- edgedf[filter_index,]

  if (!directed){
    temp_min = pmin(edges$i, edges$j)
    edges$j = pmax(edges$i, edges$j)
    edges$i = temp_min
    
    edges$absLiftScore <- abs(edges$LiftScore)
    edges <- edges[order(-edges$absLiftScore, edges$i, edges$j, decreasing = FALSE ),]
    edges <- subset(edges, !duplicated(subset(edges, select=c(i, j))))
  }
  
  if (return_weights){
    if (return_significance) {
      edges <- edges[, c('i', 'j', 'LiftScore', 'StdDist')]
    } else {
      edges <- edges[, c('i', 'j', 'LiftScore')]
    }
  } else {
    if (return_significance) {
      edges <- edges[, c('i', 'j', 'Sign', 'StdDist')]
    } else {
      edges <- edges[, c('i', 'j', 'Sign')]
    }
  }
  
  
  edges <- edges[order(edges$i, edges$j),]
  rownames(edges) <- NULL
  
  message(nrow(edges), 'edges are retained.')
  
  return(edges)
}


calculateStatistics = function(edgelist, max_iteration, precision){
  
  
  edgedf <- edgelist
  colnames(edgedf) <- c('i', 'j', 'Wij')
  
  n_self_loops <- sum(edgedf['i'] == edgedf['j'])
  
  if (n_self_loops > 0){
    edgedf <- edgedf[edgedf['i'] != edgedf['j'], ]
    message(n_self_loops, "self loops are removed.\n")
  }
  
  if (any(duplicated(edgedf[,c('i', 'j')]))){
    edgedf <- aggregate(.~i+j, edgedf ,FUN = sum)
    message("Duplicate edges are identified. The edges are merged by summing up their weights.\n")
  }
  
  nodes <- union(unique(edgedf[['i']]), unique(edgedf[['j']]))
  adj <- data.frame(matrix(nrow = length(nodes), ncol = length(nodes)), row.names = nodes)
  colnames(adj) <- nodes
  
  for (row in 1:nrow(edgedf)) {
    adj[edgedf[row, 'i'], edgedf[row, 'j']] <- edgedf[row, 'Wij']
  }
  
  adj[is.na(adj)] <- 0
  edgedf <- reshape2::melt(as.matrix(adj))
  edgedf[,c(1, 2)] <- lapply(edgedf[,c(1, 2)], as.character)
  colnames(edgedf) <- c('i', 'j', 'Wij')
  
  edgedf <- edgedf[edgedf['i'] != edgedf['j'], ]
  row.names(edgedf) <- NULL
  
  Wi. <- aggregate(.~i, edgedf[c('i', 'Wij')], sum)
  colnames(Wi.) <- c('i', 'Wi.')
  W.j <- aggregate(.~j, edgedf[c('j', 'Wij')], sum)
  colnames(W.j) <- c('j', 'W.j')
  W.i <- aggregate(.~j, edgedf[c('j', 'Wij')], sum)
  colnames(W.i) <- c('i', 'W.i')
  W.. <- sum(edgedf['Wij'])
  
  edgedf <- merge(edgedf, W.j, by.x = 'j', by.y = 'j', all.x = TRUE) 
  edgedf <- merge(edgedf, Wi., by.x = 'i', by.y = 'i', all.x = TRUE)
  edgedf <- merge(edgedf, W.i, by.x = 'i', by.y = 'i', all.x = TRUE)
  edgedf['W..'] <- W..
  
  edgedf['Pij'] <- (edgedf['Wi.'] * edgedf['W.j']) / (edgedf['W..'] - edgedf['W.i'])
  
  pri <- data.frame(matrix(nrow = length(nodes), ncol = length(nodes)), row.names = nodes)
  colnames(pri) <- nodes
  
  for (row in 1:nrow(edgedf)) {
    pri[edgedf[row, 'i'], edgedf[row, 'j']] <- edgedf[row, 'Pij']
  }
  
  pri[is.na(pri)] <- 0
  pri[pri==0] <- 0.00001

  nullmat <- getNullMatrix(adj, pri, max_iteration, precision)
  
  
  rownames(nullmat) <- rownames(adj)
  colnames(nullmat) <- colnames(adj)
  Nij <- reshape2::melt(nullmat)
  colnames(Nij) <- c('i', 'j', 'Nij')
  
  edgedf <- merge(edgedf, Nij, by.x = c('i', 'j'), by.y = c('i', 'j'), all.x = TRUE) 
  edgedf <- within(edgedf, rm('Pij'))
  
  N <- edgedf[['W..']] - edgedf[['W.i']]
  K <- ifelse(edgedf[['W.j']] == 0, precision, edgedf[['W.j']])
  n <- ifelse(edgedf[['Wi.']] == 0, precision, edgedf[['Wi.']])
  
  edgedf['Var'] <- n * (K/N) * ((N-K)/N) * ((N-n)/(N-1))
  edgedf['Std'] <- edgedf[['Var']] ** .5
  edgedf['StdDist'] <- (edgedf[['Wij']] - edgedf[['Nij']])/edgedf[['Std']]
  edgedf['Lift'] <- edgedf[['Wij']] / ifelse(edgedf[['Nij']]== 0, 1, edgedf[['Nij']])
  edgedf['LiftScore'] <- (edgedf[['Lift']] - 1) / (edgedf[['Lift']] + 1)
  edgedf['Sign'] = sign(edgedf[['LiftScore']])
  edgedf = edgedf[edgedf[['i']] != edgedf[['j']], ]

  return(edgedf)
}




getNullMatrix <- function(wmatrix, pmatrix, max_iteration, precision){
  
  n <- nrow(wmatrix)
  
  marginal_row    <- rowSums(wmatrix)
  marginal_column <- colSums(wmatrix)
  
  marginal_row[marginal_row == 0]       <- precision**2
  marginal_column[marginal_column == 0] <- precision**2
  
  prior <- pmatrix
  null <- cbind(prior)
  null[null==0] <- precision
  
  for (i in 1:n){
    null[i,i] <- 0
  }
  
  null <- as.matrix(null)
  
  iteration <- 0
  while (iteration < max_iteration){
    
    row_scaler    <- marginal_row    / rowSums(null)
    null <- diag(row_scaler) %*% null
    
    column_scaler <- marginal_column / colSums(null)
    null <- as.matrix(null) %*% diag(column_scaler) 
    
    iteration = iteration + 1
    
    MAE_row    <- mean(abs(marginal_row    - rowSums(null)))
    MAE_column <- mean(abs(marginal_column - colSums(null)))
    
    if ( MAE_row < precision & MAE_column < precision){
      message(sep='', 'Iterative fitting procedure converged at iteration ', iteration, '.\n')
      return(null)
    }
  }
  
  message(sep='', 'Iterative Fitting Procedure ended at iteration ', iteration,' with Mean Absolute Error (MAE) ', MAE_row, ' for row totals and MAE ', MAE_column,' for column totals.\n')
  return(null)
}




