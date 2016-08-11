
#' get_correlations calculate correlations for plot function 5
#'
#' @param AnnualMeans This value can be retrieved from getObservedStatistics(...)$AnnualMeans or getSimulatedStatistics(..)$AnnualMeans
#' @param vars
#' @export
#' @examples
#' obsstats$cors <- get_correlations(AnnualMeans=obssStats$AnnualMeans[-1, -1, , , drop=FALSE])
#' 
#'
getCorrelations <- function(AnnualMeans, vars= c("prcp","tmax","tmin")){
  num_trials <- dim(AnnualMeans)[1]
  
  #Intersites correlations of annual series
  intersites <- t(combn(dim(AnnualMeans)[2], 2))
  cors <- array(NA, c(1 + num_trials, nrow(intersites), length(vars)), dimnames=list(NULL, NULL, vars))
  
  for(k1 in 1:num_trials){
    for(is in 1:nrow(intersites)){
      for(iv in 1:length(vars)){
        cors[1 + k1, is, iv] <- cor(AnnualMeans[k1, intersites[is, 1], iv, ], AnnualMeans[k1, intersites[is, 2], iv, ], method="spearman")
      }
    }
  }
  cors[1, , ] <- apply(cors[-1, , , drop=FALSE], MARGIN=c(2, 3), FUN=median) #Median across trials
  intersites2 <- matrix(NA, nrow = nrow(intersites), ncol = 2 + length(vars), dimnames = list(NULL, c("SitePair1", "SitePair2", vars)))
  intersites2[, 1:2] <- intersites
  intersites2[, vars] <- cors[1, , , drop=TRUE]
  #intersites <- cbind(intersites, cors[1, , , drop=TRUE])
  
  #Lag-1 autocorrelation of annual series
  lag1auto <- array(NA, c(1 + num_trials, dim(AnnualMeans)[2], length(vars)), dimnames=list(NULL, NULL, vars))
  for(k1 in 1:num_trials){
    for(is in 1:dim(AnnualMeans)[2]){
      for(iv in 1:length(vars)){
        lag1auto[1 + k1, is, iv] <- acf(AnnualMeans[k1, is, iv, ], type="correlation", plot=FALSE)[1, ]$acf
        stopifnot(!is.na(lag1auto[1 + k1, is, iv]))
      }
    }
  }
  lag1auto[1, , ] <- apply(lag1auto[-1, , , drop=FALSE], MARGIN=c(2, 3), FUN=median) #Median across trials
  
  return(list(intersites=intersites2, lag1auto=lag1auto[1, , , drop=TRUE]))
}