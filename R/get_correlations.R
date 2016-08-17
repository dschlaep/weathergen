
#' get_correlations calculate correlations for plot function 5
#'
#' @param annual_means This value can be retrieved from get_observed_statistics(...)$AnnualMeans or get_simulated_statistics(..)$AnnualMeans
#' @param vars variables for whic a calculation is done. Default is c("prcp","tmax","tmin","wind")
#' @export
#' @examples
#' obsstats$cors <- get_correlations(annual_means=obssStats$AnnualMeans[-1, -1, , , drop=FALSE])
#' 
#'
get_correlations <- function(annual_means, vars= c("prcp","tmax","tmin","wind")){
  num_trials <- dim(annual_means)[1]
  
  #Intersites correlations of annual series
  intersites <- t(combn(dim(annual_means)[2], 2))
  cors <- array(NA, c(1 + num_trials, nrow(intersites), length(vars)), dimnames=list(NULL, NULL, vars))
  
  for(k1 in 1:num_trials){
    for(is in 1:nrow(intersites)){
      for(iv in 1:length(vars)){
        cors[1 + k1, is, iv] <- cor(annual_means[k1, intersites[is, 1], iv, ], annual_means[k1, intersites[is, 2], iv, ], method="spearman")
      }
    }
  }
  cors[1, , ] <- apply(cors[-1, , , drop=FALSE], MARGIN=c(2, 3), FUN=median) #Median across trials
  intersites2 <- matrix(NA, nrow = nrow(intersites), ncol = 2 + length(vars), dimnames = list(NULL, c("SitePair1", "SitePair2", vars)))
  intersites2[, 1:2] <- intersites
  intersites2[, vars] <- cors[1, , , drop=TRUE]
  #intersites <- cbind(intersites, cors[1, , , drop=TRUE])
  
  #Lag-1 autocorrelation of annual series
  lag1auto <- array(NA, c(1 + num_trials, dim(annual_means)[2], length(vars)), dimnames=list(NULL, NULL, vars))
  for(k1 in 1:num_trials){
    for(is in 1:dim(annual_means)[2]){
      for(iv in 1:length(vars)){
        lag1auto[1 + k1, is, iv] <- acf(annual_means[k1, is, iv, ], type="correlation", plot=FALSE)[1, ]$acf
        stopifnot(!is.na(lag1auto[1 + k1, is, iv]))
      }
    }
  }
  lag1auto[1, , ] <- apply(lag1auto[-1, , , drop=FALSE], MARGIN=c(2, 3), FUN=median) #Median across trials
  
  return(list(intersites=intersites2, lag1auto=lag1auto[1, , , drop=TRUE]))
}