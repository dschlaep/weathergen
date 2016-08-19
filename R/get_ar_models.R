#' Get ARIMA Models
#'
#' @param annual_prcp yearly prcp 
#' @param years  the years of the annual prcp data
#' @param wavelet_decomposition, can be retrieved using wavelet_components() 
#' @param use_warm_model boolean if warm modelis used
#' @export
#' @examples
#' ar_models <- get_ar_models( mysimStats$AnnualMeans[1,1,1,], names(mysimStats$AnnualMeans[1,1,1,]), wavelet_components(mysimStats$AnnualMeans[1,1,1,],wavelet_analysis(mysimStats$AnnualMeans[1,1,1,], names(mysimStats$AnnualMeans[1,1,1,]),sig.level=0.9, "white"), n.periods=1, sig.periods=c(8,9,10,11,12,13,14,15), n.comp.periods=8))
#' 
#' 
#'
get_ar_models <- function(annual_prcp, years, wavelet_decomposition, use_warm_model = TRUE){
  ar_models <- list()
  if (use_warm_model) {
    for (k in 1:dim(wavelet_decomposition)[2]) {
      ar_models[[k]] <- auto.arima(wavelet_decomposition[,k],max.p=2,max.q=2,max.P=0,max.Q=0,stationary=TRUE)
    }
  } else {
    ar_models[[1]] <- auto.arima(annual_prcp,max.p=2,max.q=2,max.P=0,max.Q=0,stationary=TRUE)
  }
  return(ar_models)
}