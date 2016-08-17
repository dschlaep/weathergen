get_ar_models <- function(annual_prcp, years, Wavelet_Decomposition, use_warm_model = TRUE){
  ar_models <- list()
  if (use_warm_model) {
    for (k in 1:dim(Wavelet_Decomposition)[2]) {
      ar_models[[k]] <- auto.arima(Wavelet_Decomposition[,k],max.p=2,max.q=2,max.P=0,max.Q=0,stationary=TRUE)
    }
  } else {
    ar_models[[1]] <- auto.arima(annual_prcp,max.p=2,max.q=2,max.P=0,max.Q=0,stationary=TRUE)
  }
  return(ar_models)
}