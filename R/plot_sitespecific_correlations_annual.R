#' Plot sitespecific annual correlations
#'
#' @param obscors observation data correlations retrieved by get_correlations()
#' @param simcors simulation output correlations retrieved by get_correlations()
#' @param fig_dir  Directory to create PDF containing the plot
#' @param filename filename of PDF to create. PDF Ending will be added, not needed in filename. Default filename is 5_Sites_MedianTrials_AnnualCorrelations.pdf
#' @param vars variables, default is c("prcp","tmax","tmin")
#' @export
#' @examples
#' plot_sitespecific_correlations_annual (obscors, simcors)
#' 
#'
plot_sitespecific_correlations_annual <- function(obscors, simcors, fig_dir=NULL, filename="5_Sites_MedianTrials_AnnualCorrelations", vars=c("prcp","tmax","tmin")){
  if(!is.null(fig_dir)){
    pdf(width=3*length(vars), height=3*2, file=file.path(fig_dir, paste0(filename, ".pdf")))
  }
  layout(matrix(1:(length(vars)*2), nrow=2, ncol=length(vars), byrow=FALSE))
  
  for(iv in seq_along(vars)){
    for(ic in 1:2){
      if(ic == 1){
        tag <- "Intersite correlation"
        y <- as.vector(simcors$intersites[, 2 + iv])
        x <- as.vector(obscors$intersites[, 2 + iv])
      } else {
        tag <- "Lag-1 autocorrelation"
        y <- as.vector(simcors$lag1auto[, iv])
        x <- as.vector(obscors$lag1auto[, iv])
      }
      
      xlim <- ylim <- c(min(x, y, na.rm = TRUE), max(x, y, na.rm = TRUE))
      plot(x, y, xlim=xlim, ylim=ylim, xlab="Observed", ylab="Simulated", main=paste(tag, vars[iv]))
      abline(0,1)
    }
  }
  
  if(!is.null(fig_dir)){
    dev.off()
  }
}