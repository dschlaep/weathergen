#' Plot sitespecific annual spells
#'
#' @param obs observation data
#' @param sim simulation output data
#' @param fig_dir  Directory to create PDF containing the plot. Can be NULL to plot on screen
#' @param filename  filename of PDF to create, containing the plot. The .pdf ending is added in the function, not needed in the filename. Default filename is 6_Sites_MedianTrials_AnnualSpells 
#' @param stats, default is c("mean","stdev","skew")
#' @param vars, variables default is c("prcp","tmax","tmin")
#' @param months, month numbers, default is c(1:12)
#' @param spells, default is c("wet","dry")  TODO: The package uses dry, wet, extreme check if we should switch defaults to c('d', 'w', 'e')
#' @export
#' @examples
#' plot_sitespecific_spells_annual (obs, sim)
#' 
#'
plot_sitespecific_spells_annual <- function(obs, sim, fig_dir=NULL, filename="6_Sites_MedianTrials_AnnualSpells", stats=c("mean","stdev","skew"), vars=c("prcp","tmax","tmin"), months=c(1:12),spells=c("wet","dry")){
  if(!is.null(fig_dir)){
    pdf(width=3*length(spells), height=3*length(stats), file=file.path(fig_dir, paste0(filename, ".pdf")))
  }
  layout(matrix(1:(length(spells)*length(stats)), ncol=length(spells), nrow=length(stats), byrow=FALSE))
  
  for(iv in 1:length(spells)){
    for(is in seq_along(stats)){
      y <- apply(sim[1, -1, (iv - 1) * 2 + is, months], MARGIN=1, FUN=mean) #mean across months
      x <- apply(obs[1, -1, (iv - 1) * 2 + is, months], MARGIN=1, FUN=mean)
      xlim <- ylim <- c(min(x, y, na.rm = TRUE), max(x, y, na.rm = TRUE))
      try(plot(x, y, xlim=xlim, ylim=ylim, xlab="Observed spells (days)", ylab="Simulated spells (days)", main=paste(spells[iv], stats[is])))
      abline(0,1)
    }
  }
  
  if(!is.null(fig_dir)){
    dev.off()
  }
}
