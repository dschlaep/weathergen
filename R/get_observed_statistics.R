
#' Generate Observation Statistics from daily observation data
#' Purpose is usage in plot functions to compare with simulation results
#'
#' @param wgen_out_array, array of outputs created by wgen_daily. Function will use the required observation data which is contained in the output
#' @param vars which data/variables are used for statistics. Default is daily precipitation,  max and min temperature c("prcp","tmax","tmin") 
#' @param timesteps used for the data, default is c("daily","monthly", "annual")
#' @param stats to be calculated over the provided data, defa c("mean","stdev","skew"). Only change if you know what you do.
#' @param spells usually keep on default c("mean","stdev","skew"). Only change if you know what you do.
#' @param return.periods 
#' 
#' @param Month contains the year of the observed data  

#' @export
#' @examples
#' get_observed_statistics(wgen_out_array)
#'  
#'
get_observed_statistics <- function(wgen_out_array, vars=c("prcp","tmax","tmin"), timesteps= c("daily","monthly", "annual"), stats=  c("mean","stdev","skew"),spells=c("wet","dry"), return.periods=c(10,30,100), mo=c(1:12)){
  
  num_site   <- dim(wgen_out_array)[1]
  num_trials <- 1 # dim(wgen_out_array)[2] # no trials because these are observational data
  wyear <-  as.matrix(format(wgen_out_array[[1,1]]$out['DATE'], "%Y"))  #All simulation outputs must have same length??
  
  #Stats = median + realizations, median + sites, statistics, variables, timesteps, months
  Stats <- array(NA, dim=c(1 + num_trials, 1 + num_site, length(stats), length(vars), length(timesteps), length(mo)),
                 dimnames=list(NULL, NULL, stats, vars, timesteps, mo))
  #AnnualMeans = median + realizations, median + sites, variables, years
  AnnualMeans <- array(0, dim=c(1 + num_trials, 1 + num_site, length(vars), dim(unique(wyear))[1]), #WATER_YEAR_SIM
                       dimnames=list(NULL, NULL, vars, sort(unique(wyear)[,1])))
  #AnnualEDV = median + realizations, median + sites, variables, return periods
  AnnualEDV <- array(NA, dim=c(1 + num_trials, 1 + num_site, length(vars), length(return.periods)),
                     dimnames=list(NULL, NULL, vars, return.periods))
  #Spells = median + realizations, median + sites, stats*{dry,wet}, months
  Spells <- array(0, dim=c(1 + num_trials, 1 + num_site, length(spells) * length(stats), length(mo)),
                  dimnames=list(NULL, NULL, paste(rep(spells, each=length(stats)), stats, sep="_"), mo))
  
  
  month <- as.matrix(obs[[1,1]]$out['MONTH'])
  
  
  for (k1 in 1:num_trials) for(iv in seq_along(vars)){
    
    # site_prcp <- wgen_out_array[[1,k1]]$out['PRCP']
    
    #SITE_OBS <- switch(iv, sites_prcp, sites_tmax, sites_tmin)
    SITE_OBS <- switch(iv,  wgen_out_array[[1,k1]]$out['PRCP'],  wgen_out_array[[1,k1]]$out['TMAX'],  wgen_out_array[[1,k1]]$out['TMIN'])
    #Monthly mean daily values
    for (k2 in (seg_along(stats))) {
      if(vars[iv] == "prcp"){
        Stats[1 + k1, -1, k2, iv, 1, mo] <- t(apply(SITE_OBS, c(2), function(x) {y <- which(x>0); aggregate(x[y], by=list(month[y]), FUN=get(stats[k2]))[,2]}))
        #Stats[1 + k1, -1, k2, iv, 1, mo] <- t(apply(SITE_OBS, c(2), function(x) {y <- which(x>0); aggregate(x[y], by=list(month[y]), FUN=mean)[,2]}))
        #Stats[1 + k1, -1, k2, iv, 1, mo] <- t(apply(SITE_OBS, c(2), function(x) {y <- which(x>0); aggregate(x[y], by=list(month[y]), FUN=sd)[,2]}))
        #Stats[1 + k1, -1, k2, iv, 1, mo] <- t(apply(SITE_OBS, c(2), function(x) {y <- which(x>0); aggregate(x[y], by=list(month[y]), FUN=psych::skew)[,2]}))
      } else {
        Stats[1 + k1, -1, k2, iv, 1, mo] <- t(apply(SITE_OBS, c(2), function(x) {y <- which(x>0); aggregate(x[y], by=list(month[y]), FUN=get(stats[k2]))[,2]}))
        # Stats[1 + k1, -1, k2, iv, 1, mo] <- t(apply(SITE_OBS, c(2), function(x) aggregate(x, by=list(month), FUN=mean)[,2]))
        # Stats[1 + k1, -1, k2, iv, 1, mo] <- t(apply(SITE_OBS, c(2), function(x) aggregate(x, by=list(month), FUN=sd)[,2]))
        # Stats[1 + k1, -1, k2, iv, 1, mo] <- t(apply(SITE_OBS, c(2), function(x) aggregate(x, by=list(month), FUN=psych::skew)[,2]))
      }
    }  
   
    # if(vars[iv] == "prcp"){
    #   Stats[1 + k1, -1, 1, iv, 1, mo] <- t(apply(SITE_OBS, c(2), function(x) {y <- which(x>0); aggregate(x[y], by=list(month[y]), FUN=mean)[,2]}))
    #   Stats[1 + k1, -1, 2, iv, 1, mo] <- t(apply(SITE_OBS, c(2), function(x) {y <- which(x>0); aggregate(x[y], by=list(month[y]), FUN=sd)[,2]}))
    #   Stats[1 + k1, -1, 3, iv, 1, mo] <- t(apply(SITE_OBS, c(2), function(x) {y <- which(x>0); aggregate(x[y], by=list(month[y]), FUN=psych::skew)[,2]}))
    # } else {
    #   Stats[1 + k1, -1, 1, iv, 1, mo] <- t(apply(SITE_OBS, c(2), function(x) aggregate(x, by=list(month), FUN=mean)[,2]))
    #   Stats[1 + k1, -1, 2, iv, 1, mo] <- t(apply(SITE_OBS, c(2), function(x) aggregate(x, by=list(month), FUN=sd)[,2]))
    #   Stats[1 + k1, -1, 3, iv, 1, mo] <- t(apply(SITE_OBS, c(2), function(x) aggregate(x, by=list(month), FUN=psych::skew)[,2]))
    # }
    # 
    #Annual mean daily values
    
    AnnualMeans[1 + k1, -1, iv, ] <- t(apply(SITE_OBS, c(2), function(x) {aggregate(x, by=list(wyear), FUN=mean)[,2]})) #years, sites
    
    # print(apply(AnnualMeans[1 + k1, -1, iv, ], 1, FUN=mean))
    for (k2 in (seg_along(stats))) {
      Stats[1 + k1, -1, 1, iv, 3, 1] <- apply(AnnualMeans[1 + k1, -1, iv, ], 1, FUN=get(stats[k2]))
    }
    # Stats[1 + k1, -1, 1, iv, 3, 1] <- apply(AnnualMeans[1 + k1, -1, iv, ], 1, FUN=mean)
    # Stats[1 + k1, -1, 2, iv, 3, 1] <- apply(AnnualMeans[1 + k1, -1, iv, ], 1, FUN=sd)
    # Stats[1 + k1, -1, 3, iv, 3, 1] <- apply(AnnualMeans[1 + k1, -1, iv, ], 1, FUN=psych::skew)
    
    
    #Annual extreme daily values  TODO: -min for tmin, what if only values > 0 are observed? 
    ext <- t(apply(SITE_OBS * ifelse(vars[iv] == "tmin", -1, 1), c(2), function(x) {aggregate(x, by=list(wyear), FUN=max)[,2]})) #dim(ext) = sites, years; find max for prcp and tmax, but -min for tmin
    AnnualEDV[1 + k1, -1, iv, ] <- t(apply(ext, 1, FUN=function(x) return.level(fevd(x, type="GEV", method="MLE", time.units="year"), return.period=return.periods)))
    AnnualEDV[1 + k1, -1, iv, ] <- ifelse(AnnualEDV[1 + k1, -1, iv, ] > 2 * max(ext), NA, AnnualEDV[1 + k1, -1, iv, ]) #set outliers to NA
    if(vars[iv] == "tmin") AnnualEDV[1 + k1, -1, iv, ] <- -AnnualEDV[1 + k1, -1, iv, ] #set -min to min for tmin
    
    #PRCP mean monthly wet/dry spells
    for(is in 1:num_site){
      # Spells[1 + k1, 1 + is, , ] <- t(aggregate(sites_prcp[, is], by=list(MONTH_D), FUN=function(x) {
      Spells[1 + k1, 1 + is, , ] <- t(aggregate(wgen_out_array[[is,k1]]$out['PRCP'], by=list(month), FUN=function(x) {
        temp <- rle(x > 0)
        res <- sapply(list(temp$lengths[temp$values], temp$lengths[!temp$values]), FUN=function(x) c(mean(x), sd(x), psych::skew(x)))
        return(res)})[, 2])
    }
  }
  AnnualMeans[, , vars == "prcp", ] <- 365 * AnnualMeans[, , vars == "prcp", ]
  
  #Medians across trials
  Stats[1, , , , , ] <- apply(Stats[-1, , , , , , drop=FALSE], MARGIN=c(2, 3, 4, 5, 6), FUN=median, na.rm=TRUE)
  AnnualMeans[1, , , ] <- apply(AnnualMeans[-1, , , , drop=FALSE], MARGIN=c(2, 3, 4), FUN=median, na.rm=TRUE)
  AnnualEDV[1, , , ] <- apply(AnnualEDV[-1, , , , drop=FALSE], MARGIN=c(2, 3, 4), FUN=median, na.rm=TRUE)
  Spells[1, , , ] <- apply(Spells[-1, , , , drop=FALSE], MARGIN=c(2, 3, 4), FUN=median, na.rm=TRUE)
  
  #Basinwide median across sites
  Stats[, 1, , , , ] <- apply(Stats[, -1, , , , , drop=FALSE], MARGIN=c(1, 3, 4, 5, 6), FUN=median, na.rm=TRUE)
  AnnualMeans[, 1, , ] <- apply(AnnualMeans[, -1, , , drop=FALSE], MARGIN=c(1, 3, 4), FUN=median, na.rm=TRUE)
  AnnualEDV[, 1, , ] <- apply(AnnualEDV[, -1, , , drop=FALSE], MARGIN=c(1, 3, 4), FUN=median, na.rm=TRUE)
  Spells[, 1, , ] <- apply(Spells[, -1, , , drop=FALSE], MARGIN=c(1, 3, 4), FUN=median, na.rm=TRUE)
  
  return(list(Stats=Stats, AnnualMeans=AnnualMeans, AnnualEDV=AnnualEDV, Spells=Spells))
}
