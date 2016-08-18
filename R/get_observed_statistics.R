
#' Generate Observation Statistics from daily observation data
#' Purpose is usage in plot functions to compare with simulation results
#'
#' @param wgen_out_array, array of outputs created by wgen_daily. Function will use the required observation data which is contained in the output
#' @param vars which data/variables are used for statistics. Default is daily precipitation,  max and min temperature and wind c("prcp","tmax","tmin", "wind") 
#' @param spells default is c("wet","dry"). TODO: Package also uses Extreme, will need implementation, currently ignored but won't lead to error
#' @param return.periods, default is c(10,30,100)
#' @export
#' @examples
#' get_observed_statistics(wgen_out_array)
#'  
#'
get_observed_statistics <- function(wgen_out_array, vars=c("prcp","tmax","tmin","wind"), spells=c("wet","dry"), return.periods=c(10,30,100)){
  
  timesteps <- c("daily","monthly", "annual")
  stats <-  c("mean","stdev","skew")
  months <- c(1:12)
  num_site   <- dim(wgen_out_array)[1]
  num_trials <- 1 # dim(wgen_out_array)[2] # no trials because these are observational data
  wyear <-  as.matrix(format(wgen_out_array[[1,1]]$obs['DATE'], "%Y"))  #All simulation outputs must have same length??
  
  #Stats = median + realizations, median + sites, statistics, variables, timesteps, months
  Stats <- array(NA, dim=c(1 + num_trials, 1 + num_site, length(stats), length(vars), length(timesteps), length(months)),
                 dimnames=list(NULL, NULL, stats, vars, timesteps, months))
  #AnnualMeans = median + realizations, median + sites, variables, years
  AnnualMeans <- array(0, dim=c(1 + num_trials, 1 + num_site, length(vars), dim(unique(wyear))[1]), #WATER_YEAR_SIM
                       dimnames=list(NULL, NULL, vars, sort(unique(wyear)[,1])))
  #AnnualEDV = median + realizations, median + sites, variables, return periods
  AnnualEDV <- array(NA, dim=c(1 + num_trials, 1 + num_site, length(vars), length(return.periods)),
                     dimnames=list(NULL, NULL, vars, return.periods))
  #Spells = median + realizations, median + sites, stats*{dry,wet}, months
  Spells <- array(0, dim=c(1 + num_trials, 1 + num_site, length(spells) * length(stats), length(months)),
                  dimnames=list(NULL, NULL, paste(rep(spells, each=length(stats)), stats, sep="_"), months))
  
  
  month <- as.matrix(wgen_out_array[[1,1]]$obs['MONTH'])
  
  
  for (k1 in 1:num_trials) for(iv in seq_along(vars)){
    SITE_OBS <- sapply(wgen_out_array[, 1], simplify = "array", FUN=function(x) unlist(x$obs[toupper(vars[iv])] ))
    
    #Monthly mean daily values
    if(toupper(vars[iv]) == "PRCP"){
      Stats[1 + k1, -1, 1, iv, 1, months] <- t(apply(SITE_OBS, c(2), function(x) {y <- which(x>0); aggregate(x[y], by=list(month[y]), FUN=mean)[,2]}))
      Stats[1 + k1, -1, 2, iv, 1, months] <- t(apply(SITE_OBS, c(2), function(x) {y <- which(x>0); aggregate(x[y], by=list(month[y]), FUN=sd)[,2]}))
      Stats[1 + k1, -1, 3, iv, 1, months] <- t(apply(SITE_OBS, c(2), function(x) {y <- which(x>0); aggregate(x[y], by=list(month[y]), FUN=psych::skew)[,2]}))
    } else {
      Stats[1 + k1, -1, 1, iv, 1, months] <- t(apply(SITE_OBS, c(2), function(x) aggregate(x, by=list(month), FUN=mean)[,2]))
      Stats[1 + k1, -1, 2, iv, 1, months] <- t(apply(SITE_OBS, c(2), function(x) aggregate(x, by=list(month), FUN=sd)[,2]))
      Stats[1 + k1, -1, 3, iv, 1, months] <- t(apply(SITE_OBS, c(2), function(x) aggregate(x, by=list(month), FUN=psych::skew)[,2]))
    }
    
    #Annual mean daily values
    
    AnnualMeans[1 + k1, -1, iv, ] <- t(apply(SITE_OBS, c(2), function(x) {aggregate(x, by=list(wyear), FUN=mean)[,2]})) #years, sites
    
    # print(apply(AnnualMeans[1 + k1, -1, iv, ], 1, FUN=mean))
    
    Stats[1 + k1, -1, 1, iv, 3, 1] <- apply(AnnualMeans[1 + k1, -1, iv, ], 1, FUN=mean)
    Stats[1 + k1, -1, 2, iv, 3, 1] <- apply(AnnualMeans[1 + k1, -1, iv, ], 1, FUN=sd)
    Stats[1 + k1, -1, 3, iv, 3, 1] <- apply(AnnualMeans[1 + k1, -1, iv, ], 1, FUN=psych::skew)
    
    
    #Annual extreme daily values  TODO: -min for tmin, what if only values > 0 are observed? 
    #ext <- t(apply(SITE_OBS * ifelse(toupper(vars[iv]) == "TMIN", -1, 1), c(2), function(x) {aggregate(x, by=list(wyear), FUN=max)[,2]})) #dim(ext) = sites, years; find max for prcp and tmax, but -min for tmin
    ext <- if(toupper(vars[iv]) == "TMIN") {
      t(apply(SITE_OBS, 2, function(x) aggregate(x, by=list(wyear), min)[,2])) #dim(ext) = sites, years; find max for prcp and tmax, but min for tmin
    } else {
      t(apply(SITE_OBS, 2, function(x) aggregate(x, by=list(wyear), max)[,2])) #dim(ext) = sites, years; find max for prcp and tmax, but -min for tmin
    }   
    AnnualEDV[1 + k1, -1, iv, ] <- t(apply(ext, 1, FUN=function(x) return.level(fevd(x, type="GEV", method="MLE", time.units="year"), return.period=return.periods)))
    #AnnualEDV[1 + k1, -1, iv, ] <- ifelse(AnnualEDV[1 + k1, -1, iv, ] > 2 * max(ext), NA, AnnualEDV[1 + k1, -1, iv, ]) #set outliers to NA
    AnnualEDV[AnnualEDV[1 + k1, -1, iv, ] > 2 * max(ext)] <- NA
    
    if(toupper(vars[iv]) == "TMIN") AnnualEDV[1 + k1, -1, iv, ] <- -AnnualEDV[1 + k1, -1, iv, ] #set -min to min for tmin
    
    #PRCP mean monthly wet/dry spells  TODO check this code, necessary to execute for each var again? should be enough for PRCP, no?
    if(toupper(vars[iv]) == "PRCP"){
      for(is in 1:num_site){
      # Spells[1 + k1, 1 + is, , ] <- t(aggregate(sites_prcp[, is], by=list(MONTH_D), FUN=function(x) {
      # Spells[1 + k1, 1 + is, , ] <- t(aggregate(wgen_out_array[[is,1]]$obs['PRCP'], by=list(month), FUN=function(x) {
      #   temp <- rle(x > 0)
      #   res <- sapply(list(temp$lengths[temp$values], temp$lengths[!temp$values]), FUN=function(x) c(mean(x), sd(x), psych::skew(x)))
      #   return(res)})[, 2])
      
        if (any(toupper(spells)=="EXTREME")) { # Currently, if extreme is in spells, ignore it. Implementation will follow
          Spells[1 + k1, 1 + is,c(-7,-8,-9) , ] <- t(aggregate(wgen_out_array[[is,1]]$obs['PRCP'], by=list(month), FUN=function(x) {
            temp <- rle(x > 0)
            res <- sapply(list(temp$lengths[temp$values], temp$lengths[!temp$values]), FUN=function(x) c(mean(x), sd(x), psych::skew(x)))
            return(res)})[, 2])    
        } else {
          # Spells[1 + k1, 1 + is, , ] <- t(aggregate(sites_prcp[, is], by=list(MONTH_D), FUN=function(x) {
          Spells[1 + k1, 1 + is, , ] <- t(aggregate(wgen_out_array[[is,1]]$obs['PRCP'], by=list(month), FUN=function(x) {
            temp <- rle(x > 0)
            res <- sapply(list(temp$lengths[temp$values], temp$lengths[!temp$values]), FUN=function(x) c(mean(x), sd(x), psych::skew(x)))
            return(res)})[, 2])
        }
      }
    }
  }
  AnnualMeans[, , toupper(vars) == "PRCP", ] <- 365 * AnnualMeans[, ,  toupper(vars) == "PRCP", ]
  
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