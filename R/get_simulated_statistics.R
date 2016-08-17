#' Generate Statistics from simulation data. Simulation data is expected to be a list of wgen_daily() results with multiple sites and trials per site.
#' Each row is a site, each column a trial.
#'
#' @param wgen_out_array simulation outputs of wgen_daily per site and trials collected in array.
#' @param vars variables of which stats are calculated. Default is daily precipitation,  max and min temperature and wind c("prcp","tmax","tmin", "wind").
#' @param return.periods default is c(10,30,100).
#' @param spells, default is c("wet","dry"). TODO: The package uses dry, wet, extreme check if we should switch defaults to c('d', 'w', 'e')
#' @export
#' @examples
#' get_simulated_statistics(wgen_out_array, stats=c("mean","stdev","skew"), vars=c("prcp","tmax","tmin"), timesteps= c("daily","monthly", "annual"), months=c(1:12), return.periods=c(10,30,100), spells=c("wet","dry"))
#'  

get_simulated_statistics <- function(wgen_out_array, vars=c("prcp","tmax","tmin", "wind"), return.periods=c(10,30,100), spells=c("wet","dry")){
  
  stats <- c("mean","stdev","skew")
  months <- c(1:12) 
  timesteps=c("daily","monthly", "annual")
  cur_data <- wgen_out_array[1,1][[1]]$out
  num_site   <- dim(wgen_out_array)[1]
  num_trials <- dim(wgen_out_array)[2]
  wyear <-  as.matrix(format(wgen_out_array[[1,1]]$out['DATE'], "%Y"))  #All simulation outputs must have same length
  month <- as.matrix(wgen_out_array[[1,1]]$out['MONTH'])
  
  #Stats = median + realizations, median + sites, statistics, variables, timesteps, months
  Stats <- array(NA, dim=c(1 + num_trials, 1 + num_site, length(stats), length(vars), length(timesteps), length(months)),
                 dimnames=list(NULL, NULL, stats, vars, timesteps, months))
  #AnnualMeans = median + realizations, median + sites, variables, years
  AnnualMeans <- array(0, dim=c(1 + num_trials, 1 + num_site, length(vars), length(unique(wyear))),
                       dimnames=list(NULL, NULL, vars, sort(unique(wyear))))
  #AnnualEDV = median + realizations, median + sites, variables, return periods
  AnnualEDV <- array(NA, dim=c(1 + num_trials, 1 + num_site, length(vars), length(return.periods)),
                     dimnames=list(NULL, NULL, vars, return.periods))
  #Spells = median + realizations, median + sites, stats*{dry,wet}, months
  Spells <- array(0, dim=c(1 + num_trials, 1 + num_site, length(spells) * length(stats), length(months)),
                  dimnames=list(NULL, NULL, paste(rep(spells, each=length(stats)), stats, sep="_"), months))
  
  for (ss in 1:num_site) {
    Sim_Variable <- array(NA,c(dim(cur_data)[1],num_trials,  length(vars)))
    
    
    for (k1 in 1:num_trials) {
      cur_data <- data.matrix(wgen_out_array[ss,k1][[1]]$out[, toupper(vars)])
      Sim_Variable[, k1, 1:length(vars)] <- cur_data
    }
    
    #######Site Statistics
    #sim_temp <- array(NA, dim=c(NUM_TRIALS, length(stats), length(vars), length(timesteps), length(months)))
    
    for (k1 in 1:num_trials) {
      for(iv in seq_along(vars)){
        if (toupper(vars[iv]) == "PRCP") {
          y <- which(Sim_Variable[, k1, iv]>0)
        } else {
          y <- 1:length(Sim_Variable[, k1, iv])
        }
        #Daily
        Stats[1 + k1, 1 + ss, 1, iv, 1, months] <- aggregate(Sim_Variable[y, k1, iv], by=list(data.matrix(wgen_out_array[ss,k1][[1]]$out[, c("MONTH")])[y]), FUN=mean)[,2]
        Stats[1 + k1, 1 + ss, 2, iv, 1, months] <- aggregate(Sim_Variable[y, k1, iv], by=list(data.matrix(wgen_out_array[ss,k1][[1]]$out[, c("MONTH")])[y]), FUN=sd)[,2]
        Stats[1 + k1, 1 + ss, 3, iv, 1, months] <- aggregate(Sim_Variable[y, k1, iv], by=list(data.matrix(wgen_out_array[ss,k1][[1]]$out[, c("MONTH")])[y]), FUN=skew)[,2]
        #Annual
        #AnnualMeans[1 + k1, 1 + ss, iv, ] <- t(aggregate(Sim_Variable[, k1, iv], by=list(WATER_YEAR_SIM), FUN=mean)[,2])
        AnnualMeans[1 + k1, 1 + ss, iv, ] <- t(aggregate(Sim_Variable[, k1, iv], by=list(wyear), FUN=mean)[,2])
        Stats[1 + k1, 1 + ss, 1, iv, 3, 1] <- mean(AnnualMeans[1 + k1, 1 + ss, iv, ])
        Stats[1 + k1, 1 + ss, 2, iv, 3, 1] <- sd(AnnualMeans[1 + k1, 1 + ss, iv, ])
        Stats[1 + k1, 1 + ss, 3, iv, 3, 1] <- skew(AnnualMeans[1 + k1, 1 + ss, iv, ])
        
        #Annual extreme daily values
        ext <- if(toupper(vars[iv]) == "TMIN") {
          aggregate(Sim_Variable[, k1, iv], by=list(cur_data[,1]), FUN=min)[,2] #dim(ext) = years; find max for prcp and tmax, but -min for tmin
        } else {
          aggregate(Sim_Variable[, k1, iv], by=list(cur_data[,1]), FUN=max)[,2] #dim(ext) = years; find max for prcp and tmax, but -min for tmin
        }           
        #ext <- aggregate(Sim_Variable[, k1, iv] * ifelse(vars[iv] == "tmin", -1, 1), by=list(cur_data[,1]), FUN=max)[,2] #dim(ext) = years; find max for prcp and tmax, but -min for tmin
        AnnualEDV[1 + k1, 1 + ss, iv, ] <- return.level(fevd(ext, type="GEV", method="MLE", time.units="year"), return.period=return.periods)
        AnnualEDV[1 + k1, 1 + ss, iv, ] <- ifelse(AnnualEDV[1 + k1, 1 + ss, iv, ] > 2 * max(ext), NA, AnnualEDV[1 + k1, 1 + ss, iv, ]) #set outliers to NA
        if(toupper(vars[iv]) == "TMIN") AnnualEDV[1 + k1, 1 + ss, iv, ] <- -AnnualEDV[1 + k1, 1 + ss, iv, ] #set -min to min for tmin
        
        #PRCP mean monthly wet/dry spells
        if(toupper(vars[iv]) == "PRCP"){
          if (any(toupper(testspells)=="EXTREME")) { # Currently, if extreme is in spells, ignore it. Implementation will follow
            Spells[1 + k1, 1 + ss, c(-7,-8,-9), ] <- t(aggregate(Sim_Variable[, k1, iv], by=list(data.matrix(wgen_out_array[ss,k1][[1]]$out[, c("MONTH")])), FUN=function(x) {
              temp <- rle(x > 0)
              res <- sapply(list(temp$lengths[temp$values], temp$lengths[!temp$values]), FUN=function(x) c(mean(x), sd(x), skew(x)))
              return(res)})[, 2])           
          } else {  
            #Spells[1 + k1, 1 + ss, , ] <- t(aggregate(Sim_Variable[, k1, iv], by=list(MONTH_SIM), FUN=function(x) {
            Spells[1 + k1, 1 + ss, , ] <- t(aggregate(Sim_Variable[, k1, iv], by=list(data.matrix(wgen_out_array[ss,k1][[1]]$out[, c("MONTH")])), FUN=function(x) {
              temp <- rle(x > 0)
              res <- sapply(list(temp$lengths[temp$values], temp$lengths[!temp$values]), FUN=function(x) c(mean(x), sd(x), skew(x)))
              return(res)})[, 2])
          }
        }
      }
    }
  }
  AnnualMeans[, , toupper(vars) == "PRCP", ] <- 365 * AnnualMeans[, , toupper(vars) == "PRCP", ]
  
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