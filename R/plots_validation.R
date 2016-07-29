#' Generate Observation Statistics from source data
#'
#' @param wgen_out_array, 2 dimensional array of dataframes as collected of returns my wgen_daily, rows are sites, cols are trials
#' @param year contains the year of the observed data  
#' @param Month contains the year of the observed data  
#' @param vars usually keep on default c("prcp","tmax","tmin"). Only change if you know what you do.
#' @param timesteps usually keep on default c("daily","monthly", "annual"). Only change if you know what you do.
#' @param stats usually keep on default c("mean","stdev","skew"). Only change if you know what you do.
#' @export
#' @examples
#' getObservedStatistics3(wgen_out_array)
#'  
#'
getObservedStatistics <- function(wgen_out_array, vars=c("prcp","tmax","tmin"), timesteps= c("daily","monthly", "annual"), stats=  c("mean","stdev","skew"),spells=c("wet","dry"), return.periods=c(10,30,100), mo=c(1:12)){
  
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
  
  
  month <- as.matrix(wgen_out_array[[1,1]]$out['MONTH'])
  
  
  for (k1 in 1:num_trials) for(iv in seq_along(vars)){
    
    # site_prcp <- wgen_out_array[[1,k1]]$out['PRCP']
    
    #SITE_OBS <- switch(iv, sites_prcp, sites_tmax, sites_tmin)
    SITE_OBS <- switch(iv,  wgen_out_array[[1,k1]]$out['PRCP'],  wgen_out_array[[1,k1]]$out['TMAX'],  wgen_out_array[[1,k1]]$out['TMIN'])
    #Monthly mean daily values
    if(vars[iv] == "prcp"){
      Stats[1 + k1, -1, 1, iv, 1, mo] <- t(apply(SITE_OBS, c(2), function(x) 
      {y <- which(x>0); 
      aggregate(x[y], by=list(month[y]), FUN=mean)[,2]}))
      Stats[1 + k1, -1, 2, iv, 1, mo] <- t(apply(SITE_OBS, c(2), function(x) {y <- which(x>0); aggregate(x[y], by=list(month[y]), FUN=sd)[,2]}))
      Stats[1 + k1, -1, 3, iv, 1, mo] <- t(apply(SITE_OBS, c(2), function(x) {y <- which(x>0); aggregate(x[y], by=list(month[y]), FUN=skewness)[,2]}))
    } else {
      Stats[1 + k1, -1, 1, iv, 1, mo] <- t(apply(SITE_OBS, c(2), function(x) aggregate(x, by=list(month), FUN=mean)[,2]))
      Stats[1 + k1, -1, 2, iv, 1, mo] <- t(apply(SITE_OBS, c(2), function(x) aggregate(x, by=list(month), FUN=sd)[,2]))
      Stats[1 + k1, -1, 3, iv, 1, mo] <- t(apply(SITE_OBS, c(2), function(x) aggregate(x, by=list(month), FUN=psych::skew)[,2]))
    }
    
    #Annual mean daily values
    
    AnnualMeans[1 + k1, -1, iv, ] <- t(apply(SITE_OBS, c(2), function(x) {aggregate(x, by=list(wyear), FUN=mean)[,2]})) #years, sites
    
    # print(apply(AnnualMeans[1 + k1, -1, iv, ], 1, FUN=mean))
    Stats[1 + k1, -1, 1, iv, 3, 1] <- apply(AnnualMeans[1 + k1, -1, iv, ], 1, FUN=mean)
    Stats[1 + k1, -1, 2, iv, 3, 1] <- apply(AnnualMeans[1 + k1, -1, iv, ], 1, FUN=sd)
    Stats[1 + k1, -1, 3, iv, 3, 1] <- apply(AnnualMeans[1 + k1, -1, iv, ], 1, FUN=psych::skew)
    
    #Annual extreme daily values
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

#' Generate Simulated Statistics from calculated data
#'
#' @param wgen_out_array outputs of wgen_daily per site and trials collected in array.
#' @param stats c("mean","stdev","skew")
#' @param vars usually keep on default c("prcp","tmax","tmin"). Only change if you know what you do.
#' @param timesteps usually keep on default c("daily","monthly", "annual"). Only change if you know what you do.
#' @param mo usually keep on default c(1:12). Only change if you know what you do.
#' @param return.periods usually keep on default c(10,30,100). Only change if you know what you do.
#' @param s??ells usually keep on default c("wet","dry"). Only change if you know what you do.
#' @export
#' @examples
#' getSimulatedStatistics(sites_prcp, sites_tmax, sites_tmin, num_trials=1, num_site=3, as.matrix(format(sim$out["DATE"], "%Y")),  as.matrix(sim$out["MONTH"]), vars=c("prcp","tmax","tmin"), timesteps= c("daily","monthly", "annual"), stats=c("mean","stdev","skew"))
#'  

getSimulatedStatistics <- function(wgen_out_array, stats=c("mean","stdev","skew"),vars=c("prcp","tmax","tmin"), timesteps=c("daily","monthly", "annual"), mo=c(1:12), return.periods=c(10,30,100), spells=c("wet","dry")){
  
  cur_data <- wgen_out_array[1,1][[1]]$out
  num_site   <- dim(wgen_out_array)[1]
  num_trials <- dim(wgen_out_array)[2]
  wyear <-  as.matrix(format(wgen_out_array[[1,1]]$out['DATE'], "%Y"))  #All simulation outputs must have same length
  month <- as.matrix(wgen_out_array[[1,1]]$out['MONTH'])
  
  #Stats = median + realizations, median + sites, statistics, variables, timesteps, months
  Stats <- array(NA, dim=c(1 + num_trials, 1 + num_site, length(stats), length(vars), length(timesteps), length(mo)),
                 dimnames=list(NULL, NULL, stats, vars, timesteps, mo))
  #AnnualMeans = median + realizations, median + sites, variables, years
  AnnualMeans <- array(0, dim=c(1 + num_trials, 1 + num_site, length(vars), length(unique(wyear))),
                       dimnames=list(NULL, NULL, vars, sort(unique(wyear))))
  #AnnualEDV = median + realizations, median + sites, variables, return periods
  AnnualEDV <- array(NA, dim=c(1 + num_trials, 1 + num_site, length(vars), length(return.periods)),
                     dimnames=list(NULL, NULL, vars, return.periods))
  #Spells = median + realizations, median + sites, stats*{dry,wet}, months
  Spells <- array(0, dim=c(1 + num_trials, 1 + num_site, length(spells) * length(stats), length(mo)),
                  dimnames=list(NULL, NULL, paste(rep(spells, each=length(stats)), stats, sep="_"), mo))
  
  for (ss in 1:num_site) {
    Sim_Variable <- array(NA,c(dim(cur_data)[1],num_trials,  length(vars)))
    
    
    for (k1 in 1:num_trials) {
      # dir_ending1 <- paste0("TRIAL_", k1)
      # dir_ending2 <- paste(dir_ending1, ending2, sep="_")
      # filename <- paste("TRIAL",k1,"SITE",ss,sep="_")
      #cur_data <- data.matrix(read.table(file.path(dat_dir, dir_ending1, dir_ending2, filename)))
      # print(file.path(dat_dir, filename))
      #cur_data <- data.matrix(read.table(file.path(dat_dir, filename)))
      cur_data <- data.matrix(wgen_out_array[ss,k1][[1]]$out[, c("PRCP","TMAX","TMIN")])
      Sim_Variable[, k1, 1:length(vars)] <- cur_data
      # Sim_Variable[, k1, 1:length(vars)] <- cur_data[, 4:(3+length(vars))]
    }
    
    #######Site Statistics
    #sim_temp <- array(NA, dim=c(NUM_TRIALS, length(stats), length(vars), length(timesteps), length(mo)))
    
    for (k1 in 1:num_trials) {
      for(iv in seq_along(vars)){
        if (vars[iv] == "prcp") {
          y <- which(Sim_Variable[, k1, iv]>0)
        } else {
          y <- 1:length(Sim_Variable[, k1, iv])
        }
        #Daily
        Stats[1 + k1, 1 + ss, 1, iv, 1, mo] <- aggregate(Sim_Variable[y, k1, iv], by=list(data.matrix(wgen_out_array[ss,k1][[1]]$out[, c("MONTH")])[y]), FUN=mean)[,2]
        Stats[1 + k1, 1 + ss, 2, iv, 1, mo] <- aggregate(Sim_Variable[y, k1, iv], by=list(data.matrix(wgen_out_array[ss,k1][[1]]$out[, c("MONTH")])[y]), FUN=sd)[,2]
        Stats[1 + k1, 1 + ss, 3, iv, 1, mo] <- aggregate(Sim_Variable[y, k1, iv], by=list(data.matrix(wgen_out_array[ss,k1][[1]]$out[, c("MONTH")])[y]), FUN=skew)[,2]
        #Annual
        #AnnualMeans[1 + k1, 1 + ss, iv, ] <- t(aggregate(Sim_Variable[, k1, iv], by=list(WATER_YEAR_SIM), FUN=mean)[,2])
        AnnualMeans[1 + k1, 1 + ss, iv, ] <- t(aggregate(Sim_Variable[, k1, iv], by=list(wyear), FUN=mean)[,2])
        Stats[1 + k1, 1 + ss, 1, iv, 3, 1] <- mean(AnnualMeans[1 + k1, 1 + ss, iv, ])
        Stats[1 + k1, 1 + ss, 2, iv, 3, 1] <- sd(AnnualMeans[1 + k1, 1 + ss, iv, ])
        Stats[1 + k1, 1 + ss, 3, iv, 3, 1] <- skew(AnnualMeans[1 + k1, 1 + ss, iv, ])
        
        #Annual extreme daily values
        ext <- aggregate(Sim_Variable[, k1, iv] * ifelse(vars[iv] == "tmin", -1, 1), by=list(cur_data[,1]), FUN=max)[,2] #dim(ext) = years; find max for prcp and tmax, but -min for tmin
        AnnualEDV[1 + k1, 1 + ss, iv, ] <- return.level(fevd(ext, type="GEV", method="MLE", time.units="year"), return.period=return.periods)
        AnnualEDV[1 + k1, 1 + ss, iv, ] <- ifelse(AnnualEDV[1 + k1, 1 + ss, iv, ] > 2 * max(ext), NA, AnnualEDV[1 + k1, 1 + ss, iv, ]) #set outliers to NA
        if(vars[iv] == "tmin") AnnualEDV[1 + k1, 1 + ss, iv, ] <- -AnnualEDV[1 + k1, 1 + ss, iv, ] #set -min to min for tmin
        
        #PRCP mean monthly wet/dry spells
        if(vars[iv] == "prcp"){
          #Spells[1 + k1, 1 + ss, , ] <- t(aggregate(Sim_Variable[, k1, iv], by=list(MONTH_SIM), FUN=function(x) {
          Spells[1 + k1, 1 + ss, , ] <- t(aggregate(Sim_Variable[, k1, iv], by=list(data.matrix(wgen_out_array[ss,k1][[1]]$out[, c("MONTH")])), FUN=function(x) {
            temp <- rle(x > 0)
            res <- sapply(list(temp$lengths[temp$values], temp$lengths[!temp$values]), FUN=function(x) c(mean(x), sd(x), skew(x)))
            return(res)})[, 2])
        }
      }
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

#' Plot Monthly variables
#'
#' @param fig_dir  Directory to create PDF containing the plot. Can be NA to plot on screen
#' @param obs observation data
#' @param sim simulation output
#' @param ftag string to appear in filename
#' @export
#' @examples
#' plot_basinwide_vars_monthly(fig_dir, obs, sim)
#' 
#'
plot_basinwide_vars_monthly <- function(obs, sim, ftag="", fig_dir=NULL, stats=c("mean","stdev","skew"), vars=c("prcp","tmax","tmin"), mo=c(1:12)){
  if(!is.null(fig_dir)){
    pdf(width=3*length(vars), height=3*length(stats), file=file.path(fig_dir, paste0("1_BasinWide_Trials_MonthlyVariables", ifelse(nchar(ftag > 0), "_"), ftag, ".pdf")))
  }
  layout(matrix(1:(length(vars)*length(stats)), ncol=length(vars), nrow=length(stats), byrow=FALSE))
  for(iv in seq_along(vars)){
    for(is in seq_along(stats)){
      x <- sim[-1, 1, is, iv, 1, mo]
      y <- obs[-1, 1, is, iv, 1, mo]
      ylim <- c(min(x, y, na.rm = TRUE), max(x, y, na.rm = TRUE))
      boxplot(x, ylim=ylim, xlab=ifelse(is == length(stats), "Months", ""), ylab=ifelse(iv == 1, stats[is], ""), main=ifelse(is == 1, vars[iv], ""))
      points(1:12, y, col="red")
    }
  }
  legend(x="bottom", legend=c("Observations", "Simulated trials"), pch=19, col=c("red", "black"))
  if(!is.null(fig_dir)){
    dev.off()
  }
}

#' Plot Monthly spells
#'
#' @param fig_dir  Directory to create PDF containing the plot
#' @param obs observation data
#' @param sim simulation output
#' @param ftag string to appear in filename
#' @export
#' @examples
#' plot_basinwide_spells_monthly(fig_dir, obs, sim)
#' 
#'
plot_basinwide_spells_monthly <- function(obs, sim, fig_dir=NULL, ftag="", stats=c("mean","stdev","skew"), vars=c("prcp","tmax","tmin"), mo=c(1:12), spells=c("wet","dry")){
  if(!is.null(fig_dir)){
    pdf(width=3*length(spells), height=3*length(stats), file=file.path(fig_dir, paste0("2_BasinWide_Trials_MonthlySpells", ifelse(nchar(ftag > 0), "_"), ftag, ".pdf")))
  }
  layout(matrix(1:(length(spells)*length(stats)), ncol=length(spells), nrow=length(stats), byrow=FALSE))
  
  for(iv in 1:length(spells)){
    for(is in seq_along(stats)){
      x <- sim[-1, 1, (iv - 1) * 2 + is, mo]
      y <- obs[-1, 1, (iv - 1) * 2 + is, mo]
      ylim <- c(min(x, y, na.rm = TRUE), max(x, y, na.rm = TRUE))
      boxplot(x, ylim=ylim, xlab=ifelse(is == length(stats), "Months", ""), ylab=ifelse(iv == 1, paste(stats[is], "(days)"), ""), main=ifelse(is == 1, spells[iv], ""))
      points(1:12, y, col="red")
    }
  }
  legend(x="bottom", legend=c("Observations", "Simulated trials"), pch=19, col=c("red", "black"))
  if(!is.null(fig_dir)){
    dev.off()
  }
}

#' Plot annual extremes
#'
#' @param obs observation data
#' @param sim simulation output
#' @param fig_dir  Directory to create PDF containing the plot
#' @param ftag string to appear in filename
#' @export
#' @examples
#' plot_basinwide_extreme_annual(obs, sim)
#' 
#'
plot_basinwide_extreme_annual <- function(obs, sim, fig_dir=NULL, ftag="", vars=c("prcp","tmax","tmin"), return.periods=c(10,30,100)){
  if(!is.null(fig_dir)){
    pdf(width=3*length(vars), height=3*1, file=file.path(fig_dir, paste0("3_BasinWide_Trials_AnnualExtremeEvents", ifelse(nchar(ftag > 0), "_"), ftag, ".pdf")))
  }
  layout(matrix(1:(length(vars)*1), ncol=length(vars), nrow=1, byrow=FALSE))
  
  for(iv in seq_along(vars)){
    x <- sim[-1, 1, iv, ]
    y <- obs[-1, 1, iv, ]
    ylim <- c(min(x, y, na.rm = TRUE), max(x, y, na.rm = TRUE))
    boxplot(x, ylim=ylim, xlab="Return period (years)", ylab=vars[iv])
    points(1:length(return.periods), y, col="red")
  }
  legend(x="bottom", legend=c("Observations", "Simulated trials"), pch=19, col=c("red", "black"))
  if(!is.null(fig_dir)){
    dev.off()
  }
}

#' Plot sitespecific variables
#'
#' @param obs observation data
#' @param sim simulation output
#' @param fig_dir  Directory to create PDF containing the plot
#' @param ftag string to appear in filename
#' @param timesteps
#' @export
#' @examples
#' plot_sitespecific_vars(obs, sim)
#' 
#'
plot_sitespecific_vars <- function(obs, sim, fig_dir=NULL,  ftag="", timestep="daily", stats=c("mean","stdev","skew"), vars=c("prcp","tmax","tmin"), mo=c(1:12)){
  timestep <- match.arg(timestep, timesteps)
  it <- switch(timestep, daily=1, monthly=2, annual=3)
  mo <- switch(timestep, daily=1:12, monthly=1:12, annual=1)
  if(!is.null(fig_dir)){
    pdf(width=3*length(vars), height=3*length(stats), file=file.path(fig_dir, paste0("4_Sites_MedianTrials_", timestep, "Variables", ifelse(nchar(ftag > 0), "_"), ftag, ".pdf")))
  }
  layout(matrix(1:(length(vars)*length(stats)),length(vars),length(stats),byrow=FALSE))
  
  for(iv in seq_along(vars)){
    for(is in seq_along(stats)){
      y <- as.vector(sim[1, -1, is, iv, it, mo])
      x <- as.vector(obs[1, -1, is, iv, it, mo])
      xlim <- ylim <- c(min(x, y, na.rm = TRUE), max(x, y, na.rm = TRUE))
      plot(x, y, xlim=xlim, ylim=ylim, xlab="Observed", ylab="Simulated", main=paste(timestep, stats[is], vars[iv]))
      abline(0,1)
    }
  }
  if(!is.null(fig_dir)){
    dev.off()
  }
}

#' Plot sitespecific annual correlations
#'
#' @param obsCors observation data
#' @param simCors simulation output
#' @param fig_dir  Directory to create PDF containing the plot
#' @param ftag string to appear in filename
#' @param vars variables
#' @export
#' @examples
#' plot_sitespecific_correlations_annual (obsCors, simCors)
#' 
#'
plot_sitespecific_correlations_annual <- function(obsCors, simCors, fig_dir=NULL, ftag="Validation", vars=c("prcp","tmax","tmin")){
  if(!is.null(fig_dir)){
    pdf(width=3*length(vars), height=3*2, file=file.path(fig_dir, paste0("5_Sites_MedianTrials_AnnualCorrelations", ifelse(nchar(ftag > 0), "_"), ftag, ".pdf")))
  }
  layout(matrix(1:(length(vars)*2), nrow=2, ncol=length(vars), byrow=FALSE))
  
  for(iv in seq_along(vars)){
    for(ic in 1:2){
      if(ic == 1){
        tag <- "Intersite correlation"
        y <- as.vector(simCors$intersites[, 2 + iv])
        x <- as.vector(obsCors$intersites[, 2 + iv])
      } else {
        tag <- "Lag-1 autocorrelation"
        y <- as.vector(simCors$lag1auto[, iv])
        x <- as.vector(obsCors$lag1auto[, iv])
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

#' Plot sitespecific annual spells
#'
#' @param obs observation data
#' @param sim simulation output data
#' @param fig_dir  Directory to create PDF containing the plot
#' @param ftag string to appear in filename
#' @export
#' @examples
#' plot_sitespecific_spells_annual (obs, sim)
#' 
#'
plot_sitespecific_spells_annual <- function(obs, sim, fig_dir=NULL, ftag="Validation", stats=c("mean","stdev","skew"), vars=c("prcp","tmax","tmin"), mo=c(1:12),spells=c("wet","dry")){
  if(!is.null(fig_dir)){
    pdf(width=3*length(spells), height=3*length(stats), file=file.path(fig_dir, paste0("6_Sites_MedianTrials_AnnualSpells", ifelse(nchar(ftag > 0), "_"), ftag, ".pdf")))
  }
  layout(matrix(1:(length(spells)*length(stats)), ncol=length(spells), nrow=length(stats), byrow=FALSE))
  
  for(iv in 1:length(spells)){
    for(is in seq_along(stats)){
      y <- apply(sim[1, -1, (iv - 1) * 2 + is, mo], MARGIN=1, FUN=mean) #mean across months
      x <- apply(obs[1, -1, (iv - 1) * 2 + is, mo], MARGIN=1, FUN=mean)
      xlim <- ylim <- c(min(x, y, na.rm = TRUE), max(x, y, na.rm = TRUE))
      try(plot(x, y, xlim=xlim, ylim=ylim, xlab="Observed spells (days)", ylab="Simulated spells (days)", main=paste(spells[iv], stats[is])))
      abline(0,1)
    }
  }
  
  if(!is.null(fig_dir)){
    dev.off()
  }
}

#' Plot PowerSpectrum
#'
#' @param PowerSpectrum, create this data this using Final_Annual_Sim_All and with its output getLowFrequencyVariability
#' @param fig_dir  Directory to create PDF containing the plot
#' @param ftag string to appear in filename
#' @export
#' @examples
#' plot7.basinwide.PowerSpectruml (fig_dir, PowerSpectrum)
#' 
#'
plot_basinwide_powerspectrum <- function(PowerSpectrum, fig_dir=NULL, ftag=""){
  if(!is.null(fig_dir)){
    pdf(width=9, height=3, file=file.path(fig_dir, paste0("7_Sites_MedianTrials_PowerSpectrum", ifelse(nchar(ftag > 0), "_"), ftag, ".pdf")))
  }
  plot(PowerSpectrum$Power_Spectrum_Period, PowerSpectrum$Power_Spectrum_Obs,xlab="Period (years)", ylab="Power (mm2)", type="l",col="black",lwd=2)
  polygon(c(PowerSpectrum$Power_Spectrum_Period,rev(PowerSpectrum$Power_Spectrum_Period)),c(PowerSpectrum$Sim_Power_Spectrum_3,rev(PowerSpectrum$Sim_Power_Spectrum_2)),col="grey")
  lines(PowerSpectrum$Power_Spectrum_Period, PowerSpectrum$Sim_Power_Spectrum_1,type="l",col="black",lwd=2)
  lines(PowerSpectrum$Power_Spectrum_Period, PowerSpectrum$Power_Spectrum_Obs,col="red",lwd=2)
  
  if(!is.null(fig_dir)){
    dev.off()
  }
}

get_new_ar_models <- function(ar_models) {
  AUTOCOR_ANNUAL_CHANGES <- c(1)

  #Change in Standard Deviation of Annual Series (multiplicative)
  # - makes the annual time series more variable, but will not change the persistence
  # - use for ???Inter-annual severity of dry/wet years???
  MAGNITUDE_NOISE_ANNUAL_CHANGES <- c(.5,1,2)
  NUM_AUTOCOR_ANNUAL_CHANGES <- length(AUTOCOR_ANNUAL_CHANGES)
  NUM_MAGNITUDE_NOISE_ANNUAL_CHANGES <-length(MAGNITUDE_NOISE_ANNUAL_CHANGES)
  NEW_AR_MODEL <- list()
  count <- 0
  for (k2 in 1:NUM_AUTOCOR_ANNUAL_CHANGES) {
    for (k3 in 1:NUM_MAGNITUDE_NOISE_ANNUAL_CHANGES) {
      count <- count + 1
      #cur_dir <- main_dir
      Lag1_Change <- AUTOCOR_ANNUAL_CHANGES[k2]
      Stdev_Change <- MAGNITUDE_NOISE_ANNUAL_CHANGES[k3]
      CUR_MODEL <- ar_models[[1]]
      e_var_org <- CUR_MODEL$sigma2
    AR_org <- as.vector(CUR_MODEL$coef)[which(names(CUR_MODEL$coef)!="intercept" & substr(names(CUR_MODEL$coef),1,2)!="ma")]
    MA_org <- as.vector(CUR_MODEL$coef)[which(names(CUR_MODEL$coef)!="intercept" & substr(names(CUR_MODEL$coef),1,2)!="ar")]		
    Other_e_var <- list()
    Other_AR <- list()
    Other_MA <- list()
    if (USE_WARM_MODEL) {
      for (l in 2:dim(Wavelet_Decomposition)[2]) {
        CUR_MODEL <- ar_models[[l]]
        Other_e_var[[l-1]] <- CUR_MODEL$sigma2
        Other_AR[[l-1]] <- as.vector(CUR_MODEL$coef)[which(names(CUR_MODEL$coef)!="intercept" & substr(names(CUR_MODEL$coef),1,2)!="ma")]
        Other_MA[[l-1]] <- as.vector(CUR_MODEL$coef)[which(names(CUR_MODEL$coef)!="intercept" & substr(names(CUR_MODEL$coef),1,2)!="ar")]	
      }
    }
    new_parameters <- Adjust_WARM_Model(Stdev_Change,Lag1_Change,e_var_org,AR_org,MA_org,Other_e_var,Other_AR,Other_MA,cur_dir)
    NEW_AR_MODEL[[count]] <- ar_models[[1]]
    NEW_AR_MODEL[[count]]$sigma2 <- new_parameters[1]
    non_intercept <- which(names(ar_models[[1]]$coef)!="intercept")
    if (length(non_intercept)>0) {
      for (cur_par in 1:length(non_intercept)) {NEW_AR_MODEL[[count]]$coef[[cur_par]] <- new_parameters[[cur_par+1]]}
    }
  }
}	

get_ar_models <- function(annual_prcp, years, use_warm_model = FALSE, sig.level =0.9, n.periods=1, sig.periods=c(8,9,10,11,12,13,14,15), n.comp.periods=8)
ar_models <- list()
if (use_warm_model) {
  #Wavelet_Decomposition <- WAVELET_DECOMPOSITION(CLIMATE_VARIABLE,NUM_FINAL_PERIODS,ALL_SIG_PERIODS,NUM_PERIODS_ALL_COMPS,plot_flag=TRUE, figname="WaveletDecomposition_AnnualPRCP.pdf")
  Wavelet_Decomposition  <- wavelet_components(annual_prcp,wavelet_analysis(ANNUAL_PRCP,years,sig.level, "white"), n.periods, sig.periods,n.comp.periods)
  for (k in 1:dim(Wavelet_Decomposition)[2]) {
    ar_models[[k]] <- auto.arima(Wavelet_Decomposition[,k],max.p=2,max.q=2,max.P=0,max.Q=0,stationary=TRUE)
  }
} else {
  ar_models[[1]] <- auto.arima(ANNUAL_PRCP,max.p=2,max.q=2,max.P=0,max.Q=0,stationary=TRUE)
}
return(ar_models)
}

#' get_final_annual_sim_all helper function to retrieve powerspectrum.
#'
#' @param fig_dir  Directory to create PDF containing the plot
#' @param PowerSpectrum, create this data this using Final_Annual_Sim_All and with its output getLowFrequencyVariability
#' @param ftag string to appear in filename
#' @export
#' @examples
#' ar_models <- list()
#' climate_variable <- annual_prcp
# if (use_warm_model) {
#   Wavelet_Decomposition  <- wavelet_components(CLIMATE_VARIABLE,wavelet_analysis(CLIMATE_VARIABLE,years,0.9, "white"), NUM_FINAL_PERIODS,ALL_SIG_PERIODS,NUM_PERIODS_ALL_COMPS)
#   for (k in 1:dim(Wavelet_Decomposition)[2]) {
#     ar_models[[k]] <- auto.arima(Wavelet_Decomposition[,k],max.p=2,max.q=2,max.P=0,max.Q=0,stationary=TRUE)
#   }
# } else {
#   ar_models[[1]] <- auto.arima(climate_variable,max.p=2,max.q=2,max.P=0,max.Q=0,stationary=TRUE)
# }
#' get_final_annual_sim_all(num_autocor_annual_changes=c(1),num_magnitude_noise_annual_changes=c(1), ar_models, new_ar_model)
#' 
#'
get_final_annual_sim_all <- function( annual_prcp, years,num_year_sim, num_trials=1, num_autocor_annual_changes=1,num_magnitude_noise_annual_changes=1){
  ar_models <- get_ar_models(annual_prcp, years)
  new_ar_model <- get_new_ar_models(ar_models)
  Final_Annual_Sim_All <- array(NA,c(num_year_sim,num_trials,num_autocor_annual_changes,num_magnitude_noise_annual_changes))
  Annual_Sim <- array(NA,c(num_year_sim,length(ar_models)))
  for (k1 in 1:num_trials) {
    count <- 0	
    for (k2 in 1:num_magnitude_noise_annual_changes) {
      for (k3 in 1:num_magnitude_noise_annual_changes) {
        count <- count + 1
        for (l in 1:length(ar_models)) {
          CUR_MODEL <- ar_models[[l]]
          if (l==1) {CUR_MODEL <- NEW_AR_MODEL[[count]]}
          AR <- as.vector(CUR_MODEL$coef)[which(names(CUR_MODEL$coef)!="intercept" & substr(names(CUR_MODEL$coef),1,2)!="ma")]
          MA <- as.vector(CUR_MODEL$coef)[which(names(CUR_MODEL$coef)!="intercept" & substr(names(CUR_MODEL$coef),1,2)!="ar")]
          INTERCEPT <- 0		#in auto.arima, the "intercept" is actually the mean (see http://www.stat.pitt.edu/stoffer/tsa2/Rissues.htm)
          if (length(which(names(CUR_MODEL$coef)=="intercept"))>0) {INTERCEPT <- as.vector(CUR_MODEL$coef)[which(names(CUR_MODEL$coef)=="intercept")]}
          Annual_Sim[,l] <- arima.sim(n = num_year_sim, list(ar = AR, ma=MA),sd = sqrt(CUR_MODEL$sigma2[[1]])) + INTERCEPT
        }
        if (dim(Annual_Sim)[2]>1) {PRCP_FINAL_ANNUAL_SIM <- apply(Annual_Sim,FUN=sum,1)} else {PRCP_FINAL_ANNUAL_SIM <- Annual_Sim[,1]}
        Final_Annual_Sim_All[,k1,k2,k3] <- PRCP_FINAL_ANNUAL_SIM
      }
    }
  }
  return( Final_Annual_Sim_All)
}

#' getCorrelations calculate correlations for plot function 5
#'
#' @param AnnualMeans This value can be retrieved from getObservedStatistics(...)$AnnualMeans or getSimulatedStatistics(..)$AnnualMeans
#' @param vars
#' @export
#' @examples
#' myobsStats$cors <- getCorrelations(AnnualMeans=myobsStats$AnnualMeans[-1, -1, , , drop=FALSE])
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

#' getCorrelations calculate correlations for plot function 5
#'
#' @param AnnualMeans This value can be retrieved from getObservedStatistics(...)$AnnualMeans or getSimulatedStatistics(..)$AnnualMeans
#' @param vars
#' @export
#' @examples
#' myobsStats$cors <- getCorrelations(AnnualMeans=myobsStats$AnnualMeans[-1, -1, , , drop=FALSE])
#' 
#'
#'
getLowFrequencyVariability <- function(ANNUAL_PRCP, years, Final_Annual_Sim_All){
  #Observed
  #Wavelet_Analysis <- WAVELET_ANALYSIS(ANNUAL_PRCP,siglvl=0.90,background_noise="white",plot_flag=FALSE)
  Wavelet_Analysis <- wavelet_analysis(ANNUAL_PRCP,years, sig.level=0.90, noise.type="white")  ## function WAVELET_ANALYSIS and wavelet analysis should be compared.
  #num_power<- length(Wavelet_Analysis[[1]])
  num_power <- length(Wavelet_Analysis$gws)                  # -1 ??
  # Power_Spectrum_Obs <- Wavelet_Analysis[[1]]    
  Power_Spectrum_Obs <- Wavelet_Analysis$gws
  #Power_Spectrum_Sig <- Wavelet_Analysis[[2]]
  Power_Spectrum_Sig <- var(ANNUAL_PRCP) * Wavelet_Analysis$gws.sig$signif
  #Power_Spectrum_Period <- Wavelet_Analysis[[3]]
  Power_Spectrum_Period <- Wavelet_Analysis$period
  
  # plot(out_WA[[2]] / var(in_annual_prcp), col = "red")
  # lines(out_wn[["gws.sig"]]$signif)
  
  #Simulated
  NUM_TRIALS <- dim(Final_Annual_Sim_All)[2]
  
  #	Annual_Area_Avg_Wavelet_Sim <- array(NA,c(num_power,NUM_TRIALS))
  #	for (k1 in 1:NUM_TRIALS) {
  #		area_prcp <- array(0,length(cur_data[,4]))
  #		for (ss in 1:num_site) {
  #			dir_ending1 <- paste("TRIAL_",k1,sep="")
  #			dir_ending2 <- paste(dir_ending1,"_A.P.COR_1_A.P.MAG_1_D.P.WETSPELL_1_D.P.DRYSPELL_1_D.P.MEAN_1_D.P.CV_1_D.T.MEAN_1",sep="")
  #			dir_name1 <- paste(dat_dir,dir_ending1,dir_ending2,sep="/")
  #			setwd(dir_name1)
  #			cur_data <- data.matrix(read.table(filenames[ss]))
  #			area_prcp <- area_prcp + cur_data[,4]
  #		}
  #		area_prcp <- area_prcp/num_site
  #		area_annual_prcp <- aggregate(area_prcp,FUN=mean,list(WATER_YEAR_SIM))[,2]*365
  #		Wavelet_Analysis <- WAVELET_ANALYSIS(area_annual_prcp,siglvl=0.90,background_noise="white",plot_flag=FALSE)
  #		Annual_Area_Avg_Wavelet_Sim[,k1] <- Wavelet_Analysis[[1]][1:num_power]
  #
  #	}
  
  Annual_Area_Avg_Wavelet_Sim <- array(NA,c(num_power,NUM_TRIALS))
  for (k1 in 1:NUM_TRIALS) {
    #Wavelet_Analysis <- WAVELET_ANALYSIS(Final_Annual_Sim_All[,k1,,],siglvl=0.90,background_noise="white",plot_flag=FALSE)
    Wavelet_Analysis <- wavelet_analysis(Final_Annual_Sim_All[,k1,,],years, sig.level=0.90,noise.type="white")
    Annual_Area_Avg_Wavelet_Sim[,k1] <- Wavelet_Analysis[[1]][1:num_power]
  }
  
  Sim_Power_Spectrum_1 <- apply(Annual_Area_Avg_Wavelet_Sim,FUN=mean,1)
  Sim_Power_Spectrum_2 <- apply(Annual_Area_Avg_Wavelet_Sim,FUN=quantile,1,.95)
  Sim_Power_Spectrum_3 <- apply(Annual_Area_Avg_Wavelet_Sim,FUN=quantile,1,.05)
  
  return(list(
    Power_Spectrum_Obs=Power_Spectrum_Obs,
    Power_Spectrum_Sig=Power_Spectrum_Sig,
    Power_Spectrum_Period=Power_Spectrum_Period,
    Sim_Power_Spectrum_1=Sim_Power_Spectrum_1,
    Sim_Power_Spectrum_2=Sim_Power_Spectrum_2,
    Sim_Power_Spectrum_3=Sim_Power_Spectrum_3))
} 
