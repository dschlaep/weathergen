#' Generate Observation Statistics from source data
#'
#' @param sites_prcp  Observed prcp data as matrix with as many columns as we have sites. Each column represents data of another site 
#' @param tmax data.frame with annual climate timeseries at selected location (YEAR, MONTH, PRCP, TMAX, TMIN, TAVG)
#' @param tmin number of iterations for arima simulation
#' @param num_trials number of runs of weathergenerator, how many datasets were created?
#' @param num_sites number of sites. This should be > 3 to use the validation plots
#' @param year contains the year of the observed data  
#' @param Month contains the year of the observed data  
#' @param vars usually keep on default c("prcp","tmax","tmin"). Only change if you know what you do.
#' @param timesteps usually keep on default c("daily","monthly", "annual"). Only change if you know what you do.
#' @param stats usually keep on default c("mean","stdev","skew"). Only change if you know what you do.
#' @export
#' @examples
#' getObservedStatistics3(sites_prcp, sites_tmax, sites_tmin, num_trials=1, num_site=3, as.matrix(format(sim$out["DATE"], "%Y")),  as.matrix(sim$out["MONTH"]), vars=c("prcp","tmax","tmin"), timesteps= c("daily","monthly", "annual"), stats=c("mean","stdev","skew"))
#'  
#'

getObservedStatistics <- function(sites_prcp, sites_tmax, sites_tmin, num_trials, num_site, wyear, month=c(1:12), vars=c("prcp","tmax","tmin"), timesteps= c("daily","monthly", "annual"), stats=  c("mean","stdev","skew")){
  mo <- (1:12)
  #Stats = median + realizations, median + sites, statistics, variables, timesteps, months
  Stats <- array(NA, dim=c(1 + num_trials, 1 + num_site, length(stats), length(vars), length(timesteps), length(mo)),
                 dimnames=list(NULL, NULL, stats, vars, timesteps, mo))
  #AnnualMeans = median + realizations, median + sites, variables, years
  AnnualMeans <- array(0, dim=c(1 + num_trials, 1 + num_site, length(vars), length(unique(wyear))), #WATER_YEAR_SIM
                       dimnames=list(NULL, NULL, vars, sort(unique(wyear))))
  #AnnualEDV = median + realizations, median + sites, variables, return periods
  AnnualEDV <- array(NA, dim=c(1 + num_trials, 1 + num_site, length(vars), length(return.periods)),
                     dimnames=list(NULL, NULL, vars, return.periods))
  #Spells = median + realizations, median + sites, stats*{dry,wet}, months
  Spells <- array(0, dim=c(1 + num_trials, 1 + num_site, length(spells) * length(stats), length(mo)),
                  dimnames=list(NULL, NULL, paste(rep(spells, each=length(stats)), stats, sep="_"), mo))
  
  for (k1 in 1:num_trials) for(iv in seq_along(vars)){
    
    SITE_OBS <- switch(iv, sites_prcp, sites_tmax, sites_tmin)
    
    #Monthly mean daily values
    if(vars[iv] == "prcp"){
      Stats[1 + k1, -1, 1, iv, 1, mo] <- t(apply(SITE_OBS, c(2), function(x) {y <- which(x>0); aggregate(x[y], by=list(month[y]), FUN=mean)[,2]}))
      Stats[1 + k1, -1, 2, iv, 1, mo] <- t(apply(SITE_OBS, c(2), function(x) {y <- which(x>0); aggregate(x[y], by=list(month[y]), FUN=sd)[,2]}))
      Stats[1 + k1, -1, 3, iv, 1, mo] <- t(apply(SITE_OBS, c(2), function(x) {y <- which(x>0); aggregate(x[y], by=list(month[y]), FUN=skew)[,2]}))
    } else {
      Stats[1 + k1, -1, 1, iv, 1, mo] <- t(apply(SITE_OBS, c(2), function(x) aggregate(x, by=list(month), FUN=mean)[,2]))
      Stats[1 + k1, -1, 2, iv, 1, mo] <- t(apply(SITE_OBS, c(2), function(x) aggregate(x, by=list(month), FUN=sd)[,2]))
      Stats[1 + k1, -1, 3, iv, 1, mo] <- t(apply(SITE_OBS, c(2), function(x) aggregate(x, by=list(month), FUN=skew)[,2]))
    }
    
    #Annual mean daily values
    
    AnnualMeans[1 + k1, -1, iv, ] <- t(apply(SITE_OBS, c(2), function(x) {aggregate(x, by=list(wyear), FUN=mean)[,2]})) #years, sites
    
    # print(apply(AnnualMeans[1 + k1, -1, iv, ], 1, FUN=mean))
    Stats[1 + k1, -1, 1, iv, 3, 1] <- apply(AnnualMeans[1 + k1, -1, iv, ], 1, FUN=mean)
    Stats[1 + k1, -1, 2, iv, 3, 1] <- apply(AnnualMeans[1 + k1, -1, iv, ], 1, FUN=sd)
    Stats[1 + k1, -1, 3, iv, 3, 1] <- apply(AnnualMeans[1 + k1, -1, iv, ], 1, FUN=skew)
    
    #Annual extreme daily values
    ext <- t(apply(SITE_OBS * ifelse(vars[iv] == "tmin", -1, 1), c(2), function(x) {aggregate(x, by=list(WATER_YEAR_D), FUN=max)[,2]})) #dim(ext) = sites, years; find max for prcp and tmax, but -min for tmin
    AnnualEDV[1 + k1, -1, iv, ] <- t(apply(ext, 1, FUN=function(x) return.level(fevd(x, type="GEV", method="MLE", time.units="year"), return.period=return.periods)))
    AnnualEDV[1 + k1, -1, iv, ] <- ifelse(AnnualEDV[1 + k1, -1, iv, ] > 2 * max(ext), NA, AnnualEDV[1 + k1, -1, iv, ]) #set outliers to NA
    if(vars[iv] == "tmin") AnnualEDV[1 + k1, -1, iv, ] <- -AnnualEDV[1 + k1, -1, iv, ] #set -min to min for tmin
    
    #PRCP mean monthly wet/dry spells
    for(is in 1:num_site){
      Spells[1 + k1, 1 + is, , ] <- t(aggregate(sites_prcp[, is], by=list(MONTH_D), FUN=function(x) {
        temp <- rle(x > 0)
        res <- sapply(list(temp$lengths[temp$values], temp$lengths[!temp$values]), FUN=function(x) c(mean(x), sd(x), skew(x)))
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
#' @param sites_prcp  Observed prcp data as matrix with as many columns as we have sites. Each column represents data of another site 
#' @param tmax data.frame with annual climate timeseries at selected location (YEAR, MONTH, PRCP, TMAX, TMIN, TAVG)
#' @param tmin number of iterations for arima simulation
#' @param num_trials number of runs of weathergenerator, how many datasets were created?
#' @param num_sites number of sites. This should be > 3 to use the validation plots
#' @param year contains the year of the observed data  
#' @param Month contains the year of the observed data  
#' @param vars usually keep on default c("prcp","tmax","tmin"). Only change if you know what you do.
#' @param timesteps usually keep on default c("daily","monthly", "annual"). Only change if you know what you do.
#' @param stats usually keep on default c("mean","stdev","skew"). Only change if you know what you do.
#' @export
#' @examples
#' getSimulatedStatistics(sites_prcp, sites_tmax, sites_tmin, num_trials=1, num_site=3, as.matrix(format(sim$out["DATE"], "%Y")),  as.matrix(sim$out["MONTH"]), vars=c("prcp","tmax","tmin"), timesteps= c("daily","monthly", "annual"), stats=c("mean","stdev","skew"))
#'  

getSimulatedStatistics <- function(num_site, num_trials, myarray, wyear, stats=c("mean","stdev","skew"),vars=c("prcp","tmax","tmin"), timesteps=c("daily","monthly", "annual"), mo=c(1:12), return.periods=c(10,30,100), spells=c("wet","dry")){

  cur_data <- myarray[1,1][[1]]$out
  #Stats = median + realizations, median + sites, statistics, variables, timesteps, months
  Stats <- array(NA, dim=c(1 + NUM_TRIALS, 1 + num_site, length(stats), length(vars), length(timesteps), length(mo)),
                 dimnames=list(NULL, NULL, stats, vars, timesteps, mo))
  #AnnualMeans = median + realizations, median + sites, variables, years
  AnnualMeans <- array(0, dim=c(1 + NUM_TRIALS, 1 + num_site, length(vars), length(unique(wyear))),
                       dimnames=list(NULL, NULL, vars, sort(unique(wyear))))
  #AnnualEDV = median + realizations, median + sites, variables, return periods
  AnnualEDV <- array(NA, dim=c(1 + NUM_TRIALS, 1 + num_site, length(vars), length(return.periods)),
                     dimnames=list(NULL, NULL, vars, return.periods))
  #Spells = median + realizations, median + sites, stats*{dry,wet}, months
  Spells <- array(0, dim=c(1 + NUM_TRIALS, 1 + num_site, length(spells) * length(stats), length(mo)),
                  dimnames=list(NULL, NULL, paste(rep(spells, each=length(stats)), stats, sep="_"), mo))
  
  for (ss in 1:num_site) {
    Sim_Variable <- array(NA,c(dim(cur_data)[1],NUM_TRIALS,  length(vars)))
    
    
    for (k1 in 1:NUM_TRIALS) {
      # dir_ending1 <- paste0("TRIAL_", k1)
      # dir_ending2 <- paste(dir_ending1, ending2, sep="_")
      # filename <- paste("TRIAL",k1,"SITE",ss,sep="_")
      #cur_data <- data.matrix(read.table(file.path(dat_dir, dir_ending1, dir_ending2, filename)))
      # print(file.path(dat_dir, filename))
      #cur_data <- data.matrix(read.table(file.path(dat_dir, filename)))
      cur_data <- data.matrix(myarray[ss,k1][[1]]$out[, c("PRCP","TMAX","TMIN")])
      Sim_Variable[, k1, 1:length(vars)] <- cur_data
      # Sim_Variable[, k1, 1:length(vars)] <- cur_data[, 4:(3+length(vars))]
    }
    
    #######Site Statistics
    #sim_temp <- array(NA, dim=c(NUM_TRIALS, length(stats), length(vars), length(timesteps), length(mo)))
    
    for (k1 in 1:NUM_TRIALS) {
      for(iv in seq_along(vars)){
        if (vars[iv] == "prcp") {
          y <- which(Sim_Variable[, k1, iv]>0)
        } else {
          y <- 1:length(Sim_Variable[, k1, iv])
        }
        #Daily
        Stats[1 + k1, 1 + ss, 1, iv, 1, mo] <- aggregate(Sim_Variable[y, k1, iv], by=list(data.matrix(myarray[ss,k1][[1]]$out[, c("MONTH")])[y]), FUN=mean)[,2]
        Stats[1 + k1, 1 + ss, 2, iv, 1, mo] <- aggregate(Sim_Variable[y, k1, iv], by=list(data.matrix(myarray[ss,k1][[1]]$out[, c("MONTH")])[y]), FUN=sd)[,2]
        Stats[1 + k1, 1 + ss, 3, iv, 1, mo] <- aggregate(Sim_Variable[y, k1, iv], by=list(data.matrix(myarray[ss,k1][[1]]$out[, c("MONTH")])[y]), FUN=skew)[,2]
        #Annual
        #AnnualMeans[1 + k1, 1 + ss, iv, ] <- t(aggregate(Sim_Variable[, k1, iv], by=list(WATER_YEAR_SIM), FUN=mean)[,2])
        AnnualMeans[1 + k1, 1 + ss, iv, ] <- t(aggregate(Sim_Variable[, k1, iv], by=list(data.matrix(myarray[ss,k1][[1]]$out[, c("SIM_YEAR")])), FUN=mean)[,2])
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
          Spells[1 + k1, 1 + ss, , ] <- t(aggregate(Sim_Variable[, k1, iv], by=list(data.matrix(myarray[ss,k1][[1]]$out[, c("MONTH")])), FUN=function(x) {
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
#' @param fig_dir  Directory to create PDF containing the plot
#' @param obs observation data
#' @param sim simulation output
#' @param ftag string to appear in filename
#' @export
#' @examples
#' plot1.basinwide.vars.monthly(fig_dir, obs, sim)
#' 
#'
plot1.basinwide.vars.monthly <- function(fig_dir, obs, sim, ftag=""){
  pdf(width=3*length(vars), height=3*length(stats), file=file.path(fig_dir, paste0("1_BasinWide_Trials_MonthlyVariables", ifelse(nchar(ftag > 0), "_"), ftag, ".pdf")))
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
  dev.off()
}
#' Plot Monthly spells
#'
#' @param fig_dir  Directory to create PDF containing the plot
#' @param obs observation data
#' @param sim simulation output
#' @param ftag string to appear in filename
#' @export
#' @examples
#' plot2.basinwide.spells.monthly(fig_dir, obs, sim)
#' 
#'
plot2.basinwide.spells.monthly <- function(fig_dir, obs, sim, ftag=""){
  pdf(width=3*length(spells), height=3*length(stats), file=file.path(fig_dir, paste0("2_BasinWide_Trials_MonthlySpells", ifelse(nchar(ftag > 0), "_"), ftag, ".pdf")))
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
  dev.off()
}

#' Plot annual extremes
#'
#' @param fig_dir  Directory to create PDF containing the plot
#' @param obs observation data
#' @param sim simulation output
#' @param ftag string to appear in filename
#' @export
#' @examples
#' plot3.basinwide.extreme.annual(fig_dir, obs, sim)
#' 
#'
plot3.basinwide.extreme.annual <- function(fig_dir, obs, sim, ftag=""){
  pdf(width=3*length(vars), height=3*1, file=file.path(fig_dir, paste0("3_BasinWide_Trials_AnnualExtremeEvents", ifelse(nchar(ftag > 0), "_"), ftag, ".pdf")))
  layout(matrix(1:(length(vars)*1), ncol=length(vars), nrow=1, byrow=FALSE))
  
  for(iv in seq_along(vars)){
    x <- sim[-1, 1, iv, ]
    y <- obs[-1, 1, iv, ]
    ylim <- c(min(x, y, na.rm = TRUE), max(x, y, na.rm = TRUE))
    boxplot(x, ylim=ylim, xlab="Return period (years)", ylab=vars[iv])
    points(1:length(return.periods), y, col="red")
  }
  legend(x="bottom", legend=c("Observations", "Simulated trials"), pch=19, col=c("red", "black"))
  dev.off()
}

#' Plot sitespecific variables
#'
#' @param fig_dir  Directory to create PDF containing the plot
#' @param obs observation data
#' @param sim simulation output
#' @param ftag string to appear in filename
#' @export
#' @examples
#' plot4.sitespecific.vars(fig_dir, obs, sim)
#' 
#'
plot4.sitespecific.vars <- function(fig_dir, timestep="daily", obs, sim, ftag=""){
  timestep <- match.arg(timestep, timesteps)
  it <- switch(timestep, daily=1, monthly=2, annual=3)
  mo <- switch(timestep, daily=1:12, monthly=1:12, annual=1)
  
  pdf(width=3*length(vars), height=3*length(stats), file=file.path(fig_dir, paste0("4_Sites_MedianTrials_", timestep, "Variables", ifelse(nchar(ftag > 0), "_"), ftag, ".pdf")))
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
  dev.off()
}

#' Plot sitespecific annual correlations
#'
#' @param fig_dir  Directory to create PDF containing the plot
#' @param obsCors observation data
#' @param simCors simulation output
#' @param ftag string to appear in filename
#' @export
#' @examples
#' plot5.sitespecific.correlations.annual (fig_dir, obsCors, simCors)
#' 
#'
plot5.sitespecific.correlations.annual <- function(fig_dir, obsCors, simCors, ftag="Validation"){
  pdf(width=3*length(vars), height=3*2, file=file.path(fig_dir, paste0("5_Sites_MedianTrials_AnnualCorrelations", ifelse(nchar(ftag > 0), "_"), ftag, ".pdf")))
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
  
  dev.off()
}

#' Plot sitespecific annual spells
#'
#' @param fig_dir  Directory to create PDF containing the plot
#' @param obs observation data
#' @param sim simulation output data
#' @param ftag string to appear in filename
#' @export
#' @examples
#' plot6.sitespecific.spells.annual (fig_dir, obs, sim)
#' 
#'
plot6.sitespecific.spells.annual <- function(fig_dir, obs, sim, ftag="Validation"){
  pdf(width=3*length(spells), height=3*length(stats), file=file.path(fig_dir, paste0("6_Sites_MedianTrials_AnnualSpells", ifelse(nchar(ftag > 0), "_"), ftag, ".pdf")))
  layout(matrix(1:(length(spells)*length(stats)), ncol=length(spells), nrow=length(stats), byrow=FALSE))
  
  for(iv in 1:length(spells)){
    for(is in seq_along(stats)){
      y <- apply(sim[1, -1, (iv - 1) * 2 + is, mo], MARGIN=1, FUN=mean) #mean across months
      x <- apply(obs[1, -1, (iv - 1) * 2 + is, mo], MARGIN=1, FUN=mean)
      print(x)
      print(y)
      xlim <- ylim <- c(min(x, y, na.rm = TRUE), max(x, y, na.rm = TRUE))
      print(xlim)
      print(ylim)
      try(plot(x, y, xlim=xlim, ylim=ylim, xlab="Observed spells (days)", ylab="Simulated spells (days)", main=paste(spells[iv], stats[is])))
      abline(0,1)
    }
  }
  
  dev.off()
}

#' Plot PowerSpectrum
#'
#' @param fig_dir  Directory to create PDF containing the plot
#' @param PowerSpectrum, create this data this using Final_Annual_Sim_All and with its output getLowFrequencyVariability
#' @param ftag string to appear in filename
#' @export
#' @examples
#' plot7.basinwide.PowerSpectruml (fig_dir, PowerSpectrum)
#' 
#'
plot7.basinwide.PowerSpectrum <- function(fig_dir, PowerSpectrum, ftag=""){
  pdf(width=9, height=3, file=file.path(fig_dir, paste0("7_Sites_MedianTrials_PowerSpectrum", ifelse(nchar(ftag > 0), "_"), ftag, ".pdf")))
  plot(PowerSpectrum$Power_Spectrum_Period, PowerSpectrum$Power_Spectrum_Obs,xlab="Period (years)", ylab="Power (mm2)", type="l",col="black",lwd=2)
  polygon(c(PowerSpectrum$Power_Spectrum_Period,rev(PowerSpectrum$Power_Spectrum_Period)),c(PowerSpectrum$Sim_Power_Spectrum_3,rev(PowerSpectrum$Sim_Power_Spectrum_2)),col="grey")
  lines(PowerSpectrum$Power_Spectrum_Period, PowerSpectrum$Sim_Power_Spectrum_1,type="l",col="black",lwd=2)
  lines(PowerSpectrum$Power_Spectrum_Period, PowerSpectrum$Power_Spectrum_Obs,col="red",lwd=2)
  
  dev.off()
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

get_final_annual_sim_all <- function( ar_models, new_ar_model, num_trials=1, num_autocor_annual_changes=1,num_magnitude_noise_annual_changes=1){
  Final_Annual_Sim_All <- array(NA,c(num_year_sim,num_trials,num_autocor_annual_changes,num_magnitude_noise_annual_changes))
  Annual_Sim <- array(NA,c(num_year_sim,length(ar_models)))
  for (k1 in 1:NUM_TRIALS) {
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
#' @param AnnualMeans 
#' @param vars
#' @param ftag string to appear in filename
#' @export
#' @examples
#' get_final_annual_sim_all(NUM_AUTOCOR_ANNUAL_CHANGES,NUM_MAGNITUDE_NOISE_ANNUAL_CHANGES,AR_MODELS,NEW_AR_MODEL)
#' 
#'
getCorrelations <- function(AnnualMeans, vars= c("prcp","tmax","tmin"), num_trials=1){
  #NUM_TRIALS <- dim(AnnualMeans)[1]
  
  #Intersites correlations of annual series
  intersites <- t(combn(dim(AnnualMeans)[2], 2))
  cors <- array(NA, c(1 + NUM_TRIALS, nrow(intersites), length(vars)), dimnames=list(NULL, NULL, vars))
  
  for(k1 in 1:num_trials){
    for(is in 1:nrow(intersites)){
      for(iv in 1:length(vars)){
        cors[1 + k1, is, iv] <- cor(AnnualMeans[k1, intersites[is, 1], iv, ], AnnualMeans[k1, intersites[is, 2], iv, ], method="spearman")
      }
    }
  }
  cors[1, , ] <- apply(cors[-1, , , drop=FALSE], MARGIN=c(2, 3), FUN=median) #Median across trials
  intersites <- cbind(intersites, cors[1, , , drop=TRUE])
  
  #Lag-1 autocorrelation of annual series
  lag1auto <- array(NA, c(1 + NUM_TRIALS, dim(AnnualMeans)[2], length(vars)), dimnames=list(NULL, NULL, vars))
  for(k1 in 1:num_trials){
    for(is in 1:dim(AnnualMeans)[2]){
      for(iv in 1:length(vars)){
        lag1auto[1 + k1, is, iv] <- acf(AnnualMeans[k1, is, iv, ], type="correlation", plot=FALSE)[1, ]$acf
        stopifnot(!is.na(lag1auto[1 + k1, is, iv]))
      }
    }
  }
  lag1auto[1, , ] <- apply(lag1auto[-1, , , drop=FALSE], MARGIN=c(2, 3), FUN=median) #Median across trials
  
  return(list(intersites=intersites, lag1auto=lag1auto[1, , , drop=TRUE]))
}
