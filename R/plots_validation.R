

#' Helper functions implementations
#'
#' @param 
#' @param 
#' @param 
#' @param 
#' @param 
#' @param 
#' @export
#' @examples
#' getObservedStatistics3(wgen_out_array)
#########################################################################################################################################
sceDefaults <- function()
  list(ncomplex = 5, ## number of complexes
       cce.iter = NA, ## number of iteration in inner loop (CCE algorithm)
       fnscale = 1, ## function scaling factor (set to -1 for maximisation)
       elitism = 1, ## controls amount of weighting in sampling towards the better parameter sets
       initsample = "latin", ## sampling scheme for initial values -- "latin" or "random"
       reltol = 1e-5, ## convergence threshold: relative improvement factor required in an SCE iteration
       tolsteps = 7, ## number of iterations within reltol to confirm convergence
       maxit = 10000, ## maximum number of iterations
       maxeval = Inf, ## maximum number of function evaluations
       maxtime = Inf, ## maximum duration of optimization in seconds
       returnpop = FALSE, ## whether to return populations from all iterations
       trace = 0, ## level of user feedback
       REPORT = 1) ## number of iterations between reports when trace >= 1


SCEoptim <- function(FUN, par, ...,
                     lower = -Inf, upper = Inf,
                     control = list())
{
  FUN <- match.fun(FUN)
  stopifnot(is.numeric(par))
  stopifnot(length(par) > 0)
  stopifnot(is.numeric(lower))
  stopifnot(is.numeric(upper))
  ## allow `lower` or `upper` to apply to all parameters
  if (length(lower) == 1)
    lower <- rep(lower, length = length(par))
  if (length(upper) == 1)
    upper <- rep(upper, length = length(par))
  stopifnot(length(lower) == length(par))
  stopifnot(length(upper) == length(par))
  
  ## determine number of variables to be optimized
  NDIM <- length(par)
  
  ## update default options with supplied options
  stopifnot(is.list(control))
  control <- modifyList(sceDefaults(), control)
  isValid <- names(control) %in% names(sceDefaults())
  if (any(!isValid))
    stop("unrecognised options: ",
         toString(names(control)[!isValid]))
  
  returnpop <- control$returnpop
  trace <- control$trace
  nCOMPLEXES <- control$ncomplex
  CCEITER <- control$cce.iter
  MAXIT <- control$maxit
  MAXEVAL <- control$maxeval
  
  ## recommended number of CCE steps in Duan et al 1994:
  if (is.na(CCEITER))
    CCEITER <- 2 * NDIM + 1
  
  if (is.finite(MAXEVAL)) {
    ## upper bound on number of iterations to reach MAXEVAL
    MAXIT <- min(MAXIT, ceiling(MAXEVAL / (nCOMPLEXES * CCEITER)))
  }
  
  ## define number of points in each complex
  nPOINTS_COMPLEX <- 2 * NDIM + 1
  
  ## define number of points in each simplex
  nPOINTS_SIMPLEX <- NDIM+1
  
  ## define total number of points
  nPOINTS <- nCOMPLEXES * nPOINTS_COMPLEX
  
  ## initialize counters
  funevals <- 0
  
  
  costFunction <- function(FUN, par, ...)
  {
    ## check lower and upper bounds
    i <- which(par < lower)
    if (any(i)) {
      i <- i[1]
      return( 1e12 + (lower[i] - par[i]) * 1e6 )
    }
    i <- which(par > upper)
    if (any(i)) {
      i <- i[1]
      return( 1e12 + (par[i] - upper[i]) * 1e6 )
    }
    funevals <<- funevals + 1
    result <- FUN(par, ...) * control$fnscale
    if (is.na(result))
      result <- 1e12
    result
  }
  
  simplexStep <- function(P, FAC)
  {
    ## Extrapolates by a factor FAC through the face of the simplex across from
    ## the highest (i.e. worst) point.
    worst <- nPOINTS_SIMPLEX
    centr <- apply(P[-worst,,drop=FALSE], 2, mean)
    newpar <- centr*(1-FAC) + P[worst,]*FAC
    newpar
  }
  
  
  ## initialize population matrix
  POPULATION <- matrix(as.numeric(NA), nrow = nPOINTS, ncol = NDIM)
  if (!is.null(names(par)))
    colnames(POPULATION) <- names(par)
  POP.FITNESS <- numeric(length = nPOINTS)
  POPULATION[1,] <- par
  
  ## generate initial parameter values by random uniform sampling
  finitelower <- ifelse(is.infinite(lower), -(abs(par)+2)*5, lower)
  finiteupper <- ifelse(is.infinite(upper), +(abs(par)+2)*5, upper)
  if (control$initsample == "latin") {
    for (i in 1:NDIM) {
      tmp <- seq(finitelower[i], finiteupper[i], length = nPOINTS-1)
      tmp <- jitter(tmp, factor = 2)
      tmp <- pmax(finitelower[i], pmin(finiteupper[i], tmp))
      POPULATION[-1,i] <- sample(tmp)
    }
  } else {
    for (i in 1:NDIM)
      POPULATION[-1,i] <- runif(nPOINTS-1, finitelower[i], finiteupper[i])
  }
  
  ## only store all iterations if requested -- could be big!
  if (!is.finite(MAXIT)) {
    MAXIT <- 10000
    warning("setting maximum iterations to 10000")
  }
  if (returnpop) {
    POP.ALL <- array(as.numeric(NA), dim = c(nPOINTS, NDIM, MAXIT))
    if (!is.null(names(par)))
      dimnames(POP.ALL)[[2]] <- names(par)
  }
  POP.FIT.ALL <- matrix(as.numeric(NA), ncol = nPOINTS, nrow = MAXIT)
  BESTMEM.ALL <- matrix(as.numeric(NA), ncol = NDIM, nrow = MAXIT)
  if (!is.null(names(par)))
    colnames(BESTMEM.ALL) <- names(par)
  
  ## the output object
  obj <- list()
  class(obj) <- c("SCEoptim", class(obj))
  obj$call <- match.call()
  obj$control <- control
  
  EXITFLAG <- NA
  EXITMSG <- NULL
  
  ## initialize timer
  tic <- as.numeric(Sys.time())
  toc <- 0
  
  ## calculate cost for each point in initial population
  for (i in 1:nPOINTS)
    POP.FITNESS[i] <- costFunction(FUN, POPULATION[i,], ...)
  
  ## sort the population in order of increasing function values
  idx <- order(POP.FITNESS)
  POP.FITNESS <- POP.FITNESS[idx]
  POPULATION <- POPULATION[idx,,drop=FALSE]
  
  ## store one previous iteration only
  POP.PREV <- POPULATION
  POP.FIT.PREV <- POP.FITNESS
  
  if (returnpop) {
    POP.ALL[,,1] <- POPULATION
  }
  POP.FIT.ALL[1,] <- POP.FITNESS
  BESTMEM.ALL[1,] <- POPULATION[1,]
  
  ## store best solution from last two iterations
  prevBestVals <- rep(Inf, control$tolsteps)
  prevBestVals[1] <- POP.FITNESS[1]
  
  ## for each iteration...
  i <- 0
  while (i < MAXIT) {
    
    i <- i + 1
    
    ## The population matrix POPULATION will now be rearranged into complexes.
    
    ## For each complex ...
    for (j in 1:nCOMPLEXES) {
      
      ## construct j-th complex from POPULATION
      
      k1 <- 1:nPOINTS_COMPLEX
      k2 <- (k1-1) * nCOMPLEXES + j
      
      COMPLEX <- POP.PREV[k2,,drop=FALSE]
      COMPLEX_FITNESS <- POP.FIT.PREV[k2]
      
      ## Each complex evolves a number of steps according to the competitive
      ## complex evolution (CCE) algorithm as described in Duan et al. (1992).
      ## Therefore, a number of 'parents' are selected from each complex which
      ## form a simplex. The selection of the parents is done so that the better
      ## points in the complex have a higher probability to be selected as a
      ## parent. The paper of Duan et al. (1992) describes how a trapezoidal
      ## probability distribution can be used for this purpose.
      
      for (k in 1:CCEITER) {
        
        ## select simplex by sampling the complex
        
        ## sample points with "trapezoidal" i.e. linear probability
        weights <- rev(ppoints(nPOINTS_COMPLEX))
        ## 'elitism' parameter can give more weight to the better results:
        weights <- weights ^ control$elitism
        LOCATION <- sample(seq(1,nPOINTS_COMPLEX), size = nPOINTS_SIMPLEX,
                           prob = weights)
        
        LOCATION <- sort(LOCATION)
        
        ## construct the simplex
        SIMPLEX <- COMPLEX[LOCATION,,drop=FALSE]
        SIMPLEX_FITNESS <- COMPLEX_FITNESS[LOCATION]
        
        worst <- nPOINTS_SIMPLEX
        
        ## generate new point for simplex
        
        ## first extrapolate by a factor -1 through the face of the simplex
        ## across from the high point,i.e.,reflect the simplex from the high point
        parRef <- simplexStep(SIMPLEX, FAC = -1)
        fitRef <- costFunction(FUN, parRef, ...)
        
        ## check the result
        if (fitRef <= SIMPLEX_FITNESS[1]) {
          ## gives a result better than the best point,so try an additional
          ## extrapolation by a factor 2
          parRefEx <- simplexStep(SIMPLEX, FAC = -2)
          fitRefEx <- costFunction(FUN, parRefEx, ...)
          if (fitRefEx < fitRef) {
            SIMPLEX[worst,] <- parRefEx
            SIMPLEX_FITNESS[worst] <- fitRefEx
            ALGOSTEP <- 'reflection and expansion'
          } else {
            SIMPLEX[worst,] <- parRef
            SIMPLEX_FITNESS[worst] <- fitRef
            ALGOSTEP <- 'reflection'
          }
        } else if (fitRef >= SIMPLEX_FITNESS[worst-1]) {
          ## the reflected point is worse than the second-highest, so look
          ## for an intermediate lower point, i.e., do a one-dimensional
          ## contraction
          parCon <- simplexStep(SIMPLEX, FAC = -0.5)
          fitCon <- costFunction(FUN, parCon, ...)
          if (fitCon < SIMPLEX_FITNESS[worst]) {
            SIMPLEX[worst,] <- parCon
            SIMPLEX_FITNESS[worst] <- fitCon
            ALGOSTEP <- 'one dimensional contraction'
          } else {
            ## can't seem to get rid of that high point, so better contract
            ## around the lowest (best) point
            SIMPLEX <- (SIMPLEX + rep(SIMPLEX[1,], each=nPOINTS_SIMPLEX)) / 2
            for (k in 2:NDIM)
              SIMPLEX_FITNESS[k] <- costFunction(FUN, SIMPLEX[k,], ...)
            ALGOSTEP <- 'multiple contraction'
          }
        } else {
          ## if better than second-highest point, use this point
          SIMPLEX[worst,] <- parRef
          SIMPLEX_FITNESS[worst] <- fitRef
          ALGOSTEP <- 'reflection'
        }
        
        if (trace >= 3) {
          message(ALGOSTEP)
        }
        
        ## replace the simplex into the complex
        COMPLEX[LOCATION,] <- SIMPLEX
        COMPLEX_FITNESS[LOCATION] <- SIMPLEX_FITNESS
        
        ## sort the complex
        idx <- order(COMPLEX_FITNESS)
        COMPLEX_FITNESS <- COMPLEX_FITNESS[idx]
        COMPLEX <- COMPLEX[idx,,drop=FALSE]
      }
      
      ## replace the complex back into the population
      POPULATION[k2,] <- COMPLEX
      POP.FITNESS[k2] <- COMPLEX_FITNESS
    }
    
    ## At this point, the population was divided in several complexes, each of which
    ## underwent a number of iteration of the simplex (Metropolis) algorithm. Now,
    ## the points in the population are sorted, the termination criteria are checked
    ## and output is given on the screen if requested.
    
    ## sort the population
    idx <- order(POP.FITNESS)
    POP.FITNESS <- POP.FITNESS[idx]
    POPULATION <- POPULATION[idx,,drop=FALSE]
    if (returnpop) {
      POP.ALL[,,i] <- POPULATION
    }
    POP.FIT.ALL[i,] <- POP.FITNESS
    BESTMEM.ALL[i,] <- POPULATION[1,]
    
    curBest <- POP.FITNESS[1]
    
    ## end the optimization if one of the stopping criteria is met
    
    prevBestVals <- c(curBest, head(prevBestVals, -1))
    reltol <- control$reltol
    if (all(abs(diff(prevBestVals)) <= reltol * (abs(curBest)+reltol))) {
      EXITMSG <- 'Change in solution over [tolsteps] less than specified tolerance (reltol).'
      EXITFLAG <- 0
    }
    
    ## give user feedback on screen if requested
    if (trace >= 1) {
      if (i == 1) {
        message(' Nr Iter Nr Fun Eval Current best function Current worst function')
      }
      if ((i %% control$REPORT == 0) || (!is.na(EXITFLAG)))
      {
        message(sprintf(' %5.0f %5.0f %12.6g %12.6g',
                        i, funevals, min(POP.FITNESS), max(POP.FITNESS)))
        if (trace >= 2)
          message("parameters: ", toString(signif(POPULATION[1,], 3)))
      }
    }
    
    if (!is.na(EXITFLAG))
      break
    
    if ((i >= control$maxit) || (funevals >= control$maxeval)) {
      EXITMSG <- 'Maximum number of function evaluations or iterations reached.'
      EXITFLAG <- 1
      break
    }
    
    toc <- as.numeric(Sys.time()) - tic
    if (toc > control$maxtime) {
      EXITMSG <- 'Exceeded maximum time.'
      EXITFLAG <- 2
      break
    }
    
    ## go to next iteration
    POP.PREV <- POPULATION
    POP.FIT.PREV <- POP.FITNESS
  }
  if (trace >= 1)
    message(EXITMSG)
  
  ## return solution
  obj$par <- POPULATION[1,]
  obj$value <- POP.FITNESS[1]
  obj$convergence <- EXITFLAG
  obj$message <- EXITMSG
  
  ## store number of function evaluations
  obj$counts <- funevals
  ## store number of iterations
  obj$iterations <- i
  ## store the amount of time taken
  obj$time <- toc
  
  if (returnpop) {
    ## store information on the population at each iteration
    obj$POP.ALL <- POP.ALL[,,1:i]
    dimnames(obj$POP.ALL)[[3]] <- paste("iteration", 1:i)
  }
  obj$POP.FIT.ALL <- POP.FIT.ALL[1:i,]
  obj$BESTMEM.ALL <- BESTMEM.ALL[1:i,]
  
  obj
}



#' Helper function for power spectrum
#'
#' @param 
#' @param 
#' @param 
#' @param 
#' @param 
#' @param 
#' @export
#' @examples
#' getObservedStatistics3(wgen_out_array)

var_cov_arma <- function(parameters,num_ar,num_ma) {
  
  e_var <- parameters[1]
  ar <-NULL
  ma <- NULL
  if (num_ar>0) {ar <- parameters[2:(2+num_ar-1)]}
  if (num_ma>0) {ma <- parameters[(num_ar+2):(1+num_ar+num_ma)]}
  
  #test for nonstationarity
  nonstationary <- FALSE
  if (length(ar)==1) {
    if(abs(ar[1])>=1) {nonstationary <- TRUE}
  } else if (length(ar==2)) {
    a1 <- ar[1]
    a2 <- ar[2]
    d <- a1^2 + 4*a2
    if(d>0) {
      test1 <- a1 + a2
      test2 <- a2-a1
      if (test1>=1 | test2>=1) {
        nonstationary<-TRUE
      }
    } else if (d==0) {
      test1 <- abs(a1)
      if (test1>=2) {
        nonstationary<-TRUE
      }		
    } else {
      test1 <- -a2
      if (test1>=1) {
        
        nonstationary<-TRUE
      }				
    }
  }
  
  a1 <- NULL
  a2 <- NULL
  b1 <- NULL
  b2 <- NULL
  #ar(0)
  if (length(ar)==0 & length(ma)==0 & nonstationary==FALSE) {
    var <- e_var
    lag1 <- 0
    
    #ar(1)
  } else if (length(ar)==1 & length(ma)==0 & nonstationary==FALSE) {
    a1 <- ar[1]
    var <- e_var/(1-a1^2)
    lag1 <- a1
    
    #ar(2)
  } else if (length(ar)==2 & length(ma)==0 & nonstationary==FALSE) {
    a1 <- ar[1]
    a2 <- ar[2]
    var <- ((1-a2)/(1+a2))*e_var/(a1+a2-1)/(a2-a1-1)
    lag1 <- a1/(1-a2)
    
    #ma(1)
  } else if (length(ar)==0 & length(ma)==1 & nonstationary==FALSE) {
    b1 <- ma[1]
    var <- (1+b1^2)*e_var
    lag1 <- b1/(1+b1^2)
    
    #ma(2)
  } else if (length(ar)==0 & length(ma)==2 & nonstationary==FALSE) {
    b1 <- ma[1]
    b2 <- ma[2]
    var <- e_var*(1+b1^2+b2^2)
    lag1 <- b1*(1+b2)/(1+b1^2+b2^2)
    
    #arma(1,1)
  } else if (length(ar)==1 & length(ma)==1 & nonstationary==FALSE) {
    a1 <- ar[1]
    b1 <- ma[1]
    var <- (1+b1^2+2*a1*b1)/(1-a1^2)*e_var
    lag1 <- (1+a1*b1)*(a1+b1)/(1+b1^2+2*a1*b1)
    
    #arma(2,1)
  } else if (length(ar)==2 & length(ma)==1 & nonstationary==FALSE) {
    a1 <- ar[1]
    a2 <- ar[2]
    b1 <- ma[1]
    var <- e_var*(a1*b1+a1*a2*b1+(1-a2)*(b1*(a1+b1)+1))/(1-a2-a1^2-a1^2*a2-a2^2*(1-a2))
    lag1 <- (a1 + (1/var)*b1*e_var)/(1-a2)
    
    #arma(1,2)
  } else if (length(ar)==1 & length(ma)==2 & nonstationary==FALSE) {
    a1 <- ar[1]
    b1 <- ma[1]
    b2 <- ma[2]
    var <- e_var*(a1*b1+b2*a1^2+b2*b1*a1+1+b1*a1+b1^2+b2^2+b2*a1^2+b2*a1*b1)/(1-a1^2)
    lag1 <- a1 + e_var*(b1+b2*a1+b2*b1)/var
    
    #arma(2,2)
  } else if (length(ar)==2 & length(ma)==2 & nonstationary==FALSE) {
    a1 <- ar[1]
    a2 <- ar[2]
    b1 <- ma[1]
    b2 <- ma[2]
    x <- 1-a2
    var <- e_var*(a1*b1 + a1^2*b2 + a1*b1*b2 + a1*a2*b1 + a1^2*a2*b2 + a1*a2*b1*b2 + x*(a2*b2 + 1 + b1*a1 + b1^2 + b2^2 + a2*b2 + b2*a1^2 + b2*a1*b1)) / (1 - a2 - a1^2 - a2*a1^2 - a2^2*(1-a2))
    lag1 <- (a1*var + b1*e_var + b2*a1*e_var + b2*b1*e_var)/(var*(1-a2))
    
  } else if (nonstationary==FALSE) {
    if(.platform$os.type == "windows"){
      windialog(type="ok",message="we only accept models with order less than or equal to arma(2,2)")
    } else {
      if(capabilities()["tcltk"]){
        library("tcltk")
        tkmessagebox(type="ok",message="we only accept models with order less than or equal to arma(2,2)")
      } else {
        temp <- readline("we only accept models with order less than or equal to arma(2,2)")
      }
    }
  } else {
    var <- 9999999999999999999 
    lag1 <- 999999999999999999
  }
  
  
  return(list(var,lag1,e_var,a1,a2,b1,b2,nonstationary))
}






#' Helper function for power spectrum
#'
#' @param 
#' @param 
#' @param 
#' @param 
#' @param 
#' @param 
#' @export
#' @examples
#' getObservedStatistics3(wgen_out_array)

#2014/07/09;	drs;	reformatted line breaks

var_cov_arma_alteration<- function(parameters,num_ar,num_ma,var_target,lag1_target,e_var_org,ar_org,ma_org,other_e_var,other_ar,other_ma,cur_dir) {
  #etwd(cur_dir)
  #source("variance_autocorrelation_arma.r")
  #get variance and lag 1 autocorrelation of noise model
  result <- var_cov_arma(parameters,num_ar,num_ma)
  var <- result[[1]]
  lag1 <- result[[2]]
  e_var <-result[[3]]
  a1 <- result[[4]]
  a2 <- result[[5]]
  b1 <- result[[6]]
  b2 <- result[[7]]
  nonstationary <- result[[8]]
  final_var <- var
  final_lag1 <- lag1
  #get variance and lag 1 autocorrelation of additional models if they exist
  num_low_freq_models <- length(other_e_var)
  if (num_low_freq_models>0) {
    var_other <- array(NA,num_low_freq_models)
    lag1_other <- array(NA,num_low_freq_models)
    lag1_other_weighted <- array(NA,num_low_freq_models)
    for (j in 1:num_low_freq_models) {
      cur_ar <- other_ar[[j]]
      cur_ma <- other_ma[[j]]
      cur_parameters <- c(other_e_var[[j]],cur_ar,cur_ma)
      cur_num_ar <- length(cur_ar)
      cur_num_ma <- length(cur_ma)
      cur_result <- var_cov_arma(cur_parameters,cur_num_ar,cur_num_ma)
      var_other[j] <- cur_result[[1]]
      lag1_other[j] <- cur_result[[2]]
      lag1_other_weighted[j] <- cur_result[[2]]*cur_result[[1]]
    }
    final_var <- sum(var_other,var)
    final_lag1 <- sum(lag1_other_weighted,(lag1*var))/final_var
  }
  e_dif <- e_var - e_var_org
  a1_dif <- 0 ; a2_dif <- 0; b1_dif <- 0; b2_dif <- 0
  if (num_ar==1 & nonstationary==FALSE) {a1_dif <- a1 - ar_org[1]}
  if (num_ar==2 & nonstationary==FALSE) {a1_dif <- a1 - ar_org[1] ; a2_dif <- a2 - ar_org[2]}
  if (num_ma==1 & nonstationary==FALSE) {b1_dif <- b1 - ma_org[1]}
  if (num_ma==2 & nonstationary==FALSE) {b1_dif <- b1 - ma_org[1] ; b2_dif <- b2 - ma_org[2]}
  penalty <- abs(e_dif) + abs(a1_dif) + 10*abs(b1_dif) + 100*abs(a2_dif) + 1000*abs(b2_dif)
  distance <- 100000000*sqrt(((var_target - final_var)/var_target)^2 + ((lag1_target-final_lag1)/lag1_target)^2)
  return(distance + penalty)
}


adjust_warm_model <- function (stdev_change,lag1_change,e_var_org,ar_org,ma_org,other_e_var,other_ar,other_ma,cur_dir) {
  #e_var_org: the variance for the noise in the arma model to be adjusted
  #ar_org: a vector of ar coefficients for the arma model to be adjusted
  #ma_org: a vector of ma coefficients for the arma model to be adjusted
  #other_e_var: a list of the variances for the noise in low frequency arma models. these are not adjusted. 
  #other_ar: ar components of low frequency arma models. this is a list of vectors. these are not adjusted.
  #other_ma: ma components of low frequency arma models. this is a list of vectors. these are not adjusted.
  #cur_dir: directory with adjustment functions
  
  #setwd(cur_dir)
  #source("variance_autocorrelation_arma_alteration.r")
  #source("variance_autocorrelation_arma.r")
  #source("sce_algorithm.r")
  
  num_ar <- length(ar_org)
  num_ma <- length(ma_org)
  e_var_min <- 0.001*e_var_org
  e_var_max <- 100*e_var_org
  parameters <- c(e_var_org,ar_org,ma_org)
  result <- var_cov_arma(parameters,num_ar,num_ma)
  var <- result[[1]]
  lag1 <- result[[2]]
  final_var <- var
  final_lag1 <- lag1
  
  #get variance and lag 1 autocorrelation of additional models if they exist
  num_low_freq_models <- length(other_e_var)
  var_other <- 0
  lag1_other_weighted <- 0
  if (num_low_freq_models>0) {
    var_other <- array(NA,num_low_freq_models)
    lag1_other <- array(NA,num_low_freq_models)
    lag1_other_weighted <- array(NA,num_low_freq_models)
    for (j in 1:num_low_freq_models) {
      cur_ar <- other_ar[[j]]
      cur_ma <- other_ma[[j]]
      cur_parameters <- c(other_e_var[[j]],cur_ar,cur_ma)
      cur_num_ar <- length(cur_ar)
      cur_num_ma <- length(cur_ma)
      cur_result <- var_cov_arma(cur_parameters,cur_num_ar,cur_num_ma)
      var_other[j] <- cur_result[[1]]
      lag1_other[j] <- cur_result[[2]]
      lag1_other_weighted[j] <- cur_result[[2]]*cur_result[[1]]		
    }
    final_var <- sum(var_other,var)
    final_lag1 <- sum(lag1_other_weighted,(lag1*var))/final_var		
  }
  
  var_target <- (stdev_change*sqrt(final_var))^2
  lag1_target <- max(min(lag1_change*final_lag1,.99),-.99)
  
  if (length(ar_org)==0 & length(ma_org)==0) {
    start_par <- c(e_var_org)	#e_var, a1
    lowerb <- c(e_var_min)
    upperb <- c(e_var_max)
  } else if (length(ar_org)==1 & length(ma_org)==0) {
    start_par <- c(e_var_org, ar_org[1])	#e_var, a1
    lowerb <- c(e_var_min,-1)
    upperb <- c(e_var_max,1)
  } else if (length(ar_org)==2 & length(ma_org)==0) {
    start_par <- c(e_var_org, ar_org[1], ar_org[2])	 #e_var, a1, a2
    lowerb <- c(e_var_min,-1,-1)
    upperb <- c(e_var_max,1,1)	
  } else if (length(ar_org)==0 & length(ma_org)==1) {
    start_par <- c(e_var_org,ma_org[1])	#e_var, b1
    lowerb <- c(e_var_min,-1)
    upperb <- c(e_var_max,1)  	
  } else if (length(ar_org)==0 & length(ma_org)==2) {
    start_par <- c(e_var_org,ma_org[1],ma_org[2])	#e_var, b1, b2
    lowerb <- c(e_var_min,-1,-1)
    upperb <- c(e_var_max,1,1)  	
  } else if (length(ar_org)==1 & length(ma_org)==1) {
    start_par <- c(e_var_org,ar_org[1],ma_org[1])	#e_var, a1, b1
    lowerb <- c(e_var_min,-1,-1)
    upperb <- c(e_var_max,1,1)  	
  } else if (length(ar_org)==2 & length(ma_org)==1) {
    start_par <- c(e_var_org,ar_org[1],ar_org[2],ma_org[1])	#e_var, a1, a2, b1
    lowerb <- c(e_var_min,-1,-1,-1)
    upperb <- c(e_var_max,1,1,1)	
  } else if (length(ar_org)==1 & length(ma_org)==2) {	
    start_par <- c(e_var_org,ar_org[1],ma_org[1],ma_org[2])	#e_var, a1, b1, b2
    lowerb <- c(e_var_min,-1,-1,-1)
    upperb <- c(e_var_max,1,1,1)
  } else if (length(ar_org)==2 & length(ma_org)==2) {
    start_par <- c(e_var_org,ar_org[1],ar_org[2],ma_org[1],ma_org[2])	#e_var, a1, a2, b1, b2
    lowerb <- c(e_var_min,-1,-1,-1,-1)
    upperb <- c(e_var_max,1,1,1,1)
  } else {
    if(.platform$os.type == "windows"){
      windialog(type="ok",message="we only accept models with order less than or equal to arma(2,2)")
    } else {
      if(capabilities()["tcltk"]){
        library("tcltk")
        tkmessagebox(type="ok",message="we only accept models with order less than or equal to arma(2,2)")
      } else {
        temp <- readline("we only accept models with order less than or equal to arma(2,2)")
      }
    }
  }
  
  
  result <- SCEoptim(var_cov_arma_alteration,start_par,num_ar,num_ma,var_target,lag1_target,e_var_org,ar_org,ma_org,other_e_var,other_ar, other_ma,cur_dir,lower=lowerb,upper=upperb)
  new_parameters <- result$par
  result <- var_cov_arma(new_parameters,num_ar,num_ma)
  var_update <- sum(var_other,result[[1]])
  lag1_update <- sum(lag1_other_weighted,(result[[2]]*result[[1]]))/var_update		
  
  percent_deviation_var <- round((var_target - var_update)*100/var_target,2)
  percent_deviation_lag1 <- round((lag1_target - lag1_update)*100/var_target,2)
  if (percent_deviation_var>5) {print("warning: adjusted annual variance is >5% different than target annual variance")}
  if (percent_deviation_lag1>5) {print("warning: adjusted annual lag-1 autocorrelation is >5% different than target annual lag-1 autocorrelation")}
  
  
  e_var_final <- new_parameters[1]
  ar_final <-NULL
  ma_final <- NULL
  if (num_ar>0) {ar_final <- new_parameters[2:(2+num_ar-1)]}
  if (num_ma>0) {ma_final <- new_parameters[(num_ar+2):(1+num_ar+num_ma)]}
  
  if (num_ar>0 & num_ma==0) {final <- list(e_var_final,ar_final)}
  if (num_ar==0 & num_ma>0) {final <- list(e_var_final,ma_final)}
  if (num_ar>0 & num_ma>0) {final <- list(e_var_final,ar_final,ma_final)}
  if (num_ar==0 & num_ma==0) {final <- list(e_var_final)}  
  
  return(final)
}



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
plot_basinwide_powerspectrum <- function(annual_prcp, years, simlength, fig_dir=NULL, ftag=""){
  Final_Annual_Sim_All <- get_final_annual_sim_all( annual_prcp, years, simlength )
  PowerSpectrum <- getLowFrequencyVariability(annual_prcp,years, Final_Annual_Sim_All)
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

#' Plot PowerSpectrum Manual, meaning you have to create powerspectrum data using get_final_annual_sim_all and getLowFrequencyVariability
#' there is plot_basinwide_powerspectrum available which will call those functions
#'
#' @param PowerSpectrum, create this data this using Final_Annual_Sim_All and with its output getLowFrequencyVariability
#' @param fig_dir  Directory to create PDF containing the plot
#' @param ftag string to appear in filename
#' @export
#' @examples
#' plot7.basinwide.PowerSpectruml (fig_dir, PowerSpectrum)
#' 
#'
plot_basinwide_powerspectrum_man <- function(PowerSpectrum, fig_dir=NULL, ftag=""){
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

get_new_ar_models <- function(ar_models, Wavelet_Decomposition,  use_warm_model = TRUE) {
  AUTOCOR_ANNUAL_CHANGES <- c(1)
  
  #Change in Standard Deviation of Annual Series (multiplicative)
  # - makes the annual time series more variable, but will not change the persistence
  # - use for ???Inter-annual severity of dry/wet years???
  MAGNITUDE_NOISE_ANNUAL_CHANGES <- c(.5,1,2)
  NUM_AUTOCOR_ANNUAL_CHANGES <- length(AUTOCOR_ANNUAL_CHANGES)
  NUM_MAGNITUDE_NOISE_ANNUAL_CHANGES <-length(MAGNITUDE_NOISE_ANNUAL_CHANGES)
  new_ar_model <- list()
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
      if (use_warm_model) {
        for (l in 2:dim(Wavelet_Decomposition)[2]) {
          CUR_MODEL <- ar_models[[l]]
          Other_e_var[[l-1]] <- CUR_MODEL$sigma2
          Other_AR[[l-1]] <- as.vector(CUR_MODEL$coef)[which(names(CUR_MODEL$coef)!="intercept" & substr(names(CUR_MODEL$coef),1,2)!="ma")]
          Other_MA[[l-1]] <- as.vector(CUR_MODEL$coef)[which(names(CUR_MODEL$coef)!="intercept" & substr(names(CUR_MODEL$coef),1,2)!="ar")]	
        }
      }
      new_parameters <- adjust_warm_model(Stdev_Change,Lag1_Change,e_var_org,AR_org,MA_org,Other_e_var,Other_AR,Other_MA,cur_dir)
      new_ar_model[[count]] <- ar_models[[1]]
      new_ar_model[[count]]$sigma2 <- new_parameters[[1]]
      non_intercept <- which(names(ar_models[[1]]$coef)!="intercept")
      if (length(non_intercept)>0) {
        for (cur_par in non_intercept) {
          new_ar_model[[count]]$coef[[cur_par]] <- new_parameters[[2]][cur_par]
        }
      }
    }
  }
  return(new_ar_model) 
}

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
get_final_annual_sim_all <- function( annual_prcp, years,num_year_sim, num_trials=1, num_autocor_annual_changes=1,num_magnitude_noise_annual_changes=1, sig.level =0.9, n.periods=1, sig.periods=c(8,9,10,11,12,13,14,15), n.comp.periods=8){
  if (any(n.periods >= length(annual_prcp)) || max(n.comp.periods) > length(sig.periods))
    stop("any(n.periods >= length(annual_prcp)) || max(n.comp.periods) > length(sig.periods) is true, simulating more years may help")
  
  Wavelet_Decomposition  <- wavelet_components(annual_prcp,wavelet_analysis(annual_prcp,years,sig.level, "white"), n.periods, sig.periods,n.comp.periods)
  ar_models <- get_ar_models(annual_prcp, years, Wavelet_Decomposition)
  new_ar_model <- get_new_ar_models(ar_models, Wavelet_Decomposition)
  Final_Annual_Sim_All <- array(NA,c(num_year_sim,num_trials,num_autocor_annual_changes,num_magnitude_noise_annual_changes))
  Annual_Sim <- array(NA,c(num_year_sim,length(ar_models)))
  for (k1 in 1:num_trials) {
    count <- 0	
    for (k2 in 1:num_magnitude_noise_annual_changes) {
      for (k3 in 1:num_magnitude_noise_annual_changes) {
        count <- count + 1
        for (l in 1:length(ar_models)) {
          CUR_MODEL <- ar_models[[l]]
          if (l==1) {CUR_MODEL <- new_ar_model[[count]]}
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
    # WAVELET_ANALYSIS returns 3 lists, wavelet_analysis returns much more data.
    # wavelet_analysis$gws == WAVELET_ANALYSIS[[1]]
    # wavelet_analysis$period == WAVELET_ANALYSIS[[3]]
    # wavelet_analysis$ == WAVELET_ANALYSIS[[2]]    NOT FOUND.... maybe some additional calculation
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
