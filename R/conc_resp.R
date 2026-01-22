library(numDeriv)

fitKd = function(conc, resp, bidirectional = FALSE, verbose = FALSE, nofit = FALSE, force_sd = F){
  
  logc = log10(conc)
  fenv <- environment()
  
  pars <- paste0(c("tp", "ga", "b", "er"))
  sds <- paste0(c("tp", "ga", "b", "er"), "_sd")
  myparams = c("success", "aic", "cov", "rme", "modl", pars, sds, "pars", "sds")
  
  #returns myparams with appropriate NAs
  if(nofit){
    out = as.list(rep(NA_real_, length(myparams)))
    names(out) = myparams
    out[["success"]] = out[["cov"]] = NA_integer_
    out[["pars"]] = pars
    out[["sds"]] = sds
    return(out)
  }
  
  rmds <- tapply(resp, logc, median)
  if(!bidirectional) mmed = rmds[which.max(rmds)] else mmed = rmds[which.max(abs(rmds))] #shortened this code
  mmed_conc <- as.numeric(names(mmed)) #fixed this bug
  
  resp_max <- max(resp)
  resp_min <- min(resp)
  logc_min <- min(logc)
  logc_max <- max(logc)
  
  er_est <- if ((rmad <- mad(resp)) > 0) log(rmad) else log(1e-32)
  
  ###------------------------ Fit the Kd Model ------------------------###
  ## Starting parameters for the Kd Model
  g <- c(mmed, # top (tp)
         mmed_conc - 0.5, # logAC50 (ga)
         0.01, # bkg (b)
         er_est) # logSigma (er)
  if (g[1] == 0) g[1] <- 0.1
  
  ## Generate the bound matrices to constrain the model.
  #                tp   ga   b   er
  Ui <- matrix(c( 1,   0,   0,   0,
                  -1,   0,   0,   0,
                  0,   1,   0,   0,
                  0,  -1,   0,   0,
                  0,   0,   1,   0,
                  0,   0,  -1,   0),
               byrow = TRUE, nrow = 6, ncol = 4)
  if(!bidirectional) {
    hbnds <- c(0, -2*resp_max, # tp bounds
               logc_min - 1, -(logc_max + 2), # ga bounds
               0, -1) # b bounds
  } else {
    val <- 2*max(abs(resp_min),abs(resp_max))
    hbnds <- c(-val,-val, # tp bounds
               logc_min - 1, -(logc_max + 2), # ga bounds
               0, -1) # b bounds
  }
  Ci <- matrix(hbnds, nrow = 6, ncol = 1)
  
  ## Optimize the hill model
  fit <- try(constrOptim(g,
                         tcplObj,
                         ui = Ui,
                         ci = Ci,
                         mu = 1e-6,
                         method = "Nelder-Mead",
                         control = list(fnscale = -1,
                                        reltol = 1e-10,
                                        maxit = 6000),
                         conc = logc,
                         resp = resp,
                         fname = "logKd"),
             silent = !verbose)
  
  ## Generate some summary statistics
  if (!is(fit, "try-error")) { # Hill model fit the data
    if(verbose) cat("Kd >>>",fit$counts[1],fit$convergence,"\n")
    
    success <- 1L
    aic <- 2*length(fit$par) - 2*fit$value # 2*length(fit$par) - 2*fit$value
    mapply(assign,
           c(pars),
           fit$par,
           MoreArgs = list(envir = fenv))
    
    ## Calculate rmse for hill
    modl = logKd(fit$par, logc)
    rme <- sqrt(mean((modl - resp)^2, na.rm = TRUE))
    
    #Set ga to regular ac50 to match other fits
    # ga = 10^(ga)
    
    ## Calculate the sd for the hill parameters
    fit$cov <- try(solve(-hessian(tcplObj,
                                  fit$par,
                                  conc = logc,
                                  resp = resp,
                                  fname = "logKd")),
                   silent = !verbose)
    
    if (!is(fit$cov, "try-error")) { # Could invert hill Hessian
      
      cov <- 1L
      hdiag_sqrt <- suppressWarnings(sqrt(diag(fit$cov)))
      if (any(is.nan(hdiag_sqrt)) & !force_sd) {
        mapply(assign,
               sds,
               NaN,
               MoreArgs = list(envir = fenv))
      } else {
        mapply(assign,
               sds,
               hdiag_sqrt,
               MoreArgs = list(envir = fenv))
        #use taylor's theorem to approximate sd's in change of units
        #(only valid when sd's are much smaller than ln(10))
        # ga_sd = ga*log(10)*ga_sd
      }
      
    } else { # Could not invert hill Hessian
      
      cov <- 0L
      mapply(assign,
             c(sds),
             NA_real_,
             MoreArgs = list(envir = fenv))
      
    }
    
  } else { # Hill model did not fit the data
    success <- 0L
    aic <- NA_real_
    cov <- NA_integer_
    rme <- NA_real_
    modl = NA_real_
    
    mapply(assign,
           c(pars, sds),
           NA_real_,
           MoreArgs = list(envir = fenv))
    
  }
  
  return(mget(myparams))
  
}

fitpoly1 = function(conc, resp, bidirectional = FALSE, verbose = FALSE, nofit = FALSE){
  
  fenv <- environment()
  #initialize myparams
  pars <- paste0(c("a", "b", "er"))
  sds <- paste0(c("a", "b", "er"), "_sd")
  myparams = c("success", "aic", "cov", "rme", "modl", pars, sds, "pars", "sds")
  
  #returns myparams with appropriate NAs
  if(nofit){
    out = as.list(rep(NA_real_, length(myparams)))
    names(out) = myparams
    out[["success"]] = out[["cov"]] = NA_integer_
    out[["pars"]] = pars
    out[["sds"]] = sds
    return(out)
  }
  
  #median at each conc, for multi-valued responses
  rmds <- tapply(resp, conc, median)
  #get max response and corresponding conc
  if(!bidirectional) mmed = rmds[which.max(rmds)] else mmed = rmds[which.max(abs(rmds))] #shortened this code
  mmed_conc <- as.numeric(names(mmed)) #fixed this bug
  
  resp_max <- max(resp)
  resp_min <- min(resp)
  conc_min <- min(conc)
  conc_max <- max(conc)
  
  er_est <- if ((rmad <- mad(resp)) > 0) log(rmad) else log(1e-16)
  
  ###--------------------- Fit the Model ----------------------###
  ## Starting parameters for the Model
  a0 = mmed/conc_max #use largest response with desired directionality
  if(a0 == 0) a0 = .01  #if 0, use a smallish number
  g <- c(a0, # y scale (a); set to run through the max resp at the max conc
         0.01, # x scale (b); set to max conc
         er_est) # logSigma (er)
  
  ## Generate the bound matrices to constrain the model.
  #                a   b    er
  Ui <- matrix(c( 1,   0,   0,
                  -1,   0,   0,
                  0,   1,   0,
                  0,  -1,   0),
               byrow = TRUE, nrow = 4, ncol = 3)
  
  if(!bidirectional){
    bnds <- c(0, -1e8*abs(a0), # a bounds (always positive)
              0, -5) # b bounds (always increasing)
  } else {
    bnds <- c(-1e8*abs(a0), -1e8*abs(a0), # a bounds (positive or negative)
              0, -5) # b bounds (always increasing or always decreasing)
  }
  
  Ci <- matrix(bnds, nrow = 4, ncol = 1)
  
  ## Optimize the model
  fit <- try(constrOptim(g,
                         tcplObj,
                         ui = Ui,
                         ci = Ci,
                         mu = 1e-6,
                         method = "Nelder-Mead",
                         control = list(fnscale = -1,
                                        reltol = 1e-10,
                                        maxit = 6000),
                         conc = conc,
                         resp = resp,
                         fname = "poly1"),
             silent = !verbose)
  
  ## Generate some summary statistics
  if (!is(fit, "try-error")) { # The model fit the data
    if(verbose) cat("poly1 >>>",fit$counts[1],fit$convergence,"\n")
    
    success <- 1L
    aic <- 2*length(fit$par) - 2*fit$value # 2*length(fit$par) - 2*fit$value
    mapply(assign,
           c(pars),
           fit$par,
           MoreArgs = list(envir = fenv))
    
    ## Calculate rmse for gnls
    modl <- poly1(fit$par,conc)
    rme <- sqrt(mean((modl - resp)^2, na.rm = TRUE))
    
    ## Calculate the sd for the gnls parameters
    fit$cov <- try(solve(-hessian(tcplObj,
                                  fit$par,
                                  conc = conc,
                                  resp = resp,
                                  fname = "poly1")),
                   silent = !verbose)
    
    if (!is(fit$cov, "try-error")) { # Could invert gnls Hessian
      
      cov <- 1L
      diag_sqrt <- suppressWarnings(sqrt(diag(fit$cov)))
      if (any(is.nan(diag_sqrt))) {
        mapply(assign,
               sds,
               NaN,
               MoreArgs = list(envir = fenv))
      } else {
        mapply(assign,
               sds,
               diag_sqrt,
               MoreArgs = list(envir = fenv))
      }
      
    } else { # Could not invert gnls Hessian
      
      cov <- 0L
      mapply(assign,
             c(sds),
             NA_real_,
             MoreArgs = list(envir = fenv))
      
    }
    
  } else { # Curve did not fit the data
    
    success <- 0L
    aic <- NA_real_
    cov <- NA_integer_
    rme <- NA_real_
    modl = NA_real_
    
    mapply(assign,
           c(pars, sds),
           NA_real_,
           MoreArgs = list(envir = fenv))
    
  }
  
  return(mget(myparams))
  
}

fitconst = function(conc, resp, bidirectional = FALSE, verbose = FALSE, nofit = FALSE, force_sd = F){
  
  fenv <- environment()
  #initialize myparams
  pars <- paste0(c("b", "er"))
  sds <- paste0(c("b", "er"), "_sd")
  myparams = c("success", "aic", "cov", "rme", "modl", pars, sds, "pars", "sds")
  
  #returns myparams with appropriate NAs
  if(nofit){
    out = as.list(rep(NA_real_, length(myparams)))
    names(out) = myparams
    out[["success"]] = out[["cov"]] = NA_integer_
    out[["pars"]] = pars
    out[["sds"]] = sds
    return(out)
  }
  
  #median at each conc, for multi-valued responses
  rmds <- tapply(resp, conc, median)
  #get max response and corresponding conc
  if(!bidirectional) mmed = rmds[which.max(rmds)] else mmed = rmds[which.max(abs(rmds))] #shortened this code
  mmed_conc <- as.numeric(names(mmed)) #fixed this bug
  
  resp_max <- max(resp)
  resp_min <- min(resp)
  conc_min <- min(conc)
  conc_max <- max(conc)
  
  er_est <- if ((rmad <- mad(resp)) > 0) log(rmad) else log(1e-16)
  
  ###--------------------- Fit the Model ----------------------###
  ## Starting parameters for the Model
  a0 = .01  #if 0, use b smallish number
  g <- c(a0, # linear coeff (b); set to run through the max resp at the max conc
         er_est) # logSigma (er)
  
  ## Generate the bound matrices to constrain the model.
  #                b   er
  Ui <- matrix(c( 1,   0,
                  -1,   0),
               byrow = TRUE, nrow = 2, ncol = 2)
  
  if(!bidirectional){
    bnds <- c(0, -5) # b bounds (always positive)
    
  } else {
    bnds <- c(5, -5) # b bounds (positive or negative)
    
  }
  
  Ci <- matrix(bnds, nrow = 2, ncol = 1)
  
  ## Optimize the model
  fit <- try(constrOptim(g,
                         tcplObj,
                         ui = Ui,
                         ci = Ci,
                         mu = 1e-6,
                         method = "Nelder-Mead",
                         control = list(fnscale = -1,
                                        reltol = 1e-10,
                                        maxit = 6000),
                         conc = conc,
                         resp = resp,
                         fname = "cnst"),
             silent = !verbose)
  
  ## Generate some summary statistics
  if (!is(fit, "try-error")) { # The model fit the data
    if(verbose) cat("cnst >>>",fit$counts[1],fit$convergence,"\n")
    
    success <- 1L
    aic <- 2*length(fit$par) - 2*fit$value # 2*length(fit$par) - 2*fit$value
    mapply(assign,
           c(pars),
           fit$par,
           MoreArgs = list(envir = fenv))
    
    ## Calculate rmse for gnls
    modl <- cnst(fit$par,conc)
    rme <- sqrt(mean((modl - resp)^2, na.rm = TRUE))
    
    ## Calculate the sd for the gnls parameters
    fit$cov <- try(solve(-hessian(tcplObj,
                                  fit$par,
                                  conc = conc,
                                  resp = resp,
                                  fname = "cnst")),
                   silent = !verbose)
    
    if (!is(fit$cov, "try-error")) { # Could invert gnls Hessian
      
      cov <- 1L
      diag_sqrt <- suppressWarnings(sqrt(diag(fit$cov)))
      if (any(is.nan(diag_sqrt)) & !force_sd) {
        mapply(assign,
               sds,
               NaN,
               MoreArgs = list(envir = fenv))
      } else {
        mapply(assign,
               sds,
               diag_sqrt,
               MoreArgs = list(envir = fenv))
      }
      
    } else { # Could not invert gnls Hessian
      
      cov <- 0L
      mapply(assign,
             c(sds),
             NA_real_,
             MoreArgs = list(envir = fenv))
      
    }
    
  } else { # Curve did not fit the data
    
    success <- 0L
    aic <- NA_real_
    cov <- NA_integer_
    rme <- NA_real_
    modl = NA_real_
    
    mapply(assign,
           c(pars, sds),
           NA_real_,
           MoreArgs = list(envir = fenv))
    
  }
  
  return(mget(myparams))
  
}

tcplObj = function(p, conc, resp, fname, errfun = "dt4", err = NULL) {
  
  mu = do.call(fname, list(ps = p, x = conc)) #get model values for each conc
  n = length(p)
  if(is.null(err)) err = exp(p[n]) #set error term
  
  #objective function is sum of log-likelihood of response given the model at each concentration
  #scaled by variance (err)
  if(errfun == "dt4") return( sum( dt((resp - mu)/err, df = 4, log = TRUE) - log(err) ) )
  if(errfun == "dnorm") return( sum( dnorm((resp - mu)/err, log = TRUE) - log(err) ) )
  
}

logKd = function(ps,x){
  #hill function with log units: x = log10(conc) and ga = log10(ac50)
  #tp = ps[1], ga = ps[2], b = ps[3]
  return(ps[1]/(1 + 10^(ps[2]-x) ) + ps[3] )
}

logKd_nobkg = function(ps,x){
  #hill function with log units: x = log10(conc) and ga = log10(ac50)
  #tp = ps[1], ga = ps[2], b = ps[3]
  return(ps[1]/(1 + 10^(ps[2]-x) ) )
}

poly1 = function(ps,x){
  #a = ps[1], b = ps[2]
  return(ps[1]*x + ps[2])
}

cnst = function(ps,x){
  #constant background, b = ps[1]
  return(rep(ps[1],length(x)))
}

Kd = function(ps, x){
  #hill function with log units: x = log10(conc) and ga = log10(ac50)
  #tp = ps[1], ga = ps[2], b = ps[3]
  
  
  return(x*ps[1]/(x + ps[2] ) + ps[3] )
}

mypredict = function(fitout, x){
  
  if(fitout$type == "Kd"){
    return(logKd(c(fitout$tp, fitout$ga, fitout$b), log10(x)))
  } else if(fitout$type == "Linear"){
    return(poly1(c(fitout$a, fitout$b), x))
  } else if(fitout$type == "Constant"){
    return(cnst(c(fitout$b), x))
  }
  
}



