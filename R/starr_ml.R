library(Matrix)
library(data.table)
library(Rfast)
library(keras)
source("R/ml_features.R")
source("R/ml_tuner.R")
source("R/ml_preprocess.R")
source("R/ml_plot.R")
source("R/ml_modeler.R")
source("R/score_funs.R")

#version 2
#datatype = "bind", "escape", "log_escape"
#wts = "none", "uncert", "fifths", "mahala"
#to use protr feats ex: feats = MoreauBroto, protr.params = list(nlag = 5)
# ex: feats = ProtFP, protr.params = list(index = 1:544, lag = 8, pc = 8)
starr_ml = function(datatype = "bind", feats = "lt", wts = "none", test = F, comb_antib = F,
                    tunefold = 1, do.plot = T, seed = 1234, ...){

  #beginning stuff
  start_time = Sys.time()
  tensorflow::set_random_seed(seed)
  
  #unpack ellipses
  inlist = list(...)
  ml = inlist[["ml"]]
  ml.params = inlist[["ml.params"]]
  nfold = inlist[["nfold"]]
  ntune = inlist[["ntune"]]
  whichtune = inlist[["whichtune"]]
  if(is.null(ntune)) ntune = 0
  if(is.null(nfold)) nfold = 5
  if(is.null(ml)) ml = "xgboost"
  
  #create output dirname based on options
  if(wts != "none") addname = paste0("_", wts) else addname = NULL
  if(!is.null(ml.params)) addname = paste0(addname, "_", 
                                           paste(names(ml.params), ml.params, sep = "_", collapse = "_"))
  if(test) addname = paste0(addname, "_test")
  if(comb_antib) addname = paste0(addname, "_comb")
  if(ntune > 0 && tunefold <= 1 && is.null(whichtune)) addname = paste0(addname, "_tuned", ntune)
  if(!is.null(whichtune)) addname = paste0(addname, "_tuned", whichtune[1], "_", whichtune[length(whichtune)])
  if(ntune > 0 && tunefold > 1) addname = paste0(addname, "_fulltuned", ntune)
  if(nfold <= 1) addname = paste0(addname, "_onefold")
  dirname = paste0("output/ml_", datatype, "_", ml, "_", paste0(feats, collapse = "_"), addname)
  if(!dir.exists(dirname)) dir.create(dirname)
  
  #do pre-processing, if necessary
  inname = paste0("input/", datatype, ".rds")
  if(!file.exists(inname)) ml_preprocess(datatype)
  indf = readRDS(inname)
  #this test gives 1000 entries per name
  if(test) indf = as.data.frame(rbindlist(lapply(unique(indf$name), function(name){indf[indf$name %in% name,][1:1000,]})))

  #maybe do data plot
  if(do.plot) data_plot(datatype)

  #get feature set
  if(feats[1] == "full") xylist = fullonehot(indf) 
  else if(feats[1] == "lt") xylist = locandtrans(indf, count = F, triads = F)
  else if(feats[1] == "aggonly") xylist = fullonehot(indf, aggonly = T) 
  else xylist = protr_feats(indf, feats = feats, ...)
  
  if(datatype %in% c("bindv2", "bindv3")) {
    source = as.numeric(as.factor(xylist$source))-1
    test = cbind(xylist$x, source)
  }
  
  #a cludge to combine antibody models into one
  if(comb_antib){
    name.df = list( Atype_ = xylist$name)
    anames = sparse.model.matrix( ~ Atype_, name.df)[,-1]
    xylist$name = rep("Antibody", length(xylist$name))
    xylist$x = cbind(xylist$x, anames)
  }
  
  #this is a cludge to be able to handle binding/expression(one model/name) or antib (ten models/names)
  allnames = unique(xylist$name)
  # if(test) allnames = allnames[1]
  
  #loop over all models to build
  no.out = lapply(allnames, function(name){
    
    usefile = paste0(dirname, "/", name) 
    #get specific x/y/uncert to use for this model
    x = xylist$x[xylist$name == name,]
    y = xylist$y[xylist$name == name]
    uncert = xylist$uncert[xylist$name == name]

    #get weights, will now be carried as a column of x
    weight = get_w(wt.type = wts, x = x, y = y, uncert = uncert)
    if(!is.null(weight)) x = cbind(x, weight)

    #actually build the model
    ################################new paradigm: 
    #tuning will save the plotting data for the best run, and save tuning stats (inside ml_tuner)
    #plots will be moved to separate function/file that reads directly from file
    #modeler/crossval moved to separate file
    #helper functions will be moved appropriately
    #compare_plot will be deprecated
    
    
    if(ntune <= 0){
      if(nfold <= 1){
        #single fold modeler - saves model internally
        no.out = modeler(xin = x, yin = y, filename = usefile, ...)
        return()
      } else {
        #cross-validated, untuned
        pred.df = crossval(x = x, y = y, ...) 
      }

    } else {
      if(tunefold <= 1){
        #tuned model - saves tuning stats internally, outputs pred.df of best one
        pred.df = ml_tuner(x = x, y = y, filename = usefile, ...)
      } else {
        #cross-validates with tuning inside each fold, for accurate estimation of final model acc
        pred.df = full_tuner(x = x, y = y, tunefold = tunefold, filename = usefile, ...)
      }

    }
    
    #add other useful data to pred.df and save
    pred.df$uncert = uncert
    pred.df$weight = weight
    pred.df$nsubs = xylist$nsubs[xylist$name == name]
    saveRDS(pred.df, file = paste0(usefile, ".rds"))
    return()
    
  })

  #maybe do plotting
  if(nfold == 1) do.plot = F #to avoid overwriting other results
  if(do.plot) ml_plot(gsub("output/", "", dirname), datatype)
  
  return(Sys.time() - start_time)
}

#############################################Weights#################################################
#Calculates weights of various kinds: inverse uncertainty, equal fifths, mahalanobis distance
get_w = function(wt.type, x = NULL, y = NULL, uncert = NULL){
  if(wt.type %in% c("uncert", "sqrt_uncert")) {
    if(wt.type == "uncert") w = 1/uncert else w = 1/sqrt(uncert)
    q95 = quantile(w, .95, na.rm = T)
    w[w > q95] = q95 #weights max out at 95% quantile to avoid dividing by very low errors
    w = w/sum(w, na.rm = T)*length(w[!is.na(w)]) # normalize weights to 1 on average
    w[is.na(w)] = 1 # missing uncertainty get a weight of one
  } else if(wt.type == "fifths") {
    #cuts endpoint in equal fifths and upweights underrepresented fifths
    breaks = seq(min(y), max(y), length.out = 6)
    cutfac = cut(y, breaks, include.lowest = T)
    wtab = length(y)/table(cutfac)/5
    w = as.numeric(wtab[cutfac])
  } else if(wt.type == "mahala") {
    #mahalnobis is slow and relies on Rfast library to speed up
    print("Mahalanobis weighting beginning.")
    covx = cova(as.matrix(x), large = T)
    choldecomp = cholesky(covx)
    cm = colMeans(x)
    d = mahala(x = as.matrix(x), mu = cm, sigma = choldecomp, ischol = T)
    w = sqrt(d)
    w = w/sum(w)*length(w)
    print("Mahalanobis weighting ended.")
  } else { 
    w = NULL
  }
  
  return(w)
}


#########################################Deprecated#################################################
# #plots two model errors on same plot
# compare_plot = function(names = c("xgboost_lt_1000", "keras_full")){
#   
#   infiles = paste0("output/ml_", names, ".rds")
#   dfs = lapply(infiles, readRDS)
#   
#   filename = paste0("output/compare_",  paste0(names, collapse = "__"),".pdf")
#   pdf(filename, width = 8.5, height = 11)
#   
#   par(mfrow = c(2,1))
#   
#   errs = sapply(dfs, function(x){x$preds - x$yout})
#   
#   #Error Dist
#   use_wmses = sapply(dfs, function(x) signif(wmse(x$preds - x$yout, 1/x$uncert), 2))
#   ds = lapply(errs, density)
#   plot(ds[[1]], main = paste0("Prediction Error Distribution"), xlab = "Error (log10(M))", 
#        ylim = c(0, max(c(ds[[1]]$y, ds[[2]]$y))) )
#   points(ds[[2]], col = "red", type = "l")
#   legend("topleft", legend = paste0(names, "; wt_rmse = ", use_wmses), col = c("black", "red"), lwd = c(1,1))
#   
#   #Error ratio dist
#   to.plots =  sapply(dfs, function(x){(x$preds - x$yout)/x$uncert})
#   ds = lapply(to.plots, density)
#   plot(ds[[1]], main = paste0("Prediction Error/Observation Error Distribution"), xlab = "Standard Errors",
#        ylim = c(0, max(c(ds[[1]]$y, ds[[2]]$y))))
#   points(ds[[2]], col = "red", type = "l")
#   legend("topleft", legend = names, col = c("black", "red"), lwd = c(1,1))
#   
#   #Error ratio zoom-in
#   levels = c(1,2,3)
#   underses = sapply(to.plots, function(to.plot){
#     paste0(sapply(levels, function(level) round(sum(abs(to.plot) < level)/length(to.plot)*100, 1)), collapse = ", ")
#   })
#   ds = lapply(to.plots, function(to.plot){density(to.plot[abs(to.plot) < 5])})
#   plot(ds[[1]], main = paste0("Prediction Error/Observation Error Distribution"), xlab =  "Standard Errors",
#        ylim = c(0, max(c(ds[[1]]$y, ds[[2]]$y))))
#   points(ds[[2]], col = "red", type = "l")
#   abline(v = c(-1, -2, -3, 1,2,3))
#   legend("topleft", legend = paste0(names, "\n% w/i 1,2,3 SEs:  ", underses), 
#          col = c("black", "red"), lwd = c(1,1))
#   
#   #Regular Error Plot
#   rmses = sapply(dfs, function(x){round(rmse(x$preds - x$yout), 2)})
#   cors = sapply(dfs, function(x){round(cor(x$yout, x$preds), 2)})
#   cols = c(rgb(0,0,0,.1), rgb(1,0,0,.1))
#   
#   maxrange = range(c(dfs[[1]]$preds, dfs[[2]]$preds))
#   plot(dfs[[1]]$yout, dfs[[1]]$preds, xlab = "Observed Delta(log10(Ka)) Binding", 
#        ylab = "Predicted Delta(log10(Ka)) Binding", pch = 16, col = rgb(0,0,0,.1), ylim = maxrange)
#   points(dfs[[2]]$yout, dfs[[2]]$preds, col = rgb(1,0,0,.1), pch = 16)
#   abline(0, 1, lwd = 3)
#   legend("bottomright", legend = paste0(names, ", RMSE: ", rmses, ", cor: ", cors), 
#          col = cols, pch = c(16,16))
#   
#   plot(dfs[[2]]$yout, dfs[[2]]$preds, xlab = "Observed Delta(log10(Ka)) Binding", 
#        ylab = "Predicted Delta(log10(Ka)) Binding", pch = 16, col = rgb(1,0,0,.1), 
#        ylim = maxrange)
#   points(dfs[[1]]$yout, dfs[[1]]$preds, col = rgb(0,0,0,.1), pch = 16)
#   abline(0, 1, lwd = 3)
#   legend("bottomright", legend = paste0(names, ", RMSE: ", rmses, ", cor: ", cors), 
#          col = cols, pch = c(16,16))
#   
#   dev.off()
#     
# }


# score_var = function(ni, nf, F1, fl = 0, ce = 1){
#   ni = ni + .5
#   nf = nf + .5
#   Ni = sum(ni)
#   Nf = sum(nf)
#   
#   # F1 = .4042*Nf/Ni
#   # print(F1)
#   E = F1*nf*Ni/ni/Nf
#   Edpre = F1*nf*(Ni + 1)/(ni + 1)/Nf
#   Edpost = F1*(nf + 1)*Ni/ni/(Nf + 1)
#   
#   E = floor_ceiling(E, fl, ce)
#   Edpre = floor_ceiling(Edpre, fl, ce)
#   Edpost = floor_ceiling(Edpost, fl, ce)
#   
#   svar = (Edpre - E)^2*ni + (Edpost - E)^2*nf
#   return(list(E = E, svar = svar))
# }
# 
# floor_ceiling = function(x, fl, ce){
#   x[x < fl] = fl
#   x[x > ce] = ce
#   return(x)
# }



