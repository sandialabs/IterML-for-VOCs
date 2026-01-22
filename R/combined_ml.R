library(Matrix)
library(data.table)
library(Rfast)
# pyconf = reticulate::py_discover_config()
# if(pyconf$python != "C:/Users/tysheff/.conda/envs/R_keras/python.exe") reticulate::use_condaenv("R_keras")
library(keras)
source("R/combined_features.R")
source("R/ml_tuner.R")
source("R/combined_preprocess.R")
source("R/combined_plot.R")
source("R/ml_modeler.R")
source("R/score_funs.R")

#version 3 of starr_ml, focused on combining with our data
#datatype = "binding", "antibody", "antibody_X" (where X is a single antibody name)
#wts = "none", "uncert", "fifths", "mahala", "rebal50" 
#(or "rebalx" where 0 <= x <= 100, x is percentage of weight given to mai source)
#to use protr feats ex: var.feats = "MoreauBroto", max.var.feats =...
# var.locs = c("conform_color", "conform_sep", "conform_paste", "conform_E", 
# "interact_color", "interact_sep", "interact_paste", "interact_E"); paste/sep/color only affect protr features
#var.locs uses only regions specified by mai, according to the setting
combined_ml = function(datatype = "antibody", var.feats = "full", antib.feats = "mixed", wts = "none", test = F,
                    tunefold = 1, seed = 1234, max.antib.feats = 16, max.var.feats = 16, do.plot = T, delta.kd = "wuhan", 
                    source = "all", var.locs = "all", get.dir.name = F, predict.folder = NULL, fix.scale = T, ...){

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
  protr.params = inlist[["protr.params"]]
  if(is.null(ntune)) ntune = 0
  if(is.null(nfold)) nfold = 5
  if(is.null(ml)) ml = "xgboost"
  if(datatype == "binding" || ml == "epistasis") antib.feats = "none"
  if(ml == "epistasis") var.feats = "epistasis"
  
  #create output dirname based on options
  if(wts != "none") addname = paste0("_", wts) else addname = ""
  if(!is.null(ml.params)) addname = paste0(addname, "_", 
                                           paste(names(ml.params), ml.params, sep = "_", collapse = "_"))
  if(source != "all") addname = paste0(addname, "_", source)
  if(test) addname = paste0(addname, "_test")
  if(max.antib.feats != 16) addname = paste0(addname, "_max_afeats_", max.antib.feats)
  if(max.var.feats != 16) addname = paste0(addname, "_max_vfeats_", max.var.feats)
  if(var.locs != "all") addname = paste0(addname, "_", var.locs)
  if(ntune > 0 && tunefold <= 1 && is.null(whichtune)) addname = paste0(addname, "_tuned", ntune)
  if(!is.null(whichtune)) addname = paste0(addname, "_tuned", whichtune[1], "_", whichtune[length(whichtune)])
  if(ntune > 0 && tunefold > 1) addname = paste0(addname, "_fulltuned", ntune)
  if(nfold <= 1) addname = paste0(addname, "_onefold")
  if(delta.kd != "none") addname = paste0(addname, "_delta", delta.kd)
  if(!is.null(predict.folder)) addname = paste0(addname, "_predicted")
  if(fix.scale) addname = paste0(addname, "_fixscale")
  dirname = paste0("output/", datatype, "_", ml, "_", paste0(c(var.feats, antib.feats), collapse = "_"), addname)
  if(get.dir.name) return(dirname)
  if(!dir.exists(dirname)) dir.create(dirname)
  
  #do pre-processing, if necessary
  #pull out single antibody if specified in datatype, then change datatype to antibody
  if(grepl("antibody_", datatype)){
    use.antib = gsub("antibody_", "", datatype)
    datatype = "antibody"
  } else use.antib = NULL
  inname = paste0("input/", datatype, "_preprocessed.rds")
  if(!file.exists(inname)) do.call(paste0(datatype, "_preprocess"))
  indf = readRDS(inname)
  if(!is.null(use.antib)) indf = indf[indf$Antibody %in% use.antib, ]
  
  #this test gives maximum 1000 entries per antibody
  if(test) {
    if(datatype == "antibody"){
      indf = as.data.frame(rbindlist(lapply(unique(indf$Antibody), function(antib){indf[indf$Antibody %in% antib,][1:1000,]})))
      indf = indf[!is.na(indf$Antibody),]
    } else indf = rbind(indf[indf$source == "mai",], indf[indf$source == "starr", ][1:10000, ])
  }
  if(source != "all") indf = indf[indf$source == source, ]
  
  #delta kd relative to wuhan
  if(delta.kd == "wuhan"){
    to.sub = indf[indf$aa_substitutions == "",]
    indf$wuhan = to.sub$endpoint[match(paste0(indf$source, "_", indf$Antibody), 
                                       paste0(to.sub$source, "_", to.sub$Antibody))]
    indf = indf[!is.na(indf$wuhan) & indf$aa_substitutions != "",] #missing wuhans will be removed
    indf$endpoint = indf$endpoint - indf$wuhan
  } #else if(delta.kd == "ace2"){
  #   #delta kd relative to ACE2
  #   to.sub = indf[indf$Antibody == "huACE2-Fc",]
  #   indf$ace2 = to.sub$endpoint[match(indf$Variant, to.sub$Variant)]
  #   indf = indf[!is.na(indf$ace2) & indf$Antibody != "huACE2-Fc",]
  #   indf$endpoint = indf$endpoint - indf$ace2
  # }
  indf = indf[!indf$Antibody %in% c("huACE2-Fc", "IgG"),] #can't use ACE2 under any circumstance

  #maybe do data plot
  # if(do.plot) data_plot(datatype)

  #get feature set
  xylist = combined_features(indf, var.feats, antib.feats, max.antib.feats = max.antib.feats, max.var.feats = max.var.feats,
                             var.locs = var.locs, datatype = datatype, use.antib = use.antib, fix.scale = fix.scale, ...)

  usefile = paste0(dirname, "/combined_", datatype)
  x = xylist$x
  y = xylist$y
  uncert = xylist$uncert
  if(ml == "epistasis") x = cbind(x, uncert)

  #get weights, will now be carried as a column of x
  weight = get_w(wt.type = wts, x = x, y = y, rebalrows = (xylist$source == "mai"))
  if(!is.null(weight)) x = cbind(x, weight)

  #actually build the model
  ################################new paradigm: 
  #tuning will save the plotting data for the best run, and save tuning stats (inside ml_tuner)
  #plots will be moved to separate function/file that reads directly from file
  #modeler/crossval moved to separate file
  
  #want to evenly distribute data from both sources among folds
  splitvec = even_splitter(xylist$source, nfold)
  
  if(!is.null(predict.folder)){
    preds = modeler(xout = x, predict.folder = predict.folder, to.scale = xylist$to.scale, ...)
    pred.df = data.frame(list(preds = preds, yout = y), stringsAsFactors = F)
    write(predict.folder, paste0(dirname, "/predict_folder_used.txt"))
  } else if(ntune <= 0){
    if(nfold <= 1){
      #single fold modeler - saves model internally
      no.out = modeler(xin = x, yin = y, filename = usefile, to.scale = xylist$to.scale, ...)
      return()
    } else {
      #cross-validated, untuned
      pred.df = crossval(x = x, y = y, splitvec = splitvec, to.scale = xylist$to.scale, ...) 
    }

  } else {
    if(tunefold <= 1){
      #tuned model - saves tuning stats internally, outputs pred.df of best one
      pred.df = ml_tuner(x = x, y = y, filename = usefile, splitvec = splitvec, to.source = "R/combined_ml.R", 
                         to.scale = xylist$to.scale, ...)
    } else {
      #cross-validates with tuning inside each fold, for accurate estimation of final model acc
      pred.df = full_tuner(x = x, y = y, tunefold = tunefold, filename = usefile, splitvec = splitvec, 
                           to.source = "R/combined_ml.R", to.scale = xylist$to.scale, ...)
    }

  }
  
  #add other useful data to pred.df and save
  pred.df$weight = weight
  pred.df$antibody = xylist$antibody
  pred.df$nsubs = xylist$nsubs
  pred.df$source = xylist$source
  saveRDS(pred.df, file = paste0(usefile, ".rds"))


  #maybe do plotting
  if(nfold == 1) do.plot = F #to avoid overwriting other results
  if(do.plot) combined_plot(gsub("output/", "", dirname), datatype)
  
  return(Sys.time() - start_time)
}

#predicts on already generated model located at predict folder
#uses new specifications to choose data, but infers ml, features, weights, delta.kd from the folder name
combined_predict = function(predict.folder, datatype = "antibody", delta.kd = "none", test = F, source = "all", 
                            var.locs = "all", get.dir.name = F){
  
  #infer model parameters
  #ml
  mls = c("rf", "keras", "xgboost")
  ml = mls[sapply(mls, grepl, x = predict.folder)]
  
  #feats
  var.feats = gsub(paste0(".*", ml, "_"), "", predict.folder)
  var.feats = gsub("_.*", "", var.feats)
  
  antib.feats = gsub(paste0(".*", var.feats, "_"), "", predict.folder)
  antib.feats = gsub("_.*", "", antib.feats)
  
  #weights
  all.wts = c("uncert", "sqrt_uncert", "fifths", "mahala", "rebal")
  wts = mls[sapply(all.wts, grepl, x = predict.folder)]
  if(length(wts) == 0) wts = "none"
  if(wts == "rebal") stop("Haven't written predict code for rebal wts")
  
  #delta.kd
  if(grepl("delta", predict.folder)) {
    delta.kd = gsub(".*_delta", "", predict.folder)
    delta.kd = gsub("_.*", "", delta.kd)
  } else delta.kd = "none"
  
  #max.feats
  if(grepl("max_vfeats", predict.folder)){
    max.var.feats = gsub(".*max_vfeats_", "", predict.folder)
    max.var.feats = as.numeric(gsub("_.*", "", max.var.feats))
  } else max.var.feats = 16
  
  if(grepl("max_afeats", predict.folder)){
    max.antib.feats = gsub(".*max_afeats_", "", predict.folder)
    max.antib.feats = as.numeric(gsub("_.*", "", max.antib.feats))
  } else max.antib.feats = 16
  
  combined_ml(datatype = datatype, ml = ml, var.feats = var.feats, antib.feats = antib.feats, wts = wts, 
              test = test, max.antib.feats = max.antib.feats, max.var.feats = max.var.feats,
              delta.kd = delta.kd, source = source, var.locs = var.locs, predict.folder = predict.folder,
              get.dir.name = get.dir.name)
  
}

#############################################Weights#################################################
#Calculates weights of various kinds: inverse uncertainty, equal fifths, mahalanobis distance
#rebalrows is a logical vector, same length as y, indicating which 
get_w = function(wt.type, x = NULL, y = NULL, uncert = NULL, rebalrows = NULL){
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
  } else if(grepl("rebal", wt.type)){
    wt.fac = as.numeric(gsub("rebal", "", wt.type))/100
    N = length(rebalrows)
    w = rep(1, N)
    w[rebalrows] = wt.fac*N/sum(rebalrows)
    w[!rebalrows] = (1 - wt.fac)*N/sum(!rebalrows)
    
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



