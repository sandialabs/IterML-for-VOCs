library(parallel)
# library(pbapply)

#version 2
#optimizes meta-parameters using random search for given ml model according to best RMSE
#writes tried parameter sets and error stats to file and returns best pred.df from crossval
ml_tuner = function(x, y, filename, which.fold = 0, ...){
  
  #unpack ellipses
  inlist = list(...)
  ml = inlist[["ml"]]
  ncores = inlist[["ncores"]]
  kcores = inlist[["kcores"]]
  ntune = inlist[["ntune"]]
  whichtune = inlist[["whichtune"]] #can can part of tuning list for long runs
  to.source = inlist[["to.source"]]
  if(is.null(to.source)) to.source = "R/starr_ml.R"
  if(is.null(ntune)) ntune = 1
  if(is.null(ml)) ml = "xgboost"
  if(is.null(ncores)) ncores = 6 
  if(is.null(kcores)) kcores = 1L
  
  
  #get tuning param list
  ml.param.list = get_tune_params(ntune, ml)
  if(!is.null(whichtune)) ml.param.list = ml.param.list[whichtune]

  start.time = Sys.time()
  if(ncores > 1){
    tensorflow::tf$config$threading$set_inter_op_parallelism_threads(kcores)
    tensorflow::tf$config$threading$set_intra_op_parallelism_threads(kcores)
    
    #set up parallel cluster, reseed
    cl = makeCluster(ncores, outfile = paste0(filename, "_tuneprog_ncore_", ncores, "_fold_", which.fold,".txt"))
    newseed = sample(1:.Machine$integer.max, 1)
    no.out = clusterCall(cl, function(newseed){
      source(to.source)
      tensorflow::set_random_seed(newseed)
      tensorflow::tf$config$threading$set_inter_op_parallelism_threads(kcores)
      tensorflow::tf$config$threading$set_intra_op_parallelism_threads(kcores)
    }, newseed = newseed)
    
    #parallel runs on crossval with timing output
    predlist = parLapply(cl = cl, X = 1:length(ml.param.list), 
                        fun = function(i, ml.param.list, myx, y, start.time, ...){
      new.start = Sys.time()
      cat("Begin tuning round", i, "with layer sizes", ml.param.list[[i]]$layer.sizes, "and relu alphas", 
          ml.param.list[[i]]$relu.alphas, "at:", capture.output(new.start - start.time),"\n")
      pred.df = crossval(x = myx, y = y, ml.params = ml.param.list[[i]], ...)
      cat("Finish tuning round", i, "with layer sizes", ml.param.list[[i]]$layer.sizes, "and relu alphas", 
          ml.param.list[[i]]$relu.alphas, "after:", capture.output(Sys.time() - new.start),"\n")
      return(pred.df)
    }, myx = x, y = y, ml.param.list = ml.param.list, start.time = start.time, ...)
    
    stopCluster(cl)
  } else {
    #non-parallel runs on crossval with timing output
    predlist = lapply(1:length(ml.param.list), function(i, ml.param.list, myx, y, start.time, ...){
      cat("Begin tuning round", i, "with layer sizes", ml.param.list[[i]]$layer.sizes, "and relu alphas",
          ml.param.list[[i]]$relu.alphas, "at:", capture.output(Sys.time() - start.time),"\n")
      pred.df = crossval(x = myx, y = y, ml.params = ml.param.list[[i]], ...)
      return(pred.df)
    }, myx = x, y = y, ml.param.list = ml.param.list,  start.time = start.time, ...)
    
    # predlist = vector(mode = "list", length = length(ml.param.list))
    # for(i in 1:length(ml.param.list)){
    #   cat("Begin tuning round", i, "with layer sizes", ml.param.list[[i]]$layer.sizes, "and relu alphas",
    #       ml.param.list[[i]]$relu.alphas, "at:", capture.output(Sys.time() - start.time),"\n")
    #   predlist[[i]] = crossval(x = x, y = y, ml.params = ml.param.list[[i]], ...)
    # }
  }

  # return(list(outlist = predlist, ml.param.list = ml.param.list))

  #get tune stats
  tune.rmses = sapply(predlist, function(x){rmse(x$yout - x$preds)})
  tune.cors = sapply(predlist, function(x){cor(x$yout, x$preds)})
  tune.q2s = sapply(predlist, function(x){q2(x$yout, x$preds)})
  
  #make tunedf and output, return best pred.df
  tunedf = as.data.frame(list(RMSE = tune.rmses, Q2 = tune.q2s, Correlation = tune.cors), stringsAsFactors = F)
  for(i in 1:4){ tunedf[, paste0("Layer.", i)] = sapply(ml.param.list, function(x){x$layer.sizes[i]}) }
  for(i in 1:4){ tunedf[, paste0("Relu.alpha.", i)] = sapply(ml.param.list, function(x){x$relu.alphas[i]}) }
  tunedf = tunedf[order(tunedf$RMSE),]
  
  write.csv(tunedf, file = paste0(filename, "_tunedf", "_fold_", which.fold, ".csv"))
  
  #return best params for full_tuner, otherwise return best predictions
  if(which.fold > 0) return(ml.param.list[[which.min(tune.rmses)]])
  return(predlist[[which.min(tune.rmses)]])
  
}

#cross-validates with tuning inside each fold, for accurate estimation of final model acc
#very similar to crossval
full_tuner = function(x, y, tunefold, filename, ...){
  
  #get randomized splitting
  splitvec = splitter(nrow(x), tunefold)
  predlist = lapply(1:tunefold, function(outfold){
    #apply splitting for given outfold as the test set
    xin = x[splitvec != outfold,]
    xout = x[splitvec == outfold,]
    yin = y[splitvec != outfold]
    yout = y[splitvec == outfold]
    # if(!is.null(w)) win = w[splitvec != outfold] else win = NULL
    
    ml.params = ml_tuner(x = xin, y = yin, filename = filename, which.fold = outfold, ...)
    preds = modeler(xin, yin, xout, ml.params = ml.params, ...)
    
    return(data.frame(list(preds = preds, yout = yout), stringsAsFactors = F))
    
  })
  
  #format output and reorder to original
  finalout = as.data.frame(rbindlist(predlist))
  splitorder = unlist(lapply(1:tunefold, function(x){which(splitvec == x)})) #tells how data was reordered
  finalout = finalout[match(1:nrow(finalout), splitorder),] #reorders back to the original
  return(finalout)
  
}

#randomly selects ntune sets of meta-parameters from preset choices
#only keras is currently supported
get_tune_params = function(ntune, ml){
  if(ml == "keras"){
    ml.param.list = unique(lapply(1:ntune, keras_tune_gen))
    #while loop could take a long time or not terminate if number of possible combinations is similar to ntune
    while(length(ml.param.list) < ntune){
      ml.param.list = unique(c(ml.param.list, list(keras_tune_gen())))
    }
    return(ml.param.list)
  }
  
  
}

#generates a single set of random keras parameters within presets
keras_tune_gen = function(x){
  nlayers = sample(2:3, 1)
  layer.sizes = 2^(sample(2:8, nlayers, replace = T)) #8 to 512
  # relu.alphas = sample(c(0,0,.01,.01), nlayers, replace = T) #0 or .01 (leaky)
  relu.alphas = rep(sample(c(0,0,.01), 1, replace = T), nlayers) #all 0 or all .01 (leaky)
  return(list(layer.sizes = layer.sizes, relu.alphas = relu.alphas))
}