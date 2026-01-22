# library(liquidSVM)
library(xgboost)
library(keras)
library(ranger)

#if filename is not null, model will be saved using appropriate method and return NULL
modeler = function(xin = NULL, yin = NULL, xout = NULL, filename = NULL, predict.folder = NULL, ...){
  
  #unpack ellipses
  inlist = list(...)
  ml = inlist[["ml"]]
  fs = inlist[["fs"]]
  ml.params = inlist[["ml.params"]]
  to.scale = inlist[["to.scale"]]
  # site.imp.quant = inlist[["site.imp.quant"]]
  # if(is.null(site.imp.quant)) site.imp.quant = -1
  if(is.null(ml)) ml = "xgboost"
  
  #extract training weight
  if(any(colnames(xin) %in% "weight")){
    weight = xin[,"weight"]
    xin = xin[ , !colnames(xin) %in% "weight"]
    xout = xout[ , !colnames(xout) %in% "weight"]
  } else weight = NULL
  
  #feature selection
  if(!is.null(fs)){
    #correlation test
    if(fs == "cor"){
      pvals = apply(xin, 2, function(x, y){cor.test(x, y)$p.value}, y = yin)
      if(sum(pvals < .05) < 1) keep.feats = colnames(xin)[which.min(pvals)] else keep.feats = colnames(xin)[pvals < .05]
      # if(sum(pvals < .05) < 1) xin = xin[, which.min(pvals), drop = F] else xin = xin[, pvals < .05, drop = F]
      xin = xin[, keep.feats, drop = F]
      xout = xout[, keep.feats, drop = F]
    }
  }
  
  #scaling
  if(!is.null(to.scale) && length(to.scale) > 0){
    means = apply(xin[, to.scale], 2, mean)
    sds = apply(xin[, to.scale], 2, sd)
    
    #it's way faster to assign rows on sparse matrices vs. columns, so transpose, assign, transpose again
    scale.cols.in = Matrix(sapply(to.scale, function(col){ (xin[, col] - means[col])/sds[col] }), sparse = T)
    txin = t(xin)
    txin[to.scale,] = t(scale.cols.in)
    xin = t(txin)
    
    scale.cols.out = Matrix(sapply(to.scale, function(col){ (xout[, col] - means[col])/sds[col] }), sparse = T)
    txout = t(xout)
    txout[to.scale,] = t(scale.cols.out)
    xout = t(txout)
  }
  


  #xgboost model
  if(ml == "xgboost"){
    #set up ml.params and build model
    if(is.null(ml.params)) ml.params = list(nrounds = 1000)
    if(!is.null(weight)) ml.params = c(ml.params, list(weight = weight))
    mod = do.call("xgboost", c(ml.params, list(print_every_n = 1000, data = xin, label = yin)))
    #save model
    if(!is.null(filename)) {
      xgb.save(mod, paste0(filename, ".xgb"))
      return()
    }
    #predict
    preds = predict(mod, newdata = xout)
    # if(site.imp.quant >= 0) xgb_site_imp(xgb.importance(model = mod), site.imp.quant = site.imp.quant)
  }
  
  #keras neural net
  if(ml == "keras"){
    if(!is.null(predict.folder)) {
      dtype = gsub("_.*|output/", "", predict.folder)
      model = load_model_hdf5(paste0(predict.folder, "/combined_", dtype, ".h5"))
      cols.orig = scan(file = paste0(predict.folder, "/combined_", dtype, "_colnames.txt"), what = "character")
      miss.cols = cols.orig[!cols.orig %in% colnames(xout)]
      miss.mat = Matrix(data = 0, ncol = length(miss.cols), nrow = nrow(xout), sparse = T, 
                        dimnames = list(NULL, miss.cols))
      xout = cbind(xout, miss.mat)
      xout = xout[, cols.orig]
    } else {
      #ml.params default
      if(is.null(ml.params)) ml.params = list(layer.sizes = c(128, 32), relu.alphas = c(0, 0))
      #model is set up layer by layer, size and relu activation only
      model <- keras_model_sequential()
      n = length(ml.params$layer.sizes)
      model %>%
        layer_dense(units = ml.params$layer.sizes[1], input_shape = c(dim(xin)[2])) %>%
        layer_activation_relu(negative_slope = ml.params$relu.alphas[1])
      # model %>%
      #   layer_dense(units = ml.params$layer.sizes[1]) %>%
      #   layer_activation_relu(negative_slope = ml.params$relu.alphas[1])
      if(n > 1){
        for(i in 2:n){
          model %>% 
            layer_dense(units = ml.params$layer.sizes[i]) %>% 
            layer_activation_relu(negative_slope = ml.params$relu.alphas[i])
        } 
      }
      model %>% layer_dense(units = 1)
      model %>% compile(optimizer = optimizer_rmsprop(), loss = 'mse')
      # model %>% compile(optimizer = optimizer_sgd(), loss = 'mse')
      # print(summary(model))
      
      #actual model fitting, saving, prediction
      model %>% fit(xin, yin, epochs=10, batch_size=32, sample_weight = weight, verbose = 0)
      if(!is.null(filename)) {
        save_model_hdf5(model, paste0(filename, ".h5"))
        write(colnames(xin), file = paste0(filename, "_colnames.txt") )
        return()
      }
    }
    
    preds = as.vector(model %>% predict(xout))
  }
  
  if(ml == "epistasis"){
    require(reticulate)
    source_python("python/epistasis_modeler.py")
    df = cbind(xin, yin)
    df = df[is.finite(df[, "uncert"]) & is.finite(df[, "yin"]), ] #won't use points with missing uncertainty/scores
    colnames(df)[colnames(df) == "yin"] = "func_score"
    colnames(df)[colnames(df) == "uncert"] = "func_score_var"
    if(is.null(xout)) xout = xin[, "aa_substitutions", drop = F]
    out = epistasis_modeler(df, xout)
    dt = out[[1]]
    if(!is.null(filename)) {
      actual = yin
      write.csv(cbind(dt, actual), file = paste0(filename, ".csv"), row.names = F)
      lat = out[[2]]
      obs = out[[3]]
      lat$observed.effect = obs$effect[match(lat$mutation, obs$mutation)]
      colnames(lat)[colnames(lat) == "effect"] = "latent.effect"
      lat = lat[order(lat$site, lat$mutant),]
      write.csv(lat, file = paste0(filename, "_single_mut_effects.csv"), row.names = F)
      return()
    }
    
    preds = dt$preds
    preds[is.na(preds)] = median(preds[!is.na(preds)])
  }
  
  #our model is too big for this to work well
  if(ml == "svm"){
    # mod = lsSVM(x = as.matrix(xin) , y = yin, nthread = 0)
    mod = do.call("lsSVM", c(ml.params, list(x = xin, y = yin)))
    preds = predict(mod, xout)
  }
  
  #ranger random forest - deprecated
  if(ml == "rf"){
    if(is.null(ml.params)) ml.params = list(case.weights = weight) else ml.params = c(ml.params, list(case.weights = weight))
    mod = do.call("ranger", c(ml.params, list(x = xin, y = yin, num.threads = 7, importance = "impurity")))
    if(!is.null(filename)) {
      saveRDS(mod, file = paste0(filename, ".rds"))
      return()
    }
    # importance(mod)
    # mod = ranger(x = xin, y = yin, num.trees = 100, mtry = 20)
    preds = predict(mod, data = xout)$predictions
  }
  
  #Linear Regression
  if(ml == "lr"){
    if(is.null(ml.params)) ml.params = list(weights = weight) else ml.params = c(ml.params, list(weights = weight))
    newx = as.data.frame(cbind(xin, yin))
    mod = do.call("lm", c(ml.params, list(formula = yin ~ ., data = newx)))
    # mod = lm(yin ~ ., newx)
    preds = predict(mod, as.data.frame(xout))
  }

  
  return(preds)
  
}

#crossvalidates model over nfold folds
crossval = function(x, y, splitvec = NULL, ...){
  
  #unpack ellipses
  inlist = list(...)
  # w = inlist[["w"]]
  nfold = inlist[["nfold"]]
  if(is.null(nfold)) nfold = 5
  if(is.null(splitvec)) splitvec = inlist[["splitvec"]] #allows splitvec to be passed directly or through ...
  
  #get randomized splitting
  if(!is.null(splitvec)) nfold = max(splitvec) else splitvec = splitter(nrow(x), nfold)
  predlist = lapply(1:nfold, function(outfold){
    #apply splitting for given outfold as the test set
    xin = x[splitvec != outfold,]
    xout = x[splitvec == outfold, , drop = F]
    yin = y[splitvec != outfold]
    yout = y[splitvec == outfold]
    # if(!is.null(w)) win = w[splitvec != outfold] else win = NULL
    preds = modeler(xin, yin, xout, ...)

    return(data.frame(list(preds = preds, yout = yout), stringsAsFactors = F))
    
  })
  
  #format output and reorder to original
  finalout = as.data.frame(rbindlist(predlist))
  splitorder = unlist(lapply(1:nfold, function(x){which(splitvec == x)})) #tells how data was reordered
  finalout = finalout[match(1:nrow(finalout), splitorder),] #reorders back to the original
  return(finalout)
  
}

#randomly assigns a fold number to each row of data
splitter = function(rows, nfold) {
  folds = 1:rows%%nfold + 1
  folds = sample(folds)
  return(folds)
}

#randomly assigns a fold number to each category in cat_col separately, so each is evenly spread among folds
even_splitter = function(cat_col, nfold){
  out = rep(NA_integer_, length(cat_col))
  for(cat in unique(cat_col)){
    out[cat_col == cat] = splitter(sum(cat_col == cat), nfold)
  }
  if(sum(is.na(out)) != 0) stop("Error in even_splitter")
  return(out)
}