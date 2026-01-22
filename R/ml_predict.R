source("R/starr_ml.R")
library(seqinr)
library(data.table)
library(gplots)
library(RColorBrewer)
library(openxlsx)

#This assumes you've already made the model with starr_ml
#aa_subs is a vector of variants, each entry is a string with mutations separated by a space
#binding: inpath = "ml_bind_keras_full_layer.sizes_c(64, 4, 16)_relu.alphas_c(0.01, 0.01, 0.01)_onefold0/Binding.h5"
#expression: inpath = "ml_express_keras_full_layer.sizes_c(16, 128)_relu.alphas_c(0.01, 0.01)_onefold/Expression.h5"
#antibody: inpath = "ml_log_escape_keras_full_layer.sizes_c(256, 32, 16)_relu.alphas_c(0, 0, 0)_comb_onefold/Antibody.h5"
ml_predict = function(aa_subs, inpath){
  
  #get the premade model first (so it fails fast if it wasn't made or was misspelled)
  #this new version is a workaround for some janky keras version checking that doesn't work
  # filepath = paste0("output/", inpath)
  # custom_objects <- keras:::objects_with_py_function_names(NULL)
  # args <- list(filepath = keras:::normalize_path(filepath), custom_objects = custom_objects)
  # args$compile <- TRUE
  # mymodel = do.call(keras$models$load_model, args)
  mymodel = load_model_hdf5(paste0("output/", inpath))
  
  #extract parameters from the path name
  datatype = gsub("ml_", "", gsub("_keras.*", "", inpath))
  feats = gsub(".*keras_", "", gsub("_layer.*", "", inpath))
  
  #get the values for experimental/starr along with the variant mutations
  preddummy = as.data.frame(list(aa_substitutions = aa_subs, endpoint = rep(0, length(aa_subs)), 
                                 uncert = rep(0, length(aa_subs))))

  #still need to match the original feature encoding
  inname = paste0("input/", datatype, ".rds")
  indf = readRDS(inname)
  if(feats[1] == "full") {
    predlist = fullonehot(preddummy, noagg = T)
    xylist = fullonehot(indf)
  } else if(feats[1] == "lt"){
    predlist = locandtrans(preddummy, count = F, triads = F, noagg = T)
    xylist = locandtrans(indf, count = F, triads = F)
  } 
  
  # else {
  #   xylist = protr_feats(indf, feats = feats, ...)
  #   predlist = protr_feats(preddummy, feats = feats, ...)
  # }
  
  #fix predlist x to have same features as original
  misscols = colnames(xylist$x)[!colnames(xylist$x) %in% colnames(predlist$x)]
  newx = cbind(predlist$x, Matrix(data = 0, nrow = nrow(predlist$x), 
                                  ncol = length(misscols), dimnames = list(NULL, misscols)))
  finalx = newx[, match(colnames(xylist$x), colnames(newx)), drop = F]
  
  #need all the antibody columns, too
  if(grepl("comb", inpath)){
    all_antib = unique(xylist$name)
    name.df = list( Atype_ = rep(all_antib, each = nrow(finalx)))
    anames = sparse.model.matrix( ~ Atype_, name.df)[,-1]
    
    bigx = finalx
    for(i in 2:length(all_antib)){
      bigx = rbind(bigx, finalx)
    }
    finalx = cbind(bigx, anames)
    
    # name.df = list( Atype_ = xylist$name)
    # anames = sparse.model.matrix( ~ Atype_, name.df)[,-1]
    # xylist$name = rep("Antibody", length(xylist$name))
    # xylist$x = cbind(xylist$x, anames)
  }
  
  #get the prediction
  #for lots of predictions, do it in batches to avoid 'problem too large' error
  batch.size = 200000
  N = nrow(finalx)
  if(N > batch.size){
    ind.list = split(1:N, cut(1:N, c(seq(from = 0, to = N, by = batch.size), N) ))
    preds = unlist(lapply(ind.list, function(i){
      return(as.vector(mymodel %>% keras:::predict.keras.engine.training.Model(finalx[i, ])))
    }))
  } else {
    preds = as.vector(mymodel %>% keras:::predict.keras.engine.training.Model(finalx))
  }
  
  if(grepl("comb", inpath)){
    return(as.data.frame(list(aa_subs = rep(aa_subs, times = length(all_antib)), 
                         antibs = name.df$Atype_, preds = preds)))
  }
  return(preds)
  
}

#driver that runs binding, expression, and antibody predictions in succession, removing variants with
#poor expression or binding before the next step, then saves the results for plotting
#mutmax = c(1, 2)
#target = c("Omicron", "BA4_BA5")
# inpaths = c("ml_bind_keras_full_layer.sizes_c(64, 4, 16)_relu.alphas_c(0.01, 0.01, 0.01)_onefold0/Binding.h5",
#             "ml_express_keras_full_layer.sizes_c(16, 128)_relu.alphas_c(0.01, 0.01)_onefold/Expression.h5",
#             "ml_log_escape_keras_full_layer.sizes_c(256, 32, 16)_relu.alphas_c(0, 0, 0)_comb_onefold/Antibody.h5")
pred_multi_driver = function(inpaths, target = "Omicron", mutmax = 1){
  
  #starr
  # bind_min = -2.35
  # expr_min = -1.5
  bind_min = -1.5
  expr_min = -1
  #relax minimums, so that omicron/Ba5 predictions are caught
  # if(target %in% c("Omicron", "BA4_BA5")){
  #   #omicron/Ba5 predictions 
  #   bind_min = -1.5
  #   expr_min = -1
  # } else {
  #   #need higher minimums because ML algs
  #   bind_min = -2.35
  #   expr_min = -1.5
  # }

  
  var_muts = read.csv("output/var_muts.csv")
  targ_subs = var_muts$aa_subs[var_muts$Variant == target]
  
  #start with binding
  mut.df = pred_muts(target, mutmax, inpaths[1])
  mut.df$bind.pred = ml_predict(mut.df$new_subs, inpaths[1])
  
  #convert values to relative to target
  targ.pred = mut.df$bind.pred[mut.df$new_subs == targ_subs]
  mut.df$bind.pred = mut.df$bind.pred - targ.pred
  
  #apply cutoff
  mut.df$bind.pass = mut.df$bind.pred >= bind_min
  cat("Removing", sum(!mut.df$bind.pass), "of", nrow(mut.df), "variants with binding under", bind_min, 
      ", totaling", round(sum(!mut.df$bind.pass)/nrow(mut.df),3)*100, "%\n")
  
  #now expression for rest
  mut.df$expr.pred = NA
  mut.df$expr.pred[mut.df$bind.pass] = ml_predict(mut.df$new_subs[mut.df$bind.pass], inpaths[2])
  targ.pred = mut.df$expr.pred[mut.df$new_subs == targ_subs]
  mut.df$expr.pred = mut.df$expr.pred - targ.pred
  
  #apply cutoff
  mut.df$expr.pass = (mut.df$bind.pass & mut.df$expr.pred >= expr_min)
  cat("Removing", sum(mut.df$bind.pass & !mut.df$expr.pass), "of", sum(mut.df$bind.pass), 
      "variants with expression under", expr_min, ", totaling", 
      round(sum(mut.df$bind.pass & !mut.df$expr.pass)/sum(mut.df$bind.pass),3)*100, "%\n")
  
  #now antibodies
  antib.preds = ml_predict(mut.df$new_subs[mut.df$expr.pass], inpaths[3])
  antib.targ = antib.preds[antib.preds$aa_subs == targ_subs,]
  antib.preds$antib.pred = antib.preds$preds - antib.targ$preds[match(antib.preds$antibs, antib.targ$antibs)]
  
  #merge
  colnames(mut.df)[colnames(mut.df) %in% "new_subs"] = "aa_subs"
  outdf = merge(mut.df, antib.preds, by = "aa_subs", all = T, sort = F)
  
  #add inpaths
  outdf$inpaths = NA
  outdf$inpaths[1:3] = inpaths
  
  #write out
  outpath = paste0("output/multipreds/multipreds_", target, "_mrad", mutmax, ".csv")
  write.csv(outdf, file = outpath, row.names = F)

}

#generates all variants in mutmax number of mutations distance from target variant
#targets from "output/var_muts.csv
#target = "Omicron", "BA4_BA5", etc
pred_muts = function(target = "Omicron", mutmax = 1, inpath){
  
  var_muts = read.csv("output/var_muts.csv")
  targ_subs = var_muts$aa_subs[var_muts$Variant == target]
  
  allseqs = read.fasta("Datasets/var_to_align.fasta", seqtype = "AA")
  wuhan = allseqs$Wuhan
  targseq = allseqs[[target]]
  
  allaa = c('*','A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y')
  datatype = gsub("ml_", "", gsub("_keras.*", "", inpath))
  if(datatype == "bind") allaa = allaa[!allaa %in% "*"]
  targ_subs = strsplit(targ_subs, " ")[[1]]
  targ_locs = as.numeric(gsub("[[:alpha:]]|\\*", "", targ_subs))
  targ_muts = gsub(".*[[:digit:]]", "", targ_subs)
  
  #get all possible aa configurations
  m = length(wuhan)
  n = length(allaa)
  
  #remove configs in target
  allmuts = as.data.frame(list(locs = rep(1:m, each = n), aas = rep(allaa, times = m)))
  allmuts = as.data.frame(rbindlist(lapply(1:m, 
                       function(i){ return( allmuts[allmuts$locs == i & allmuts$aas != targseq[i],] ) } )))
  allmuts$orig = wuhan[allmuts$locs]
  outmuts = c(NA, paste0(allmuts$orig, allmuts$locs, allmuts$aas)) #this is the set of single mutations
  
  #formula for final # of mutations for given mutmax, m, n
  # N = 0
  # for(i in 1:mutmax){
  #   N = N + prod((m:1)[1:i])/prod(1:i)*(n-1)^i
  # }
  if(mutmax >= 1) {
    mcols = paste0("mut", 1:mutmax)
    lcols = paste0("loc", 1:mutmax)
  } else {
    mcols = "mut1"
    lcols = "loc1"
  }
  
  #expand.grid style: generate all permutations and then reduce to unique only
  #get all double, triple, etc mutations (output length is ~ 4000^mutmax)
  if(mutmax > 1) {
    uselist = list(mut1 = outmuts)
    for(nmut in 2:mutmax){
      uselist = c(uselist, list(outmuts))
    }
    mutdf = expand.grid(uselist, stringsAsFactors = F)
    colnames(mutdf) = mcols
    
    #add locations
    # for(i in 1:mutmax){
    #   mutdf[,paste0("loc", i)] = as.numeric(gsub("[[:alpha:]]|\\*", "", mutdf[,paste0("mut", i)]))
    # }
    locdf = mutdf
    locdf[] = lapply(mutdf, function(col){as.numeric(gsub("[[:alpha:]]|\\*", "", col))})
    colnames(locdf) = lcols
    mutdf = cbind(mutdf, locdf)

    #remove mutations at locations <= any previous locations, to remove duplicates
    to.remove = rep(F, nrow(locdf))
    for(i in 2:ncol(locdf)){
      #maximum of previous columns (= NA if previous are all NA)
      if(i == 2) prevmax = locdf[, i-1] else prevmax = do.call("pmax.int", c(as.list(locdf[, 1:i-1]), list(na.rm = T)))
      #remove if number <= non-NA previous max or if is missing and a previous value is not
      to.remove = to.remove | (!is.na(locdf[, i]) & !is.na(prevmax) & (locdf[, i] <= prevmax)) | 
        (is.na(locdf[, i]) & !is.na(prevmax))
    }
    finalout = mutdf[!to.remove, ]
      
  } else if(mutmax == 1){
    finalout = as.data.frame(list(mut1 = outmuts, loc1 = as.numeric(gsub("[[:alpha:]]|\\*", "", outmuts))))
  } else finalout = as.data.frame(list(mut1 = NA_character_, loc1 = NA_real_))
  
  #convert to aa_subs
  #all-in-one apply
  finalout$new_subs = apply(finalout, 1, function(x){
    #combine locs with target, overwriting target mutations
    new_locs = as.numeric(x[lcols])
    new_locs = new_locs[!is.na(new_locs)]
    all_locs = unique(c(new_locs, targ_locs))
    
    #combine muts with target, overwriting target mutations
    new_muts = x[mcols]
    new_muts = new_muts[!is.na(new_muts)]
    new_muts = gsub(".*[[:digit:]]", "", new_muts)
    all_muts = c(new_muts, targ_muts[!targ_locs %in% new_locs]) 
    
    #don't keep mutations back to wuhan
    start_aa = wuhan[all_locs]
    to.ignore = (all_muts == start_aa) 
    all_locs = all_locs[!to.ignore]
    all_muts = all_muts[!to.ignore]
    start_aa = start_aa[!to.ignore]
    
    #reorder by site and paste
    ord = order(all_locs)
    aa_subs = paste0(start_aa[ord], all_locs[ord], all_muts[ord])
    return(paste0(aa_subs, collapse = " "))
    
  })
  
  return(finalout)
  
  # if(grepl("comb", inpath)){
  #   output = ml_predict(finalout$new_subs, inpath)
  #   output$pred_escape_fraction = 10^output$preds
  #   bigout = finalout
  #   for(i in 2:(nrow(output)/nrow(finalout))){
  #     bigout = rbind(bigout, finalout)
  #   }
  #   if(any(bigout$new_subs != output$aa_subs)) stop("non-matching subs in combined antibody prediction")
  #   to.write = cbind(bigout[,c(mcols, "new_subs")], output[,c("antibs", "preds", "pred_escape_fraction")])
  #   
  # } else {
  #   finalout$preds = ml_predict(finalout$new_subs, inpath)
  #   to.write = finalout[, c(mcols, "new_subs", "preds")]
  # }
  # 
  # outpath = gsub("/.*", paste0("/preds_", target, "_mrad", mutmax, ".csv"), inpath)
  # 
  # write.csv(to.write, file = paste0("output/", outpath), row.names = F)
  
}

#out.type = c("efrac", "1mefrac", "logefrac)
multi_plots = function(target = "Omicron", mutmax = 2, out.type = "efrac"){
  
  n = 10
  min.bind = 0
  min.expr = 0
  inpath = paste0("output/multipreds/multipreds_", target, "_mrad", mutmax)
  
  indf = fread(paste0(inpath, ".csv"))
  
  pdf(file = paste0(inpath, "_", out.type, ".pdf"), width = 8.5, height = 11)
  par(mfrow = c(2,1))
  
  usedf = indf[bind.pass & expr.pass,]
  usedf = usedf[usedf$bind.pred >= min.bind & usedf$expr.pred >= min.expr,]
  
  #correct preds > 0 back to 0
  usedf$preds[usedf$preds > 0] = 0
  targdf = usedf[is.na(usedf$mut1) & is.na(usedf$mut2),]
  usedf$antib.pred = usedf$preds - targdf$preds[match(usedf$antibs, targdf$antibs)]
  
  usedf$newmut1 = gsub(".*[[:digit:]]", "", usedf$mut1)
  usedf$newmut2 = gsub(".*[[:digit:]]", "", usedf$mut2)
  usedf$allmuts = paste0(usedf$loc1 + 330, usedf$newmut1, " ", usedf$loc2 + 330, usedf$newmut2)
  # usedf$allmuts = paste0(usedf$mut1, " ", usedf$mut2)
  
  #Better in binding/expression
  mean.escape = as.data.frame(usedf[, .(mean.escape = mean(preds)), by = .(allmuts)])
  mean.escape = mean.escape[order(mean.escape$mean.escape, decreasing = T),]
  
  to.cast = usedf[usedf$allmuts %in% c(mean.escape$allmuts[1:n], "NANA NANA"), c("allmuts", "antibs", "preds")]
  to.plot = dcast(to.cast, allmuts ~ antibs, value.var = "preds")
  to.plot = to.plot[match(c(mean.escape$allmuts[1:n], "NANA NANA"), to.plot$allmuts),]
  keepnames = to.plot$allmuts
  keepnames[keepnames == "NANA NANA"] = target
  to.plot$allmuts = NULL
  
  cols = c("black", brewer.pal(9, "Set1"))
  if(out.type == "efrac"){
    final.plot = 10^t(to.plot)
    ylab = "Escape Fraction"
    use.ylim = c(0,1.5)
  } else if(out.type == "1mefrac"){
    final.plot = 1 - 10^t(to.plot)
    ylab = "1 - Escape Fraction"
    use.ylim = c(0,1.5)
  } else if(out.type == "logefrac"){
    final.plot = apply(t(to.plot), 2, function(x){x - t(to.plot)[, 11]})
    ylab = "log10(Escape Fraction)"
  }
  colnames(final.plot) = keepnames
  barplot(final.plot, beside = T, xlab = "Variant", ylab = ylab, col = cols, ylim = use.ylim, cex.names = .5,
          main = paste0("Top 10 ", target," Variants in Mean Log Antibody Escape\nBinding/Expression >= ", target,
                        "\nMaximum of ", mutmax, " Mutations"))
  legend("topright", legend = rownames(final.plot), fill = cols, ncol = 2, cex = .8)
  
  #better in every way
  bestdf = usedf[usedf$antib.pred >= 0,]
  besttab = table(bestdf$allmuts)
  keepmuts = names(besttab[besttab == 10])
  bestdf = bestdf[bestdf$allmuts %in% keepmuts,]
  
  mean.escape = as.data.frame(bestdf[, .(mean.escape = mean(preds)), by = .(allmuts)])
  mean.escape = mean.escape[order(mean.escape$mean.escape, decreasing = T),]
  
  to.cast = usedf[usedf$allmuts %in% c(mean.escape$allmuts[1:n], "NANA NANA"), c("allmuts", "antibs", "preds")]
  to.plot = dcast(to.cast, allmuts ~ antibs, value.var = "preds")
  to.plot = to.plot[match(c(mean.escape$allmuts[1:n], "NANA NANA"), to.plot$allmuts),]
  keepnames = to.plot$allmuts
  keepnames[keepnames == "NANA NANA"] = target
  to.plot$allmuts = NULL
  
  cols = c("black", brewer.pal(9, "Set1"))
  if(out.type == "efrac"){
    final.plot2 = 10^t(to.plot)
    ylab = "Escape Fraction"
    use.ylim = c(0,1.5)
  } else if(out.type == "1mefrac"){
    final.plot2 = 1 - 10^t(to.plot)
    ylab = "1 - Escape Fraction"
    use.ylim = c(0,1.5)
  } else if(out.type == "logefrac"){
    final.plot2 = t(to.plot)
    ylab = "log10(Escape Fraction)"
  }
  colnames(final.plot2) = keepnames
  barplot(final.plot2, beside = T, xlab = "Variant", ylab = ylab, col = cols, names.arg = keepnames, 
          ylim = use.ylim, cex.names = .5, main = paste0("Top 10 ", target, 
                        " Variants in Mean Log Antibody Escape\nBinding/Expression/All Antibody Escape >= ", 
                        target, "\nMaximum of ", mutmax, " Mutations"))

  legend("topright", legend = rownames(final.plot2), fill = cols, ncol = 2, cex = .8)

  dev.off()
  
  write.xlsx(x = list(bind_exp_ge_targ = final.plot, bind_exp_esc_ge_targ = final.plot2), 
             file = paste0(inpath, "_", out.type, ".xlsx"), rowNames = T)
}

#makes both barplot version
antib_var_bar_driver = function(infolder, out.type = "efrac"){
  
  
  pdf(file = paste0("output/", infolder, "/antib_bars_", out.type,".pdf"), width = 8.5, height = 5.5)
  # par(mfrow = c(2,1))
  for(uselog in c(F)){
    out = antib_var_bar(infolder, uselog = uselog, out.type = out.type)
    write.csv(out, paste0("output/", infolder, "/antib_bars_", out.type,".csv"))
  }
  
  dev.off()
}

#plots antib escape fractions for all possible variants 
# infolder = "ml_log_escape_keras_full_layer.sizes_c(256, 32, 16)_relu.alphas_c(0, 0, 0)_comb_onefold"
antib_var_bar = function(infolder, uselog = F, out.type = "efrac"){
  
  if(uselog) endpoint = "preds" else endpoint = "pred_escape_fraction"
  
  
  wuhan.in = read.csv(paste0("output/", infolder, "/preds_Wuhan_mrad1.csv"))
  wuhan.in = wuhan.in[wuhan.in$new_subs %in% "",]
  
  var.vec = c("Alpha_V1", "Beta", "Gamma_V3", "Delta", "Kappa", "Omicron", "BA4_BA5")
  
  varmat = sapply(var.vec, function(variant){
    inmat = read.csv(paste0("output/", infolder,"/preds_", variant, "_mrad0.csv"))
    out = inmat[, endpoint]
    names(out) = inmat$antibs
    return(out)
  })
  
  Wuhan = wuhan.in[, endpoint]
  names(Wuhan) = wuhan.in$antibs
  varmat = cbind(Wuhan, varmat)
  
  if(uselog) use.ylab = "Log Escape Fraction" else use.ylab = "Escape Fraction"
  if(uselog) use.ylim = c(-5,0) else use.ylim = c(0,1)
  if(out.type == "1mefrac") {
    use.ylab = paste0("1 - ", use.ylab)
    varmat = 1 - varmat
    use.ylim = c(0, 1.5)
  }
  cols = c("black", brewer.pal(9, "Set1"))
  barplot(varmat, beside = T, xlab = "Variant", ylab = use.ylab, col = cols, ylim = use.ylim)
  legend("topleft", legend = rownames(varmat), fill = cols, ncol = 2)
  
  return(varmat)
}

#take the predictions for two mutations and display as heatmap
#filename = "ml_bind_keras_full_layer.sizes_c(64, 4, 16)_relu.alphas_c(0.01, 0.01, 0.01)_onefold0/preds_Omicron_mrad2.csv"
ml_heatmap = function(filename){
  
  #get data and cast to two dimensional df
  indata = read.csv(paste0("output/", filename))
  indata$new_subs = NULL
  cdata = as.data.frame(dcast(as.data.table(indata), mut1 ~ mut2, value.var = "preds", fill = NA, drop = F))
  
  #reformat/reorder
  rownames(cdata) = cdata$mut1
  cdata$mut1 = NULL
  cdata = cdata[reord(rownames(cdata)),]
  cdata = cdata[, c(1 + reord(colnames(cdata)[-1]), 1)]
  heatmap(as.matrix(cdata), Rowv = NA, Colv = NA)
  heatmap.2(as.matrix(cdata[1:300, 1:300]), Rowv = F, Colv = F, col = hcl.colors(20,"viridis"), trace = "none",
            dendrogram = "none", symbreaks = F)
  
}

ml_pred_sort = function(filename){
  
  #get data and cast to two dimensional df
  indata = read.csv(paste0("output/", filename))
  indata$new_subs = NULL
  
  secdat =  indata[,c("mut2", "preds")]
  colnames(secdat)[1] = "mut1"
  sing.data = rbind(as.data.table(indata[,c("mut1", "preds")]), as.data.table(secdat))
  sing.data = sing.data[!is.na(sing.data$mut1),]
  sing.data$site = as.numeric(gsub("[[:alpha:]]|\\*", "", sing.data$mut1))
  
  sing.agg  = as.data.frame(sing.data[, .(mpreds = mean(preds)), by = mut1])
  sing.agg = sing.agg[order(sing.agg$mpreds, decreasing = T),]
  sing.data = as.data.frame(sing.data)
  
  top.data = sing.data[sing.data$mut1 %in% sing.agg$mut1[1:20],]
  top.data$mut1 = factor(top.data$mut1, levels = sing.agg$mut1[1:20])
  # top.data$mpreds = sing.agg$mpreds[match(top.data$mut1, sing.agg$mut1)]
  # top.data = top.data[order(top.data$mpreds, decreasing = T),]
  boxplot(preds ~ mut1, top.data)

  # umuts = unique(indata$mut1)
  # sing.preds = sapply(umuts, function(mut){
  #   return(mean(indata$preds[indata$mut1 %in% mut | indata$mut2 %in% mut]))
  # })
  # newdf = as.data.frame(list(umuts = umuts, sing.preds = sing.preds))
  # newdf = newdf[order(newdf$sing.preds, decreasing = T),]
  
  #by site
  site.agg =  as.data.frame(sing.data[, .(mpreds = mean(preds)), by = site])
  site.agg = site.agg[order(site.agg$mpreds, decreasing = T),]
  top.data = sing.data[sing.data$site %in% site.agg$site[1:20],]
  top.data$site = factor(top.data$site, levels = site.agg$site[1:20])
  
  boxplot(preds ~ site, top.data)
  
}

reord = function(cnames){
  locs = as.numeric(gsub("[[:alpha:]]|\\*", "", cnames))
  muts = gsub(".*[[:digit:]]", "", cnames)
  return(order(locs, muts))
}


