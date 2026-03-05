library(data.table)
library(openxlsx)
library(protr)
source("R/combined_ml.R")
source("R/score_funs.R")

#model driver, save plots to file

same_holdout_plotter = function(plot.type = "scatter", ml = "rf", feats = "mixed"){
  
  set.seed(1234)
  
  #first get variants only in version 3
  indata = read.csv("output/mai_antib_kds_12_20_2022.csv")
  #v2
  newdata = read.csv("output/mai_antib_kds_02_09_2023.csv")
  indata = rbind(indata, newdata)
  v12_vars = unique(indata$Variant)
  #v3
  newdata = read.csv("output/mai_antib_kds_09_18_2023.csv")
  newdata$Variant[newdata$Variant == "Wuhan_control"] = "Wuhan" #just average it into other measurements
  indata = rbind(indata, newdata)
  new_vars = setdiff(unique(indata$Variant), v12_vars)
  holdout.vars = sample(new_vars, floor(length(new_vars)/5))

  pdf(file = paste0("output/same_holdouts_", plot.type, "_", ml, "_", feats, ".pdf"))
  par(mfrow = c(3,2))
  
  cols = c("blue", "darkgreen", "purple4")
  for(version in c(1,2,3)){
        if(version == 3) keep.var.feats = 10 else keep.var.feats = 2
        mai_antib_model_v2(ml = ml, feats = feats, delta.kd = "wuhan", version = version, holdouts = holdout.vars,
                        plot.type = plot.type, keep.var.feats = keep.var.feats, use.col = cols[version])
  }
  
  
  dev.off()
  
}

get_indata = function(version = 3, delta.kd = "none"){
  indata = read.csv("output/mai_antib_kds_12_20_2022.csv")
  if(version >= 2){
    newdata = read.csv("output/mai_antib_kds_02_09_2023.csv")
    indata = rbind(indata, newdata)
  }
  #add in 08_18_2023 data
  if(version >= 3){
    newdata = read.csv("output/mai_antib_kds_09_18_2023.csv")
    newdata$Variant[newdata$Variant == "Wuhan_control"] = "Wuhan" #just average it into other measurements
    indata = rbind(indata, newdata)
  }
  
  indata = indata[!indata$Log_Kd %in% "Constant",]
  indata$endpoint = as.numeric(indata$Log_Kd)
  indata$uncert = as.numeric(indata$Log_Kd_SE)
  
  if(version >= 2){
    #avg the endpoint
    endagg = aggregate(endpoint ~ Variant + Antibody, indata, mean)
    uncertagg = aggregate(uncert ~ Variant + Antibody, indata, function(x){sqrt(sum(x^2))/length(x)}, na.action = na.pass)
    endagg$uncert = uncertagg$uncert[match(paste0(endagg$Variant, "_", endagg$Antibody), 
                                           paste0(uncertagg$Variant, "_", uncertagg$Antibody))]
    indata = endagg
    
  }
  
  #delta kd relative to wuhan
  if(delta.kd == "wuhan"){
    to.sub = indata[indata$Variant == "Wuhan",]
    indata$wuhan = to.sub$endpoint[match(indata$Antibody, to.sub$Antibody)]
    indata = indata[!is.na(indata$wuhan) & indata$Variant != "Wuhan",]
    indata$endpoint = indata$endpoint - indata$wuhan
  } else if(delta.kd == "ace2"){
    #delta kd relative to ACE2
    to.sub = indata[indata$Antibody == "huACE2-Fc",]
    indata$ace2 = to.sub$endpoint[match(indata$Variant, to.sub$Variant)]
    indata = indata[!is.na(indata$ace2) & indata$Antibody != "huACE2-Fc",]
    indata$endpoint = indata$endpoint - indata$ace2
  }
  indata = indata[!indata$Antibody %in% c("huACE2-Fc", "IgG"),] #can't use ACE2 under any circumstance
  
  seqdata = read.xlsx("Datasets/Library of 59 antibodies_ Binding and sequences_MP 12202022.xlsx", 
                      sheet = "Antibody sequences", cols = 16:19)
  colnames(seqdata)[colnames(seqdata) == "Name"] = "Antibody"
  indata = merge(indata, seqdata, by = "Antibody") #removes IgG and ACE2 
  
  return(indata)
}

#This version adapted specifically for holdouts and simplified
mai_antib_model_revision = function(ml = "rf", feats = "mixed", keep.var.feats = 2, wts = "none", max.antib.feats = 24, split.cdr = F, variant.fold = F, 
                           seed = 1234, filename = NULL, cdr.len.feat = F, delta.kd = "none", version = 3, 
                           plot.type = "scatter", new.preds = "holdouts", holdouts = NULL, use.col = "black", 
                           vars = "213_variants", do.plot = T, holdout.version = 3, ...){
  
  tensorflow::set_random_seed(seed)
  inlist = list(...)
  ml.params = inlist[["ml.params"]]
  if(ml %in% "xgboost") ml.params = c(ml.params, list(nrounds = 20))
  
  indata = get_indata(version = version, delta.kd = delta.kd)

  if(!is.null(holdouts) ){
    holddata = get_indata(version = holdout.version, delta.kd = delta.kd)
    indata = indata[!indata$Variant %in% holdouts,]
    holddata = holddata[holddata$Variant %in% holdouts,]
  }
  
  #features
  xylist = mai_antib_features(indata, feats, keep.var.feats, wts, max.antib.feats, split.cdr, cdr.len.feat, vars = vars)

  if(!is.null(holdouts) ){
    xylist.holdout = mai_antib_features(holddata, feats, keep.var.feats, wts, max.antib.feats, split.cdr, cdr.len.feat, vars = vars)
    preds = modeler(xin = xylist$x, yin = xylist$y, xout = xylist.holdout$x, ml.params = ml.params, ml = ml, filename = filename, ...)
    pred.df = data.frame(preds = preds, yout = xylist.holdout$y)
    } else {
    if(variant.fold) splitvec = as.numeric(as.factor(xylist$variant)) else splitvec = NULL
    pred.df = crossval(x = xylist$x, y = xylist$y, splitvec = splitvec, ml.params = ml.params, ml = ml, ...)
  }

  #Plotting
  
  # pred.df$uncert = xylist$uncert
  # pred.df$weight = xylist$weight
  # pred.df$variant = xylist$variant
  
  if(!do.plot) return(pred.df)
  
  if(delta.kd == "none"){
    xlab = "Experimental Log(Kd)"
    ylab = "Predicted Log(Kd)"
    end = "None"
  } else if(delta.kd == "wuhan"){
    xlab = "Experimental Delta(Log(Kd))"
    ylab = "Predicted Delta(Log(Kd))"
    end = "Wuhan"
  } else if(delta.kd == "ace2"){
    xlab = "Experimental Delta(Log(Kd))"
    ylab = "Predicted Delta(Log(Kd))"
    end = "ACE2"
  }
  if(variant.fold) add = ", Var CV" else add = ""
  
  
  if(plot.type == "scatter"){
    plot(pred.df$yout, pred.df$preds, xlab = xlab, ylab = ylab, col = use.col,
         main = paste0("Kd Delta: ", end, ", Version: ", version, add,
                       "\nRMSE = ", round(rmse(pred.df$yout - pred.df$preds), 2), 
                       ", cor = ", round(cor(pred.df$yout, pred.df$preds), 2), 
                       ", Q^2 = ", round(q2(pred.df$yout, pred.df$preds), 2) ))
    
    abline(0, 1, lwd = 3)
  }
  
  
}

mai_antib_features = function(indata, feats, keep.var.feats = 4, wts = "none", max.antib.feats = 24, split.cdr = F,
                              cdr.len.feat = F, vars = "213_variants"){

  if(vars == "std") pdata = read.csv("output/mai_antib_model/var_preds.csv") #this is no longer necessary, since 213 has all these too
  if(vars == "213_variants") pdata = read.csv("output/mai_antib_model/213_variant_preds.csv")
  
  pdata = pdata[pdata$Variant %in% unique(indata$Variant),]
  pdata$aa_subs = NULL
  pdata$antibs = gsub("Antibody COV2-|Antibody CR", "A", pdata$antibs)
  pdata$antibs = gsub("_400", "", pdata$antibs)
  #based on looking at ml_log_escape_keras_full_layer.sizes_c(256, 32, 16)_relu.alphas_c(0, 0, 0)_comb_onefold/antib_bars.pdf
  keep.ord = c("A2050", "A2082", "A2094", "A2096", "A2499", "A2677", "A2832", "A2479", "A2165", "A3022") 
  pdata = pdata[pdata$antibs %in% keep.ord[1:keep.var.feats],]
  
  add.data = as.data.frame(dcast(as.data.table(pdata), Variant ~ antibs, value.var = "preds"))
  vfeats = colnames(add.data)[!colnames(add.data) %in% c("Variant")]
  fulldata = merge(indata, add.data, by = "Variant") 
  
  ################Next, get protr features for antibody cdrs, keep highly customizable

  if(split.cdr){
    featlist = lapply(c("CDR1", "CDR2", "CDR3"), function(colname){
      featset = t(sapply(fulldata[,colname], get_feats, feats = feats, max.antib.feats = floor(max.antib.feats/3),
                         min.len = min(sapply(fulldata[,colname], nchar)), cdr.len.feat = cdr.len.feat))
      colnames(featset) = paste0(colname, ".", colnames(featset))
      return(featset)
    })
    featset = do.call(cbind, featlist)
  } else {
    fulldata$ALLCDR = paste0(fulldata$CDR1, fulldata$CDR2, fulldata$CDR3)
    featset = t(sapply(fulldata$ALLCDR, get_feats, feats = feats, max.antib.feats = max.antib.feats, 
                       min.len = min(sapply(fulldata$ALLCDR, nchar)), cdr.len.feat = cdr.len.feat))
  }
  
  afeats = colnames(featset)
  fulldata = cbind(fulldata, featset)
  
  ###############weights
  x = as.matrix(fulldata[ ,c(vfeats, afeats)])
  weight = get_w(wt.type = wts, uncert = fulldata$uncert)
  if(!is.null(weight)) x = cbind(x, weight)

  return(list(x = x, y = fulldata$endpoint, uncert = fulldata$uncert, weight = weight, antibody = fulldata$Antibody, 
              variant = fulldata$Variant))
}


#get feature vector for give naa string, feats, and maximum
get_feats = function(x, feats, max.antib.feats = NULL, min.len = NULL, cdr.len.feat = F) {
  
  if(feats %in% "mixed"){
    # f1 = extractMoreauBroto(x, nlag = 1)
    # names(f1) = paste0("MB.", names(f1))
    # f2 = extractCTDT(x)
    # names(f2) = paste0("CTDT.", names(f2))
    f3 = extractSOCN(x, nlag = 4)
    names(f3) = paste0("SOCN.", names(f3))
    f4 = extractDescScales(x, propmat = "AATopo", index = c(37:41, 43:47), lag = 2, pc = 2)
    names(f4) = paste0("DScales.", names(f4))
    # f5 = extractBLOSUM(x, lag = 2, k = 2)
    # names(f5) = paste0("BLOSUM.", names(f5))
    # f6 = extractProtFP(x, index = c(160:165, 258:296), pc = 2, lag = 2)
    # names(f6) = paste0("ProtFP.", names(f6))
    out = c(f3, f4)
  } else {
    if(feats %in% c("MoreauBroto", "Moran", "Geary")) {
      args = list(nlag = floor(max.antib.feats/8))
    } else if(feats %in% c("CTDC", "CTDT")){
      args = NULL
    } else if(feats == "SOCN"){
      args = list(nlag = min(c(floor(max.antib.feats/2), min.len - 1)))
    } else if(feats == "DescScales"){
      args = c(list(propmat = "AATopo", index = c(37:41, 43:47)), get_pars(max.antib.feats))
    } else if(feats == "BLOSUM"){
      args = get_pars(max.antib.feats, usenames = c("lag", "k"))
    } else if(feats == "ProtFP"){
      args = c(list(index = c(160:165, 258:296)), get_pars(max.antib.feats))
    }
    out = do.call(paste0("extract", feats), c(list(x = x), args))
  }
  
  if(cdr.len.feat) {
    addout = nchar(x)
    names(addout) = "CDR.length"
    out = c(out, addout)
  }

  
  return(out)
}

#get highest balanced parameters with given ceiling 
get_pars = function(ceil, usenames = c("lag", "pc")){
  usex = 1
  usey = 1
  while(use_fun(usex, usey) < ceil){
    tryx = usex + 1
    if(use_fun(tryx, usey) > ceil) break else usex = tryx
    tryy = usey + 1
    if(use_fun(usex, tryy) > ceil) next else usey = tryy
  }
  
  to.return = list(usex, usey)
  names(to.return) = usenames
  return(to.return)
}

use_fun = function(x, y){ return(x*y^2)}

#given a single aa string target and a vector of variant aa strings, get the space separated targetlocationmutant aa substitutions
#assume same length
seq_diffs = function(targ, vars){
  
  targsep = unlist(strsplit(targ, ""))
  varsep = strsplit(vars, "")
  outvec = sapply(varsep, function(var, targsep){
    locs = which(var != targsep)
    out = paste0(targsep[locs], locs, var[locs], collapse = " ")
    return(out)
  }, targsep = targsep)
  
  return(outvec)
}



