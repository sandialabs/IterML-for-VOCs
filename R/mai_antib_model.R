library(data.table)
library(openxlsx)
library(protr)
source("R/combined_ml.R")
source("R/score_funs.R")
source("R/ml_predict.R")

#model driver, save plots to file

mai_antib_model_plotter = function(plot.type = "scatter", ml = "rf", feats = "mixed"){
  

  pdf(file = paste0("output/antib_model_", plot.type, "_", ml, "_", feats, "_new.pdf"))
  par(mfrow = c(3,2))

  cols = c("blue", "darkgreen", "purple4")
  for(version in c(1,2,3)){
    for(delta.kd in c("none", "wuhan", "ace2")){
      for(variant.fold in c(F, T)){
        if(version == 3) keep.var.feats = 10 else keep.var.feats = 2
        mai_antib_model(ml = ml, feats = feats, delta.kd = delta.kd, variant.fold = variant.fold, version = version,
                    plot.type = plot.type, keep.var.feats = keep.var.feats, use.col = cols[version])
      }
    }
  }

  
  dev.off()
  
}

#does the four plots from the paper A,B,C,D
#(A) CV with both Omicron/BA.5 held out
#(B) Omicron/BA5 held out and Omicron tested
#(C) CV with BA.5 held out
#(D) BA5 held out and BA5 tested
holdout_plots = function(){
  
  use.cols = c("blue", "firebrick", "darkgreen", "purple4")
  
  #A
  adata = mai_antib_model(holdouts = c("Omicron", "BA4_BA5"), version = 2, do.plot = F)
  
  #B
  bdata = read.csv("output/mai_antib_model/rf_mixed_holdout_omi_ba5.csv")
  bdata$yout = bdata$endpoint
  bdata$preds = bdata$log_Prediction
  
  #C
  cdata = mai_antib_model(holdouts = c("BA4_BA5"), version = 2, do.plot = F)
  
  #D
  ddata = read.csv("output/mai_antib_model/rf_mixed_holdout_ba5.csv")
  ddata$yout = ddata$endpoint
  ddata$preds = ddata$log_Prediction
  
  alldata = c("adata", "bdata", "cdata", "ddata")
  
  #do it as a four part figure, and separately in pdf
  # pdf(file = "output/mai_antib_model/manuscript_holdout_plots.pdf", width = 8.5, height = 11)
  tiff(filename = "output/mai_antib_model/manuscript_holdout_plots_all4.tif", width = 8, height = 8,
       units = "in", res = 300)
  par(mfrow = c(2,2))
  par(mar = c(5.1, 5.1, 5.1, 1.1))
  
  for(i in 1:4){
    usedata = get(alldata[i])
    plot(usedata$yout, usedata$preds, xlab = expression("Experimental Log"[10]*"(K"[D]*")"),
           ylab = expression("Predicted Log"[10]*"(K"[D]*")"), col = use.cols[i],
           main = paste0("Kd Delta: None, Version: 2",
                         "\nRMSE = ", round(rmse(usedata$yout - usedata$preds), 2), 
                         ", cor = ", round(cor(usedata$yout, usedata$preds), 2), 
                         ", Q^2 = ", round(q2(usedata$yout, usedata$preds), 2) ))
    box(lwd = 2)
    abline(0, 1, lwd = 3) 
  }
  
  dev.off()
  
  for(i in 1:4){
    use.label = c("a", "b", "c", "d")
    
    tiff(filename = paste0("output/mai_antib_model/manuscript_holdout_plot_", use.label[i], ".tif"), 
         width = 4, height = 4, units = "in", res = 300)
    par(mfrow = c(1,1))
    
    usedata = get(alldata[i])
    plot(usedata$yout, usedata$preds, xlab = expression("Experimental Log"[10]*"(K"[D]*")"),
         ylab = expression("Predicted Log"[10]*"(K"[D]*")"), col = use.cols[i],
         main = paste0("Kd Delta: None, Version: 2",
                       "\nRMSE = ", round(rmse(usedata$yout - usedata$preds), 2), 
                       ", cor = ", round(cor(usedata$yout, usedata$preds), 2), 
                       ", Q^2 = ", round(q2(usedata$yout, usedata$preds), 2) ))
    box(lwd = 2)
    abline(0, 1, lwd = 3) 
    
    dev.off()
  }
  
  
  
  #do it as separate .png or .tif files
  
}

#for the reviewer response: do a plot as in holdout d, except use BA.1 data to predict BA.5 data
review_response = function(){
  
  ddata = read.csv("output/mai_antib_model/rf_mixed_holdout_ba5.csv")
  ddata$yout = ddata$endpoint
  ddata$preds = ddata$log_Prediction
  
  indata = read.csv("output/mai_antib_kds_12_20_2022.csv")
  newdata = read.csv("output/mai_antib_kds_02_09_2023.csv")
  indata = rbind(indata, newdata)
  newdata = read.csv("output/mai_antib_kds_09_18_2023.csv")
  newdata$Variant[newdata$Variant == "Wuhan_control"] = "Wuhan" #just average it into other measurements
  indata = rbind(indata, newdata)
  indata = indata[!indata$Log_Kd %in% "Constant",]
  indata$endpoint = as.numeric(indata$Log_Kd)
  indata$uncert = as.numeric(indata$Log_Kd_SE)
  indata = aggregate(endpoint ~ Variant + Antibody, indata, mean)

  omidata = indata[indata$Variant == "Omicron", ]
  
  ddata$omi.pred = omidata$endpoint[match(ddata$Antibody, omidata$Antibody)] 
  
  
  tiff(filename = paste0("output/mai_antib_model/review_response.tif"), 
       width = 4, height = 4, units = "in", res = 300)
  
  plot(ddata$yout, ddata$omi.pred, xlab = expression("Omicron BA.4/BA.5 Experimental Log"[10]*"(K"[D]*")"),
       ylab = expression("Omicron BA.1 Experimental Log"[10]*"(K"[D]*")"), col = "purple4",
       main = paste0("BA.4/BA.5 Predicted using BA.1 Data",
                     "\nRMSE = ", round(rmse(ddata$yout - ddata$omi.pred), 2), 
                     ", cor = ", round(cor(ddata$yout, ddata$omi.pred), 2), 
                     ", Q^2 = ", round(q2(ddata$yout, ddata$omi.pred), 2) ))
  box(lwd = 2)
  abline(0, 1, lwd = 3) 
  
  dev.off()
}


#feats = c("MoreauBroto", "Moran", "Geary", "CTDC", "CTDT", "SOCN", "DescScales", "BLOSUM", "ProtFP", "mixed")
#delta.kd = c("none", "wuhan", "ace2")
#plot.type = c("none", "scatter", "bar")
# new.preds = c("miss.omi", "monoclonal")
#
mai_antib_model = function(ml = "rf", feats = "mixed", keep.var.feats = 2, wts = "none", max.antib.feats = 24, split.cdr = F, variant.fold = F, 
                           seed = 1234, filename = NULL, cdr.len.feat = F, delta.kd = "none", version = 3, 
                           plot.type = "scatter", new.preds = "holdouts", holdouts = NULL, use.col = "black", 
                           vars = "213_variants", do.plot = T, ...){
  
  tensorflow::set_random_seed(seed)
  inlist = list(...)
  ml.params = inlist[["ml.params"]]
  if(ml %in% "xgboost") ml.params = c(ml.params, list(nrounds = 20))
  
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
  
  #database data, assumes no repeats of existing antibodies in data
  # if(version >= 3){
  #   #col 3 is ab type, col 23 is cdr3 length
  #   newdata = read.xlsx("Datasets/Database for SAR2- full length antibody and HCAbs _ MP 03092023.xlsx",
  #                       cols = c(2, 6, 14, 20:22))
  #   newdata2 = read.xlsx("Datasets/Database for SAR2- full length antibody and HCAbs _ MP 03092023.xlsx",
  #                        cols = c(2, 6, 14, 20:22), sheet = "HC only IgG")
  #   newdata = rbind(newdata, newdata2)
  #   
  #   newdata$endpoint = log10(as.numeric(newdata$`Binding.Affinity.(KD).(nM)`))
  #   newdata = newdata[!is.na(newdata$endpoint), ]
  #   colnames(newdata)[colnames(newdata) == "Antibody.name"] = "Antibody"
  #   colnames(newdata)[colnames(newdata) == "Variant.subtypes"] = "Variant"
  #   
  #   #All data is from WT, AKA Wuhan, except 5 with multiple tags. The Hong- antibs acutally represent WT, n3113.1 not so clear.
  #   newdata[grepl("Hong", newdata$Antibody), "Variant"] = "WT"
  #   newdata = newdata[newdata$Variant == "WT", ]
  #   newdata$Variant = "Wuhan"
  #   newdata$`Binding.Affinity.(KD).(nM)` = NULL
  #   
  #   #Two CDR2s have a blank or lowercase l; discard them
  #   cdr2s = strsplit(newdata$CDR2, "")
  #   to.disc = sapply(cdr2s, function(x){any(x %in% c(" ", "l"))})
  #   newdata = newdata[!to.disc, ]
  #   
  #   #resolve same cdr strings (they are not necessarily the same)
  #   if(length(intersect(paste(indata$CDR1, indata$CDR2, indata$CDR3, sep = "_"), 
  #                       paste(newdata$CDR1, newdata$CDR2, newdata$CDR3, sep = "_"))) > 0) stop("Duplicate CDRs in different datasets need to be resolved")
  #   ends = aggregate(endpoint ~ CDR1 + CDR2 + CDR3 + Variant, newdata, mean)
  #   abs = aggregate(Antibody ~ CDR1 + CDR2 + CDR3 + Variant, newdata, paste, collapse = "__")
  #   ends$Antibody = abs$Antibody[match(paste(ends$CDR1, ends$CDR2, ends$CDR3, sep = "_"), 
  #                                      paste(abs$CDR1, abs$CDR2, abs$CDR3, sep = "_"))]
  #   
  #   ends$uncert = NA
  #   
  #   ends = ends[, match(colnames(indata), colnames(ends))]
  #   
  #   indata = rbind(indata, ends)
  # }
  
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
  
  if(!is.null(holdouts) ){
    holddata = indata[indata$Variant %in% holdouts,]
    indata = indata[!indata$Variant %in% holdouts,]
  }
  
  #features
  xylist = mai_antib_features(indata, feats, keep.var.feats, wts, max.antib.feats, split.cdr, cdr.len.feat, vars = vars)
  
  # if(!is.null(filename)){
  #   no.out = modeler(xin = xylist$x, yin = xylist$y, ml.params = ml.params, 
  #                    filename = paste0("output/mai_antib_model/", filename), ...)
  #   return()
  # }

  if(!is.null(filename)){
    filename = paste0("output/mai_antib_model/", filename)
    
    # vars = "std"
    if(new.preds == "miss.omi"){
      all_antibs = unique(indata$Antibody)
      miss.omicron = setdiff(all_antibs, indata$Antibody[indata$Variant == "Omicron"])
      
      newdata1 = as.data.frame(list(Variant = "Omicron", Antibody = miss.omicron, endpoint = 0, uncert = 0))
      newdata2 = as.data.frame(list(Variant = "BA4_BA5", Antibody = all_antibs, endpoint = 0, uncert = 0))
      usedata = rbind(newdata1, newdata2)
      
      seqdata = read.xlsx("Datasets/Library of 59 antibodies_ Binding and sequences_MP 12202022.xlsx", 
                          sheet = "Antibody sequences", cols = 16:19)
      colnames(seqdata)[colnames(seqdata) == "Name"] = "Antibody"
      usedata = merge(usedata, seqdata, by = "Antibody")
    } else if(new.preds == "monoclonal"){
      antibs = unlist(read.xlsx("Datasets/Sequencing and monoclonal ELISA data from biopanning agianst Omicron BA1 _MP 03032023.xlsx", 
                                cols = 10, startRow = 2, colNames = F))
      vars = c("Wuhan", "Alpha", "Beta", "Delta", "Omicron")
      usedata = expand.grid(vars, antibs, stringsAsFactors = F)
      colnames(usedata) = c("Variant", "Antibody")
      usedata$endpoint = usedata$uncert = 0

      monodata = read.xlsx("Datasets/Sequencing and monoclonal ELISA data from biopanning agianst Omicron BA1 _MP 03032023.xlsx",
                           cols = c(10, 18:20))
      colnames(monodata)[colnames(monodata) == "Name"] = "Antibody"
      usedata = merge(usedata, monodata, by = "Antibody")
    } else if(new.preds == "holdouts"){
      usedata = holddata[holddata$Variant %in% holdouts[1],]
    } else if(new.preds == "213.variants"){
      varmat = read.xlsx("output/Variant_Selection_2023-04-26.xlsx")
      varnames = gsub("_NA", "", paste0(varmat$Related.Variant, "_", gsub(" ", "_", varmat$Mutations.In.Related.Variant)))
      
      all_antibs = unique(indata$Antibody)
      usedata = expand.grid(varnames, all_antibs, stringsAsFactors = F)
      colnames(usedata) = c("Variant", "Antibody")
      usedata$endpoint = usedata$uncert = 0
      
      seqdata = read.xlsx("Datasets/Library of 59 antibodies_ Binding and sequences_MP 12202022.xlsx", 
                          sheet = "Antibody sequences", cols = 16:19)
      colnames(seqdata)[colnames(seqdata) == "Name"] = "Antibody"
      usedata = merge(usedata, seqdata, by = "Antibody")
      
      # vars = "213_variants"
    }
    
    newxylist = mai_antib_features(usedata, feats, keep.var.feats, wts, max.antib.feats, split.cdr, cdr.len.feat)
    
    # preds = modeler(xin = xylist$x, yin = xylist$y, xout = newxylist$x, ml.params = ml.params, ...)
    
    no.out = modeler(xin = xylist$x, yin = xylist$y, ml.params = ml.params, ml = ml, filename = filename, ...)

    mod = readRDS(paste0(filename, ".rds"))

    preds = predict(mod, data = newxylist$x )$predictions
    usedata$log_Prediction = preds
    usedata$Prediction = 10^preds
    
    write.csv(usedata, row.names = F, file = paste0(filename, ".csv"))
  }

  if(variant.fold) splitvec = as.numeric(as.factor(xylist$variant)) else splitvec = NULL
  pred.df = crossval(x = xylist$x, y = xylist$y, splitvec = splitvec, ml.params = ml.params, ml = ml, ...)
  
  #Plotting
  # if(save.plot){
  #   filename = paste0("output/mai_antib_model/", plot.type, "_", ml, "_", feat, )
  # }
  
  pred.df$uncert = xylist$uncert
  pred.df$weight = xylist$weight
  pred.df$variant = xylist$variant
  
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

  #Error by Variant
  
  if(plot.type == "bar"){
    var.dfs = split(pred.df, pred.df$variant)
    rmses = sapply(var.dfs, function(x){rmse(x$yout - x$preds)})
    cors = sapply(var.dfs, function(x){cor(x$yout, x$preds)})
    q2s = sapply(var.dfs, function(x){q2(x$yout, x$preds)})
    
    # to.plot = rbind(cors, q2s)
    # cols = c("black", "slateblue")
    # barplot(to.plot, beside = T, ylim = c(min(c(cors, q2s, 0)),1.2), col = cols, xlab = "Variant")
    # legend("topleft", legend = c("Pearson Correlation", "Q^2"), fill = cols)
    
    all.q2 = q2(pred.df$yout, pred.df$preds)
    barplot(q2s, beside = T, ylim = c(-1,1), col = "slateblue", xlab = "Variant", ylab = "Q^2", cex.names = .6,
            main = paste0("Kd Delta: ", end, ", Version: ", version, add,
                          "\nOverall Q^2 = ", round(all.q2, 2) ))
    abline(h = all.q2, lty = 2)
    

  }
  
  # omi.df = pred.df[xylist$variant == "Omicron",]
  # plot(omi.df$yout, omi.df$preds, xlab = "Experimental Log(Kd)", ylab = "Predicted Log(Kd)",
  #      main = paste0("Omicron Only; ML type: ", ml, ", Features: ", feats, 
  #                    "\nRMSE = ", round(rmse(omi.df$yout - omi.df$preds), 2), 
  #                    ", cor = ", round(cor(omi.df$yout, omi.df$preds), 2), 
  #                    ", Q^2 = ", round(q2(omi.df$yout, omi.df$preds), 2) ))
  # abline(0, 1, lwd = 3)
  
  #Prediction error relative to uncertainty
  # to.plot = (pred.df$preds - pred.df$yout)/pred.df$uncert
  # levels = c(1,2,3)
  # underse = sapply(levels, function(level) round(sum(abs(to.plot) < level, na.rm = T)/length(to.plot[!is.na(to.plot)])*100, 1) )
  # plot(density(to.plot[abs(to.plot) < 5], na.rm = T), main = paste0("\nPrediction Error/Observation Error Distribution
  #             Percent Predictions within ", paste0(levels, collapse = ", "), 
  #                                                        " Standard Errors: ", paste0(underse, collapse = ", ")), xlab =  "Standard Errors")
  # abline(v = c(-3, -2, -1, 1, 2, 3))
  
  # useform = as.formula( paste0("Log_Kd ~ ", paste0(c(vfeats, afeats), collapse = " + ")) )
  # res = lm(useform, fulldata)
  # print(summary(res))
  
}

mai_antib_features = function(indata, feats, keep.var.feats = 4, wts = "none", max.antib.feats = 24, split.cdr = F,
                              cdr.len.feat = F, vars = "213_variants"){
  ################First, get starr model output for known sequences + BA4/BA5, this hopefully covers variant effects
  #run starr model first time only and save
  if(!file.exists("output/mai_antib_model/var_preds.csv")){
    # allvars = unique(indata$Variant)
    var.muts = read.csv("output/var_muts.csv")
    var.muts$Variant = gsub("_V.*", "", var.muts$Variant)
    # use.muts = var.muts[var.muts$Variant %in% allvars,]
    use.muts = var.muts
    #Will just ignore deletions to avoid completely ignoring this data
    use.muts$aa_subs = gsub("[[:alnum:]]*\\- ", "", use.muts$aa_subs)
    use.muts$aa_subs = gsub(" [[:alnum:]]*\\-$", "", use.muts$aa_subs)
    
    preds = ml_predict(aa_subs = use.muts$aa_subs, inpath = 
                         "ml_log_escape_keras_full_layer.sizes_c(256, 32, 16)_relu.alphas_c(0, 0, 0)_comb_onefold/Antibody.h5")
    preds$Variant = use.muts$Variant[match(preds$aa_subs, use.muts$aa_subs)]
    
    write.csv(preds, file = "output/mai_antib_model/var_preds.csv", row.names =  F)
  }
  
  if(!file.exists("output/mai_antib_model/213_variant_preds.csv")){
    
    #Doing this to deal with deletion-having wilds in the 213 variants without doing an alignment again
    wilds = read.csv("output/var_muts.csv")
    wilds$Variant = gsub("_V.*", "", wilds$Variant)
    #Will just ignore deletions to avoid completely ignoring this data
    wilds$aa_subs = gsub("[[:alnum:]]*\\- ", "", wilds$aa_subs)
    wilds$aa_subs = gsub(" [[:alnum:]]*\\-$", "", wilds$aa_subs)
    
    var.muts = read.xlsx("output/Variant_Selection_2023-04-26.xlsx")
    var.muts$var.name = gsub("_NA", "", paste0(var.muts$Related.Variant, "_", gsub(" ", "_", var.muts$Mutations.In.Related.Variant)))
    #Need AA sub relative to wuhan
    wuhan = unique(var.muts$Related.Variant.RBD.Sequence[var.muts$Related.Variant == "Wuhan"])
    #separate wildtypes and use aa_subs from wilds
    var.muts$wildtype = is.na(var.muts$Mutations.In.Related.Variant)
    use.wilds = wilds[wilds$Variant %in% var.muts[var.muts$wildtype, "Related.Variant"], ]
    var.muts = var.muts[!var.muts$wildtype,]
    
    var.muts$aa_subs = seq_diffs(wuhan, var.muts$`Desired.RBD.Sequence.(Length.201)`)
    
    # preds = ml_predict(aa_subs = c(use.wilds$aa_subs, var.muts$aa_subs), inpath = 
    #                      "ml_bind_keras_full_layer.sizes_c(64, 4, 16)_relu.alphas_c(0.01, 0.01, 0.01)_onefold0/Binding.h5")
    # 
    preds = ml_predict(aa_subs = c(use.wilds$aa_subs, var.muts$aa_subs), inpath = 
                         "ml_log_escape_keras_full_layer.sizes_c(256, 32, 16)_relu.alphas_c(0, 0, 0)_comb_onefold/Antibody.h5")
    preds$Variant = rep(c(use.wilds$Variant, var.muts$var.name), 10)
    
    write.csv(preds, file = "output/mai_antib_model/213_variant_preds.csv", row.names =  F)
    # outdf = data.frame(list(Variant = c(use.wilds$Variant, var.muts$var.name), delta_log_Ka_predictions = preds), 
    #                    stringsAsFactors = F)
    # write.csv(outdf, file = "output/mai_antib_model/213_variant_binding_preds.csv", row.names =  F)
  }

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
  # seqdata = read.xlsx("Datasets/Library of 59 antibodies_ Binding and sequences_MP 12202022.xlsx", 
  #                     sheet = "Antibody sequences", cols = 16:19)
  # colnames(seqdata)[colnames(seqdata) == "Name"] = "Antibody"
  # fulldata = merge(fulldata, seqdata, by = "Antibody")

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

mai_dangerous_variants_response_09112023 = function(){
  #first transform
  indata = read.xlsx("Datasets/Top dangerous variants _MP 09112023.xlsx", sheet = 1)
  locs = strsplit(gsub("[[:alpha:]]*", "", indata$`aa_subs.(in.sequence.201.aa)`), " ")
  locs = lapply(locs, function(x){as.numeric(x) + 327})
  
  aa1s = strsplit(gsub("[[:digit:]]+[[:alpha:]]", "", indata$`aa_subs.(in.sequence.201.aa)`), " ")
  aa2s = strsplit(gsub("[[:alpha:]][[:digit:]]+", "", indata$`aa_subs.(in.sequence.201.aa)`), " ")
  
  indata$`aa_subs.(in.sequence.259.aa)` = mapply(paste0, aa1s, locs, aa2s, collapse = " ")
  indata$`aa_subs.(in.sequence.259.aa)`[indata$`aa_subs.(in.sequence.259.aa)` == "NANANA"] = NA
  
  #now predict sheet 2, don't need to do again
  # mai_antib_model(ml = "rf", feats = "mixed", keep.var.feats = 10, seed = 1234, 
  #                 filename = "preds_213_v3_rf_mixed_10antibs_ace2", delta.kd = "ace2", 
  #                 version = 3, new.preds = "213.variants")
  
  ddata = read.xlsx("Datasets/Top dangerous variants _MP 09112023.xlsx", sheet = 2, cols = 2:4, colNames = F)
  colnames(ddata) = c("TS Number", "Variant_201aa", "Variant_259aa")
  pdata = read.csv("output/mai_antib_model/preds_213_v3_rf_mixed_10antibs_ace2.csv")
  pdata = pdata[pdata$Variant %in% ddata$Variant_201aa, ]
  pdata$log_Prediction = -pdata$log_Prediction
  
  castdata = as.data.frame(dcast(as.data.table(pdata), Variant ~ Antibody, value.var = "log_Prediction"))
  castdata = castdata[match(ddata$Variant_201aa, castdata$Variant), -1]
  
  cmean = colMeans(castdata)
  castdata = rbind(castdata, cmean)
  castdata = castdata[, order(cmean)]
  outdata = cbind(rbind(ddata, rep("Mean", 3)), castdata)
  
  outlist = list(Sheet1 = indata, Sheet2 = outdata)
  write.xlsx(outlist, "output/Top dangerous variants _MP 09112023_response.xlsx")
  
}

mai_dangerous_variants_response_09212023 = function(){
  
  ddata = read.xlsx("Datasets/Top dangerous variants  - V2.xlsx", sheet = 1, cols = 1:3, startRow = 2, colNames = T)
  colnames(ddata) = c("TS Number", "Variant_201aa", "Variant_259aa")
  pdata = read.csv("output/mai_antib_model/preds_213_v3_rf_mixed_10antibs_ace2.csv")
  ddata$Variant_201aa = gsub(" ", "_", ddata$Variant_201aa)
  ddata$Variant_201aa = gsub("_$", "", ddata$Variant_201aa)
  pdata = pdata[pdata$Variant %in% ddata$Variant_201aa, ]
  pdata$log_Prediction = -pdata$log_Prediction
  
  castdata = as.data.frame(dcast(as.data.table(pdata), Variant ~ Antibody, value.var = "log_Prediction"))
  castdata = castdata[match(ddata$Variant_201aa, castdata$Variant), -1]
  
  cmean = colMeans(castdata)
  castdata = rbind(castdata, cmean)
  castdata = castdata[, order(cmean)]
  outdata = cbind(rbind(ddata, rep("Mean", 3)), castdata)
  
  write.xlsx(outdata, "output/Top dangerous variants  - V2_response.xlsx")
  
}

# mai_antib_predict = function(filename = "output/mai_antib_model/rf_mixed_nosplitcdr_2antib.rds"){
#   
#   indata = read.csv("output/mai_antib_kds.csv")
#   
#   all_antibs = unique(indata$Antibody)
#   miss.omicron = setdiff(all_antibs, indata$Antibody[indata$Variant == "Omicron"])
#   
#   newdata1 = as.data.frame(list(Variant = "Omicron", Antibody = miss.omicron, endpoint = 0, uncert = 0))
#   newdata2 = as.data.frame(list(Variant = "BA4_BA5", Antibody = all_antibs, endpoint = 0, uncert = 0))
#   usedata = rbind(newdata1, newdata2)
#   
#   xylist = mai_antib_features(indata = usedata, feats = "mixed", keep.var.feats = 2, wts = "none", max.antib.feats = 24, split.cdr = F)
# 
#   mod = readRDS(filename)
#   
#   preds = predict(mod, data = xylist$x )$predictions
#   usedata$log_Prediction = preds
#   usedata$Prediction = 10^preds
#   usedata$endpoint = usedata$uncert = NULL
#   
#   write.csv(usedata, row.names = F, file = gsub(".rds", ".csv", filename))
#   
# }

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



