#input mai 59 antibody data, clean it up, do the fits and get binding constants 
library(openxlsx)
library(data.table)

source("R/conc_resp.R")

#data_fun is a bespoke function for processing a given data set, 
#data_date in c("12_20_2022", "02_09_2023", "08_18_2023", "09_18_2023")
mai_antib_bind = function(data_date, to.plot = T, force_sd = F){
  
  finaldata = do.call(paste0("mai_data_", data_date), args = list())
  #deal with zeroes the easy way: ignore them
  finaldata = finaldata[finaldata$Conc != 0,]
  
  #Now for fitting
  var_ants = unique(finaldata$var_ant)
  
  fit.names = c("Kd", "Constant")
  fit.tries = c("fitKd", "fitconst")
  
  fit.list = lapply(var_ants, function(var_ant){
    usedata = finaldata[finaldata$var_ant == var_ant, ]
    x = usedata$Conc
    y = usedata$Binding
    
    fits = lapply(fit.tries, do.call, args = list(conc = x, resp = y, force_sd = force_sd))
    for(i in 1:length(fit.tries)){
      fits[[i]]$type = fit.names[i]
    }
    
    aics = sapply(fits, function(x){x$aic})
    out = fits[[which.min(aics)]]
    out$max_conc = max(x)
    
    return(out)
    
    # start.list = list(a = 3, b = 0, Kd = .1)
    # lower.list = list(a = 0, b = 0, Kd = .01*min(x))
    # upper.list = list(a = 4, b = 1, Kd = 100*max(x))
    # 
    # return(try(nls(y ~ a*(x/(x+Kd))+b, start = start.list, lower = lower.list,
    #               upper = upper.list, algorithm = "port")))
    # if(class(out) != "try-error"){
    #   # coefs = summary(out)$coefficients
    #   # return(unlist(list(Kd = coefs["Kd", "Estimate"], Kd_se = coefs["Kd", "Std. Error"],
    #   #                    a = coefs["a", "Estimate"], b = coefs["b", "Estimate"])))
    # } else {
    #   # return(unlist(list(Kd = NA_real_, Kd_se = NA_real_, a = NA_real_, b = NA_real_)))
    # }
    
    
  })
  
  max_concs = sapply(fit.list, function(x){x$max_conc})
  Kds = sapply(fit.list, function(x){if(!is.null(x$ga)) return(x$ga) else return(NA_real_)})
  Kd_ses = sapply(fit.list, function(x){if(!is.null(x$ga_sd)) return(x$ga_sd) else return(NA_real_)})
  tps = sapply(fit.list, function(x){if(!is.null(x$tp)) return(x$tp) else return(NA_real_)})
  as = sapply(fit.list, function(x){if(!is.null(x$a)) return(x$a) else return(NA_real_)})
  bs = sapply(fit.list, function(x){if(!is.null(x$b)) return(x$b) else return(NA_real_)})
  errs = sapply(fit.list, function(x){exp(x$er)})
  types = sapply(fit.list, function(x){x$type})
  
  # Kds = sapply(fit.list, function(x){tryCatch({summary(x)$coefficients["Kd", "Estimate"]}, error = function(e){NA_real_})})
  # Kd_ses = sapply(fit.list, function(x){tryCatch({summary(x)$coefficients["Kd", "Std. Error"]}, error = function(e){NA_real_})})
  # as = sapply(fit.list, function(x){tryCatch({summary(x)$coefficients["a", "Estimate"]}, error = function(e){NA_real_})})
  # bs = sapply(fit.list, function(x){tryCatch({summary(x)$coefficients["b", "Estimate"]}, error = function(e){NA_real_})})
  
  if(force_sd) addtag = "_forcesd" else addtag = ""
  #plotting
  
  if(to.plot){
    pdf(file = paste0("output/mai_antib_cr_plots_", data_date, addtag, ".pdf"))
    par(mfrow = c(3,2))
    fitx = 10^(seq(from = log10(min(finaldata$Conc[finaldata$Conc != 0])), 
                   to = log10(max(finaldata$Conc)), length.out = 1000))
    for(i in 1:length(var_ants)){
  
      usedata = finaldata[finaldata$var_ant %in% var_ants[i],]
      
      # fity = try(predict(fit.list[[i]], list(x = fitx)))
      fl = fit.list[[i]]
      fity = try(mypredict(fl, fitx))
      # max.err =  try(Kd_predict(list(tp = fl$tp + 2*fl$tp_sd, ga = max(c(fl$ga - 2*fl$ga_sd, 0)), b = fl$b + 2*fl$b_sd), fitx) + 2*errs[i])
      # min.err =  try(Kd_predict(list(tp = max(c(fl$tp - 2*fl$tp_sd, 0)), ga = fl$ga + 2*fl$ga_sd, b = fl$b - 2*fl$b_sd), fitx) - 2*errs[i])

      plot(usedata$Conc, usedata$Binding, bg = "black", pch = 21, lwd = 2,
           xlab = "Conc (ug/mL)", ylab = paste0(var_ants[i], " Binding"), log = "x", 
           main = paste0(var_ants[i], ", Kd = ", signif(10^Kds[i], 3), "\nKd_SE = ", signif(10^Kd_ses[i], 3),
                         ", Max Resp = ", signif(tps[i], 3), ", Fit Type = ", types[i]),
           cex.main = .95, cex.lab = .9)
      if(class(fity) != "try-error") {
        points(fitx, fity, type = "l", lwd = 2)
        # points(fitx, fity + 2*errs[i], type = "l", lwd = 2, lty = 2)
        # points(fitx, fity - 2*errs[i], type = "l", lwd = 2, lty = 2)
        polygon(c(fitx, rev(fitx)), c(fity - 2*errs[i], rev(fity + 2*errs[i])), col = rgb(0,0,0, .3), border = NA)
        # polygon(c(fitx, rev(fitx)), c(min.err, rev(max.err)), col = rgb(0,0,0, .3), border = NA)
      }
  
      abline(v = 10^Kds[i], lwd = 2)
      
      box(lwd=2)
    }
    dev.off()
    
    #plots of constant fits only
    pdf(file = paste0("output/mai_antib_cr_plots_constant_", data_date, addtag, ".pdf"))
    par(mfrow = c(3,2))
    fitx = 10^(seq(from = log10(min(finaldata$Conc[finaldata$Conc != 0])), 
                   to = log10(max(finaldata$Conc)), length.out = 1000))
    for(i in 1:length(var_ants)){
      
      if(types[i] != "Constant") next
      usedata = finaldata[finaldata$var_ant %in% var_ants[i],]
      
      # fity = try(predict(fit.list[[i]], list(x = fitx)))
      fl = fit.list[[i]]
      fity = try(mypredict(fl, fitx))
      # max.err =  try(Kd_predict(list(tp = fl$tp + 2*fl$tp_sd, ga = max(c(fl$ga - 2*fl$ga_sd, 0)), b = fl$b + 2*fl$b_sd), fitx) + 2*errs[i])
      # min.err =  try(Kd_predict(list(tp = max(c(fl$tp - 2*fl$tp_sd, 0)), ga = fl$ga + 2*fl$ga_sd, b = fl$b - 2*fl$b_sd), fitx) - 2*errs[i])
      
      plot(usedata$Conc, usedata$Binding, bg = "black", pch = 21, lwd = 2,
           xlab = "Conc (ug/mL)", ylab = paste0(var_ants[i], " Binding"), log = "x", 
           main = paste0(var_ants[i], ", Kd = ", signif(10^Kds[i], 3), "\nKd_SE = ", signif(10^Kd_ses[i], 3),
                         ", Max Resp = ", signif(tps[i], 3), ", Fit Type = ", types[i]),
           cex.main = .95, cex.lab = .9)
      if(class(fity) != "try-error") {
        points(fitx, fity, type = "l", lwd = 2)
        # points(fitx, fity + 2*errs[i], type = "l", lwd = 2, lty = 2)
        # points(fitx, fity - 2*errs[i], type = "l", lwd = 2, lty = 2)
        polygon(c(fitx, rev(fitx)), c(fity - 2*errs[i], rev(fity + 2*errs[i])), col = rgb(0,0,0, .3), border = NA)
        # polygon(c(fitx, rev(fitx)), c(min.err, rev(max.err)), col = rgb(0,0,0, .3), border = NA)
      }
      
      abline(v = 10^Kds[i], lwd = 2)
      
      box(lwd=2)
    }
    
    
    dev.off()
  }
  
  #Make final output
  finalkd = Kds
  finalkd[types != "Kd"] = types[types != "Kd"]
  
  outdf = as.data.frame(list(Variant = gsub("_[^_]+$", "", var_ants), Antibody = gsub(".*_", "", var_ants), Log_Kd = finalkd,
                             Log_Kd_SE = Kd_ses, Max_Conc = max_concs))
  
  write.csv(outdf, file = paste0("output/mai_antib_kds_", data_date, addtag, ".csv"), row.names = F)
  
  
}

#output data frame with variant, Conc, antibody, Binding, var_ant columns
mai_data_12_20_2022 = function(){
  
  #want to automatically parse excel file format without too many assumptions
  rawdata = read.xlsx("Datasets/Library of 59 antibodies_ Binding and sequences_MP 12202022.xlsx", sheet = "Binding data ",
                      skipEmptyRows = F)
  startrows = which(rawdata$X2 %in% "[HCAb] (ug/ml)") + 1
  
  datalist = lapply(1:length(startrows), function(i){
    if(i == length(startrows)){
      read.xlsx("Datasets/Library of 59 antibodies_ Binding and sequences_MP 12202022.xlsx", sheet = "Binding data ",
                startRow = startrows[i])
    } else {
      read.xlsx("Datasets/Library of 59 antibodies_ Binding and sequences_MP 12202022.xlsx", sheet = "Binding data ",
                rows = startrows[i]:(startrows[i+1]-1))
    }
  })
  
  if(!all(sapply(datalist, function(x){all.equal(colnames(x), colnames(datalist[[1]]))}))){ stop("column names mismatch") }
  
  indata = as.data.frame(rbindlist(datalist))
  spacings = indata[, colnames(indata) %in% "[HCAb].(ug/ml)"]
  if(!all(sapply(spacings, function(x){all.equal(x, spacings[,1])}))){ stop("inconsistent spacing columns") }
  indata = indata[, c(colnames(indata)[!colnames(indata) %in% c("[HCAb].(ug/ml)")], "[HCAb].(ug/ml)")]
  
  meltdata = as.data.frame(melt(as.data.table(indata), id.vars = c("X1", "[HCAb].(ug/ml)"), variable.name = "antibody", 
                                value.name = "Binding", variable.factor = F, na.rm = T))
  colnames(meltdata)[colnames(meltdata) == "X1"] = "variant"
  colnames(meltdata)[colnames(meltdata) == "[HCAb].(ug/ml)"] = "Conc"
  
  #deal with multiple run scales
  #first find out which scales are which using a tolerance
  allscales = unique(meltdata$Conc)
  scale10 = allscales[sapply(allscales, function(x){
    any(sapply(10/3^(0:7), function(y,x){isTRUE(all.equal(y,x))}, x = x))
  })]
  scale25 = allscales[sapply(allscales, function(x){
    any(sapply(25/3^(0:7), function(y,x){isTRUE(all.equal(y,x))}, x = x))
  })]
  if(any(!allscales %in% c(scale10, scale25))) stop("Not all scales accounted for")
  
  #now find out which var_ants have both scales, keep the 10 scale if both
  meltdata$scale10 = meltdata$Conc %in% scale10
  meltdata$scale25 = meltdata$Conc %in% scale25
  
  meltdata$variant = trimws(meltdata$variant)
  meltdata$var_ant = paste0(meltdata$variant, "_", meltdata$antibody)
  vaboth = intersect(unique(meltdata$var_ant[meltdata$scale10]), unique(meltdata$var_ant[meltdata$scale25]))
  finaldata = meltdata[!(meltdata$var_ant %in% vaboth) | !meltdata$scale25,]
  
  #now need to do ace2 data
  acedata = read.xlsx("Datasets/Library of 59 antibodies_ Binding and sequences_MP 12202022.xlsx", sheet = "ACE2-hFc binding ",
                      skipEmptyRows = F)
  colnames(acedata)[colnames(acedata) == "Ace2_Conc"] = "Conc"
  colnames(acedata)[colnames(acedata) == "Variant"] = "variant"
  acedata$variant = trimws(acedata$variant)
  acedata$antibody = "huACE2-Fc"
  acedata$scale10 = F
  acedata$scale25 = T
  acedata$var_ant = paste0(acedata$variant, "_", acedata$antibody)
  acedata = acedata[, match(colnames(finaldata), colnames(acedata))]
  
  finaldata = rbind(finaldata, acedata)
  
  return(finaldata)
}

#output data frame with variant, Conc, antibody, Binding, var_ant columns
mai_data_02_09_2023 = function(){
  
  files = c("Datasets/2_08_23_omicron_ba5_elisa.xlsx", "Datasets/2_09_23_omicron_ba5_elisa.xlsx")
  sheets = c("Wuhan", "Omicron BA1", "Omicron BA4+5")
  sheet_files = expand.grid(sheets, files, stringsAsFactors = F)
  colnames(sheet_files) = c("sheet", "file")
  sheet_files$variant = rep(c("Wuhan", "Omicron", "BA4_BA5"), 2)
  
  inlist = lapply(split(sheet_files, 1:6), function(sfv){
    indata1 = read.xlsx(sfv$file, sfv$sheet, rows = 34:42)
    indata2 = read.xlsx(sfv$file, sfv$sheet, rows = 43:51)
    if(!isTRUE(all.equal(indata1[, 1], indata2[, 1]))) stop("Data in same sheets with different concentrations")
    
    indata = cbind(indata1, indata2[,-1])
    colnames(indata)[1] = "Conc"
    
    #fix column names - assume 3 columns per antibody
    n = ncol(indata)
    anames = colnames(indata)[seq(2, n-2, 3)]
    colnames(indata)[seq(3, n-1, 3)] = anames
    colnames(indata)[seq(4, n, 3)] = anames
    
    #melt
    meltdata = as.data.frame(melt(as.data.table(indata), id.vars = c("Conc"), variable.name = "antibody", 
                                  value.name = "Binding", variable.factor = F, na.rm = T))
    meltdata$variant = sfv$variant
    meltdata$var_ant = paste0(meltdata$variant, "_", meltdata$antibody)
    return(meltdata)
  })
  
  finaldata = as.data.frame(rbindlist(inlist))
  return(finaldata)
  
}

# This data version is deprecated
# mai_data_07_13_2023 = function(){
#   
#   sheets = c("AP1B6", "AP1C3", "BP2A6", "BP3D5", "KP1C9", "HuACE2")
#   rowsets = list(group1 = 1:9, group2 = 10:18)
#   
#   inlist = lapply(sheets, function(sheet){
#     insheet = lapply(rowsets, function(rowset){
#       indata = read.xlsx("Datasets/Variant selection and binding assay _MP 07132023.xlsx", sheet = sheet, rows = rowset)
#       cols = lapply(2:4, function(startcol){ return(indata[, c(1, seq(startcol, ncol(indata), by = 3))]) })
#       colnames(cols[[2]]) = colnames(cols[[3]]) = colnames(cols[[1]])
#       return(rbind(cols[[1]], cols[[2]], cols[[3]]))
#     })
#     if(!all(insheet[[1]][,1] == insheet[[2]][,1])) stop("Concentrations within sheet are not equal")
#     
#     outdf = cbind(insheet[[1]], insheet[[2]][,-1])
#     colnames(outdf)[1] = "Conc"
#     finaldf = as.data.frame(melt(as.data.table(outdf), id.vars = "Conc", variable.name = "variant", value.name = "Binding"))
#     
#     
#     antib = gsub("HuACE2", "huACE2-Fc", sheet)
#     finaldf$antibody = antib
#     return(finaldf)
#     
#   })
#   
#   fulldf = as.data.frame(rbindlist(inlist))
#   
#   varkey = read.xlsx("Datasets/Variant selection and binding assay _MP 07132023.xlsx", sheet = "Variant lib ", cols = 1:2)
#   varkey$`Name.(201.aa)` = gsub(" ", "_", varkey$`Name.(201.aa)`)
#   varkey$`Name.(201.aa)` = gsub("_$", "", varkey$`Name.(201.aa)`)
#   varkey = rbind(varkey, c("Wuhan", "Wuhan_control"))
#   fulldf$variant = varkey$`Name.(201.aa)`[match(fulldf$variant, varkey$Name.for.Twist.library)]
#   
#   fulldf$var_ant = paste0(fulldf$variant, "_", fulldf$antibody)
#   
#   return(fulldf)
# }

mai_data_08_18_2023 = function(){
  
  sheets = c("AP1B6", "AP1C3", "BP2A6", "BP3D5", "KP1C9", "HuACE2")
  
  inlist = lapply(sheets, function(sheet){
    #autodetect data rows based on second column
    col2 = read.xlsx("Datasets/Variant selection and binding assay _MP 08182023.xlsx", sheet = sheet, cols = 2, colNames = F,
                     skipEmptyRows = F)
    rowstarts = grep("TS", col2[, 1])
    keeprows = which(!is.na(col2[, 1]))
    rowsets = lapply(1:length(rowstarts), function(i){
      if(i == length(rowstarts)) return(keeprows[keeprows >= rowstarts[i]]) else
        return(keeprows[keeprows >= rowstarts[i] & keeprows < rowstarts[i+1]])
    })
    
    insheet = lapply(rowsets, function(rowset){
      indata = read.xlsx("Datasets/Variant selection and binding assay _MP 08182023.xlsx", sheet = sheet, rows = rowset, cols = 1:25)
      cols = lapply(2:4, function(startcol){ return(indata[, c(1, seq(startcol, ncol(indata), by = 3))]) })
      colnames(cols[[2]]) = colnames(cols[[3]]) = colnames(cols[[1]])
      return(rbind(cols[[1]], cols[[2]], cols[[3]]))
    })
    allconcs = sapply(insheet, function(x){x[,1]})
    if(!is.matrix(allconcs) || !all(allconcs == allconcs[,1])) stop("Concentrations within sheet are not equal")
    
    keepconcs = insheet[[1]][, 1]
    
    outdf = cbind(keepconcs, do.call(cbind, lapply(insheet, function(x){x[, -1]})))
    outdf = outdf[, !duplicated(colnames(outdf))] #remove accidental doubles
    colnames(outdf)[1] = "Conc"
    finaldf = as.data.frame(melt(as.data.table(outdf), id.vars = "Conc", variable.name = "variant", value.name = "Binding",
                                 variable.factor = F))
    
    
    antib = gsub("HuACE2", "huACE2-Fc", sheet)
    finaldf$antibody = antib
    return(finaldf)
    
  })
  
  fulldf = as.data.frame(rbindlist(inlist))
  fulldf$variant[fulldf$variant == "TS188_A4"] = "TS188"
  fulldf$ts_variant = fulldf$variant
  
  
  varkey = read.xlsx("Datasets/Variant selection and binding assay _MP 08162023.xlsx", sheet = "Variant lib ", cols = 1:2)
  varkey$`Name.(201.aa)` = gsub(" ", "_", varkey$`Name.(201.aa)`)
  varkey$`Name.(201.aa)` = gsub("_$", "", varkey$`Name.(201.aa)`)
  varkey = rbind(varkey, c("Wuhan", "Wuhan_control"))
  fulldf$variant = varkey$`Name.(201.aa)`[match(fulldf$ts_variant, varkey$Name.for.Twist.library)]

  fulldf$var_ant = paste0(fulldf$variant, "_", fulldf$antibody)
  
  return(fulldf)
}

mai_data_09_18_2023 = function(){
  
  sheets = c("AP1B6", "AP1C3", "BP2A6", "BP3D5", "KP1C9", "HuACE2")
  
  inlist = lapply(sheets, function(sheet){
    #autodetect data rows based on second column
    col2 = read.xlsx("Datasets/Variant selection and binding assay _MP 09182023.xlsx", sheet = sheet, cols = 2, colNames = F,
                     skipEmptyRows = F)
    rowstarts = grep("TS", col2[, 1])
    keeprows = which(!is.na(col2[, 1]))
    rowsets = lapply(1:length(rowstarts), function(i){
      if(i == length(rowstarts)) return(keeprows[keeprows >= rowstarts[i]]) else
        return(keeprows[keeprows >= rowstarts[i] & keeprows < rowstarts[i+1]])
    })
    
    insheet = lapply(rowsets, function(rowset){
      indata = read.xlsx("Datasets/Variant selection and binding assay _MP 09182023.xlsx", sheet = sheet, rows = rowset, cols = 1:25)
      cols = lapply(2:4, function(startcol){ return(indata[, c(1, seq(startcol, ncol(indata), by = 3))]) })
      colnames(cols[[2]]) = colnames(cols[[3]]) = colnames(cols[[1]])
      return(rbind(cols[[1]], cols[[2]], cols[[3]]))
    })
    allconcs = sapply(insheet, function(x){x[,1]})
    if(!is.matrix(allconcs) || !all(allconcs == allconcs[,1])) stop("Concentrations within sheet are not equal")
    
    keepconcs = insheet[[1]][, 1]
    
    outdf = cbind(keepconcs, do.call(cbind, lapply(insheet, function(x){x[, -1]})))
    outdf = outdf[, !duplicated(colnames(outdf))] #remove accidental doubles
    colnames(outdf)[1] = "Conc"
    finaldf = as.data.frame(melt(as.data.table(outdf), id.vars = "Conc", variable.name = "variant", value.name = "Binding",
                                 variable.factor = F))
    
    
    antib = gsub("HuACE2", "huACE2-Fc", sheet)
    finaldf$antibody = antib
    return(finaldf)
    
  })
  
  fulldf = as.data.frame(rbindlist(inlist))
  fulldf$variant[fulldf$variant == "TS188_A4"] = "TS188"
  fulldf$ts_variant = fulldf$variant
  
  
  varkey = read.xlsx("Datasets/Variant selection and binding assay _MP 08162023.xlsx", sheet = "Variant lib ", cols = 1:2)
  varkey$`Name.(201.aa)` = gsub(" ", "_", varkey$`Name.(201.aa)`)
  varkey$`Name.(201.aa)` = gsub("_$", "", varkey$`Name.(201.aa)`)
  varkey = rbind(varkey, c("Wuhan", "Wuhan_control"))
  fulldf$variant = varkey$`Name.(201.aa)`[match(fulldf$ts_variant, varkey$Name.for.Twist.library)]
  
  fulldf$var_ant = paste0(fulldf$variant, "_", fulldf$antibody)
  
  return(fulldf)
}

sharepoint_format = function(data_date = "08_18_2023"){
  indata = read.csv(paste0("output/mai_antib_kds_", data_date,".csv"))
  indata$Kd = 10^(as.numeric(indata$Log_Kd))
  indata$Kd[indata$Log_Kd == "Constant"] = "Constant"
  
  varkey = read.xlsx("Datasets/Variant selection and binding assay _MP 08162023.xlsx", sheet = "Variant lib ", cols = 1:2)
  varkey$`Name.(201.aa)` = gsub(" ", "_", varkey$`Name.(201.aa)`)
  varkey$`Name.(201.aa)` = gsub("_$", "", varkey$`Name.(201.aa)`)
  varkey = rbind(varkey, c("Wuhan", "Wuhan_control"))
  indata$Variant = varkey$Name.for.Twist.library[match(indata$Variant, varkey$`Name.(201.aa)`)]
  
  outdata = as.data.frame(dcast(as.data.table(indata), Variant ~ Antibody, value.var = "Kd"))
  
  wantcols = c("AP1B6", "AP1C3", "BP2A6", "BP3D5", "KP1C9", "huACE2-Fc")
  wantrows = paste0("TS", 1:213)
  outdata = outdata[outdata$Variant %in% wantrows,]
  rownames(outdata) = outdata$Variant
  
  finaldata = outdata[wantrows, wantcols]
  write.csv(finaldata, paste0("output/mai_antib_kds_", data_date,"_sharepoint.csv"))
}


