
#makes diag plots for all .rds files in dirname and outputs to same dir
ml_plot = function(runname, datatype){
  
  #label names based on data type
  if(datatype %in% c("bind", "bindv2", "bindv3")){
    xlabs = c("Error (log10(M))", "Observed Delta(log10(Ka)) Binding")
    ylabs = c("Predicted Delta(log10(Ka)) Binding")
  } else if(datatype == "escape"){
    xlabs = c("Error", "Observed Escape Fraction")
    ylabs = c("Predicted Escape Fraction")
  } else if(datatype == "log_escape"){
    xlabs = c("Error", "Observed log10(Escape Fraction)")
    ylabs = c("Predicted log10(Escape Fraction)")
  } else if(datatype == "express"){
    xlabs = c("Error", "Observed Delta Mean Fluorescence")
    ylabs = c("Predicted Delta Mean Fluorescence")
  }
  
  #get antib names from rds files in directory
  dirname = paste0("output/", runname)
  allnames = gsub(".rds", "", list.files(path = dirname, pattern = ".rds"))
  
  #begin pdf
  pdf(paste0(dirname, "/", runname,".pdf"), width = 8.5, height = 11)
  par(mfrow = c(2,1))
  
  for(i in 1:length(allnames)){
    #get data
    name = allnames[i]
    pred.df = readRDS(paste0(paste0(dirname, "/", name) , ".rds"))
    
    #Error density plot with WMSE
    use_wmse = signif(wmse(pred.df$preds - pred.df$yout, get_w("uncert", uncert = pred.df$uncert)), 2)
    plot(density(pred.df$preds - pred.df$yout), main = paste0(name, "\nPrediction Error Distribution
       Weighted RMSE: ", use_wmse), xlab = xlabs[1])
    
    #Prediction error relative to uncertainty
    to.plot = (pred.df$preds - pred.df$yout)/pred.df$uncert
    levels = c(1,2,3)
    underse = sapply(levels, function(level) round(sum(abs(to.plot) < level)/length(to.plot)*100, 1) )
    plot(density(to.plot), main = paste0(name, "\nPrediction Error/Observation Error Distribution"), 
         xlab = "Standard Errors")
    plot(density(to.plot[abs(to.plot) < 5]), main = paste0(name, "\nPrediction Error/Observation Error Distribution
              Percent Predictions within ", paste0(levels, collapse = ", "), 
              " Standard Errors: ", paste0(underse, collapse = ", ")), xlab =  "Standard Errors")
    abline(v = c(-3, -2, -1, 1, 2, 3))
    
    #Main error plot with stats
    plot(pred.df$yout, pred.df$preds, xlab = xlabs[2], ylab = ylabs[1], 
         main = paste0(name, 
                       "\nRMSE = ", round(rmse(pred.df$yout - pred.df$preds), 2), 
                       ", cor = ", round(cor(pred.df$yout, pred.df$preds), 2), 
                       ", Q^2 = ", round(q2(pred.df$yout, pred.df$preds), 2) ))
    abline(0, 1, lwd = 3)
    
    #Error by # mutations
    pred.df$err = pred.df$yout - pred.df$preds
    pred.df$nsub.comb = pred.df$nsubs
    pred.df$nsub.comb[pred.df$nsubs < 2] = "0-1"
    pred.df$nsub.comb[pred.df$nsubs > 8] = "9+"
    cn = length(table(pred.df$nsub.comb))
    
    #RMSES by mutations
    rmses = aggregate(err ~ nsub.comb, pred.df, rmse)
    lengs = aggregate(err ~ nsub.comb, pred.df, length)
    alpha_level = 0.05
    n = lengs$err
    c_l <- sqrt((n - 1)/qchisq(alpha_level/2, n-1, lower.tail = FALSE))*rmses$err
    c_u <- sqrt((n - 1)/qchisq(alpha_level/2, n-1, lower.tail = TRUE))*rmses$err
    
    stripchart(as.list(rmses$err), group.names = rmses$nsub.comb, vertical = T, 
               ylim = range(c(c_l, c_u, rmses$err), na.rm = T), pch = 16,
               xlab = "Number of Mutations", ylab = "RMSE", main = "RMSE with 95% Confidence Interval")
    arrows(x0 = 1:cn, y0 = rmses$err, x1 = 1:cn, y1 = c_u, angle = 90, length = .1)
    arrows(x0 = 1:cn, y0 = rmses$err, x1 = 1:cn, y1 = c_l,  angle = 90, length = .1)
    
    #Q2/Cor by mutations
    q2s = by(pred.df, pred.df$nsub.comb, FUN = function(x){
      return(q2(x$yout, x$preds))
    })
    cor.tests = by(pred.df, pred.df$nsub.comb, FUN = function(x){
      if(length(x$yout) < 4) return(list(estimate = NA, conf.int = c(NA, NA)))
      return(cor.test(x$yout, x$preds))
    })
    cors = sapply(cor.tests, function(x){x$estimate})
    cis = sapply(cor.tests, function(x){x$conf.int})
    
    stripchart(as.list(cors), group.names = names(cor.tests), vertical = T, ylim = range(as.vector(cis), na.rm = T), 
               pch = 16, xlab = "Number of Mutations", ylab = "Pearson Correlation", 
               main = "Correlations with 95% Confidence Interval")
    arrows(x0 = 1:cn, y0 = cors, x1 = 1:cn, y1 = cis[2,], angle = 90, length = .1)
    arrows(x0 = 1:cn, y0 = cors, x1 = 1:cn, y1 = cis[1,],  angle = 90, length = .1)
    
    stripchart(as.list(q2s), group.names = names(q2s), vertical = T, ylim = range(q2s[is.finite(q2s)]), 
               pch = 16, xlab = "Number of Mutations", ylab = "Q^2", 
               main = "Q^2 vs. Mutation Number")
    
    vars = aggregate(yout ~ nsub.comb, pred.df, function(x){mean((x - mean(x))^2)})
    mses = (rmses$err)^2
    stripchart(as.list(vars$yout), group.names = vars$nsub.comb, vertical = T, ylim = range(c(vars$yout, mses)), 
               pch = 16, xlab = "Number of Mutations", ylab = "Mean Squared Error", col = "red",
               main = "Variance vs. MSE by Mutation Number")
    stripchart(as.list(mses), vertical = T, pch = 16, add = T, col = "black")
    legend("topright", legend = c("Prediction MSE", "Data Variance"), pch = c(16, 16), col = c("black", "red"))
    
    #obs by number mutations
    boxplot(yout ~ nsub.comb, pred.df, xlab = "Number of Mutations", ylab = xlabs[2])
    
  }
  
  dev.off()
  
}

#plots just from data
data_plot = function(datatype){
  #get input data
  inname = paste0("input/", datatype, ".rds")
  indf = readRDS(inname)
  
  #begin pdf
  pdf(paste0("output/", datatype, "_dataplot.pdf"), width = 8.5, height = 11)
  par(mfrow = c(2,1))

  to.plot = rbind(table(indf$n_aa_substitutions), table(indf$n_aa_substitutions[!duplicated(indf$aa_substitutions)]))
  barplot(to.plot, beside = T, xlab = "Number of Mutations", ylab = "Number Observations", col = c("black", "blue"))
  legend("topright", legend = c("Total", "Unique Variant     "), fill = c("black", "blue"))
  
  dev.off()
}

