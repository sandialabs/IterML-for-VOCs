library(Matrix)

#use epistasis modeling results to test for interactions on the residue
binding_interactions = function(){
  
  usefolder = "output/binding_epistasis_epistasis_none_onefold_deltawuhan_fixscale/"
  indf = read.csv(paste0(usefolder, "combined_binding.csv"))
  residue = indf$preds - indf$actual
  dists = loc_dists(indf$aa_substitutions)
  all_dists = sort(unique(unlist(dists)))
  all_dists = all_dists[all_dists != 0]
  
  #construct data frame with counts of distances in each row
  # distmat = matrix(data = 0, nrow = nrow(indf), ncol = length(all_dists), 
  #                  dimnames = list(indf$aa_substitutions, all_dists))
  # 
  # for(i in 1:length(dists)){
  #   dist = dists[[i]]
  #   if(length(dist) == 1 && dist == 0) next
  #   tab = table(dist)
  #   distmat[i, match(names(tab), colnames(distmat))] = tab
  # }
  # 
  # 
  # fitdf = as.data.frame(cbind(residue, distmat))
  # res = lm(residue ~ . - 1 , fitdf)
  # 
  # coefs = summary(res)$coefficients
  # # write.csv(coefs, paste0(usefolder, "distance_fit.csv"))
  # plot(abs(coefs[,3]), xlab = "Distance Between Interacting Mutations", ylab = "Fit Estimate Standard Errors from 0")
  # abline(h = 2)
  # plot(-log10(coefs[,4]))
  # abline(h = -log10(.05))
  
  #Test interactions between regions
  # regs = c(rep("reg1", 50), rep("reg2", 50), rep("reg3", 50), rep("reg4", 51)) #original version
  regs = c(rep("G1", 70), rep("G2", 25), rep("G3", 10), rep("G4", 20), rep("G5", 40), rep("G6", 15), rep("G7", 21)) #mai
  names(regs) = 1:201
  reg.ints = loc_dists(indf$aa_substitutions, regs = regs)
  
  all_regints = NULL
  allregs = unique(regs)
  n = length(allregs)
  for(i in 1:n){
    all_regints = c(all_regints, paste(allregs[i], allregs[i:n], sep = "_"))
  }
  
  #construct data frame with counts of distances in each row
  intmat = matrix(data = 0, nrow = nrow(indf), ncol = length(all_regints), 
                   dimnames = list(indf$aa_substitutions, all_regints))
  
  for(i in 1:length(reg.ints)){
    reg.int = reg.ints[[i]]
    if(length(reg.int) == 1 && reg.int == 0) next
    tab = table(reg.int)
    intmat[i, match(names(tab), colnames(intmat))] = tab
  }
  
  res.int = lm(residue ~ . - 1 , as.data.frame(cbind(residue, intmat)))
  
  coefs = as.data.frame(summary(res.int)$coefficients)
  coefs$Fit_Estimate_Standard_Errors_From_0 = abs(coefs[,"t value"])
  write.csv(coefs, paste0(usefolder, "reg_fit_final.csv"))
  
  
  pdf(file = paste0(usefolder, "reg_int_plot_final.pdf"), width = 14, height = 7)
  par(mar = c(5.1, 4.1, 2.1, 2.1))
  
  cols = rep(rainbow(7), 7:1)
  bplot = barplot(abs(coefs[,3]), ylab = "Fit Estimate Standard Errors from 0", xlab = "Regions Interacting", cex.names = .8,
          col = cols, names.arg = NA)
  which.sig = abs(coefs[,3]) >= 2
  axis(1, at = bplot[which.sig], labels = rownames(coefs)[which.sig], font = 2, cex.axis = .7, las = 2)
  axis(1, at = bplot[!which.sig], labels = rownames(coefs)[!which.sig], cex.axis = .7, las = 2)
  rnames = unique(regs)
  rsites= sapply(rnames, function(rname){return(paste(range(as.numeric(names(regs)[regs == rname]) + 330), collapse = "-"))})
  
  legend("topright", legend = paste0(rnames, ": Sites ", rsites), fill = rainbow(7))
  abline(h = 2, col = "black", lwd = 2)
  
  dev.off()
  
  # #test interactions within best distances
  # sigdists = rownames(coefs)[coefs[, 4] < .001]
  # sigdists = as.numeric(gsub("`", "", sigdists))
  # 
  # ints = loc_dists(indf$aa_substitutions, keep_dists = sigdists)
  # all_ints = unique(unlist(ints))
  # all_ints = all_ints[all_ints != "0"]
  # 
  # all_ints = all_ints[order(as.numeric(gsub("_.*", "", all_ints)), as.numeric(gsub(".*_", "", all_ints)))]
  # 
  # #construct data frame with counts of distances in each row
  # intmat = matrix(data = 0, nrow = nrow(indf), ncol = length(all_ints), 
  #                   dimnames = list(indf$aa_substitutions, all_ints))
  # 
  # for(i in 1:length(ints)){
  #   int = ints[[i]]
  #   if(length(int) == 1 && int == 0) next
  #   intmat[i, match(unique(int), colnames(intmat))] = 1
  # }
  # 
  # residue = indf$preds - indf$actual
  # fitdf = as.data.frame(cbind(residue, intmat))
  # 
  # res.int = lm(residue ~ . - 1 , fitdf)
  # coefs = summary(res.int)$coefficients
  # best.ints = rownames(coefs[coefs[,4] < .001,])
  # best.ints = gsub("`", "", best.ints)
  # 
  # intmat2 = intmat[, best.ints]
  # residue = indf$preds - indf$actual
  # fitdf = as.data.frame(cbind(residue, intmat2))
  # 
  # res.int2 = lm(residue ~ . - 1 , fitdf)
  
  #now without counting instances of distances inside rows, about the same as above
  #construct data frame with counts of distances in each row
  # distmat2 = matrix(data = 0, nrow = nrow(indf), ncol = length(all_dists), 
  #                   dimnames = list(indf$aa_substitutions, all_dists))
  # 
  # for(i in 1:length(dists)){
  #   dist = dists[[i]]
  #   if(length(dist) == 1 && dist == 0) next
  #   distmat2[i, match(unique(dist), colnames(distmat))] = 1
  # }
  # 
  # residue = indf$preds - indf$actual
  # fitdf2 = as.data.frame(cbind(residue, distmat2))
  # 
  # res2 = lm(residue ~ ., fitdf2)

}

#turn aa_substitutions into sparse matrix of all present location interactions
loc_dists = function(aa_subs, keep_dists = NULL, regs = NULL){
  locs = gsub("[[:alpha:]]|\\*", "", aa_subs)
  out = sapply(locs, single_dist, keep_dists = keep_dists, regs = regs)
  return(out)
}

#does interactions within one data point, ordered and once per combo
#keep_dists will return only interactions with given distances
#regs is a list of regions defined by vectors of sites, and will return interactions between regions
single_dist = function(inlocs, keep_dists = NULL, regs = NULL){
  invec = unlist(strsplit(inlocs, " "))
  n = length(invec)
  if(n < 2) return(0)
  out = NULL
  interacts = NULL
  if(!is.null(regs)){
    inregs = regs[invec]
    for(i in 1:(n-1)){
      out = c(out, paste(inregs[i], inregs[(i+1):n], sep = "_"))
    }
    return(out)
  }
  
  
  for(i in 1:(n-1)){
    out = c(out, abs(as.numeric(invec[i]) - as.numeric(invec[(i+1):n])))
    if(!is.null(keep_dists)) interacts = c(interacts, paste(invec[i], invec[(i+1):n], sep = "_"))
  }
  
  if(!is.null(keep_dists)) return(interacts[out %in% keep_dists])
  return(sort(out))
}

#too large in memory to be workable
# binding_interactions = function(){
#   
#   indf = read.csv("output/binding_epistasis_epistasis_none_onefold_onefold/combined_binding.csv")
#   interacts = loc_interactions(indf$aa_substitutions)
#   all_ints = unique(unlist(interacts))
#   all_ints = all_ints[all_ints != ""]
#   
#   all_ints = all_ints[order(as.numeric(gsub("_.*", "", all_ints)), as.numeric(gsub(".*_", "", all_ints)))]
#   
#   #onehots = matrix(data = rep(0,nrow(data)*length(unsubs)), nrow = nrow(data), ncol = length(unsubs))
#   #do one hot encoding as sparse matrix
#   un_index = seq_along(all_ints)
#   names(un_index) = all_ints
#   
#   #gives column indices of subs for each row
#   sub_index = lapply(interacts, function(interact){
#     if(length(interact) == 1 && interact == "") return(NULL)
#     out = un_index[interact]
#     # regpat = gsub(" ","|",aa_sub)
#     # out = grep(regpat, unsubs)
#     return(out)
#   })
#   
#   rowinds = rep(1:nrow(indf), sapply(sub_index, length))
#   colinds = unlist(sub_index)
#   
#   onehots = sparseMatrix(i = rowinds, j = colinds, x = rep(1, length(rowinds)))
#   colnames(onehots) = all_ints
#   rownames(onehots) = indf$aa_substitutions
#   
#   residue = indf$preds - indf$actual
#   onehots = cbind(residue, onehots)
#   fitdf = as.data.frame(as.matrix(onehots))
#   
#   res = lm(residue ~ ., fitdf)
#   
# }
# #turn aa_substitutions into sparse matrix of all present location interactions
# loc_interactions = function(aa_subs){
#   locs = gsub("[[:alpha:]]|\\*", "", aa_subs)
#   out = sapply(locs, single.interaction)
#   return(out)
# }
# 
# #does interactions within one data point, ordered and once per combo
# single.interaction = function(inlocs){
#   invec = unlist(strsplit(inlocs, " "))
#   n = length(invec)
#   if(n < 2) return("")
#   out = NULL
#   for(i in 1:(n-1)){
#     out = c(out, paste(invec[i], invec[(i+1):n], sep = "_"))
#   }
#   return(out)
# }