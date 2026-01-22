library(Matrix)
library(data.table)
library(protr)
library(parallel)

#indf needs at least three columns: "endpoint", "uncert", "aa_substitutions"
#default: use number of instances of location and transition types.
#toadd: existence of mutation types. conjoint triad category transitions
#noagg ignores the aggregation to retain original order in special cases
locandtrans = function(indf, count = T, triads = F, noagg = F){
  indf$locs = gsub("[[:alpha:]]|\\*", "", indf$aa_substitutions)
  indf$trans = gsub("[[:digit:]]", "", indf$aa_substitutions)
  
  #conjoint categories
  if(triads){
    indf$trans = gsub("A|G|V",   "X1", indf$trans)
    indf$trans = gsub("I|L|F|P", "X2", indf$trans)
    indf$trans = gsub("Y|M|T|S", "X3", indf$trans)
    indf$trans = gsub("H|N|Q|W", "X4", indf$trans)
    indf$trans = gsub("R|K",     "X5", indf$trans)
    indf$trans = gsub("D|E",     "X6", indf$trans)
    indf$trans = gsub("C",       "X7", indf$trans)
    indf$trans = gsub("\\*",     "X8", indf$trans)
  }
  
  #ensure locs and trans are sorted the same way to ensure rows are identical
  #this isn't really necessary since locs are already orderd
  #locs are guaranteed to occur only once each per row, so doesn't need to be in the if statement
  indf$locs = sapply(strsplit(indf$locs, " "), 
                        function(x){paste0(sort(as.numeric(x)), collapse = " ")})
  if(count){
    indf$trans = sapply(strsplit(indf$trans, " "), function(x){paste0(sort(x), collapse = " ")})
  } else {
    indf$trans = sapply(strsplit(indf$trans, " "), function(x){paste0(sort(unique(x)), collapse = " ")})
  }
  
  if(noagg){
    data = indf[, c("locs", "trans", "endpoint", "uncert")]
  } else {
    data = aggregate(endpoint ~ name + locs + trans + n_aa_substitutions, indf, mean)
    uncerts = aggregate(uncert ~ name + locs + trans + n_aa_substitutions, indf, mean_uncert)
    data$uncert = uncerts$uncert
  }

  
  unlocs = sort(as.numeric(unique(unlist(strsplit(data$locs, " ")))))
  untrans = sort(unique(unlist(strsplit(data$trans, " "))))
  allcols = c(unlocs, untrans)
  
  #do one hot encoding as sparse matrix
  un_index = 1:length(allcols)
  names(un_index) = allcols
  
  loc_tabs = mapply(function(aa_loc, aa_trans, count){
    if(aa_loc == "") return(NULL)
    indivs = unlist(strsplit(aa_loc, " "))
    out1 = un_index[indivs]
    
    if(aa_trans == "") return(NULL)
    indivs = unlist(strsplit(aa_trans, " "))
    out2 = un_index[indivs]
    
    #table is for counting, if not counting, should all come out to one
    if(count) return(table(c(out1, out2))) else return(c(out1, out2))
  }, aa_loc = data$locs, aa_trans = data$trans, count = count)
  
  if(count){
    loc_index = lapply(loc_tabs, function(x){as.numeric(names(x))})
    loc_nums = unlist(sapply(loc_tabs, function(x){as.numeric(x)}))
  } else {
    loc_index = loc_tabs
    loc_nums = rep(1, length(unlist(loc_index)))
  }
  
  rowinds = rep(1:nrow(data), sapply(loc_index, length))
  colinds = unlist(loc_index)
  
  onehots = sparseMatrix(i = rowinds, j = colinds, x = loc_nums)
  colnames(onehots) = allcols
  rownames(onehots) = paste(data$name, data$locs, data$trans, sep = "_")
  
  return(list(x = onehots, y = data$endpoint, uncert = data$uncert, name = data$name, 
              nsubs = data$n_aa_substitutions))
}

fullonehot = function(indf, noagg = F, aggonly = F){
  #get mean endpoints and combined uncertainties
  if(noagg) {
    data = indf[, c("aa_substitutions", "endpoint", "uncert")]
  } else {
    if(is.null(indf$source)){
      data = aggregate(endpoint ~ name + aa_substitutions + n_aa_substitutions, indf, mean)
      uncerts = aggregate(uncert ~ name + aa_substitutions + n_aa_substitutions, indf, mean_uncert)
      #missing uncerts removed
      data$uncert = uncerts$uncert
    } else {
      data = aggregate(endpoint ~ name + aa_substitutions + n_aa_substitutions + source, indf, mean)
      uncerts = aggregate(uncert ~ name + aa_substitutions + n_aa_substitutions + source, indf, mean_uncert)
      #missing uncerts removed
      data$uncert = uncerts$uncert[match(paste(data$name, data$aa_substitutions, data$n_aa_substitutions, data$source),
                                         paste(uncerts$name, uncerts$aa_substitutions, uncerts$n_aa_substitutions, uncerts$source))]
      data = data[!is.na(data$uncert), ]
    }

  }
  
  if(aggonly) {
    return(list(x = data[,c("aa_substitutions", "uncert")], y = data$endpoint, uncert = data$uncert, name = data$name,
                nsubs = data$n_aa_substitutions))
  }
  
  #get all mutations
  unsubs = unique(unlist(strsplit(data$aa_substitutions, " ")))
  
  #full one hot encoding
  #order mutations
  pos = gsub("[[:alpha:]]|\\*","",unsubs)
  subs = gsub(".*[[:digit:]]","",unsubs)
  unsubs = unsubs[order(as.numeric(pos), subs)]
  
  #onehots = matrix(data = rep(0,nrow(data)*length(unsubs)), nrow = nrow(data), ncol = length(unsubs))
  #do one hot encoding as sparse matrix
  un_index = seq_along(unsubs)
  names(un_index) = unsubs
  
  #gives column indices of subs for each row
  sub_index = lapply(data$aa_substitutions, function(aa_sub){
    if(aa_sub == "") return(NULL)
    indivs = unlist(strsplit(aa_sub, " "))
    out = un_index[indivs]
    # regpat = gsub(" ","|",aa_sub)
    # out = grep(regpat, unsubs)
    return(out)
  })
  
  rowinds = rep(1:nrow(data), sapply(sub_index, length))
  colinds = unlist(sub_index)
  
  onehots = sparseMatrix(i = rowinds, j = colinds, x = rep(1, length(rowinds)))
  colnames(onehots) = unsubs
  rownames(onehots) = paste0(data$name, "_", data$aa_substitutions)
  return(list(x = onehots, y = data$endpoint, uncert = data$uncert, name = data$name,
              nsubs = data$n_aa_substitutions, source = data$source))
}

protr_feats = function(indf, feats, aa_num = 201, ...){
  inlist = list(...)
  protr.params = inlist[["protr.params"]]
  
  #get mean endpoints and combined uncertainties
  #most methods give unique features for unique sequences, except triad based ones
  data = aggregate(endpoint ~ aa_substitutions, indf, mean)
  uncerts = aggregate(uncert ~ aa_substitutions, indf, mean_uncert)
  data$uncert = uncerts$uncert
  
  #need to reconstruct actual sequences, starting with wildtype
  unsubs = unique(unlist(strsplit(data$aa_substitutions, " ")))
  unstarts = unique(gsub("^(.*[[:digit:]]).*", "\\1", unsubs))
  pos = as.numeric(gsub("[[:alpha:]]", "", unstarts))
  starts = unstarts[order(pos)]
  wuhan_chars = gsub("[[:digit:]]", "", starts)
  if(length(wuhan_chars) != aa_num) stop("Error in protr_feats: 
                                 original sequence could not be obtained (aa number not equal to aa_num)")

  #now can get full sequence of each row of data
  data$fullseqs = sapply(data$aa_substitutions, seq_getter, wildtype = wuhan_chars)
  
  #need to parallelize it
  #moreaubroto nlag 5 takes about 7-8 minutes
  cl = makeCluster(8)
  clusterEvalQ(cl, {
    library(protr)
  })
  for(i in 1:length(feats)){
    feat = feats[i]
    if(feat == "ProtFP" && is.null(protr.params)) temp.params = list(index = 1:544, pc = 5, lag = 7) else temp.params = NULL
    call.list = c(list(cl = cl, X = data$fullseqs, FUN = paste0("extract", feat)), protr.params, temp.params)
  
    # x = parSapply(cl, data$fullseqs, extractMoreauBroto, nlag = 5)
    x = do.call(parSapply, call.list)

    if(i == 1) xout = t(x) else xout = cbind(xout, t(x))
  }
  stopCluster(cl)

  return(list(x = xout, y = data$endpoint, uncert = data$uncert))
}

seq_getter = function(aa_subs, wildtype){
  
  aa_sub_vec = unlist(strsplit(aa_subs, " "))
  aa_pos_vec = as.numeric(gsub("[[:alpha:]]", "", aa_sub_vec))
  aa_mut_vec = gsub(".*[[:digit:]]", "", aa_sub_vec)
  
  out = wildtype
  out[aa_pos_vec] = aa_mut_vec
  
  return(paste(out, collapse = ""))
}

#combines uncertainties when taking the mean
#using sqrt(sum(delta^2))/n for vector of uncertainties
mean_uncert = function(uncerts){
  return(sqrt(sum(uncerts^2, na.rm = T))/length(uncerts))
  
}