library(Matrix)
library(data.table)
library(protr)
library(parallel)
library(openxlsx)

combined_features = function(indf, var.feats, antib.feats, max.antib.feats, max.var.feats, var.locs = "all", 
                             datatype = NULL, use.antib = NULL, fix.scale = T, ...){
  
  inlist = list(...)
  var.cores =  inlist[["var.cores"]]
  # antib.cores =  inlist[["antib.cores"]]
  if(is.null(var.cores)) var.cores = 8
  # if(is.null(antib.cores)) antib.cores = 1
  
  #Use only desired variant sequence locations
  if(var.locs != "all") {
    out = use_locs(indf$aa_substitutions, var.locs, datatype = datatype, use.antib = use.antib)
    # indf$aa_substitutions = out$newsubs #old way, change substitutions
    indf = indf[!out$remvec, ] #new way, remove data with muts outside epitope
    deftable = out$deftable
  } else deftable = NULL
  
  source = as.numeric(as.factor(indf$source))
  
  #first variant features
  to.scale = c()
  if(var.feats[1] == "full") x1 = fullonehot(indf)
  # else if(grepl("partial", var.feats[1])) x1 = fullonehot(indf, as.numeric(gsub("partial", "", var.feats)))
  #not yet implemented
  else if(var.feats[1] == "lt") {
    lt.count = inlist[["lt.count"]]
    lt.triads = inlist[["lt.triads"]]
    if(is.null(lt.count)) lt.count = T
    if(is.null(lt.triads)) lt.triads = F
    x1 = locandtrans(indf, count = lt.count, triads = lt.triads)
  } else if(var.feats[1] %in% c("epistasis")){
    x1 = indf[ , "aa_substitutions", drop = F]
  } else {
    #This is for PROTR
    #first get whole sequences for all variants
    varseqdf = get_var_seqs(indf$aa_substitutions, deftable)
    
    x1.list = lapply(1:ncol(varseqdf), function(i){
      #feed into get_feats each region separately
      var.seqs = varseqdf[, i]
      un.seq = unique(var.seqs)
      cl = makeCluster(var.cores)
      clusterEvalQ(cl, {
        library(protr)
        source("R/combined_features.R")
      })
      
      featset = t(parSapply(cl = cl, X = un.seq, FUN = get_feats, feats = var.feats, 
                            max.antib.feats = floor(max.var.feats/ncol(varseqdf)), 
                            max.lag = min(sapply(un.seq, nchar)) - 1, chunk.size = ceiling(length(un.seq)/var.cores/10)))
      stopCluster(cl)
      
      # featset = t(sapply(un.seq, get_feats, feats = var.feats, max.antib.feats = max.var.feats))
      if(!fix.scale) featset = apply(featset, 2, function(x){return((x-mean(x))/sd(x))}) #zscore scaling
      x1 = Matrix(featset[match(var.seqs, un.seq), ], sparse = T)
      colnames(x1) = paste0(colnames(x1), ".reg", i)
      if(fix.scale) to.scale = c(to.scale, colnames(x1))
      return(x1)
      
    })

    x1 = do.call(cbind, x1.list)
   
  }
  
  #then antibody features
  if(antib.feats != "none"){
    un.antib = unique(indf$HC_seq)
    featset = t(sapply(un.antib, get_feats, feats = antib.feats, max.antib.feats = max.antib.feats))
    if(!fix.scale) featset = apply(featset, 2, function(x){return((x-mean(x))/sd(x))}) #zscore scaling
    x2 = Matrix(featset[match(indf$HC_seq, un.antib), ], sparse = T)
    if(fix.scale) to.scale = c(to.scale, colnames(x2))
    if(var.feats != "none") x = cbind(x1, x2, source) else x = cbind(x2, source)
  } else x = cbind(x1, source)

  rownames(x) = NULL
  
  return(list(x = x, y = indf$endpoint, antibody = indf$Antibody,
              nsubs = indf$n_aa_substitutions, source = indf$source, uncert = indf$uncert, to.scale = to.scale))
}

#removes mutations or locations from aa_substitutions that we don't want to use
#also returns region defs
#var.locs = c("conform/interact_E/sep/color/paste", "singlemut_X"), i.e. conform_E or interact_sep
#where X is the percentage of highest + lowest single mutant effect sites to keep
#var.locs = "control", "frequent" remove 117 (avg of all E and color methods) random or least frequent locs
#conform/interact use different epitope definitions (columns in input)
#color splits epitopes by color (default as in file), sep renumbers them so each separate region is its own group,
#paste renumbers them all as group 1, E uses the direct E based definition, which has no inherent groups, so is pasted
use_locs = function(aa.subs, var.locs, datatype = NULL, use.antib = NULL){
  
  #remove 117 random locations
  control.num = 117
  aa.list = strsplit(aa.subs, " ")
  if(grepl("reverse", var.locs)) {
    reverse = T
    var.locs = gsub("reverse", "", var.locs)
  } else reverse = F
  
  if(var.locs == "control"){
    all.locs = unique(unlist(lapply(aa.list, function(aavec){ as.numeric(gsub("[[:alpha:]]|\\*", "", aavec)) })))
    remlocs = sort(sample(all.locs, control.num, replace = F))
    deftable = as.data.frame(list(loc = sort(setdiff(all.locs, remlocs)), feat = 1))
  } else if(var.locs == "frequent"){
    loc.tab = table(unlist(lapply(aa.list, function(aavec){ as.numeric(gsub("[[:alpha:]]|\\*", "", aavec)) })))
    all.locs = sort(as.numeric(names(loc.tab)))
    remlocs = sort(as.numeric(names(sort(loc.tab))[1:control.num]))
    deftable = as.data.frame(list(loc = sort(setdiff(all.locs, remlocs)), feat = 1))
  } else if(grepl("singlemut", var.locs)){
    if(is.null(use.antib)) add.label = datatype else add.label = paste0(datatype, "_", use.antib)
    sitedf = read.csv(paste0("output/", add.label,"_combined_site_regions.csv"), row.names = NULL, check.names = F)
    colnames(sitedf)[c(2,3)] = c("Sum Binding >0", "Sum Binding <0")
    #number of highest + lowest to keep
    keep.num = floor(as.numeric(gsub("singlemut_", "", var.locs))*.01*nrow(sitedf))
    high.sites = sitedf$site[order(sitedf$`Sum Binding >0`, decreasing = T)]
    low.sites = sitedf$site[order(sitedf$`Sum Binding <0`)]
    #alternate taking a high and low, discarding dupes until keep.num is filled
    keep.sites = rep(0, keep.num)
    for(i in 1:keep.num){
      if(i %% 2 == 1) keep.sites[i] = high.sites[1] else keep.sites[i] = low.sites[1]
      high.sites = high.sites[!high.sites %in% keep.sites[i]]
      low.sites = low.sites[!low.sites %in% keep.sites[i]]
    }
    remlocs = setdiff(sitedf$site, keep.sites)
    deftable = as.data.frame(list(loc = sort(keep.sites), feat = 1))
  } else {
    usecols = list(conform_color = c(1, 2), conform_sep = c(1, 2), conform_paste = c(1, 2), conform_E = c(1, 6),
                   interact_color = c(8, 9), interact_sep = c(8, 9), interact_paste = c(8, 9), interact_E = c(8, 12))
    
    edata = read.xlsx("Datasets/Epitope prediction for SARS2 RBD _MP 11152023_TYS.xlsx", rows = 15:215, 
                      cols = usecols[[var.locs]], rowNames = F, colNames = F)
    colnames(edata) = c("loc", "feat")
    
    if(grepl("color", var.locs)){
      remlocs = edata$loc[edata$feat == 0]
      deftable = edata[edata$feat != 0, ]
    } else if(grepl("sep", var.locs)){
      remlocs = edata$loc[edata$feat == 0]
      grp = 1
      for(i in 1:nrow(edata)){
        if(edata$feat[i] != 0) edata$feat[i] = grp
        if(i == nrow(edata)) break
        if(edata$feat[i] != 0 && edata$feat[i + 1] == 0) grp = grp + 1
      }
      deftable = edata[edata$feat != 0, ]
    } else if(grepl("E", var.locs)){
      remlocs = edata$loc[!edata$feat %in% "E"]
      deftable = edata[edata$feat == "E", ]
      deftable$feat = 1
    } else if(grepl("paste", var.locs)){
      remlocs = edata$loc[edata$feat == 0]
      deftable = edata[edata$feat != 0, ]
      deftable$feat = 1
    } 
  } 
  
  #remove unwanted locations
  # newsubs = unlist(lapply(aa.list, function(aavec){
  #   aaloc = as.numeric(gsub("[[:alpha:]]|\\*", "", aavec))
  #   return(paste0(aavec[!aaloc %in% remlocs], collapse = " "))
  # }))
  
  #remove data with mutations outside epitope region
  remvec = unlist(lapply(aa.list, function(aavec){
    aaloc = as.numeric(gsub("[[:alpha:]]|\\*", "", aavec))
    return(any(aaloc %in% remlocs))
  }))
  if(reverse) remvec = !remvec

  return(list(remvec = remvec, deftable = deftable))
}

#indf needs at least three columns: "endpoint", "uncert", "aa_substitutions"
#default: use number of instances of location and transition types.
#toadd: existence of mutation types. conjoint triad category transitions
#noagg ignores the aggregation to retain original order in special cases
locandtrans = function(indf, count = T, triads = F){
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

  data = indf
  # if(noagg){
  #   data = indf[, c("locs", "trans", "endpoint", "uncert")]
  # } else {
  #   data = aggregate(endpoint ~ name + locs + trans + n_aa_substitutions, indf, mean)
  #   uncerts = aggregate(uncert ~ name + locs + trans + n_aa_substitutions, indf, mean_uncert)
  #   data$uncert = uncerts$uncert
  # }

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

  return(onehots)
}

fullonehot = function(data){

  #get all mutations
  unsubs = unique(unlist(strsplit(data$aa_substitutions, " ")))
  # if(numsubs > 0){
  #   unsubs = names(sort(table(unlist(strsplit(data$aa_substitutions, " "))), decreasing = T))[1:numsubs]
  # }
  
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
    out = out[!is.na(out)]
    if(length(out) == 0) return(NULL)
    # regpat = gsub(" ","|",aa_sub)
    # out = grep(regpat, unsubs)
    return(out)
  })
  
  rowinds = rep(1:nrow(data), sapply(sub_index, length))
  colinds = unlist(sub_index)
  
  onehots = sparseMatrix(i = rowinds, j = colinds, x = rep(1, length(rowinds)), dims = c(nrow(data), length(unsubs)))
  colnames(onehots) = unsubs
  rownames(onehots) = paste0(data$source, "_", data$aa_substitutions)
  return(onehots)
  
}

# protr_feats = function(indf, feats, aa_num = 201, ...){
#   inlist = list(...)
#   protr.params = inlist[["protr.params"]]
# 
#   browser()
#   #get mean endpoints and combined uncertainties
#   #most methods give unique features for unique sequences, except triad based ones
#   data = aggregate(endpoint ~ aa_substitutions, indf, mean)
#   uncerts = aggregate(uncert ~ aa_substitutions, indf, mean_uncert)
#   data$uncert = uncerts$uncert
# 
#   #need to reconstruct actual sequences, starting with wildtype
#   unsubs = unique(unlist(strsplit(data$aa_substitutions, " ")))
#   unstarts = unique(gsub("^(.*[[:digit:]]).*", "\\1", unsubs))
#   pos = as.numeric(gsub("[[:alpha:]]", "", unstarts))
#   starts = unstarts[order(pos)]
#   wuhan_chars = gsub("[[:digit:]]", "", starts)
#   if(length(wuhan_chars) != aa_num) stop("Error in protr_feats:
#                                  original sequence could not be obtained (aa number not equal to aa_num)")
# 
#   #now can get full sequence of each row of data
#   data$fullseqs = sapply(data$aa_substitutions, seq_getter, wildtype = wuhan_chars)
# 
#   #need to parallelize it
#   #moreaubroto nlag 5 takes about 7-8 minutes
#   cl = makeCluster(8)
#   clusterEvalQ(cl, {
#     library(protr)
#   })
#   for(i in 1:length(feats)){
#     feat = feats[i]
#     if(feat == "ProtFP" && is.null(protr.params)) temp.params = list(index = 1:544, pc = 5, lag = 7) else temp.params = NULL
#     call.list = c(list(cl = cl, X = data$fullseqs, FUN = paste0("extract", feat)), protr.params, temp.params)
# 
#     # x = parSapply(cl, data$fullseqs, extractMoreauBroto, nlag = 5)
#     x = do.call(parSapply, call.list)
# 
#     if(i == 1) xout = t(x) else xout = cbind(xout, t(x))
#   }
#   stopCluster(cl)
# 
#   return(list(x = xout, y = data$endpoint, uncert = data$uncert))
# }

#optimize to get antib feats only once per unique antibody
antib_feats = function(x, feats, max.antib.feats = NULL){
  
  un.antib = unique(x)
  featset = t(sapply(un.antib, get_feats, feats = feats, max.antib.feats = max.antib.feats))
  out = Matrix(featset[match(x, un.antib), ], sparse = T)
  
}

#get feature vector for given aa string, feats, and maximum
get_feats = function(x, feats, max.antib.feats = NULL, max.lag = Inf) {
  
  if(feats %in% "mixed"){
    # f1 = extractMoreauBroto(x, nlag = 1)
    # names(f1) = paste0("MB.", names(f1))
    # f2 = extractCTDT(x)
    # names(f2) = paste0("CTDT.", names(f2))
    f3 = extractSOCN(x, nlag = min(c(floor(max.antib.feats/4), max.lag)))
    names(f3) = paste0("SOCN.", names(f3))
    args = c(list(propmat = "AATopo", index = c(37:41, 43:47)), get_pars(max.antib.feats/2, max.x = max.lag))
    f4 =  do.call("extractDescScales", c(list(x = x), args))
    # f3 = extractSOCN(x, nlag = 4)
    # f4 = extractDescScales(x, propmat = "AATopo", index = c(37:41, 43:47), lag = 2, pc = 2)
    names(f4) = paste0("DScales.", names(f4))
    # f5 = extractBLOSUM(x, lag = 2, k = 2)
    # names(f5) = paste0("BLOSUM.", names(f5))
    # f6 = extractProtFP(x, index = c(160:165, 258:296), pc = 2, lag = 2)
    # names(f6) = paste0("ProtFP.", names(f6))
    out = c(f3, f4)
  } else {
    if(feats %in% c("MoreauBroto", "Moran", "Geary")) {
      args = list(nlag = min(c(floor(max.antib.feats/8), max.lag)))
    } else if(feats %in% c("CTDC", "CTDT")){
      args = NULL
    } else if(feats == "SOCN"){
      args = list(nlag = min(c(floor(max.antib.feats/2), max.lag)) )
    } else if(feats == "DescScales"){
      args = c(list(propmat = "AATopo", index = c(37:41, 43:47)), get_pars(max.antib.feats, max.x = max.lag))
    } else if(feats == "BLOSUM"){
      args = get_pars(max.antib.feats, usenames = c("lag", "k"), max.x = max.lag)
    } else if(feats == "ProtFP"){
      args = c(list(index = c(160:165, 258:296)), get_pars(max.antib.feats, max.x = max.lag))
    }

    out = do.call(paste0("extract", feats), c(list(x = x), args))


  }
  
  return(out)
}

#get highest balanced parameters with given ceiling 
get_pars = function(ceil, usenames = c("lag", "pc"), max.x = Inf){
  usex = 1
  usey = 1
  #assumes that y grows faster
  while(par_fun(usex, usey) < ceil){
    tryx = usex + 1
    if(tryx <= max.x) {if(par_fun(tryx, usey) > ceil) break else usex = tryx}
    tryy = usey + 1
    if(par_fun(usex, tryy) > ceil){
      if(tryx > max.x) break else next
    } else usey = tryy
  }
  
  to.return = list(usex, usey)
  names(to.return) = usenames
  return(to.return)
}

par_fun = function(x, y){ return(x*y^2)}

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
  return(sqrt(sum(uncerts^2))/length(uncerts))
  
}

#start with char vector of aa.subs, convert to whole variant sequences based on wuhan sequence
#return matrix with as many columns as regions in deftable, 1 column if deftable is null
get_var_seqs = function(aa.subs, deftable = NULL){
  require(seqinr)
  
  un.subs = unique(aa.subs)
  
  var.seqs = read.fasta("Datasets/var_to_align.fasta", seqtype = "AA", forceDNAtolower = F)
  wuhan = var.seqs$Wuhan
  
  mut.seqs = sapply(un.subs, function(aa.sub, wuhan, deftable){
    locs = as.numeric(gsub("[[:alpha:]]|\\*", "", unlist(strsplit(aa.sub, " "))))
    muts = gsub(".*[[:digit:]]", "", unlist(strsplit(aa.sub, " ")))
    muts[muts == "*"] = ""
    out = wuhan
    out[locs] = muts
    if(!is.null(deftable)){
      regs = unique(deftable$feat)
      vout = sapply(regs, function(reg){ paste0(out[deftable$loc[deftable$feat == reg]], collapse = "")})
    } else {
      vout = paste0(out, collapse = "")
    }
    return(vout)
    
  }, wuhan = wuhan, deftable = deftable)
  if(is.matrix(mut.seqs)){
    finalout = t(mut.seqs)[match(aa.subs, un.subs),]
  } else {
    finalout = matrix(mut.seqs[match(aa.subs, un.subs)], ncol = 1)
  }
  
  return(finalout)
}



