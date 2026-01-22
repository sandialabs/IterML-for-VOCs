library(openxlsx)
library(seqinr)
library(data.table)

#4/10/23
#selects two hundred variants to experiment on based on various criteria
#selects n.top most dangerous specific variant pairs, looks for n.sing top most dangerous single mutations
#and selects all combinations (single or paired) without dupes, for each variant
#So, (n.top + (n.sing^2 + n.sing)/2)*n.variants so far. Uses the remaining of n.tot allotment for wild variants and replicates
variant_selection = function(n.top = 25, n.sing = 10, n.over = 5, n.tot = 200, surf.type = "psa"){
  
  set.seed(1234)
  #Most dangerous variants around omicron, ba.5 based on results of starr model pred_multi_driver
  #don't want variants with worse binding/expr than originals, since binding already low
  min.bind = 0
  min.expr = 0
  vseqs = read.csv("Datasets/variant_sequences.csv")
  variants = c("Omicron", "BA4_BA5")

  danglist = lapply(variants, function(variant){
    targseq = seq_getter(vseqs, variant)
    # targseq = vseqs$sequence[vseqs$name == variant]
    # targseq = unlist(strsplit(targseq, ""))[13:213]
    
    mpreds = fread(file = paste0("output/multipreds/multipreds_", variant, "_mrad2.csv"))
    
    #Make sure variants pass binding/expression minimums
    usedf = mpreds[bind.pass & expr.pass,]
    usedf = usedf[usedf$bind.pred >= min.bind & usedf$expr.pred >= min.expr,]
    
    usedf$preds[usedf$preds > 0] = 0
    targdf = usedf[is.na(usedf$mut1) & is.na(usedf$mut2),]
    usedf$antib.pred = usedf$preds - targdf$preds[match(usedf$antibs, targdf$antibs)]
    usedf$mut1 = sub("[[:alpha:]]", "", usedf$mut1) #remove original Wuhan AA, since these are based on variants
    usedf$mut2 = sub("[[:alpha:]]", "", usedf$mut2) 
    usedf$allmuts = paste0(usedf$mut1, " ", usedf$mut2)
    
    # mean.escape = as.data.frame(usedf[, .(mean.escape = mean(antib.pred)), by = .(allmuts)])
    # mean.escape = mean.escape[order(mean.escape$mean.escape, decreasing = T),]
    
    #only mutations with better escape for all antibodies
    bestdf = usedf[usedf$antib.pred >= 0,]
    besttab = table(bestdf$allmuts)
    keepmuts = names(besttab[besttab == 10])
    bestdf = bestdf[bestdf$allmuts %in% keepmuts,]
    #add solvent accessibility
    # bestdf$surf1 = surf.df[match(bestdf$loc1, surf.df$loc), paste0(surf.type, ".", variant)]
    # bestdf$surf2 = surf.df[match(bestdf$loc2, surf.df$loc), paste0(surf.type, ".", variant)]
    
    mean.escape = as.data.frame(bestdf[, .(mean.escape = mean(antib.pred)), by = .(allmuts)])
    mean.escape = mean.escape[order(mean.escape$mean.escape, decreasing = T),]
    
    #format top n mutations for final
    topdf = mean.escape[1:n.top,]
    seqs = sapply(topdf$allmuts, seq_change_vec, targseq = targseq)
    topn.df = as.data.frame(list(seq = seqs, reason = paste0("Top ", n.top, " Most Dangerous Variants Within 2 Mutations of ", variant),
                                  variant = variant, targseq = targseq, aa_subs = topdf$allmuts, 
                                 met.type = "Mean Log Escape of All Antibodies", metric = topdf$mean.escape))
    
    # finaldf$mut1 = gsub(" .*", "", finaldf$allmuts)
    # finaldf$mut2 = gsub(".* ", "", finaldf$allmuts)
    # seqs1 = lapply(finaldf$mut1, seq_changer, targseq = targseq) #do first mutation
    # seqs2 = mapply(seq_changer, mut = finaldf$mut2, targseq = seqs1, SIMPLIFY = F) #second mutation after that
    # finaldf$seq = sapply(seqs2, paste0, collapse = "")
    # finaldf$reason = paste0("Top ", n.top, " Most Dangerous Variants Within 2 Mutations of ", variant)
    
    #to get importance of each single mutation, sum mean.escape for each appearance on list, so both quantity and quality are accounted
    # mean.escape$height = nrow(mean.escape):1
    mean.escape$mut = gsub(" .*", "", mean.escape$allmuts)
    mean.escape2 = mean.escape
    mean.escape2$mut = gsub(".* ", "", mean.escape2$allmuts)
    mut.esc = aggregate(mean.escape ~ mut, rbind(mean.escape, mean.escape2), sum)
    mut.esc = mut.esc[order(mut.esc$mean.escape, decreasing = T),]
    top.sing = mut.esc$mut[1:n.sing]
    seqs = sapply(top.sing, seq_change_vec, targseq = targseq)
    sing.df = as.data.frame(list(seq = seqs, reason = paste0(n.sing, " Most Dangerous Single Mutations Of ", variant),
                                 variant = variant, targseq = targseq, aa_subs = top.sing, 
                                 met.type = "Sum of Mean Log Escapes of All 'Dangerous' Variants This Mutation Appears in", 
                                 metric = mut.esc$mean.escape[1:n.sing]))
    
    #get solvent accessibility
    # mut.esc$loc = gsub("[[:alpha:]]|\\*", "", mut.esc$mut)
    # mut.esc$surf = surf.df[match(mut.esc$loc, surf.df$loc), paste0(surf.type, ".", variant)]
    
    # mut.ht = aggregate(height ~ mut, rbind(mean.escape, mean.escape2), sum)
    # mut.ht = mut.ht[order(mut.ht$height, decreasing = T),]
    
    #combinations and singles of top mutants

    mut.exp = expand.grid(c(top.sing), c(top.sing)) #singles and combos of top ten
    #remove pairs with same loc
    colnames(mut.exp) = c("mut1", "mut2")
    mut.exp$loc1 = as.numeric(gsub("[[:alpha:]]|\\*", "", mut.exp$mut1))
    mut.exp$loc2 = as.numeric(gsub("[[:alpha:]]|\\*", "", mut.exp$mut2))
    mut.exp = mut.exp[is.na(mut.exp$loc1) | mut.exp$loc1 < mut.exp$loc2,]
    mut.exp = mut.exp[order(mut.exp$loc1, mut.exp$loc2),]
    #add sum of metrics
    mut.exp$met1 = mut.esc$mean.escape[match(mut.exp$mut1, mut.esc$mut)]
    mut.exp$met2 = mut.esc$mean.escape[match(mut.exp$mut2, mut.esc$mut)]
    mut.exp$met.sum = mut.exp$met1 + mut.exp$met2
    mut.exp = mut.exp[order(mut.exp$met.sum, decreasing = T),]
    
    #formatting
    # finaldf2 = mut.exp
    # finaldf2$loc1 = finaldf2$loc2 = NULL
    # #get aa subs from original df
    # # finaldf2$aa_subs = mpreds$aa_subs[match(paste0(finaldf2$mut1, "_", finaldf2$mut2), paste0(mpreds$mut1, "_", mpreds$mut2))]
    # seqs1 = lapply(finaldf2$mut1, seq_changer, targseq = targseq) #do first mutation
    # seqs2 = mapply(seq_changer, mut = finaldf2$mut2, targseq = seqs1, SIMPLIFY = F) #second mutation after that
    # finaldf2$seq = sapply(seqs2, paste0, collapse = "")
    # finaldf2$reason = paste0("Combinations of ", n.sing, " Most Dangerous Single Mutations Near ", variant, " (1 or 2 Mutations)")

    
    seqs = sapply(paste0(mut.exp$mut1, " ", mut.exp$mut2), seq_change_vec, targseq = targseq)
    scomb.df = as.data.frame(list(seq = seqs, 
                                  reason = paste0("Combinations of ", n.sing, " Most Dangerous Single Mutations Of ", variant),
                                 variant = variant, targseq = targseq, aa_subs = paste0(mut.exp$mut1, " ", mut.exp$mut2), 
                                 met.type = "Sum of Single Mutation Inclusion Metrics", 
                                 metric = mut.exp$met.sum))
    
    #combine and remove dupes
    finaldf = rbind(topn.df, sing.df, scomb.df)
    # finaldf = finaldf[!duplicated(finaldf$seq),]
    # finaldf$variant = variant
    # finaldf$targseq = paste0(targseq, collapse = "")
    return(finaldf)
    
  })
  dangdf = as.data.frame(rbindlist(danglist))
  # dangdf = dangdf[!duplicated(dangdf$seq),]
  # dangdf$aa_subs = paste0(dangdf$mut1, " ", dangdf$mut2)
  # dangdf$aa_subs = gsub("NA | NA|NA", "", dangdf$aa_subs)
  # dangdf$mut1 = dangdf$mut2 = NULL
  
  #Mai's additions to the dangerous mutations
  #Addition 1: Add all dangerous single mutations in relation to Wuhan
  wsing.list = lapply(variants, function(variant){
    sing.muts = dangdf[dangdf$reason == paste0(n.sing, " Most Dangerous Single Mutations Of ", variant), c("aa_subs", "metric")]
    # metdf = aggregate(metric ~ aa_subs, dangdf[dangdf$aa_subs %in% sing.muts,], max)
    # metdf = metdf[order(metdf$metric, decreasing = T),]
    wuhan = seq_getter(vseqs, "Wuhan")
    # wuhan = vseqs$sequence[vseqs$name == "Wuhan"]
    # wuhan = unlist(strsplit(wuhan, ""))[13:213]
    newseqs = sapply(sing.muts$aa_subs, seq_change_vec, targseq = wuhan)
    # newseqs = unlist(lapply(newseqs, paste0, collapse = ""))
    wsing.df = as.data.frame(list(seq = newseqs, 
                                  reason = paste0("Most Dangerous Single Mutations in ", variant," Applied to Wuhan"),
                                  variant = "Wuhan", targseq = wuhan, aa_subs = sing.muts$aa_subs,
                                  met.type = paste0("Same as Inclusion Metric for Single Mutation in ", variant), 
                                  metric = sing.muts$metric))
    return(wsing.df)
  })
  wsing.df = as.data.frame(rbindlist(wsing.list))
  
  #Addition 2:Want At least ten overlapping double mutations between omicron and BA5, will take 5 from each top 25
  over.list = mapply(FUN = function(var1, var2){
    dmuts = dangdf[grepl(paste0("Within 2 Mutations of ", var1), dangdf$reason), ]
    dmuts = dmuts[dmuts$variant == var1, ][1:n.over, ]
    targseq = seq_getter(vseqs, var2)
    seqs = sapply(dmuts$aa_subs, seq_change_vec, targseq = targseq)
    outdf = as.data.frame(list(seq = seqs, 
                               reason = paste0(n.over, " Most Dangerous Mutation Combinations in ", var1," Applied to ", var2),
                               variant = var2, targseq = targseq, aa_subs = dmuts$aa_subs,
                               met.type = paste0("Same as Inclusion Metric for Mutation Combination in ", var1),
                               metric = dmuts$metric))
    return(outdf)
  }, var1 = variants, var2 = rev(variants), SIMPLIFY = F)
  over.df = as.data.frame(rbindlist(over.list))
  
  # dmuts = dangdf[grepl("Within 2 Mutations",dangdf$reason), ]
  # dmut.omi = dmuts[dmuts$variant == "Omicron",][1:n.over,]
  # ba5 = seq_getter(vseqs, "BA4_BA5")
  # omiseqs = sapply(dmut.omi$aa_subs, seq_change_vec, targseq = ba5)
  # omi.df = as.data.frame(list(seq = omiseqs, 
  #                             reason = paste0(n.over, " Most Dangerous Mutation Combinations in Omicron Applied to BA4_BA5",
  #                               variant = "BA4_BA5", targseq = ba5, aa_subs = omi_subs))
  # 
  # ba5_subs = dmuts$aa_subs[dmuts$variant == "BA4_BA5"][1:n.over]
  # omicron = seq_getter(vseqs, "Omicron")
  # ba5seqs = unlist(lapply(ba5_subs, seq_change_vec, targseq = omicron))
  # ba5.df = as.data.frame(list(seq = ba5seqs, reason = "Five Most Dangerous Mutation Combinations in BA4_BA5 Applied to Omicron",
  #                             variant = "Omicron", targseq = omicron, aa_subs = ba5_subs))
  # add.df = rbind(wsing.df, omi.df, ba5.df)
  
  #combinations of important mutations not covered in starr
  #get score regions
  # reg.df = read.xlsx("scores.xlsx")
  # reg.sites = reg.df$site[!is.na(reg.df$region) & reg.df$region <= reg.max]
  # 
  # #Starr binding and expression effects
  # eff.df =  read.csv("Datasets/SARS-CoV-2-RBD_DMS/results/single_mut_effects/single_mut_effects.csv")
  # #remove wildtype (equals zer oby definition) and deletion mutations
  # eff.df = eff.df[eff.df$wildtype != eff.df$mutant & eff.df$mutant != "*",]
  # eff.df = eff.df[eff.df$site_SARS2 %in% reg.sites, ]
  # 
  # top_muts = eff.df[!is.na(eff.df$expr_avg) & eff.df$expr_avg > 0 & eff.df$bind_avg > 0,]
  # #to get top combined, add percentiles in binding and expression, where smallest rank is highest expr/binding
  # top_muts$bind_perc = ecdf(top_muts$bind_avg)(top_muts$bind_avg)
  # top_muts$expr_perc = ecdf(top_muts$expr_avg)(top_muts$expr_avg)
  # top_muts$perc_add = top_muts$bind_perc + top_muts$expr_perc
  # top_muts = top_muts[order(top_muts$perc_add, decreasing = T),]
  
  #Wild Variants
  wilds = read.fasta(file = "Datasets/var_aligned.fasta", seqtype = "AA", strip.desc = T)
  wilds = lapply(wilds, function(x){
    out = x[1:201]
    out = out[out != "-"]
    paste0(out, collapse = "")
    })
  names(wilds) = gsub("_V.*", "", names(wilds))
  wild.df = as.data.frame(list(seq = unlist(wilds), reason = "Wild Variant", variant = names(wilds), 
                               targseq = unlist(wilds), aa_subs = "", met.type = "No Inclusion Metric", metric = NA))
  # outdf = rbind(dangdf, wild.df)
  # wilds = vseqs
  # wilds$name = gsub("_V.*", "", wilds$name)
  # wilds$variant = wilds$spike_pdb = NULL
  # colnames(wilds) = c("variant", "seq")
  # wilds = wilds[!wilds$variant %in% c("Lambda", "Mu", "Epsilon")]
  # wilds$seq = sapply(wilds$seq, function(x){
  #   paste0(unlist(strsplit(x, ""))[13:213], collapse = "")
  # })
  # wilds$mut1 = wilds$mut2 = NA
  
  
  
  #Replicates of Starr data
  #starr data
  sdata = readRDS("input/log_escape.rds")
  # usedata = sdata[as.logical(sdata$pass_pre_count_filter) & as.logical(sdata$pass_ACE2bind_expr_filter), ]
  usedata = sdata
  aa_subs = unique(usedata$aa_substitutions)
  aa_subs = aa_subs[!aa_subs %in% c("", "N171Y")] #remove known overlap of wuhan, alpha with wilds
  
  #shorten subs to match format of others
  short_subs = sapply(aa_subs, function(aa_sub){
    subvec = unlist(strsplit(aa_sub, " "))
    subvec = sub("[[:alpha:]]", "", subvec)
    out = paste0(subvec, collapse = " ")
    return(out)
  }) 
  
  #Starr subs of wuhan also in dangerous subs
  wuhan = seq_getter(vseqs, "Wuhan")
  dlist = lapply(variants, function(variant){
    dang.sh = dangdf[dangdf$aa_subs %in% short_subs & dangdf$variant == variant, ]
    # dsubs = short_subs[short_subs %in% dangdf$aa_subs]

    dseqs = sapply(dang.sh$aa_subs, seq_change_vec, targseq = wuhan)
    
    ddf = as.data.frame(list(seq = dseqs, reason = paste0("'Dangerous' Mutations of ", variant," Also Found in Starr Data"), 
                             variant = "Wuhan", targseq = wuhan, aa_subs = dang.sh$aa_subs,
                             met.type = paste0("Same Inclusion Metric as Mutation Combination in ", variant), 
                             metric = dang.sh$metric))
    return(ddf)
  })
  ddf = as.data.frame(rbindlist(dlist))
  
  #Starr subs of Wuhan also in omicron/BA5
  varmuts = read.csv("output/var_muts.csv")
  vlist = lapply(variants, function(variant){
    vsubs = varmuts$aa_subs[varmuts$Variant %in% variant]
    vsubs = unique(unlist(strsplit(vsubs, " ")))
    vsubs = sub("[[:alpha:]]", "", vsubs)

    vsubs = intersect(short_subs, vsubs)
    seqs = sapply(vsubs, seq_change_vec, targseq = wuhan)
    vdf = as.data.frame(list(seq = seqs, reason = paste0(variant, " Mutations from Wuhan Also Found in Starr Data"), 
                             variant = "Wuhan", targseq = wuhan, aa_subs = vsubs,
                             met.type = "No Inclusion Metric", metric = NA))
    return(vdf)
  })
  vdf = as.data.frame(rbindlist(vlist))
  
  # sdf = rbind(ddf, vdf)
  # samp = sort(sample(1:nrow(sdf),n.tot - nrow(outdf)))
  # print(samp)
  # sdf = sdf[samp, ]
  
  #put it all together
  outdf = rbind(dangdf, wsing.df, over.df, wild.df, ddf, vdf)
  # if(sum(duplicated(outdf$seq)) > 0) stop("Duplicated Sequences")
  
  #get surface percent solvent accessibility to focus locations
  s.variants = c("Wuhan", "BA1", "BA4")
  surf.list = lapply(s.variants, function(variant){
    sheet = paste0("SASA - ", variant)
    indf = read.xlsx("Datasets/Variant selection - solvent accessibility for varaint _ MP 04202023.xlsx",
                     sheet = sheet, cols = c(4, 6, 8))
    colnames(indf) = c("loc", "rsa", "psa")
    indf = indf[,c("loc", surf.type)]
    indf = indf[indf$loc %in% 1:201,]
    #deal with possible structure solving failures by setting them to max
    indf[, surf.type] = as.numeric(indf[, surf.type])
    indf[is.na(indf[, surf.type]), surf.type] = max(indf[, surf.type], na.rm = T)

    colnames(indf)[colnames(indf) == surf.type] = variant
    return(indf)

  })
  surf.df = cbind(surf.list[[1]], surf.list[[2]], surf.list[[3]])
  surf.df = surf.df[, c("loc", s.variants)]
  colnames(surf.df) = gsub("BA1", "Omicron", colnames(surf.df))
  colnames(surf.df) = gsub("BA4", "BA4_BA5", colnames(surf.df))
  rownames(surf.df) = surf.df$loc
  
  
  outdf$loc1 = gsub(" .*|[[:alpha:]]", "", outdf$aa_subs)
  outdf$loc2 = gsub("[[:alpha:]]|.* ", "", outdf$aa_subs)
  outdf$loc2[!grepl(" ", outdf$aa_subs)] = NA
  
  outdf$surf2 = outdf$surf1 = NA
  to.change = outdf$variant %in% c("Wuhan", "Omicron", "BA4_BA5")
  outdf$surf1[to.change] = mapply(function(row, col){
    surf.df[row, col]
  }, row = outdf$loc1[to.change], col = outdf$variant[to.change])
  outdf$surf2[to.change] = mapply(function(row, col){
    surf.df[row, col]
  }, row = outdf$loc2[to.change], col = outdf$variant[to.change])
  outdf$loc1 = outdf$loc2 = NULL
  
  #final output
  colnames(outdf) = c("Desired RBD Sequence (Length 201)", "Reason for Inclusion", "Related Variant", 
                      "Related Variant RBD Sequence", "Mutations In Related Variant", "Inclusion Metric Type",
                      "Inclusion Metric", "Percent Solvent Accessibility at First Mutation Location", 
                      "Percent Solvent Accessibility at Second Mutation Location")
  outdf$Duplicated = duplicated(outdf$`Desired RBD Sequence (Length 201)`) | 
                             duplicated(outdf$`Desired RBD Sequence (Length 201)`, fromLast = T)
  outlist = list(Unique = outdf[!duplicated(outdf$`Desired RBD Sequence (Length 201)`), ], 
                 Duplicates = outdf[duplicated(outdf$`Desired RBD Sequence (Length 201)`), ])
  
  write.xlsx(outlist, file = paste0("output/Variant_Selection_", Sys.Date(), ".xlsx"))
  
  
  # samp = sample(aa_subs, n.tot - nrow(outdf))
  # 
  # dang_aas = unique(unlist(strsplit(dangdf$aa_subs, " ")))
  # test = sapply(aa_subs, seq_change_vec, targseq = wilds[["Wuhan"]])
  
  

  
}

#cuts the sequences from vseqs down to size and returns as string
seq_getter = function(vseqs, name){
  out = vseqs$sequence[vseqs$name == name]
  out = unlist(strsplit(out, ""))[13:213]
  return(paste0(out, collapse = ""))
}

#wrapper for seq_changer with multiple muts as space separated string, targseq is string of aas, so is output
seq_change_vec = function(aa_sub, targseq){
  mutvec = strsplit(aa_sub, " ")[[1]]
  out = strsplit(targseq, "")[[1]]
  for(mut in mutvec){
    out = seq_changer(mut, out)
  }
  return(paste0(out, collapse = "") )
  
}

#for a given target sequence vector and mutation (as aminonumberamino), returns the sequence with mutations applied
seq_changer = function(mut, targseq){
  mloc = as.numeric(gsub("[[:alpha:]]|\\*", "", mut))
  if(is.na(mloc)) return(targseq)
  newmut = gsub(".*[[:digit:]]", "", mut)
  out = targseq
  out[mloc] = newmut
  
  return(out)
  
}

#scoring function type thing for our epistasis models
#type = c("combined", "starr", "mai")
new_starr_regions = function(high.cutoff = 0, low.cutoff = 0, datatype = "binding", type = "combined"){
  
  if(type == "combined") folder.type = "" else folder.type = paste0("_", type)
  if(grepl("antibody_", datatype)) main.datatype = "antibody" else main.datatype = datatype
  eff.df = read.csv(paste0("output/", datatype, "_epistasis_epistasis_none", folder.type, 
                           "_onefold_deltawuhan_fixscale/", type, "_", main.datatype, "_single_mut_effects.csv"))
  
  sum.names = paste0("Sum", rep(" Single Mutation Effects ", 2), c(">", "<"), c(high.cutoff, low.cutoff))
  scale.names = paste0("Scaled Sum", rep(" Single Mutation Effects ", 2), c(">", "<"), c(high.cutoff, low.cutoff))
  
  agg.df = aggregate(latent.effect ~ site, eff.df, 
                     function(x){return(c(sum(x[x > high.cutoff]), c(sum(x[x < low.cutoff]))))})
  
  new.df = as.data.frame(cbind(agg.df[, 1], agg.df[[2]]))
  colnames(new.df) = c("site", sum.names)
  
  new.df[, scale.names] = lapply( new.df[, sum.names], function(x){x/max(abs(x))})

  write.csv(new.df, paste0("output/", datatype, "_", type, "_site_regions.csv"), row.names = F)
  
  pdf(paste0("output/", datatype, "_", type, "_site_regions.pdf"), width = 8.5, height = 11)
  par(mfrow = c(2,1))
  
  #added for figure 5 consistency
  ylims = c(-1.1, 1.1)
  use.cols = c("darkorange3", "palegreen4")
  barplot(new.df[, scale.names[1]], names.arg = c(), xaxt = "n", xlab = "Residue position on SARS-CoV-2 spike protein",
          cex.names = .4, col = use.cols[1], ylim = ylims, las = 3, ylab = "Site importance score", space = 0)
  barplot(new.df[, scale.names[2]],  col = use.cols[2], add = T, space = 0)
  legend("topleft", legend = c("Scaled Sum Binding >0", "Scaled Sum Binding <0"), fill = use.cols[1:2], bg = "white", 
         bty= "n", border = "white", cex = 1.2)
  axis(1, at = seq(0, 200, 20), labels = seq(330, 530, 20))
  axis(1, at = seq(10, 190, 20), labels = F, tck = -.01)
  
  ylims = c(-2, 7.5)
  cols = rainbow(7, start = 0, end = .7, rev = T)
  barplot(new.df[, scale.names[1]], names.arg = new.df$site + 330, 
          cex.names = .4, col = cols[2], ylim = ylims, las = 3, ylab = "Site Score", space = 0)
  barplot(new.df[, scale.names[2]],  col = cols[6], add = T, space = 0)
  legend("topleft", legend = scale.names, fill = cols[c(2,6)], cex = .8, bg = "white")
  
  ylims = c(-1, 1.5)
  cols = rainbow(7, start = 0, end = .7, rev = T)
  barplot(new.df[, scale.names[1]], names.arg = new.df$site + 330, 
          cex.names = .4, col = cols[2], ylim = ylims, las = 3, ylab = "Site Score", space = 0)
  barplot(new.df[, scale.names[2]],  col = cols[6], add = T, space = 0)
  legend("topleft", legend = scale.names, fill = cols[c(2,6)], cex = .8, bg = "white")
  
  ylims = c(-70, 20)
  cols = rainbow(7, start = 0, end = .7, rev = T)
  barplot(new.df[, sum.names[1]], names.arg = new.df$site + 330, 
          cex.names = .4, col = cols[2], ylim = ylims, las = 3, ylab = "Site Score", space = 0)
  barplot(new.df[, sum.names[2]],  col = cols[6], add = T, space = 0)
  legend("topleft", legend = sum.names, fill = cols[c(2,6)], cex = .8, bg = "white")
  
  dev.off()
  
  
}

#identifies regions with maximum sites of interest, minimum conserved sites
starr_regions = function(high.cutoff = 0, low.cutoff = 0){
  eff.df =  read.csv("Datasets/SARS-CoV-2-RBD_DMS/results/single_mut_effects/single_mut_effects.csv")
  #remove wildtype (equals zer oby definition) and deletion mutations
  eff.df = eff.df[eff.df$wildtype != eff.df$mutant & eff.df$mutant != "*",]
  
  #get sites of interest based on variance
  # bind.var = aggregate(bind_avg ~ site_SARS2, eff.df, sd)
  
  #abuse aggregate to get all four summary stats
  agg.df = aggregate(cbind(bind_avg, expr_avg) ~ site_SARS2, eff.df, 
                   function(x){return(c(sum(x[x > high.cutoff]), c(sum(x[x < low.cutoff]))))})
  #rescale as fraction of absolute maximum
  
  
  
  #top ten percent of sum of binding/expr above the cutoff
  # highbind.df = aggregate(bind_avg ~ site_SARS2, eff.df, function(x){sum(x[x > high.cutoff])})
  # highbind.sites = highbind.df$site_SARS2[highbind.df$bind_avg > quantile(highbind.df$bind_avg, high.quant)]
  # 
  # highexpr.df = aggregate(expr_avg ~ site_SARS2, eff.df, function(x){sum(x[x > high.cutoff])})
  # highexpr.sites = highexpr.df$site_SARS2[highexpr.df$expr_avg > quantile(highexpr.df$expr_avg, high.quant)]
  # 
  # #top ten percent of number of bind/expr below the cutoff
  # lowbind.df = aggregate(bind_avg ~ site_SARS2, eff.df, function(x){sum(x[x < low.cutoff])})
  # lowbind.sites = lowbind.df$site_SARS2[lowbind.df$bind_avg < quantile(lowbind.df$bind_avg, low.quant)]
  # 
  # lowexpr.df = aggregate(expr_avg ~ site_SARS2, eff.df, function(x){sum(x[x < low.cutoff])})
  # lowexpr.sites = lowexpr.df$site_SARS2[lowexpr.df$expr_avg < quantile(lowexpr.df$expr_avg, low.quant)]
  
  #congeal it all into a unified score: +1 for appearance in interesting cat, -1 for low bind/expr appearances
  site.df = read.xlsx("sites.xlsx")
  in.names = colnames(site.df)
  
  score.df = as.data.frame(cbind(331:531, matrix(0, 201, ncol(site.df) + 5)))
  #high bind/expr then low bind/expr
  # newnames = paste0(rep(c("Top ", "Bottom "), each = 2), rep(c(high.quant, low.quant), each = 2)*100, 
  #                   "% of Sum of Site Mutation", rep(c(" Expression ", " Binding "), 2),
  #                   rep(c("Above ", "Below "), each = 2), rep(c(high.cutoff, low.cutoff), each = 2))
  # newnames = paste0("Scaled Sum of Site Mutation", rep(c(" Binding ", " Expression "), 2),
  #                   rep(c("Above ", "Below "), each = 2), rep(c(high.cutoff, low.cutoff), each = 2))
  newnames = paste0("Scaled Sum", rep(c(" Binding ", " Expression "), 2),
                    rep(c(">", "<"), each = 2), rep(c(high.cutoff, low.cutoff), each = 2))
  colnames(score.df) = c("site", in.names, newnames, "score")
  
  #input sites.xlsx columns as positives
  for(in.name in in.names){
    score.df[score.df$site %in% site.df[, in.name],in.name] = 1
  }
  
  #new pos/negs
  score.df[, newnames] = cbind(agg.df$bind_avg[, 1], agg.df$expr_avg[, 1], 
                               agg.df$bind_avg[, 2], agg.df$expr_avg[, 2])
  score.df[, newnames] = lapply( score.df[, newnames], function(x){x/max(abs(x))})

  score.df$score = rowSums(score.df[,c(in.names, newnames)])
  
  #can just find every possible region: about 20,000 possible
  site.range = range(score.df$site)
  i = 1
  reg.score = rep(0, nrow(score.df)*(nrow(score.df)+1)/2)
  for(start in site.range[1]:site.range[2]){
    for(end in start:site.range[2]){
      reg.score[i] = sum(score.df$score[score.df$site %in% start:end])
      names(reg.score)[i] = paste0(start, "_", end)
      i = i + 1
    }
  }
  
  best.ranges = vector(mode = "numeric")
  rem.score = reg.score
  i = 1
  while(length(rem.score) > 0){
    best.range = sort(rem.score, decreasing = T)[1]
    if(best.range < 1 ) break
    best.name = strsplit(names(best.range), "_")[[1]]
    all.range.names = strsplit(names(rem.score), "_")
    keep.name = sapply(all.range.names, function(x){x[1] > best.name[2] | x[2] < best.name[1]})
    if(sum(keep.name) == 0) break
    rem.score = rem.score[keep.name]
    best.ranges = c(best.ranges, best.range)
    names(best.ranges)[length(best.ranges)] = names(best.range) 
    i = i + 1
  }
  rangesplit = sapply(strsplit(names(best.ranges), "_"), as.numeric)
  
  #Add regions to score.df
  score.df$region = NA
  for(i in 1:ncol(rangesplit)){
    score.df$region[rangesplit[1,i] <= score.df$site & rangesplit[2,i] >= score.df$site] = i
  }
  
  write.xlsx(score.df, "scores.xlsx")
  
  #plot the scores/regions
  pdf("output/site_regions.pdf", width = 8.5, height = 11)
  par(mfrow = c(2,1))
  
  nsplit = 3
  ylims = c(-2,7.5)
  splits = sort(0:(nrow(score.df) - 1)%%nsplit + 1) 
  cols = rainbow(7, start = 0, end = .7, rev = T)
  cols[1:5] = cols[5:1]
  
  #All together plot
  barplot(t(score.df[, c(in.names, newnames[1:2])]), names.arg = score.df$site[], 
          cex.names = .4, col = cols[1:5], ylim = ylims, las = 3, ylab = "Site Score", space = 0)
  barplot(t(score.df[ ,c(newnames[3:4])]),  col = cols[6:7], add = T, 
          names.arg = rep("", length(splits)), space = 0)
  
  
  # plot best regions
  usesplit = rangesplit - 330
  rect(usesplit[1,] - 1, ylims[1], usesplit[2,], ylims[2], col = rgb(0,0,0,.1), border = NA)
  
  legend("topleft", legend = c(in.names, newnames)[c(5:1, 6:7)], fill = cols[c(5:1, 6:7)], ncol = 2, cex = .8,
         bg = "white")
  
  text(colMeans(usesplit), 5, adj = .5, cex = .5,
       labels = paste0("R", 1:length(best.ranges), "\nS", round(best.ranges, 1)))
  
  #Split plot
  for(i in 1:nsplit){
    barplot(t(score.df[splits == i, c(in.names, newnames[1:2])]), names.arg = score.df$site[splits == i], 
            cex.names = .4, col = cols[1:5], ylim = ylims, las = 3, ylab = "Site Score", space = 0)
    barplot(t(score.df[splits == i ,c(newnames[3:4])]),  col = cols[6:7], add = T, 
            names.arg = rep("", sum(splits == i)), space = 0)

    
    # plot best regions
    usesplit = rangesplit - 330 - which(splits == i)[1] + 1
    rect(usesplit[1,] - 1, ylims[1], usesplit[2,], ylims[2], col = rgb(0,0,0,.1), border = NA)
    
    legend("topleft", legend = c(in.names, newnames)[c(5:1, 6:7)], fill = cols[c(5:1, 6:7)], ncol = 2, cex = .8,
           bg = "white")

    text(colMeans(usesplit), 5, adj = .5, cex = .7,
         labels = paste0("Region ", 1:length(best.ranges), "\nScore: ", round(best.ranges, 1)))
  }
  


  # barplot(score.df$score, names.arg = score.df$site, cex.names = .5)
  # barplot(score.df$score[score.df$score != 0], names.arg = score.df$site[score.df$score != 0])
  

  
  dev.off()
  
}

#finds sites in Starr data that can't be changed without virus failing
#specifically, those that change expression/binding by min.change for every non-wild mutation at a site
#with up to keep.exceptions mutations that don't have the negative effect
#takes in sites.xlsx and outputs sites_conserved.xlsx with added columns for each bind/expr/exception number combo
# starr_site_cons = function(min.change = -1, keep.exceptions = 4){
#   eff.df =  read.csv("Datasets/SARS-CoV-2-RBD_DMS/results/single_mut_effects/single_mut_effects.csv")
#   eff.df = eff.df[eff.df$wildtype != eff.df$mutant,]
#   
#   #-1 removes deletion-related missing values that are present for all binding sites
#   bind.df = aggregate(bind_avg ~ site_SARS2, eff.df, 
#                       function(x){ sum(x[!is.na(x)] > min.change) + sum(is.na(x)) - 1 }, na.action = na.pass)
#   exp.df = aggregate(expr_avg ~ site_SARS2, eff.df, 
#                      function(x){ sum(x[!is.na(x)] > min.change) + sum(is.na(x)) }, na.action = na.pass)
#   
#   site.df = read.xlsx("sites.xlsx")
#   rowlen = max(c(table(bind.df$bind_avg[bind.df$bind_avg <= keep.exceptions]), 
#                  table(exp.df$expr_avg[exp.df$expr_avg <= keep.exceptions]),
#                  nrow(site.df)))
#   
#   bind.exc.df = as.data.frame(sapply(0:keep.exceptions, function(i){
#     bind.df$site_SARS2[bind.df$bind_avg == i][1:rowlen]
#   }))
#   colnames(bind.exc.df) = paste0(0:keep.exceptions, " Mutants Binding Change More Than ", min.change)
#   
#   exp.exc.df = as.data.frame(sapply(0:keep.exceptions, function(i){
#     exp.df$site_SARS2[exp.df$expr_avg == i][1:rowlen]
#   }))
#   colnames(exp.exc.df) = paste0(0:keep.exceptions, " Mutants Expression Change More Than ", min.change)
#   
#   outdf = cbind(site.df[1:rowlen,], bind.exc.df, exp.exc.df)
#   
#   write.xlsx(outdf, "sites_conserved.xlsx")
# }


xgb_site_imp = function(imp.df, site.imp.quant = .9){
  
  imp.df = imp.df[!is.na(as.numeric(imp.df$Feature)),]
  site.cut = quantile(imp.df$Gain, site.imp.quant)
  top.sites = imp.df$Feature[imp.df$Gain > site.cut]
  top.sites = as.numeric(top.sites) + 330
  
  site.df =  read.xlsx("sites_conserved.xlsx")
  rowmax = max(nrow(site.df), length(top.sites))
  
  out.df = cbind(site.df[1:rowmax,], top.sites[1:rowmax])
  colnames(out.df)[ncol(out.df)] = paste0((1-site.imp.quant)*100, "% Most Important Model Sites (In Order)")
  
  write.xlsx(out.df, "sites_cons_plus_xgbimp.xlsx")
}

ncbi_sites = function(top.perc = .1){
  inlist = unique(read.fasta(file = "Datasets/new_spike_sequences/aligned_rbd_reg.fasta", seqtype = "AA"))
  
  names(inlist) = sapply(inlist, function(x){attributes(x)$name})
  wuhan = inlist[["Wuhan"]]
  
  wuhan.start = 62
  sitenames = rep(NA_integer_, length(wuhan))
  ct = 1
  for(i in wuhan.start:length(wuhan)){
    if(wuhan[i] == "-"){
      next
    } else {
      sitenames[i] = ct
      ct = ct + 1
    }
  }
  rbdlocs = which(sitenames %in% 1:201)
  wuhanrbd = wuhan[rbdlocs]
  
  rbdseqs = unique( lapply(inlist, function(x){x[rbdlocs]}) )
  names(rbdseqs) = paste0("Seq", 1:length(rbdseqs))
  rbd.df = as.data.frame(rbdseqs)
  
  mut.per.site = apply(rbd.df, 1, function(x){
    typetab = table(unlist(x))
    return(sum(typetab) - max(typetab))
  })
  
  mut.cutoff = length(rbdseqs)*top.perc
  var.sites = which(mut.per.site > mut.cutoff) + 330
  
  site.df = read.xlsx("sites_cons_plus_xgbimp.xlsx")
  maxrow = max(nrow(site.df), length(var.sites))
  
  out.df = cbind(site.df[1:maxrow,], var.sites[1:maxrow])
  colnames(out.df)[ncol(out.df)] = paste0("Sites With Mutations in ", top.perc*100, 
                                          "% or More of Unique NCBI RBD Sequences")
  
  write.xlsx(out.df, "sites_cons_xgbimp_seq.xlsx")
}



