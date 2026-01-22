library(data.table)
library(openxlsx)

#preprocess data for antibodies, both Starr and ours, aggregate duplicates, get heavy chains, etc.
antibody_preprocess = function(){
  
  #first Starr data
  indf = as.data.frame(fread("Datasets/SARS-CoV-2-RBD_MAP_Crowe_antibodies/results/escape_scores/scores.csv"))
  
  starrdf = indf[indf$pass_pre_count_filter,]
  starrdf$endpoint = log10(starrdf$score)
  starrdf = starrdf[!is.na(starrdf$endpoint),]
  starrdf$score_sd = sqrt(starrdf$score_var)
  starrdf$score_sd[starrdf$score_sd <= 0] = min(starrdf$score_sd[starrdf$score_sd > 0])/2 #add sd floor half of nonzero minimum
  starrdf$uncert = starrdf$score_sd/starrdf$score/log(10)
  
  #aggregate endpoint/uncert
  uncerts = aggregate(uncert ~ name + aa_substitutions + n_aa_substitutions, starrdf, mean_uncert)
  starrdf = aggregate(endpoint ~ name + aa_substitutions + n_aa_substitutions, starrdf, mean)
  starrdf = merge(starrdf, uncerts, sort = F)
  
  #add heavy chain seqs
  starrdf$Antibody = gsub("_400", "", starrdf$name)
  seqdata = read.xlsx("greaney antibodies short.xlsx")
  starrdf$HC_seq = seqdata$HC.Sequence[match(starrdf$Antibody, seqdata$Antibody)]
  starrdf$source = "starr"
  
  #now our data
  data1 = read.csv("output/mai_antib_kds_12_20_2022.csv")
  data2 = read.csv("output/mai_antib_kds_02_09_2023.csv")
  data3 = read.csv("output/mai_antib_kds_09_18_2023.csv")
  maidata = rbind(data1, data2, data3)
  
  maidata$Variant[maidata$Variant == "Wuhan_control"] = "Wuhan" #just average it into other measurements
  maidata = maidata[!maidata$Log_Kd %in% "Constant",]
  maidata$endpoint = as.numeric(maidata$Log_Kd)
  maidata$uncert = as.numeric(maidata$Log_Kd_SE)
  uncerts = aggregate(uncert ~ Variant + Antibody, maidata, mean_uncert, na.action = na.pass)
  maidata = aggregate(endpoint ~ Variant + Antibody, maidata, mean)
  maidata = merge(maidata, uncerts, sort = F)
  
  maidata = maidata[!maidata$Antibody %in% c("IgG"),] #can't use IgG, but keep ACE2 for now
  
  #get aa_subs
  maidata$new_subs = gsub("BA4_BA5_|BA4_BA5", "", maidata$Variant)
  maidata$new_subs = gsub("^[[:alpha:]]*$|^[[:alpha:]]*_", "", maidata$new_subs)
  maidata$new_subs = gsub("_", " ", maidata$new_subs)
  maidata$base_variant = gsub("_.*", "", maidata$Variant)
  maidata$base_variant[maidata$base_variant == "BA4"] = "BA4_BA5"
  
  #need to add base variant subs
  wilds = read.csv("output/var_muts.csv")
  wilds$Variant = gsub("_V.*", "", wilds$Variant)
  wilds$aa_subs = gsub("-", "*", wilds$aa_subs)
  maidata$base_subs = wilds$aa_subs[match(maidata$base_variant, wilds$Variant)]
  
  var.muts = read.xlsx("output/Variant_Selection_2023-04-26.xlsx")
  wuhan = unique(var.muts$Related.Variant.RBD.Sequence[var.muts$Related.Variant == "Wuhan"])
  wuhan = unlist(strsplit(wuhan, ""))
  
  #combine subs
  maidata$aa_substitutions = sub_adder(maidata$base_subs, maidata$new_subs, wuhan)
  maidata$n_aa_substitutions = sapply(strsplit(maidata$aa_substitutions, " "), length)
  
  #add HC sequences
  seqdata = read.xlsx("Datasets/Library of 59 antibodies_ Binding and sequences_MP 12202022.xlsx", 
                      sheet = "Antibody sequences", cols = c(16,20))
  maidata$HC_seq = seqdata$Nanobody.Full.Sequence[match(maidata$Antibody, seqdata$Name)]
  maidata$source = "mai"
  
  keepcols = c("Antibody", "aa_substitutions", "n_aa_substitutions", "HC_seq", "source", "endpoint", "uncert")
  
  outdata = rbind(maidata[, keepcols], starrdf[, keepcols])
  
  
  saveRDS(outdata, file = paste0("input/antibody_preprocessed.rds"))
  
}

binding_preprocess = function(){
  #starr
  starrdf = read.csv("Datasets/SARS-CoV-2-RBD_DMS/results/binding_Kds/binding_Kds.csv")
  starrdf$endpoint = starrdf$delta_log10Ka
  starrdf$uncert = starrdf$log10SE
  
  uncerts = aggregate(uncert ~ aa_substitutions + n_aa_substitutions, starrdf, mean_uncert)
  starrdf = aggregate(endpoint ~ aa_substitutions + n_aa_substitutions, starrdf, mean)
  starrdf = merge(starrdf, uncerts, sort = F)
  starrdf$source = "starr"
  
  #now our data
  data1 = read.csv("output/mai_antib_kds_12_20_2022.csv")
  data2 = read.csv("output/mai_antib_kds_02_09_2023.csv")
  data3 = read.csv("output/mai_antib_kds_09_18_2023.csv")
  maidata = rbind(data1, data2, data3)
  
  maidata$Variant[maidata$Variant == "Wuhan_control"] = "Wuhan" #just average it into other measurements
  maidata = maidata[!maidata$Log_Kd %in% "Constant",]
  maidata = maidata[maidata$Antibody %in% c("huACE2-Fc"),] #ACE2 only
  
  maidata$endpoint = as.numeric(maidata$Log_Kd)
  maidata$uncert = as.numeric(maidata$Log_Kd_SE)
  uncerts = aggregate(uncert ~ Variant, maidata, mean_uncert, na.action = na.pass)
  maidata = aggregate(endpoint ~ Variant, maidata, mean)
  maidata = merge(maidata, uncerts, sort = F)
  #convert to log10(Ka) to match Starr
  maidata$wuhan = maidata$endpoint[maidata$Variant == "Wuhan"]
  maidata$wuhan_uncert = maidata$uncert[maidata$Variant == "Wuhan"]
  maidata$endpoint = maidata$wuhan - maidata$endpoint 
  maidata$uncert = maidata$uncert + maidata$wuhan_uncert
  
  #get aa_subs
  maidata$new_subs = gsub("BA4_BA5_|BA4_BA5", "", maidata$Variant)
  maidata$new_subs = gsub("^[[:alpha:]]*$|^[[:alpha:]]*_", "", maidata$new_subs)
  maidata$new_subs = gsub("_", " ", maidata$new_subs)
  maidata$base_variant = gsub("_.*", "", maidata$Variant)
  maidata$base_variant[maidata$base_variant == "BA4"] = "BA4_BA5"
  
  #need to add base variant subs
  wilds = read.csv("output/var_muts.csv")
  wilds$Variant = gsub("_V.*", "", wilds$Variant)
  wilds$aa_subs = gsub("-", "*", wilds$aa_subs)
  maidata$base_subs = wilds$aa_subs[match(maidata$base_variant, wilds$Variant)]
  
  var.muts = read.xlsx("output/Variant_Selection_2023-04-26.xlsx")
  wuhan = unique(var.muts$Related.Variant.RBD.Sequence[var.muts$Related.Variant == "Wuhan"])
  wuhan = unlist(strsplit(wuhan, ""))
  
  #combine subs
  maidata$aa_substitutions = sub_adder(maidata$base_subs, maidata$new_subs, wuhan)
  maidata$n_aa_substitutions = sapply(strsplit(maidata$aa_substitutions, " "), length)
  
  maidata$source = "mai"
  
  keepcols = c("aa_substitutions", "n_aa_substitutions", "source", "endpoint", "uncert")
  
  outdata = rbind(maidata[, keepcols], starrdf[, keepcols])
  outdata$Antibody = "ACE2"
  
  saveRDS(outdata, file = paste0("input/binding_preprocessed.rds"))
  
  
}

#combines, base subs and new subs, over-riding base subs when necessary and adds original aa to new subs
#input and output are vectors
sub_adder = function(base_subs, new_subs, wuhan){
  
  comb_subs = base_subs
  
  for(i in 1:length(base_subs)){
    bsub = unlist(strsplit(base_subs[i], " "))
    nsub = unlist(strsplit(new_subs[i], " "))
    if(length(nsub) == 0) next
    
    blocs = as.numeric(gsub("[[:alpha:]]", "", bsub))
    nlocs = as.numeric(gsub("[[:alpha:]]", "", nsub))
    
    nsub = paste0(wuhan[nlocs], nsub) #add orig aa to new subs
    bsub = bsub[!blocs %in% nlocs]
    blocs = as.numeric(gsub("[[:alpha:]]", "", bsub))
    
    csub = c(bsub, nsub)[order(c(blocs, nlocs))]
    comb_subs[i] = paste0(csub, collapse = " ")
  }
  
  return(comb_subs)
  
}

#combines uncertainties when taking the mean
#using sqrt(sum(delta^2))/n for vector of uncertainties
mean_uncert = function(uncerts){
  if(all(is.na(uncerts))) return(NA_real_)
  return(sqrt(sum(uncerts^2, na.rm = T))/length(uncerts))
  
}