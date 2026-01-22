library(openxlsx)

#pre-process raw data and save to file
ml_preprocess = function(datatype){
  
  #get data (binding, anti escape, expression), maybe do pre-processing, and set endpoint/uncert names
  if(datatype == "bind"){
    indf = read.csv("Datasets/SARS-CoV-2-RBD_DMS/results/binding_Kds/binding_Kds.csv")
    indf$name = "Binding"
    endpoint = "delta_log10Ka"
    uncert = "log10SE"
  } else if(datatype %in% c("escape", "log_escape")){
    indf = read.csv("Datasets/SARS-CoV-2-RBD_MAP_Crowe_antibodies/results/escape_scores/scores.csv")
    indf = indf[as.logical(indf$pass_pre_count_filter),]
    #indf = indf[as.logical(indf$pass_ACE2bind_expr_filter),] #optional binding filter
    indf$score_sd = sqrt(indf$score_var)
    indf$score_sd[indf$score_sd <= 0] = min(indf$score_sd[indf$score_sd > 0])/2 #add sd floor half of nonzero minimum
    indf$log_score = log10(indf$score)
    indf$log_score_sd = indf$score_sd/indf$score/log(10)
    indf$name = paste0("Antibody ", indf$name)
    # indf = indf[indf$post_count > 0,]
    if(datatype == "escape") endpoint = "score" else endpoint = "log_score"
    if(datatype == "escape") uncert = "score_sd" else uncert = "log_score_sd"
  } else if(datatype == "express"){
    indf = read.csv("Datasets/SARS-CoV-2-RBD_DMS/results/expression_meanFs/expression_meanFs.csv")
    indf$name = "Expression"
    endpoint = "delta_ML_meanF"
    uncert = "var_ML_meanF"
  } else if(datatype %in% c("bindv2", "bindv3")){
    starrdf = read.csv("Datasets/SARS-CoV-2-RBD_DMS/results/binding_Kds/binding_Kds.csv")
    starrdf$endpoint = starrdf$delta_log10Ka
    starrdf$uncert = starrdf$log10SE
    starrdf$source = "starr"
    
    #now our data
    data1 = read.csv("output/mai_antib_kds_12_20_2022.csv")
    data2 = read.csv("output/mai_antib_kds_02_09_2023.csv")
    maidata = rbind(data1, data2)
    if(datatype == "bindv3"){
      data3 = read.csv("output/mai_antib_kds_09_18_2023.csv")
      maidata = rbind(maidata, data3)
    }

    maidata$Variant[maidata$Variant == "Wuhan_control"] = "Wuhan" #just average it into other measurements
    maidata = maidata[!maidata$Log_Kd %in% "Constant",]
    maidata = maidata[maidata$Antibody %in% c("huACE2-Fc"),] #can't use IgG, but keep ACE2 for now

    maidata$endpoint = as.numeric(maidata$Log_Kd)
    maidata$uncert = as.numeric(maidata$Log_Kd_SE)
    wuhan.end = mean(maidata$endpoint[maidata$Variant == "Wuhan"], na.rm = T)
    maidata$endpoint = wuhan.end - maidata$endpoint #delta_logKa
    # maidata = aggregate(endpoint ~ Variant + Antibody, maidata, mean)
    
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
    endpoint = "endpoint"
    uncert = "uncert"
    
    indf = rbind(maidata[, keepcols], starrdf[, keepcols])
    
    indf$name = "Binding"
    
  }
  
  #use generic endpoint/uncert names in data, remove missing value endpoint data
  colnames(indf)[colnames(indf) %in% endpoint] = "endpoint"
  colnames(indf)[colnames(indf) %in% uncert] = "uncert"
  indf = indf[!is.na(indf$endpoint),]
  
  saveRDS(indf, file = paste0("input/", datatype, ".rds"))
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