library(seqinr)


var_muts = function(){
  var_aligned = read.fasta("Datasets/var_aligned.fasta", seqtype = "AA", forceDNAtolower = F)
  #need to match variants up, the only length change is a pair of deletions, none have insertions
  var_aligned = lapply(var_aligned, function(x){x[1:201]})
  
  refseq = var_aligned$Wuhan
  
  muts = lapply(var_aligned, function(varseq){
    diffs = varseq != refseq
    muts = paste0(refseq[diffs], (1:length(varseq))[diffs], varseq[diffs])
    return(paste0(muts, collapse = " "))
  })
  
  outdf = data.frame(Variant = names(muts), aa_subs = unlist(muts), stringsAsFactors = F)
  
  write.csv(outdf, file = "output/var_muts.csv")
}
