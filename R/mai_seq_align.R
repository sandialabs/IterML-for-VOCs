library(seqinr)

mai_seq_align = function(){
  
  var_seqs = read.csv("Datasets/variant_sequences.csv")
  var_seqs$sequence = trimws(var_seqs$sequence)
  short_seqs = strsplit(var_seqs$sequence, "")
  short_seqs = lapply(short_seqs, function(x){x[13:213]}) #Mai's sequences are R319-N577. Starr's are N331-T531
  #note, we need to be more careful if there are variant rbds longer than the original (insertions)
  
  write.fasta(short_seqs, var_seqs$name, "Datasets/var_to_align.fasta")
  
  #paste this into command line to run muscle
  # C:\Users\tysheff\covid_variants\muscle-5.1\muscle.exe -align C:\Users\tysheff\covid_variants\Datasets\var_to_align.fasta -output C:\Users\tysheff\covid_variants\Datasets\var_aligned.fasta


}