library(Biostrings)
library(protr)

mySequences <- readAAStringSet("data/dedup4.OR.fa", format = "fasta") #OR protein sequence downloaded from NCBI
mySequences = removeGaps(mySequences,pattern = " ")
d2 = parSeqSim(mySequences, cores = 4, type = "global")
d2[is.na(d2)]=0
rownames(d2) = names(mySequences)
colnames(d2) = names(mySequences)
saveRDS(d2,"OR_similarity_score_default.rds")
