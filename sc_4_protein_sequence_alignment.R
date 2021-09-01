library(Biostrings)
library(protr)

mySequences <- readAAStringSet("data/dedup4.OR.fa", format = "fasta") # the protein sequence downloaded from https://www.ensembl.org/Mus_musculus/Info/Index
mySequences = removeGaps(mySequences,pattern = " ")
d2 = parSeqSim(mySequences, cores = 4, type = "global")
d2[is.na(d2)]=0
rownames(d2) = names(mySequences)
colnames(d2) = names(mySequences)
saveRDS(d2,"data/OR_similarity_score_defult.rds")
