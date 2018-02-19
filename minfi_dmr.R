require("minfi", quietly = TRUE)

options(warn = -1)
options("download.file.method"="wget")

args <- commandArgs(trailingOnly = TRUE)
input1 = args[1]
input2 = args[2]
cutoff = 0.02
B = 0 
pickCutoffQ = 0.02 
output1 = args[3]


GRset <- get(load(input1))

pheno <- read.csv(input2)

designMatrix <- model.matrix(~ pheno$Phenotype)

dmrs <- bumphunter(GRset, design = designMatrix,
                   cutoff =cutoff, B=B, type="Beta",  pickCutoff=TRUE,
                   pickCutoffQ=pickCutoffQ)

DMRTable <- dmrs$table

write.table(DMRTable, output1)

