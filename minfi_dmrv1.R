require("minfi", quietly = TRUE)

options(warn = -1)
options("download.file.method"="wget")

args <- commandArgs(trailingOnly = TRUE)
input1 = args[1]
input2 = args[2]
input3 = args[3]
input4 = args[4]
input5 = args[5]
input6 = args[6]
output1 = args[7]


GRset <- get(load(input1))

pheno <- read.csv(input2)

designMatrix <- model.matrix(~ pheno$Phenotype)

cutoff <- as.numeric(input3)

B <- as.numeric(input4)

pickCutoff <- input5

pickCutoffQ <- as.numeric(input6)

dmrs <- bumphunter(GRset, design = designMatrix,cluster = cluster, maxGap = maxGap, cutoff = cutoff, B = B, pickCutoff = pickCutoff, pickCutoffQ = pickCutoffQ, smooth = smooth, smoothFunction = smoothFunction)

DMRTable <- dmrs$table

write.table(DMRTable, output1)

