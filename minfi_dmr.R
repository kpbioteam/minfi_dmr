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
input7 = args[7]
input8 = args[8]
input9 = args[9]
input10 = args[10]
output1 = args[11]


GRSet <- get(load(input1))

pheno <- read.csv(input2)

designMatrix <- model.matrix(~ pheno$Phenotype)

maxGap <- as.numeric(input3)

if(is.null(GRSet$cluster)){
  cluster = NULL
  maxGap = maxGap
} else {
  cluster = GRSet$cluster
  maxGap = NULL
}

coef <- as.numeric(input4)
cutoff <- as.numeric(input5)
nullMethod <- input6 
B <- as.numeric(input7)
verbose <- input8
smooth <- input9
smoothFunction <- input10

if(smooth == FALSE){
  smoothFunction = NULL
} else {
  smoothFunction = smoothFunction
}
  
bumps <- bumphunter(GRSet,
                    design = designMatrix, 
                    cluster = cluster,
                    maxGap = maxGap,
                    coef = coef,
                    cutoff = cutoff, 
                    nullMethod = nullMethod,
                    B = B, 
                    verbose = verbose, 
                    smooth = smooth,
                    smoothFunction = smoothFunction)
                


DMRTable <- dmrs$table

write.table(DMRTable, output1)        
