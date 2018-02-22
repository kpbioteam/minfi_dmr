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
input11 = args[11]
output1 = args[12]


GRSet <- get(load(input1))

pheno <- read.csv(input2)

designMatrix <- model.matrix(~ pheno$Phenotype)

cluster <- input3
maxGap <- as.numeric(input4)

if(is.null(GRSet$cluster)){
  cluster = NULL
  maxGap = maxGap
} else {
  cluster = GRSet$cluster
  maxGap = NULL
}

coef <- as.numeric(input5)
cutoff <- as.numeric(input6)
nullMethod <- input7 
B <- as.numeric(input8)
verbose <- input9
smooth <- input10
smoothFunction <- input11

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
