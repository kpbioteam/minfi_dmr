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
output1 = args[8]
output2 = args[9]


GRSet <- get(load(input1))

pheno <- read.table(input2,  skip = 1)

designMatrix <- model.matrix(~ pheno$V2)

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
B <- 0 #default
verbose <- input7
  
dmrs <- bumphunter(GRSet,
                    design = designMatrix, 
                    cluster = cluster,
                    maxGap = maxGap,
                    coef = coef,
                    cutoff = cutoff, 
                    nullMethod = nullMethod,
                    B = B, 
                    verbose = verbose)

DMRTable <- dmrs$table

write.table(DMRTable, output1)

sign=sign(DMRTable$value)

sign[sign==-1]="-"
sign[sign==1]="+"

dmr_track<-cbind(as.character(DMRTable$chr),DMRTable$start,DMRTable$end,sign)

write.table(dmr_track,output2,quote = FALSE, sep = "\t\t\t\t",row.names = FALSE, col.names = FALSE,append=TRUE)        
