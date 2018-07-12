require("minfi", quietly = TRUE)

options(warn = -1)
options("download.file.method"="wget")

loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

args <- commandArgs(trailingOnly = TRUE)
input1 = args[1]
input2 = args[2]
input3 = args[3]
input4 = args[4]
input5 = args[5]
input6 = args[6]
input7 = args[7]
output1 = args[8]

GRSet <- get(load(input1))

pheno <- read.table(input2,skip = 1)

group <- pheno$V2

pair <- factor(pheno$V3)

design.matrix <- model.matrix(~ group + pair)

maxGap <- as.numeric(input3)

if(is.null(GRSet$cluster)){
  cluster = NULL
  maxGap = maxGap
} else {
  cluster = GRSet$cluster
  maxGap = NULL
}


cutoff <- as.numeric(input4)
B <- as.numeric(input5)
nullMethod <- input6
coef <- 2 #default
verbose <- input7

dmrs <- bumphunter(GRSet,
                    design = design.matrix, 
                    cluster = cluster,
                    maxGap = maxGap,
                    cutoff = cutoff, 
                    nullMethod = nullMethod,
                    B = B, 
                    verbose = verbose)


dmrGR <- with(dmrs$table,GRanges(chr,IRanges(start,end),area=area,value=value))

dmrGR <- as.data.frame(dmrGR)

colnames(dmrGR) <- c("seqnames","start","end","width","strand","area","value")

dmrGR$strand <- NULL

dmrGR$area <- NULL

write.table(dmrGR, file= output1, quote = FALSE,col.names = FALSE, row.names = FALSE, sep = "\t")
