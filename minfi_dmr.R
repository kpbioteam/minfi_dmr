require("minfi", quietly = TRUE)
require("rtracklayer", quietly = TRUE)

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
input8 = args[8]
input9 = args[9]
output1 = args[10]

GRSet <- get(load(input1))

pheno <- read.table(input2,skip = 1)

group <- pheno$V2

group <- factor(group, levels = c(input3, input4))

design.matrix <- model.matrix(~ 0 + group)

colnames(design.matrix) <- levels(group)

maxGap <- as.numeric(input5)

if(is.null(GRSet$cluster)){
  cluster = NULL
  maxGap = maxGap
} else {
  cluster = GRSet$cluster
  maxGap = NULL
}


cutoff <- as.numeric(input6)
B <- as.numeric(input7)
nullMethod <- input8
coef <- 2 #default
verbose <- input9
  
dmrs <- bumphunter(GRSet,
                    design = design.matrix, 
                    cluster = cluster,
                    maxGap = maxGap,
                    coef = coef,
                    cutoff = cutoff, 
                    nullMethod = nullMethod,
                    B = B, 
                    verbose = verbose)

DMRTable <- dmrs$table

meth  <- GRanges(seqnames=DMRTable$chr,
                  ranges=IRanges
                  (start=DMRTable$start,
                  end=DMRTable$end),
                  value_pos=DMRTable$value)

export.bed(meth,output1)
    
