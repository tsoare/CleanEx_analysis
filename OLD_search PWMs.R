### 6/12/17 TFBS motif searching

# install packages:
#source("https://bioconductor.org/biocLite.R")
#biocLite("TFBSTools")
#biocLite("Biostrings")
#biocLite("JASPAR2016")
#install.packages("seqinr")
#biocLite("ggbio")

library(TFBSTools)
library(Biostrings)
library(seqinr)
library(JASPAR2016)
#library(GenomicRanges)
#library(ggbio)

# set wd:
datadir <- ("/Users/tsoare/Documents/Insight/data")

# import fasta sequences
#seq_records <- read.fasta(file = file.path(datadir, 
#                                           "promoter_seq_clean.fa"), 
#                          as.string=T, forceDNAtolower = FALSE)
# convert each to DNAstrings class
#seq_strings <- sapply(seq_records, DNAString)

# better to import as DNAstring initially:
set = readDNAStringSet(file.path(datadir, 
                                 "promoter_seq_clean.fa"), 
                       format="fasta")




# import all PWMs from JASPAR somehow
# (http://bioconductor.org/packages/devel/bioc/vignettes/TFBSTools/inst/doc/TFBSTools.html)

opts <- list()
opts[["species"]] <- 9606
opts[["matrixtype"]] <- "PFM"
#opts[["length"]] <- "7"
#opts[["type"]] <- "SELEX"
#opts[["all_versions"]] <- TRUE
PFMatrixList <- getMatrixSet(JASPAR2016, opts)
PFMatrixList
#PFMatrixList of length 386
#names(386): MA0025.1 MA0030.1 MA0031.1 ... MA0909.1 MA0914.1
## 386 human TFBSs



# search each DNAstring for each PFM and export key metrics
siteset <- searchSeq(PFMatrixList[1], set[1], min.score="60%")

sitesetList <- searchSeq(PFMatrixList[1], set, min.score="60%")
# not working properly.


# try to export .pfm files for python:
for (i in 1:length(names(PFMatrixList))){
  n = names(PFMatrixList)[i]
  m = paste0("PFMatrixList@listData$", n, "@profileMatrix")
  mat <- as.matrix(mget(m))
  write.table(mat, file=file.path(datadir, 
                                  paste0("/matrices/", n, ".pfm")), 
              sep="\t", row.names=F, col.names=F)
}
# also not working properly.
write.table(paste0(names(PFMatrixList), collapse=", "),
            file=file.path(datadir, "human_tfbs.txt"), 
            sep="", quote=F, row.names=F, col.names=F)

# save these for downstream analysis