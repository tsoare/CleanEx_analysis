## 6/1/17 get all gene accession numbers:
install.packages('plyr')
install.packages('dplyr')
install.packages('tidyr')

## install bioc:
#source("https://bioconductor.org/biocLite.R")
#biocLite()
#biocLite("org.Hs.eg.db")
library(org.Hs.eg.db)

# get #s:
## Bimap interface:
x <- org.Hs.egACCNUM
# Get the entrez gene identifiers that are mapped to an ACCNUM
mapped_genes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_genes])
if(length(xx) > 0) {
  # Get the ACCNUM for the first five genes
  xx[1:5]
  # Get the first one
  #xx[[1]]
}
# too many.  I think I want NM here, but proceed with gene symbols
rm(x)

# Convert the object to a list
xx <- as.list(org.Hs.egALIAS2EG)
# Remove pathway identifiers that do not map to any entrez gene id
xx <- xx[!is.na(xx)]
names(xx)[1:10]

rm(list=ls())


# downloaded from web:
#

library(plyr)
library(dplyr)
library(tidyr)

datadir = "/Users/tsoare/Documents/Insight/data"

genes <- read.delim(file=file.path(datadir, 
                                   "protein-coding_gene.txt"), 
                    sep="", header=T)

list <- unique(genes$symbol)
write.table(list, file=file.path(datadir, "gene list.txt"), 
            sep=",", row.names = F, col.names = F, quote=F)

