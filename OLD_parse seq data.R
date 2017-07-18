# 6/5/17 parse file already:
library(plyr)
library(dplyr)
library(tidyr)


datadir = "/Users/tsoare/Documents/Insight/data"

dat <- read.delim(file=file.path(datadir, 
                                   "promoter_exp.fa"), 
                    sep="", header=F, 
                  stringsAsFactors = F)[c(TRUE, FALSE, FALSE), ]

length(unique(dat$V2))
#[1] 17785

# split string
g <- sapply(strsplit(dat$V2, "_"), `[`, 1)
length(unique(g))


# get expression:
dat <- read.delim(file=file.path(datadir, 
                                 "promoter_exp.fa"), 
                  sep="", header=F, 
                  stringsAsFactors = F)[c(FALSE, TRUE, FALSE), ]
expr <- sapply(strsplit(dat$V1, ">"), `[`, 2)
e <- strsplit(expr, ",")
e2 <- do.call(rbind, e)
e3 <- as.data.frame(e2)
e3[, ] <- sapply(e3[, ], function(x) as.numeric(levels(x))[x])

tis.name <- c("adipose tissue","adrenal gland","appendix","bone marrow","breast","cerebral cortex","cervix"," uterine","colon","duodenum","endometrium","epididymis","esophagus","fallopian tube","gallbladder","heart muscle","kidney","liver","lung","lymph node","ovary","pancreas","parathyroid gland","placenta","prostate","rectum","salivary gland","seminal vesicle","skeletal muscle","skin","small intestine","smooth muscle","spleen","stomach","testis","thyroid gland","tonsil","urinary bladder")


