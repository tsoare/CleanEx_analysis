## 6/6/17 This R script is for parsing the expression data.
## Here I am averaging expression data for different transcripts within
## genes as well as cleaning up the tissue names for better downstream processing.

library(plyr)
library(dplyr)
library(tidyr)

# load data
datadir = "/Users/tsoare/Documents/Insight/data"

dat <- read.delim(file=file.path(datadir, 
                                 "rna_tissue.csv"), 
                  sep=",", header=T, 
                  stringsAsFactors = F)

dat <- dat[, c(2:4)]

# how many different genes are represented?
length(unique(dat$Gene.name))
#[1] 19610

# summarize expression across transcripts within a gene
df <- dat %>% group_by(Gene.name, Sample) %>%
  summarize(mean.expr = mean(Value)) %>%
  spread(Sample, mean.expr)

# re-format and fix tissue names:
df <- as.data.frame(df)
grep("[,]", names(df), value=T)
#[1] "cervix, uterine"
names(df)[grep("[,]", names(df), value=F)] <- "cervix_uterine"

# export cleaned data file of expression (in wide format)
write.table(df, file=file.path(datadir, "rna_tissue_wide.csv"), 
          sep=",", row.names = F, quote=F)

