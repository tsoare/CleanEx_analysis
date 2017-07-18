### 6/13/17 This script is for examining the results of the MOODS analysis.
### The MOODS analysis searched each promoter sequence for known binding sites
### (using position weight matrices) and now I need to derive relevant features
### from the output: number of matches to each TFBS, maximum score of each TFBS,
### and location in base pairs (0-1000) of the best match for each TFBS.


library(plyr)
library(dplyr)
library(tidyr)

# load data:
datadir = "/Users/tsoare/Documents/Insight/data"

r <- read.delim(file=file.path(datadir, "moods_out.txt"), 
                sep=",", header=F, stringsAsFactors = F)
names(r) <- c("Seq_info", "TFBS_ID", "start_pos", "strand", 
              "score", "seq", "extra_col")

# check distribution of scores
hist(r$score)


# need to de-duplicate complementary sequences
# get gene name
Seq_name = sapply(strsplit(r$Seq_info, " "), "[[", 2)
r$gene.name = sapply(strsplit(Seq_name, "_"), "[[", 1)
rm(Seq_name)

# clean up TFBS names
tf <- sapply(strsplit(r$TFBS_ID, ".pfm"), "[[", 1)
r$TFBS_ID <- tf
rm(tf)

# de-duplicate
r <- r %>% distinct(gene.name, TFBS_ID, start_pos, .keep_all=T)
# subset to relevant columns
r <- r[, c(8, 2:5)]

# how many putative TFBS per gene?
j <- as.data.frame(table(r$gene.name))
table(j$Freq)
hist(j$Freq)
summary(j$Freq)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   24.0    95.0   117.0   124.5   146.0   672.0

length(unique(r$gene.name))
#[1] 17780
## so each gene has at least 33 TFBS matches (at input settings)
length(unique(r$TFBS_ID))
#[1] 386
## and each TFBS has a putative match with at least one gene

j <- as.data.frame(table(r$TFBS_ID))
table(j$Freq)
hist(j$Freq)
summary(j$Freq)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#    542    2701    3858    5733    5607   82112 
## and each TFBS is present multiple times...


# create matrix of genes x TFBSs: counts of TFBSs
gene.level.matrix <- r %>% group_by(gene.name, TFBS_ID) %>%
  summarize(count = n()) %>% 
  spread(TFBS_ID, count)

g <- as.data.frame(gene.level.matrix)
g[is.na(g)] <- 0
g[1:5, 1:5] # sanity check

write.table(g, file=file.path(datadir, "TFBS_counts.csv"), 
            sep=",", row.names=F, col.names=T, quote=F)
rm(g, gene.level.matrix)


# get max match score for each TFBS within a gene
gmat2 <- r %>% group_by(gene.name, TFBS_ID) %>%
  summarize(max.score = max(score)) %>% 
  spread(TFBS_ID, max.score)
g2 <- as.data.frame(gmat2)
g2[is.na(g2)] <- 0
g2[1:5, 1:5] # sanity check
colnames(g2)[2:ncol(g2)] <- paste0(colnames(g2)[2:ncol(g2)], "_max.score")
write.table(g2, file=file.path(datadir, "TFBS_max.score.csv"), 
            sep=",", row.names=F, col.names=T, quote=F)
rm(g2, gmat2)


# get max position (i.e. closest to TSS)
gmat3 <- r %>% group_by(gene.name, TFBS_ID) %>%
  summarize(max.position = max(start_pos)) %>% 
  spread(TFBS_ID, max.position)
g3 <- as.data.frame(gmat3)
g3[is.na(g3)] <- 0
g3[1:5, 1:5] # sanity check
colnames(g3)[2:ncol(g3)] <- paste0(colnames(g3)[2:ncol(g3)], 
                                   "_max.start.pos")
write.table(g3, file=file.path(datadir, "TFBS_max.start.position.csv"), 
            sep=",", row.names=F, col.names=T, quote=F)
rm(g3, gmat3, r)




