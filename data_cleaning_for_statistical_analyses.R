### 6/13/17 This R script is for preparing the input files for statistical analyses.
### Specifically, I am calulating the liver:mean expression ratio as well as 
### combining all predictors (N-grams[1:6] + TFBS_counts) into a single dataset.  
### I am also taking the top 10% and bottom 10% of all genes on this ratio,
### making these the two groups for classification (N~1700 genes in each group).

library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(moments)

datadir = "/Users/tsoare/Documents/Insight/data"


# load N-gram data
ngrams <- read.delim(file=file.path(datadir, 
                                    "ngram_data.csv"), 
                     sep=",", header=F, stringsAsFactors = F)
# fix column names
ngrams[1,1] <- "gene.name"
ng.names = ngrams[1, ]
ng.names = gsub("[']", "", ng.names)
ng.names = gsub("[,]", "", ng.names)
ng.names = gsub(" ", "", ng.names)
ng.names = gsub("[(]", "", ng.names)
ng.names = gsub("[)]", "", ng.names)
# re-load data without header
ngrams <- read.delim(file=file.path(datadir, 
                                    "ngram_data.csv"), 
                     sep=",", header=F, skip=1)
names(ngrams) <- ng.names
rm(ng.names)


# calc background for MOODS
bases <- c("A", "C", "G", "T")
base.counts <- apply(ngrams[, colnames(ngrams) %in% bases], 
                     2, sum)
base.counts
#      A       C       G       T 
#4191139 4786108 4741373 4096940 
sum(base.counts)
#[1] 17815560
round(base.counts/sum(base.counts), 3)
#    A     C     G     T 
#0.235 0.269 0.266 0.230 

############
## run PWM search in each sequence in MOODS
## compile results in MOODS results.R
############


# load TFBS motif data
tfbs1 <- read.delim(file=file.path(datadir, "TFBS_counts.csv"), 
                   sep=",", header=T, stringsAsFactors = F)

tfbs2 <- read.delim(file=file.path(datadir, "TFBS_max.score.csv"), 
                    sep=",", header=T, stringsAsFactors = F)

tfbs3 <- read.delim(file=file.path(datadir, "TFBS_max.start.position.csv"), 
                    sep=",", header=T, stringsAsFactors = F)

# merge into a single dataframe
tfbs <- join_all(list(tfbs1, tfbs2, tfbs3), by="gene.name")
rm(tfbs1, tfbs2, tfbs3)
X_preds <- inner_join(ngrams, tfbs, by="gene.name")
rm(ngrams, tfbs)


# load raw (not log-transformed) expression data:
raw.expr <- read.delim(file=file.path(datadir, "rna_tissue_wide.csv"), 
                       sep=",", header=T, stringsAsFactors = F)
names(raw.expr)[1] <- "gene.name"



# normalize to mean
raw.expr$mean <- rowMeans(raw.expr[, -1])
# calculate liver:mean expression ratio and check distribution
raw.expr$liver.ratio <- (raw.expr$liver + 1) / (raw.expr$mean + 1)
hist(raw.expr$liver.ratio)
# check distribution of log-transformed ratio
hist(log2(raw.expr$liver.ratio))



# join outcome with predictors
expr = inner_join(X_preds, raw.expr[, c(1, ncol(raw.expr))], 
                 by="gene.name")

# rank genes and then take top and bottom 10%
expr$log2.liver.ratio <- log2(expr$liver.ratio)


expr$decile <- cut(expr$log2.liver.ratio, 
                   breaks=c(quantile(expr$log2.liver.ratio, 
                                     probs = seq(0, 1, by = 0.1))),
                   labels=c(1:10), 
                   include.lowest=TRUE)

expr$highly.liver.expressed <- ifelse(expr$decile == 1, 0, 
                                       ifelse(expr$decile == 10, 1, NA))

expr2 <- expr[!is.na(expr$highly.liver.expressed), ]

# export dataset for analysis
write.table(expr2, file=file.path(datadir, "input_data_top-bottom_classify.csv"), 
            sep=",", row.names=F, col.names=T, quote=F)



# plot liver vs. mean expression with labeled datapoints:
expr3 = inner_join(expr[, c("gene.name", "decile")], 
                   raw.expr[, c("gene.name", "liver", "mean")], 
                   by="gene.name")
expr3$liver.specific <- ifelse(expr3$decile == 10, 1, #highly liver-specific
                               ifelse(expr3$decile == 1, 3, 2))
expr3$liver.specific <- factor(expr3$liver.specific, 
                               levels = c(1:3), 
                               labels = c("On", "Neutral", "Off"))


p <- ggplot(expr3, aes(x=log2(mean + 1), y=log2(liver + 1)), 
            group=liver.specific) + 
  geom_point(aes(color=liver.specific), alpha=0.5) + 
  scale_color_manual("Liver specificity", 
                     values = c("springgreen2", "grey40", "tomato")) + 
  labs(x="Mean expression (log2)", 
       y="Liver expression (log2)") + 
  #geom_abline(slope=1, intercept=0, color="red") + 
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=16), 
        legend.position = c(0.27, 0.8), 
        legend.title = element_text(size=14), 
        legend.text = element_text(size=12), 
        panel.grid = element_blank())
p

savedir = "/Users/tsoare/Documents/Insight/pres"
ggsave(p, file=file.path(savedir, "Liver vs. mean expression.pdf"), 
       h=5, w=7)


