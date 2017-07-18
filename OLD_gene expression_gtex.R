## 6/2/17 get some descriptives of expression

library(plyr)
library(dplyr)
library(tidyr)

datadir = "/Users/tsoare/Documents/Insight/data"

expr = read.delim(file=file.path(datadir, 
                                 "GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct"), 
                  sep="\t", header=T, stringsAsFactors = F, skip=2)

names(expr)
expr$Description[1:10]
#1] "DDX11L1"      "WASH7P"       "MIR1302-11"   "FAM138A"     
#[5] "OR4G4P"       "OR4G11P"      "OR4F5"        "RP11-34P13.7"
#[9] "CICP27"       "AL627309.1"  

length(unique(expr$Description))
#[1] 54301
## too many genes.
## need to truncate the .X off the end

expr2 = expr %>% select(Description, Liver, Brain...Cerebellum) %>%
  group_by(Description) %>%
  summarize(mean.liver = mean(log10(Liver + 1)), 
            mean.cereb = mean(log10(Brain...Cerebellum + 1)))

plot(mean.liver ~ mean.cereb, data=expr2)
# need nice density plot


# what about liver vs. mean of all tissues?
expr2b = expr
expr2b$all = rowMeans(expr2b[, 3:ncol(expr)])
expr2b = expr2b %>% group_by(Description) %>%
  summarize(mean.liver = mean(log2(Liver + 1)), 
            mean.all = mean(log2(all + 1)))
plot(mean.liver ~ mean.all, data=expr2b)
# plot density plot with 1,1 line

#install.packages('ggplot2')
#install.packages('gridExtra')
library(ggplot2)



# PCA of log expression values?
expr3 = expr
expr3[, 3:ncol(expr3)] <- lapply(expr3[, 3:ncol(expr3)], 
                                 function(x) log2(x+1))
pcs <- prcomp(~ ., data=expr3[, 3:ncol(expr3)])

p <- pcs$rotation
plot(PC1 ~ PC2, data=p)
# plot with labels: which tissue(s) are most different?

#x <- pcs$x
#plot(PC1 ~ PC2, data=x)
pp = as.data.frame(p[, 1:2])
pp$label = row.names(pp)
pp$label2 = ifelse(pp$label == "Liver", "Liver", NA)
pp$group = factor(ifelse(pp$label == "Liver", 1, 0), 
                  levels= c("1", "0"))

library(ggrepel)
g <- ggplot(pp, aes(x=PC2, y=PC1, group=group, label=label2)) + 
  geom_point(aes(color=group), size=3, alpha=0.9) + 
  scale_color_manual(NULL, values=c("chartreuse", "grey40")) + 
  geom_label(fill = "chartreuse", position="jitter") + 
  #geom_text_repel(aes(label=label2, size=6), show_guide = FALSE) + 
  theme_bw() + 
  theme(axis.title = element_text(size=14), 
        axis.text = element_text(size=12), 
        legend.position = "none")
g
# liver is reasonably separated from many other tissues
savedir = "/Users/tsoare/Documents/Insight/pres"
#ggsave(g, file=file.path(savedir, "Tissue expression PCA.pdf"), 
#       h=5, w=7)
ggsave(g, file=file.path(savedir, "Tissue expression PCA_label.pdf"), 
       h=5, w=7)

# what about liver vs. mean of all tissues except liver?
expr2c = expr
other.tissues = setdiff(colnames(expr2c[, 3:ncol(expr2c)]), 
                        "Liver")
expr2c$other = rowMeans(expr2c[, names(expr2c) %in% other.tissues])
expr2c = expr2c %>% group_by(Description) %>%
  summarize(mean.liver = mean(log10(Liver + 1)), 
            mean.other = mean(log10(other + 1)))
plot(mean.liver ~ mean.other, data=expr2c)
# +1 due to issue with NA values from log-transformation


# take residuals from linear regression:
fit <- lm(mean.liver ~ mean.other, data=expr2c)
summary(fit)
expr2c$resid.other <- fit$resid
# plot distribution of residuals:
hist(expr2c$resid.other)

# rank genes by +1 SD positive residuals:
m = mean(expr2c$resid.other)
s = sd(expr2c$resid.other)
expr2c$liver.specific = ifelse(expr2c$resid.other > (m + s), 1, 0)
table(expr2c$liver.specific, useNA="ifany")
#
#    0     1 
#50972  3329 
## need to get down to gene-level.

