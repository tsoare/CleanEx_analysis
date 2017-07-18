### 6/23/17 This Rscript is for making figures for the presentation:
### 1. Plot of feature importances from random forests analysis
### 2. Plot of (directional) effect estimates from logistic regression
### 3. Bargraph of mean accuracy for all models across all numbers of features
### 4. Correlation matrix among features

library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(corrplot)


datadir = "/Users/tsoare/Documents/Insight/data"
savedir = "/Users/tsoare/Documents/Insight/pres"

# load data
rf <- read.delim(file=file.path(datadir, "RF_top20features.csv"), 
                 sep=",", header=T, stringsAsFactors = F)
rf <- rf[, -1]
names(rf) <- c("Var", "Importance")


coefs <- read.delim(file=file.path(datadir, "top20features_LogReg_coefficients.csv"), 
                    sep=",", header=T, stringsAsFactors = F)
coefs <- coefs[, -1]
names(coefs) <- c("Var", "Coefficient")

# rank by RF feature importance and merge:
rf <- rf[order(rf$Importance, decreasing=T), ]
df <- left_join(rf, coefs)


# plot:
df$Var <- factor(df$Var, levels = rev(df$Var), 
                 labels = rev(c("CG", "MA0528 - maximum score", 
                 "CC", "MA0056 - location", "TCT", "CT", "TCA", 
                 "CTG", "TGG", "AC", "MA0149 - location", "AGG", 
                 "AGGG", "MA0508 - location", "MA0528 - count", 
                 "MA0685 - maximum score", "MA0478 - location", 
                 "CGAA", "MA0597 - location", "MA0759 - location")))


# plot combination of #1 and #2
g1 <- ggplot(df, aes(x=Var, y=Importance)) + 
  geom_bar(stat="identity", fill="firebrick2") + 
  labs(x="Variable", 
       y="Variable Importance") + 
  coord_flip() + 
  theme_bw() + 
  theme(axis.text = element_text(size=12), 
        axis.title.y = element_blank())
g1

g2 <- ggplot(df, aes(x=Var, y=Coefficient)) + 
  geom_bar(stat="identity", fill="purple") + 
  geom_hline(yintercept=0, color="black") + 
  labs(y="Effect Size (Standardized)") + 
  coord_flip() + 
  theme_bw()+ 
  theme(axis.text.y = element_blank(), 
        axis.title.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.x = element_text(size=12)) # rm x-axis text
g2

g <- grid.arrange(g1, g2, 
                  layout_matrix=rbind(c(1,1,1,1,1,2,2,2), 
                                      c(1,1,1,1,1,2,2,2)))
ggsave(g, file=file.path(savedir, "Top 20 features.pdf"), 
       h=4, w=8)



## also make just single #2
g2 <- ggplot(df, aes(x=Var, y=Coefficient)) + 
  geom_bar(stat="identity", fill="purple", alpha=0.8) + 
  geom_hline(yintercept=0, color="black", lty=3) + 
  labs(y="Effect Size (Standardized)") + 
  coord_flip() + 
  theme_bw() + 
  theme(axis.text.y = element_text(size=12), 
        axis.title.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.x = element_text(size=12))
g2
ggsave(g2, file=file.path(savedir, "Top 20 features_logreg.effect.sizes.pdf"), 
       h=4, w=8)


## and single version of #1:
g1 <- ggplot(df, aes(x=Var, y=Importance)) + 
  geom_bar(stat="identity", fill="firebrick2", alpha=0.8) + 
  labs(x="Variable", 
       y="Variable Importance") + 
  coord_flip() + 
  theme_bw() + 
  theme(axis.text.y = element_text(size=12), 
        axis.title.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.x = element_text(size=12))
g1
ggsave(g1, file=file.path(savedir, "Top 20 features_RF.importance.pdf"), 
       h=4, w=8)




# here, plot #2 but in order of logistic regression effect size (abs)
df <- df[order(abs(df$Coefficient), decreasing=T), ]
df$Var <- factor(df$Var, levels=rev(df$Var))

g2 <- ggplot(df, aes(x=Var, y=Coefficient)) + 
  geom_bar(stat="identity", fill="purple", alpha=0.8) + 
  geom_hline(yintercept=0, color="black", lty=3) + 
  labs(y="Effect Size (Standardized)") + 
  coord_flip() + 
  theme_bw() + 
  theme(axis.text.y = element_text(size=12), 
        axis.title.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.x = element_text(size=12))
g2
ggsave(g2, file=file.path(savedir, "Top 20 features_logreg.effect.sizes_re-order.pdf"), 
       h=4, w=8)




# now plot #3 the mean accuracy across all models and feature numbers
rm(list=setdiff(ls(), c("datadir", "savedir")))
m <- read.delim(file=file.path(datadir, "RFE_model_output.csv"), 
                 sep=",", header=T, stringsAsFactors = F)
m$model <- factor(m$model, 
                  levels=c("LogReg", "SVM", "RF"), 
                  labels=c("Logistic Regression", 
                           "Support Vector Machine", 
                           "Random Forests"))
m$no.features <- factor(m$no.features, 
                        levels=c("6618", "200", "20", "10"))
m$se = m$sd.accuracy / sqrt(5)

p <- ggplot(m, aes(x=no.features, y=mean.accuracy, 
                   ymin=mean.accuracy-se, ymax=mean.accuracy+se, 
                   group=model)) + 
  geom_bar(aes(fill=model), stat="identity", 
           position="dodge", alpha=0.7) + 
  scale_fill_manual("Model", 
                    values=c("cadetblue1", 
                             "cadetblue3", 
                             "cadetblue4")) + 
  geom_errorbar(position="dodge") + 
  scale_y_continuous(limits=c(0,0.65)) + 
  geom_hline(yintercept=0.5, lty=2) + 
  labs(x="Number of features", 
       y="Mean accuracy") + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.text = element_text(size=14), 
        axis.title = element_text(size=16), 
        legend.title = element_text(size=14), 
        legend.text = element_text(size=12))
p
ggsave(p, file=file.path(savedir, "Mean accuracy for RFE.pdf"), 
       h=5, w=8)



# now calculate correlations among all features and plot
rm(list=setdiff(ls(), c("datadir", "savedir")))
X <- read.delim(file=file.path(datadir, "top20features_X_train.csv"),
                sep=",", header=T, stringsAsFactors = F)
X <- X[, c("gene.name", 
           "CG", "MA0528.1_max.score", 
           "CC", "MA0056.1_max.start.pos", 
           "TCT", "CT", 
           "TCA", "CTG", 
           "TGG", "AC", 
           "MA0149.1_max.start.pos", "AGG", 
           "AGGG", "MA0508.1_max.start.pos", 
           "MA0528.1", "MA0685.1_max.score", 
           "MA0478.1_max.start.pos", "CGAA", 
           "MA0597.1_max.start.pos", "MA0759.1_max.start.pos")]


names(X) <- c("gene.name", "CG", "MA0528 - maximum score", 
              "CC", "MA0056 - location", "TCT", "CT", "TCA", 
              "CTG", "TGG", "AC", "MA0149 - location", "AGG", 
              "AGGG", "MA0508 - location", "MA0528 - count", 
              "MA0685 - maximum score", "MA0478 - location", 
              "CGAA", "MA0597 - location", "MA0759 - location")

# correlation matrix
cor.mat = cor(X[, -1])
cor.coefs =  cor.mat[upper.tri(cor.mat, diag=F)]
max(abs(cor.coefs))
#[1] 0.7590927
## at least one pair highly correlated...
hist(cor.coefs)
## seems reasonable...


# make correlation matrix plot
corrplot(cor.mat, method = "square", #type="lower", 
         #order="hclust", addrect=2, 
         tl.cex=1, tl.col="black", 
         cl.ratio=c(0.3), cl.align.text="l")


