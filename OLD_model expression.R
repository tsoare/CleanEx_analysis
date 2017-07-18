### 6/8/17 modelling exression

library(MASS)
library(plyr)
library(dplyr)
library(tidyr)


# load raw (not log-transformed) expression data:
datadir = "/Users/tsoare/Documents/Insight/data"
raw.expr <- read.delim(file=file.path(datadir, "rna_tissue_wide.csv"), 
                       sep=",", header=T)

cleanex.dir <- "/Users/tsoare/Documents/Insight/flask/clean_ex"
expr <- read.delim(file=file.path(cleanex.dir, "input_data.csv"), 
          sep=",", header=T, stringsAsFactors=F)
names(expr)[1] <- "Gene.name"

expr = expr[, 1:17] # subset to only bigrams
expr = left_join(expr, raw.expr, by="Gene.name")
rm(raw.expr)

expr = expr[, -55] # bladder got chopped off for some reason


# plot:
library(ggplot2)
# 1. PCAs of tissues:
pcs <- prcomp(~ ., data=expr[, 18:54])

p <- pcs$rotation
plot(PC1 ~ PC2, data=p)
plot(PC2 ~ PC3, data=p)

pp = as.data.frame(p[, 1:2])
pp$label = row.names(pp)
g <- ggplot(pp, aes(x=PC2, y=PC1, label=label)) + 
  geom_label()
g
## not promising

# remove ovaries and testes:
expr2 = expr[, setdiff(colnames(expr), c("ovary", "testis"))]
pcs <- prcomp(~ ., data=expr2[, 18:ncol(expr2)])

p <- pcs$rotation
plot(PC1 ~ PC2, data=p)
plot(PC1 ~ PC3, data=p)
plot(PC2 ~ PC3, data=p)

pp = as.data.frame(p[, 1:2])
pp$label = row.names(pp)
g <- ggplot(pp, aes(x=PC2, y=PC1, label=label)) + 
  geom_label()
g
rm(expr2, pcs, p, pp, g)




# normalize to mean:
expr$mean <- rowMeans(expr[, 18:ncol(expr)])
expr[, 18:(ncol(expr) - 1)] <- lapply(expr[, 18:(ncol(expr) - 1)], 
                                      function(x) (x + 0.1)/(expr$mean + 0.1))

expr_train = expr[1:10070, ]

dat_X_train = expr_train[, 2:17]
dat_Y_train = expr_train[, "liver"]
dat = cbind(dat_X_train, "liver"=dat_Y_train)

fit1 <- lm(liver ~ ., data=dat)
summary(fit1)
AIC(fit1)
#[1] 34963.53

fit2 <- glm(liver ~ ., data=dat, family=poisson(link = "log2"))
summary(fit2)
AIC(fit2)
#[1] Inf

fit3 <- glm.nb(liver ~ ., data = dat)
summary(fit3)
AIC(fit3)
#[1] 25680


#pchisq(2 * (logLik(fit2) - logLik(fit3)), df = 1, lower.tail = FALSE)
# -Inf b/c fit2 so bad



# lasso:
X = as.matrix(dat[, 1:16])
y = dat[, 17]
lasso <- lars(X, y)
covTest(lasso, X, y)$results
# 13th predictor or... TA !
max(lasso$R2)
#[1] 0.005177598


# RF:
#set.seed()
rf <- randomForest(liver ~ ., data=dat, ntree=5000)
mean(rf$rsq)
#[1] -0.05090132
round(importance(rf), 2)
## AA most important, TA not near top of list... 


save(lasso, rf, 
     file=file.path(datadir, "bigrams_lasso.and.rf.rda"))


# plot:
# 1. liver vs. mean expression with 1,1 line:
# first log2 expr:
expr[, 18:ncol(expr)] <- lapply(expr[, 18:ncol(expr)], 
                                log2)
g1 <- ggplot(expr, aes(x=mean, y=liver)) + 
  geom_point(color="dodgerblue2", alpha=0.5) + 
  labs(x="Mean expression (log2)", 
       y="Liver expression (log2)") + 
  geom_abline(slope=1, intercept=0, color="red") + 
  theme_bw()
g1
ggsave(g1, file=file.path(datadir, "Liver vs. mean expression.jpg"), 
       h=5, w=7)

# 2. bigram violin plot
bg <- expr[, 1:17] %>% group_by(Gene.name) %>%
  gather(bigram, count, 2:17)
bg$bigram <- factor(bg$bigram, 
                    labels = c("AA", "AC", "AG", "AT", 
                               "CA", "CC", "CG", "CT", 
                               "GA", "GC", "GG", "GT",
                               "TA", "TC", "TG", "TT"))

g2 <- ggplot(bg, aes(y=count, x=bigram, group=bigram)) + 
  geom_violin(fill="firebrick", alpha=0.7) + 
  labs(x="Bigram", y="Frequency") + 
  theme_bw()
g2
ggsave(g2, file=file.path(datadir, "Bigram frequency.jpg"), 
       h=5, w=7)


# 3. r2 for bigram, lasso, rf
# 3. elbow plot for lasso
# 3. rf variable importance


