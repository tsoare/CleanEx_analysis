### 6/19/17 This is a script for double-checking the model and 
### creating the files for the web app input.

library(plyr)
library(dplyr)
library(tidyr)

datadir = "/Users/tsoare/Documents/Insight/data"


X_select <- read.delim(file=file.path(datadir, 
                                      "top20features_X_train.csv"), 
                       sep=",", header=T, stringsAsFactors = F)
X_test_select <- read.delim(file=file.path(datadir, 
                                           "top20features_X_test.csv"), 
                            sep=",", header=T, stringsAsFactors = F)

y_train <- read.delim(file=file.path(datadir, 
                                     "y_train.csv"), 
                      sep=",", header=T, stringsAsFactors = F)
y_test <- read.delim(file=file.path(datadir, 
                                    "y_test.csv"), 
                     sep=",", header=F, stringsAsFactors = F)

table(y_train$label)
table(y_test$V2)
## looks ok...


r <- read.delim(file=file.path(datadir, 
                               "model_output.csv"), 
                sep=",", header=T, stringsAsFactors = F)
hist(r$Probability_liver_specific)

r2 <- read.delim(file=file.path(datadir, 
                               "model_output_labels.csv"), 
                sep=",", header=T, stringsAsFactors = F)
table(r2$Class_prediction)


r2b <- r2[r2$Class_prediction == 0, ]
r2b$gene <- sapply(strsplit(r2b$Fasta_header, "_"), `[`, 2)
r2c <- merge(r2c, y_test, by.x="gene", by.y="V1")

r2c$Fasta_header[1:10]

y_test2 = y_test[y_test$V2 == 1, ]
y_test2$V1[1:10]
