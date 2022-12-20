library(devtools)
library(lionessR)
library(igraph)
library(reshape2)
library(limma)
library(data.table)
library(dplyr)
library(stringr)
library(SNFtool)
library(Hmisc)
library(glmnet)
setwd("C:/Users/durouxd/Documents/Data/April/BrainCancer/1_rawData")
source("C:/Users/durouxd/Documents/Data/April/cancerTypes/functions/svm_TraintTest.R")

###################
# gene expression #
###################
data=fread("C:/Users/durouxd/Documents/Data/April/BrainCancer/0_Data/dataTrainRNA.csv")
rownames(data)=data$TCGA_ID
yTrain=data.frame(rownames(data), data$Cancer.Type)
colnames(yTrain)=c("sample", "outcome")
yTrain$outcome=as.factor(yTrain$outcome)
exp=data.frame(data[,5:ncol(data)])
rownames(exp)=data$TCGA_ID
exp=standardNormalization(exp)
Test=fread("C:/Users/durouxd/Documents/Data/April/BrainCancer/0_Data/dataTestRNA.csv")
rownames(Test)=paste(Test$TCGA_ID)
yTest=data.frame(rownames(Test), Test$Cancer.Type)
colnames(yTest)=c("sample", "outcome")
yTest$outcome=as.factor(yTest$outcome)
exptest=data.frame(Test[,5:ncol(Test)])
rownames(exptest)=Test$TCGA_ID
exptest=standardNormalization(exptest)
exp=rbind(exp, exptest)
y=rbind(yTrain, yTest)

#Spearman correlation
Kall=rcorr(as.matrix(t(exp)), type  = "spearman")
Kall=Kall$r
fwrite(Kall, "Kgene_rcorr.txt")

##########
# images #
##########
data=fread("C:/Users/durouxd/Documents/Data/April/BrainCancer/0_Data/dataTrainImages.csv")
rownames(data)=data$TCGA_ID
yTrain=data.frame(rownames(data), data$Cancer.Type)
colnames(yTrain)=c("sample", "outcome")
yTrain$outcome=as.factor(yTrain$outcome)
exp=data.frame(data[,5:ncol(data)])
rownames(exp)=data$TCGA_ID
exp=standardNormalization(exp)
Test=fread("C:/Users/durouxd/Documents/Data/April/BrainCancer/0_Data/dataTestImages.csv")
rownames(Test)=paste(Test$TCGA_ID)
yTest=data.frame(rownames(Test), Test$Cancer.Type)
colnames(yTest)=c("sample", "outcome")
yTest$outcome=as.factor(yTest$outcome)
exptest=data.frame(Test[,5:ncol(Test)])
rownames(exptest)=Test$TCGA_ID
exptest=standardNormalization(exptest)
exp=rbind(exp, exptest)
y=rbind(yTrain, yTest)

#Spearman correlation
Kall=rcorr(as.matrix(t(exp)), type  = "spearman")
Kall=Kall$r
fwrite(Kall, "Kimage_rcorr.txt")

#######################
# Early concatenation #
#######################
rna=fread("C:/Users/durouxd/Documents/Data/April/BrainCancer/0_Data/dataTrainRNA.csv")
images=fread("C:/Users/durouxd/Documents/Data/April/BrainCancer/0_Data/dataTrainImages.csv")
data=merge(rna, images, by=c("TCGA_ID", "slide_id", "Gender", "Cancer.Type"))
rownames(data)=data$TCGA_ID
yTrain=data.frame(rownames(data), data$Cancer.Type)
colnames(yTrain)=c("sample", "outcome")
yTrain$outcome=as.factor(yTrain$outcome)
exp=data.frame(data[,5:ncol(data)])
colnames(exp)=data$TCGA_ID
exp=standardNormalization(exp)
rna=fread("C:/Users/durouxd/Documents/Data/April/BrainCancer/0_Data/dataTestRNA.csv")
images=fread("C:/Users/durouxd/Documents/Data/April/BrainCancer/0_Data/dataTestImages.csv")
Test=merge(rna, images, by=c("TCGA_ID", "slide_id","Gender", "Cancer.Type"))
rownames(Test)=paste(Test$TCGA_ID)
yTest=data.frame(rownames(Test), Test$Cancer.Type)
colnames(yTest)=c("sample", "outcome")
yTest$outcome=as.factor(yTest$outcome)
exptest=data.frame(Test[,5:ncol(Test)])
colnames(exptest)=Test$TCGA_ID
exptest=standardNormalization(exptest)
exp=rbind(exp, exptest)
y=rbind(yTrain, yTest)

#Spearman correlation
Kall=rcorr(as.matrix(t(exp)), type  = "spearman")
Kall=Kall$r
fwrite(Kall, "Kearly_rcorr.txt")

#######
# SNF #
#######
KTrain_gene=as.matrix(fread("Kgene_rcorr.txt"))
KTrain_images=as.matrix(fread("Kimage_rcorr.txt"))
Ktrain=SNF(list(KTrain_gene, KTrain_images), K=5, t=5)
fwrite(Ktrain, "K_SNF_rcorr.txt")

KTrain=Ktrain[1:nrow(data),1:nrow(data)]
KTest=Ktrain[(nrow(data)+1):nrow(Ktrain),1:nrow(data)]

###########
# Average #
###########
K_GeneExpression=as.matrix(fread("Kgene_rcorr.txt"))
KTrain_GeneExpression=K_GeneExpression[1:nrow(data),1:nrow(data)]
KTest_GeneExpression=K_GeneExpression[(nrow(data)+1):nrow(K_GeneExpression),1:nrow(data)]

K_images=as.matrix(fread("Kimage_rcorr.txt"))
KTrain_images=K_images[1:nrow(data),1:nrow(data)]
KTest_images=K_images[(nrow(data)+1):nrow(K_images),1:nrow(data)]

X <- list(KTrain_GeneExpression, KTrain_images)
Y <- do.call(cbind, X)
Y <- array(Y, dim=c(dim(X[[1]]), length(X)))
KTrain=apply(Y, c(1, 2), mean, na.rm = TRUE)

X <- list(KTest_GeneExpression, KTest_images)
Y <- do.call(cbind, X)
Y <- array(Y, dim=c(dim(X[[1]]), length(X)))
KTest=apply(Y, c(1, 2), mean, na.rm = TRUE)

fwrite(KTrain, "KTrain_average_rcorr.txt")
fwrite(KTest, "KTest_average_rcorr.txt")

#######
# SVM #
#######
KTrain=Kall[1:nrow(data),1:nrow(data)]
KTest=Kall[(nrow(data)+1):nrow(Kall),1:nrow(data)]

yTrain=as.factor(yTrain$outcome)
yTest=as.factor(yTest$outcome)

#Parameter C tunning
k=seq(1, 5, 1)
C=c(10^k)

macroF1train=c()
macroF1test=c()
CMtrain=list()
CMtest=c()
for(i in 1:length(C)){
  k=C[i]
  perf=runsvm(KTrain, yTrain, KTest, yTest, C = k)
  macroF1train <- c(macroF1train, perf[[1]])
  macroF1test <- c(macroF1test, perf[[2]])
  CMtrain[[i]]=perf[[3]]
  CMtest[[i]]=perf[[4]]
}

print("top macro-F1 test score")
print(max(macroF1train,na.rm=TRUE))
print("corresponding top macro-F1 test score")
print(macroF1test[[which.max(macroF1train)]])


