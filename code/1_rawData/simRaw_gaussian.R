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
library(KRLS)
setwd("C:/Users/durouxd/Documents/Data/April/BrainCancer/1_rawData")
source("C:/Users/durouxd/Documents/Data/April/cancerTypes/functions/svm_TraintTest.R")

#normalizes a kernel matrix by dividing through the square root product of the corresponding diagonal entries
#this is not a linear operation, so it should be trated as a hyperparameter
normalisation=function(K){
  D = diag(1/sqrt(diag(K)))
  #ensures that only non zero entries will be subjected to the normalization procedure
  #remaining entries are kept to 0. It prevents NaN values from cropping up
  D = ifelse(D==Inf, 0, D)
  K_norm = D %*% K %*% D
  return(K_norm)
}

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

#Gaussian kernel
Kall=gausskernel(X=exp, sigma=1000)
Kall=normalisation(Kall)
fwrite(Kall, "Kgene_gaussian.txt")

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

#Gaussian kernel
Kall=gausskernel(X=exp, sigma=1000)
Kall=normalisation(Kall)
fwrite(Kall, "Kimage_gaussian.txt")

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

#Gaussian kernel
Kall=gausskernel(X=exp, sigma=1000)
Kall=normalisation(Kall)
fwrite(Kall, "Kearly_gaussian.txt")

#######
# SNF #
#######
KTrain_gene=as.matrix(fread("Kgene_gaussian.txt"))
KTrain_images=as.matrix(fread("Kimage_gaussian.txt"))
Ktrain=SNF(list(KTrain_gene, KTrain_images), K=5, t=5)
Ktrain=normalisation(Ktrain)
#fwrite(Ktrain, "K_SNF_gaussian.txt")

KTrain=Ktrain[1:nrow(data),1:nrow(data)]
KTest=Ktrain[(nrow(data)+1):nrow(Ktrain),1:nrow(data)]

###########
# Average #
###########
K_GeneExpression=as.matrix(fread("Kgene_gaussian.txt"))
KTrain_GeneExpression=K_GeneExpression[1:nrow(data),1:nrow(data)]
KTest_GeneExpression=K_GeneExpression[(nrow(data)+1):nrow(K_GeneExpression),1:nrow(data)]

K_images=as.matrix(fread("Kimage_gaussian.txt"))
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

fwrite(KTrain, "KTrain_average_gaussian.txt")
fwrite(KTest, "KTest_average_gaussian.txt")




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
print(max(macroF1test,na.rm=TRUE))
print("corresponding top macro-F1 test score")
print(macroF1test[[which.max(macroF1train)]])

