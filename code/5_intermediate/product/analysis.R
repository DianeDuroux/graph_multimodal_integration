library(devtools)
library(lionessR)
library(igraph)
library(reshape2)
library(limma)
library(data.table)
library(dplyr)
library(stringr)
library(SNFtool)
setwd("/home/durouxd/April/brainCancer/5_intermediate/produit")
source("/home/durouxd/April/functions/svm_TraintTest.R")

nsel_gene=600
nsel_image=300
nsel=nsel_gene+nsel_image
cor_edges_image=0.75
cor_edges_gene=0.5
meth="edd" #Similarity measure between networks
scaled="yes"

#############
# Load data #
#############
#Load outcome
rna=fread("/home/durouxd/April/brainCancer/0_Data/dataTrainRNA.csv")
images=fread("/home/durouxd/April/brainCancer/0_Data/dataTrainImages.csv")
data=merge(rna, images, by=c("TCGA_ID", "slide_id", "Gender", "Cancer.Type"))
rownames(data)=data$TCGA_ID
yTrain=data.frame(rownames(data), data$Cancer.Type)
colnames(yTrain)=c("sample", "outcome")
yTrain$outcome=as.factor(yTrain$outcome)
exp=data.frame(t(data[,5:ncol(data)]))
colnames(exp)=data$TCGA_ID

rna=fread("/home/durouxd/April/brainCancer/0_Data/dataTestRNA.csv")
images=fread("/home/durouxd/April/brainCancer/0_Data/dataTestImages.csv")
Test=merge(rna, images, by=c("TCGA_ID", "slide_id","Gender", "Cancer.Type"))
rownames(Test)=paste(Test$TCGA_ID)
yTest=data.frame(rownames(Test), Test$Cancer.Type)
colnames(yTest)=c("sample", "outcome")
yTest$outcome=as.factor(yTest$outcome)
exptest=data.frame(t(Test[,5:ncol(Test)]))
colnames(exptest)=Test$TCGA_ID

#Gene expression
K_GeneExpression=as.matrix(fread(paste("simMatrix_", nsel_gene, "topGenes_", cor_edges_gene, "corEdges_", meth, "_scaled", scaled, "_all.txt", sep="")))
KTrain_GeneExpression=K_GeneExpression[1:nrow(data),1:nrow(data)]
KTest_GeneExpression=K_GeneExpression[(nrow(data)+1):nrow(K_GeneExpression),1:nrow(data)]

#Images
K_images=as.matrix(fread(paste("simMatrix_", nsel_image, "topGenes_", cor_edges_image, "corEdges_", meth, "_scaled", scaled, "_all.txt", sep="")))
KTrain_images=K_images[1:nrow(data),1:nrow(data)]
KTest_images=K_images[(nrow(data)+1):nrow(K_images),1:nrow(data)]

#SNF
K=SNF(list(K_GeneExpression, K_images), K=20, t=20)
KTrain_SNF=K[1:nrow(data),1:nrow(data)]
KTest_SNF=K[(nrow(data)+1):nrow(K),1:nrow(data)]

#Fusion
X <- list(KTrain_GeneExpression, KTrain_images)
Y <- do.call(cbind, X)
Y <- array(Y, dim=c(dim(X[[1]]), length(X)))
Ktrain_average=apply(Y, c(1, 2), mean, na.rm = TRUE)

X <- list(KTest_GeneExpression, KTest_images)
Y <- do.call(cbind, X)
Y <- array(Y, dim=c(dim(X[[1]]), length(X)))
Ktest_average=apply(Y, c(1, 2), mean, na.rm = TRUE)

KTrain=Ktrain_average
KTest=Ktest_average
fwrite(KTrain, "KTrain_average.txt")
fwrite(KTest, "KTest_average.txt")


fusion_function=c("average", "SNF")
Ktrain=list(Ktrain_average, KTrain_SNF)
Ktest=list(Ktest_average, KTest_SNF)

#######
# SVM #
#######
yTrain=as.factor(yTrain$outcome)
yTest=as.factor(yTest$outcome)

#Parameter C tunning
k=seq(1, 5, 1)
C=c(10^k)

for(type in 1:2){
  fusion=fusion_function[type]
  KTrain=Ktrain[[type]]
  KTest=Ktest[[type]]
  
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
print("top macro-F1 train score")
print(max(macroF1train,na.rm=TRUE))
print("corresponding top macro-F1 test score")
print(macroF1test[[which.max(macroF1train)]])
print("Best C")
print(10^which.max(macroF1train))
  
results=data.frame(c(nsel, max(macroF1train,na.rm=TRUE), macroF1test[[which.max(macroF1train)]], 10^which.max(macroF1train)))
colnames(results)=c("results")
rownames(results)=c("nbNodes",  "macroF1train", "corresponding_macroF1test", "Cparam_SVM")

  write.csv(results, paste("results_", nsel, "topGenes_", meth, "_scaled", scaled, fusion, ".txt", sep=""), row.names=T)

}




