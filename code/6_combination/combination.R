library(devtools)
library(lionessR)
library(igraph)
library(reshape2)
library(limma)
library(data.table)
library(dplyr)
library(stringr)
library(SNFtool)
setwd("C:/Users/durouxd/Documents/Data/April/BrainCancer/6_combination")
source("C:/Users/durouxd/Documents/Data/April/cancerTypes/functions/svm_TraintTest.R")

#############
# Load data #
#############
data=fread("C:/Users/durouxd/Documents/Data/April/BrainCancer/0_Data/dataTrainRNA.csv")
rownames(data)=data$TCGA_ID
yTrain=data.frame(rownames(data), data$Cancer.Type)
colnames(yTrain)=c("sample", "outcome")
yTrain$outcome=as.factor(yTrain$outcome)

Test=fread("C:/Users/durouxd/Documents/Data/April/BrainCancer/0_Data/dataTestRNA.csv")
rownames(Test)=paste(Test$TCGA_ID)
yTest=data.frame(rownames(Test), Test$Cancer.Type)
colnames(yTest)=c("sample", "outcome")
yTest$outcome=as.factor(yTest$outcome)

#Raw rcorr
setwd("C:/Users/durouxd/Documents/Data/April/BrainCancer/1_rawData")
rcorr_gene=as.matrix(fread("Kgene_rcorr.txt"))
rcorr_geneTrain=as.matrix(rcorr_gene[1:nrow(data),1:nrow(data)])
rcorr_geneTest=as.matrix(rcorr_gene[(nrow(data)+1):nrow(rcorr_gene),1:nrow(data)])
rcorr_image=as.matrix(fread("Kimage_rcorr.txt"))
rcorr_imageTrain=as.matrix(rcorr_image[1:nrow(data),1:nrow(data)])
rcorr_imageTest=as.matrix(rcorr_image[(nrow(data)+1):nrow(rcorr_image),1:nrow(data)])
rcorr_early=as.matrix(fread("Kearly_rcorr.txt"))
rcorr_earlyTrain=as.matrix(rcorr_early[1:nrow(data),1:nrow(data)])
rcorr_earlyTest=as.matrix(rcorr_early[(nrow(data)+1):nrow(rcorr_early),1:nrow(data)])
rcorr_average_Train=as.matrix(fread("KTrain_average_rcorr.txt"))
rcorr_average_Test=as.matrix(fread("KTest_average_rcorr.txt"))

#Product
setwd("C:/Users/durouxd/Documents/Data/April/BrainCancer")
productGene=as.matrix(fread("2_singleData/gene/product/simMatrix_600topGenes_0.5corEdges_edd_scaledyes_all.txt"))
product_geneTrain=as.matrix(productGene[1:nrow(data),1:nrow(data)])
product_geneTest=as.matrix(productGene[(nrow(data)+1):nrow(productGene),1:nrow(data)])
productImage=as.matrix(fread("2_singleData/image/product/simMatrix_300topGenes_0.75corEdges_edd_scaledyes_all.txt"))
product_imageTrain=as.matrix(productImage[1:nrow(data),1:nrow(data)])
product_imageTest=as.matrix(productImage[(nrow(data)+1):nrow(productImage),1:nrow(data)])
product_early=as.matrix(fread("4_early/product/simMatrix_900topGenes_0.5corEdges_edd_scaledyes_all.txt"))
product_earlyTrain=as.matrix(product_early[1:nrow(data),1:nrow(data)])
product_earlyTest=as.matrix(product_early[(nrow(data)+1):nrow(product_early),1:nrow(data)])
product_average_train=as.matrix(fread("5_intermediate/product/KTrain_average.txt"))
product_average_test=as.matrix(fread("5_intermediate/product/KTest_average.txt"))

#Lioness
lioness_geneTrain=as.matrix(fread("2_singleData/gene/lioness/simMatrix_600topGenes_0.5corEdges_edd_scaledyes_train.txt"))
lioness_geneTest=as.matrix(fread("2_singleData/gene/lioness/simMatrix_600topGenes_0.5corEdges_edd_scaledyes_test.txt"))
lioness_geneTest=as.matrix(lioness_geneTest[-1,])
lioness_imageTrain=as.matrix(fread("2_singleData/image/lioness/simMatrix_300topGenes_0.75corEdges_edd_scaledyes_train.txt"))
lioness_imageTest=as.matrix(fread("2_singleData/image/lioness/simMatrix_300topGenes_0.75corEdges_edd_scaledyes_test.txt"))
lioness_imageTest=as.matrix(lioness_imageTest[-1,])
lioness_earlyTrain=as.matrix(fread("4_early/lioness/simMatrix_900topGenes_0.75corEdges_edd_scaledyes_train.txt"))
lioness_earlyTest=as.matrix(fread("4_early/lioness/simMatrix_900topGenes_0.75corEdges_edd_scaledyes_test.txt"))
lioness_earlyTest=as.matrix(lioness_earlyTest[-1,])
lioness_average_trainrain=as.matrix(fread("5_intermediate/lioness/KTrain_average.txt"))
lioness_average_test=as.matrix(fread("5_intermediate/lioness/KTest_average.txt"))

setwd("C:/Users/durouxd/Documents/Data/April/BrainCancer/6_combination")
yTrain=as.factor(yTrain$outcome)
yTest=as.factor(yTest$outcome)

###########
# Average #
###########
X <- list(rcorr_geneTrain, product_geneTrain) 
#X <- list(rcorr_imageTrain, product_imageTrain) 
#X <- list(rcorr_earlyTrain, product_earlyTrain) 
#X <- list(rcorr_average_Train, product_average_train) 
#X <- list(rcorr_geneTrain, lioness_geneTrain) 
#X <- list(rcorr_imageTrain, lioness_imageTrain) 
X <- list(rcorr_earlyTrain, lioness_earlyTrain) 
#X <- list(rcorr_average_Train, lioness_average_train) 
Y <- do.call(cbind, X)
Y <- array(Y, dim=c(dim(X[[1]]), length(X)))
KTrain=apply(Y, c(1, 2), mean, na.rm = TRUE)

X <- list(rcorr_geneTest, product_geneTest)
#X <- list(rcorr_imageTest, product_imageTest)
#X <- list(rcorr_earlyTest, product_earlyTest)
#X <- list(rcorr_average_Test, product_average_test)
#X <- list(rcorr_geneTest, lioness_geneTest)
#X <- list(rcorr_imageTest, lioness_imageTest)
X <- list(rcorr_earlyTest, lioness_earlyTest)
#X <- list(rcorr_average_Test, lioness_average_test)
Y <- do.call(cbind, X)
Y <- array(Y, dim=c(dim(X[[1]]), length(X)))
KTest=apply(Y, c(1, 2), mean, na.rm = TRUE)

#######
# SNF #
#######
#K_SNF=SNF(list(productImage, rcorr_image), K=5, t=5)
#K_SNF=SNF(list(product_early, rcorr_early), K=5, t=5)

#K1=productImage
#K2=rcorr_image
#colnames(K2)=colnames(K1)
#K_SNF=SNF(list(K1, K2), K=20, t=20)
#KTrain=K_SNF[1:nrow(data),1:nrow(data)]
#KTest=K_SNF[(nrow(data)+1):nrow(K_SNF),1:nrow(data)]

#######
# SVM #
#######
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
  
  print("top macro-F1 train score")
  print(max(macroF1train,na.rm=TRUE))
  print("corresponding top macro-F1 train score")
  print(macroF1test[[which.max(macroF1train)]])
  
  results=data.frame(c(nsel, cor_edges, max(macroF1test,na.rm=TRUE), macroF1train[[which.max(macroF1test)]], 10^which.max(macroF1test)))
  colnames(results)=c("results")
  rownames(results)=c("nbNodes", "cor_edges", "macroF1test", "corresponding_macroF1train", "Cparam_SVM")
  
  write.csv(results, paste("results_", nsel, "topGenes_", cor_edges, "corEdges_", meth, "_scaled", scaled, "_", fusion, ".txt", sep=""), row.names=T)




