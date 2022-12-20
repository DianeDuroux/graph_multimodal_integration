library(devtools)
library(lionessR)
library(igraph)
library(reshape2)
library(limma)
library(data.table)
library(dplyr)
library(stringr)
library(splitTools)
setwd("/home/durouxd/April/brainCancer/2_FeatureSelection/genes")
source("/home/durouxd/April/functions/netANOVA.R")
source("/home/durouxd/April/functions/svm_TraintTest.R")
set.seed(2022)

###########
# Options #
###########
nsel=400 #Number of most variable features selected
cor_edges=0.5
meth="edd" #Similarity measure between networks
scaled="yes"

print("Options")
nsel
cor_edges
meth
scaled

folds=5

#########################################
# Create 5 stratified test train for cv #
#########################################
#data=fread("/home/durouxd/April/brainCancer/0_Data/dataTrainRNA.csv")

#Stratified splitting into 5 folds
#inds <- partition(data$Cancer.Type, p = c(fold1 = 0.2, fold2 = 0.2, fold3 = 0.2, fold4=0.2, fold5=0.2))

#for(i in 1:folds){
#  holdout=inds[[i]]
#  test=data[holdout, ]
#  train=data[-holdout, ]
#  fwrite(train, paste("train", i, ".txt", sep=""))
#  fwrite(test, paste("test", i, ".txt", sep=""))
#}

##################################
#  Apply model on each splitting #
##################################

final_output=matrix(nrow=6)
for(fold in 1:folds){

##########################################
# Load gene expression data and reformat #
##########################################
data=fread(paste("train", fold, ".txt", sep=""))
rownames(data)=data$TCGA_ID
yTrain=data.frame(rownames(data), data$Cancer.Type)
colnames(yTrain)=c("sample", "outcome")
yTrain$outcome=as.factor(yTrain$outcome)
exp=data.frame(t(data[,5:ncol(data)]))
colnames(exp)=data$TCGA_ID

Test=fread(paste("test", fold, ".txt", sep=""))
rownames(Test)=paste(Test$TCGA_ID)
yTest=data.frame(rownames(Test), Test$Cancer.Type)
colnames(yTest)=c("sample", "outcome")
yTest$outcome=as.factor(yTest$outcome)
exptest=data.frame(t(Test[,5:ncol(Test)]))
colnames(exptest)=Test$TCGA_ID

exp=cbind(exp, exptest)
y=rbind(yTrain, yTest)

#####################
# Feature selection #
#####################

#select the k most variably expressed genes 
cvar <- apply(as.array(as.matrix(exp)), 1, sd)
dat <- cbind(cvar, exp)
dat <- dat[order(dat[,1], decreasing=T),]
genesToKeep=rownames(dat[1:nsel,])
#fwrite(data.frame(genesToKeep), paste("variablesImages", nsel, ".txt", sep=""))
dat <- as.matrix(dat[genesToKeep, -1])

#outcome-specific networks to select edges that have large absolute differences in their co-expression levels between each pair of scores
group=list()
for(i in 1:length(levels(yTrain$outcome))){
  group[[i]]=which(as.character(yTrain$outcome)==levels(yTrain$outcome)[i])
}
names(group)=levels(yTrain$outcome)
pairs=combn(levels(yTrain$outcome), 2)

if(is.null(dim(pairs))){ # if we compare 2 groups
  len=1
} else { len=ncol(pairs) } # if more than 2 groups

tosel=c()
for(i in 1:len){
  if(len==1){ # if we compare 2 groups
    group1=pairs[1]
    group2=pairs[2]
  } else { # if more than 2 groups
    group1=pairs[1,i]
    group2=pairs[2,i]
  }
  group1=group[group1]
  group2=group[group2]
  netyes <- cor(t(dat[,group1[[1]]]))
  netno <- cor(t(dat[,group2[[1]]]))
  netdiff <- netyes-netno
  cormat2 <- rep(1:nsel, each=nsel) #convert adjacency matrices to edgelists
  cormat1 <- rep(1:nsel,nsel)
  el <- cbind(cormat1, cormat2, c(netdiff))
  melted <- reshape2::melt(upper.tri(netdiff)) #As this is a symmetric adjacency matrix, we takeconvert the upper triangle of the co-expression adjacency matrix into an edge list.
  melted <- melted[which(melted$value),]
  values <- netdiff[which(upper.tri(netdiff))]
  melted <- cbind(melted[,1:2], values)
  genes <- row.names(netdiff)
  melted[,1] <- genes[melted[,1]]
  melted[,2] <- genes[melted[,2]]
  row.names(melted) <- paste(melted[,1], melted[,2], sep=";")
  tosub <- melted
  tosel <- c(tosel, row.names(tosub[which(abs(tosub[,3])>cor_edges),])) #elect edges that have a difference in Pearson R correlation coefficient of at least 0.5
}
tosel=unique(tosel)

print("Number od edges selected")
length(tosel)

################################
# Individual specific networks #
################################
#Min Max values scaling
if(scaled=="yes"){
  all_values=unlist(exp)
  min_all=min(all_values)
  max_all=max(all_values)
  exp=as.data.frame((exp-min_all)/(max_all-min_all))
}

dat=t(dat)
data=matrix(nrow=nrow(dat))
for(i in 1:length(tosel)){ #for each selected variable pair
  edge=str_split_fixed(tosel[i], ";",2)
  tmp=dat[,which(colnames(dat) %in% edge)]
  data=cbind(data, tmp[,1]*tmp[,2])
}
data=data[,-1]
colnames(data)=tosel

#Min Max values scaling
if(scaled=="yes"){
  all_values=unlist(data)
  min_all=min(all_values)
  max_all=max(all_values)
  data=as.data.frame((data-min_all)/(max_all-min_all))
}

#Create 2 separated columns for the variable names
edges=as.data.frame(str_split_fixed(colnames(data), ";",2))
colnames(edges)=c("gene1","gene2")
data=t(data)
data=cbind(edges, data)

#Convert input to the correct format for ANOVA
data_all=list()
el=as.matrix(data[,c(1,2,3)]) 
lab=names(table(el[,1:2])) #extract the existing node IDs
print("Number of nodes")
print(length(unique(lab)))
for(i in 3:ncol(data)){
  el=as.matrix(data[,c(1,2,i)])
  mat=matrix(0,nrow=length(lab),ncol=length(lab),dimnames=list(lab,lab)) #create a matrix of 0s with the node IDs as rows and columns
  mat[el[,1:2]] <-as.numeric(el[,3])
  mat[el[,2:1]] <- as.numeric(el[,3])
  data_all[[i-2]]=mat
}

print("Number of individuals")
print(length(data_all))
print("Number of nodes/genes")
print(dim(data_all[[1]]))

##Similarity matrix
if(meth=="edd"){
  d=nd.edd(data_all)$D
}
d=as.matrix(d)
sim=as.matrix(1/(1+d))

colnames(sim)=y$sample
rownames(sim)=y$sample

KTrain=sim[1:nrow(yTrain),1:nrow(yTrain)]
KTest=sim[(nrow(yTrain)+1):nrow(y),1:nrow(yTrain)]

#fwrite(data.table(sim), paste("simMatrix_", nsel, "topGenes_", cor_edges, "corEdges_", meth, "_scaled", scaled, "_all.txt", sep=""), col.names=T)
#fwrite(data.table(KTrain), paste("simMatrix_", nsel, "topGenes_", cor_edges, "corEdges_", meth, "_scaled", scaled, "_train.txt", sep=""), col.names=T)
#fwrite(data.table(KTest), paste("simMatrix_", nsel, "topGenes_", cor_edges, "corEdges_", meth, "_scaled", scaled, "_test.txt", sep=""), col.names=T)

#######
# SVM #
#######
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
print("top macro-F1 train score")
print(max(macroF1train,na.rm=TRUE))
print("corresponding top macro-F1 test score")
print(macroF1test[[which.max(macroF1train)]])
print("Best C")
print(10^which.max(macroF1train))

results=data.frame(c(nsel, cor_edges, length(tosel), max(macroF1train,na.rm=TRUE), macroF1test[[which.max(macroF1train)]], 10^which.max(macroF1train)))
colnames(results)=c("results")
rownames(results)=c("nbNodes", "cor_edges","nb_edges", "macroF1train", "corresponding_macroF1test", "Cparam_SVM")

write.csv(results, paste("results_", nsel, "topGenes_", cor_edges, "corEdges_", meth, "_scaled", scaled,"_", fold, ".txt", sep=""), row.names=T)
  final_output=cbind(final_output, results )
}

final_output=final_output[,-1]
final_output$mean=rowMeans(final_output)
colnames(final_output)=c(colnames(final_output[,1:(ncol(final_output)-1)]), "mean")
write.csv(final_output, paste("results_", nsel, "topGenes_", cor_edges, "corEdges_", meth, "_scaled", scaled, "_all.txt", sep=""), row.names=T)
