library(devtools)
library(lionessR)
library(igraph)
library(reshape2)
library(limma)
library(data.table)
library(dplyr)
library(stringr)
setwd("/home/durouxd/April/brainCancer/3_singleSource/gene/lioness")
source("/home/durouxd/April/functions/netANOVA.R")
source("/home/durouxd/April/functions/svm_TraintTest.R")

###########
# Options #
###########
nsel=600 #Number of most variably expressed genes selected
cor_edges=0.5
meth="edd" #Similarity measure between networks
scaled="yes"

print("Options")
nsel
cor_edges
meth
scaled

##########################################
# Load gene expression data and reformat #
##########################################
data=fread("/home/durouxd/April/brainCancer/0_Data/dataTrainRNA.csv")
rownames(data)=data$TCGA_ID
yTrain=data.frame(rownames(data), data$Cancer.Type)
colnames(yTrain)=c("sample", "outcome")
yTrain$outcome=as.factor(yTrain$outcome)
exp=data.frame(t(data[,5:ncol(data)]))
colnames(exp)=data$TCGA_ID

Test=fread("/home/durouxd/April/brainCancer/0_Data/dataTestRNA.csv")
rownames(Test)=paste(Test$TCGA_ID)
yTest=data.frame(rownames(Test), Test$Cancer.Type)
colnames(yTest)=c("sample", "outcome")
yTest$outcome=as.factor(yTest$outcome)
exptest=data.frame(t(Test[,5:ncol(Test)]))
colnames(exptest)=Test$TCGA_ID

#####################
# Feature selection #
#####################

#select the k most variably expressed genes 
cvar <- apply(as.array(as.matrix(exp)), 1, sd)
dat <- cbind(cvar, exp)
dat <- dat[order(dat[,1], decreasing=T),]
genesToKeep=rownames(dat[1:nsel,])
fwrite(data.frame(genesToKeep), paste("variablesImages", nsel, ".txt", sep=""))
dat <- as.matrix(dat[genesToKeep, -1])
Test=as.matrix(exptest[genesToKeep, ])

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

#################################
# Train set ISNs and similarity #
#################################
#Train set
#model the single-sample networks based on co-expression using lionessR
cormat <- lioness(dat, netFun)
row.names(cormat) <- paste(cormat[,1], cormat[,2], sep=";")
corsub <- cormat[which(row.names(cormat) %in% tosel),3:ncol(cormat)] #subset networks to the selection of edges which we had defined above
corsub <- as.matrix(corsub)
print("Corsub dimension")
dim(corsub)
#write.table(corsub ,file=paste("ISNs_", nsel, "topGenes", cor_edges, "corEdgesTrain.txt", sep="")) # keeps the rownames

#prepare for netANOVA
data=as.data.frame(t(corsub))
data <- mutate_all(data, function(x) as.numeric(as.character(x)))
sample=rownames(data)

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
#init=initialization(data_all, meth=meth)
#sim=init[[2]]
if(meth=="edd"){
  d=nd.edd(data_all)$D
}
d=as.matrix(d)
sim=1/(1+d)

colnames(sim)=sample
rownames(sim)=sample
fwrite(data.table(sim), paste("simMatrix_", nsel, "topGenes_", cor_edges, "corEdges_", meth, "_scaled", scaled, "_train.txt", sep=""), col.names=T)
KTrain=as.matrix(sim)

#################
# Test set ISNs #
#################
KTest=matrix(ncol=ncol(KTrain))
colnames(KTest)=colnames(KTrain)

for(ind in 1:ncol(Test)){
  print(paste(((ind/ncol(Test))*100), "%", sep=""))
  tmpdat=cbind(dat,Test[,ind, drop=F])
  Testcormat=lioness(tmpdat, netFun)
  row.names(Testcormat) <- paste(Testcormat[,1], Testcormat[,2], sep=";")
  Testcorsub <- Testcormat[which(row.names(Testcormat) %in% tosel),3:ncol(Testcormat)] #subset networks to the selection of edges which we had defined above
  
  ###############################
  # Similarity between networks #
  ###############################
  Testcorsub <- as.matrix(Testcorsub)
  data=as.data.frame(t(Testcorsub))
  data <- mutate_all(data, function(x) as.numeric(as.character(x)))
  sample=rownames(data)
  
  #Min Max values scaling
  if(scaled=="yes"){
    all_values=unlist(data)
    min_all=min(all_values)
    max_all=max(all_values)
    data=as.data.frame((data-min_all)/(max_all-min_all))
  }
  
  #Create 2 separated columns for the gene names
  edges=as.data.frame(str_split_fixed(colnames(data), ";",2))
  colnames(edges)=c("gene1","gene2")
  data=t(data)
  data=cbind(edges, data)
  
  #Convert input to the correct format for ANOVA
  data_all=list()
  el=as.matrix(data[,c(1,2,3)]) 
  lab=names(table(el[,1:2])) #extract the existing node IDs
  for(i in 3:ncol(data)){
    el=as.matrix(data[,c(1,2,i)])
    mat=matrix(0,nrow=length(lab),ncol=length(lab),dimnames=list(lab,lab)) #create a matrix of 0s with the node IDs as rows and columns
    mat[el[,1:2]] <-as.numeric(el[,3])
    mat[el[,2:1]] <- as.numeric(el[,3])
    data_all[[i-2]]=mat
  }
  
  #Similarity matrix
  indNames=colnames(data[,3:ncol(data)])
  maxNodes=ncol(data_all[[1]]) #Derive the maximum number of nodes
  sim=c()
  for(i in 1:(length(indNames)-1)){
    tmp=dist_sim(data_all[[length(indNames)]], data_all[[i]], meth=meth, maxNodes)
    sim=c(sim, tmp[[2]])
  }
  sim=data.frame(t(sim))
  colnames(sim)=indNames[1:(length(indNames)-1)]
  rownames(sim)=indNames[length(indNames)]
  KTest=rbind(KTest, sim)
  fwrite(data.table(KTest), paste("simMatrix_", nsel, "topGenes_", cor_edges, "corEdges_", meth, "_scaled", scaled, "_test.txt", sep=""), col.names=T)
}
KTest=as.matrix(KTest[-1,])

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

write.csv(results, paste("results_", nsel, "topGenes_", cor_edges, "corEdges_", meth, "_scaled", scaled, ".txt", sep=""), row.names=T)
