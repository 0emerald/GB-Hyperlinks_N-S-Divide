## 02/12/2022

## Looking to see what difference Winsorizing before or after concatenating the LAD x LAD A matrices for the UK does.
# Measure this with some idea of 'variance'

library(tidyverse)
library(Matrix)
library(DescTools)
library(expm)
setwd("~/Documents/MiniProject/Network_1_pc/15-1WinsorizationPoint")

## Get the data
namedf=read.csv("LADS_list_UK.csv",header=TRUE) # names of all the LADs in the UK, 380 total
namedf <- as.data.frame(namedf[,-1])
rownames(namedf) <- 0:(nrow(namedf)-1)

# read in A matrices
years=2005:2010
files=paste0("A_LAD_", years, "_UK_v5.mtx")
allf=lapply(files,readMM)

## DO SOME DATA PREPROCESSING STUFF

## Get the sum of all A matrices, called sumAm
sumA=allf[[1]]
for(i in 2:length(allf))sumA=sumA+allf[[i]]
n=dim(sumA)[1]
sumAm=as.matrix(sumA)


## Make an order for heatmaps using sumAm and hierarchical clustering
mydata=(log10(1+sumAm))
tdata=(mydata + t(mydata))/2
thc=hclust(dist(tdata))
thd=as.dendrogram(thc)
thd=reorder(thd,colSums(tdata))
to=labels(thd)
## Get a "nice" version to plot
plotdata=tdata[to,to]
rownames(plotdata)=colnames(plotdata)=namedf[to,]

SVD_sumAm = svd(mydata)

## Similarly, get a nice version of the individual matrices to plot
allfto=lapply(allf,function(x){
  tx=as.matrix(x[to,to]) # reorder the rows and columns to be in the order the HC gave
  rownames(tx)=colnames(tx)=namedf[to,]
  tx
})
alldata=do.call("cbind",allfto) # all years column binded together



### WINSORIZING (then logging) THEN CONCATENATING

AWinsLogConc <- vector(mode="list", length=length((years)))
for (i in 1:length(years)){
  A <- as.matrix(allfto[[i]])
  nrows = nrow(A)
  A_wins <- Winsorize(A, minval = 0, probs=c(0,0.99))
  A_wins <- as.matrix(A_wins, nrow=nrows)
  A <- log10(A_wins+1)
  AWinsLogConc[[i]] = A # list of all winsorised and logged A matrices
}

# for the AWinsLog data:
allWinsLogConc <- do.call("cbind", AWinsLogConc) # make one rectangular matrix
SVD_WinsLogConc <- svd(allWinsLogConc)


### (logging then) THEN CONCATENATING THEN WINSORIZING 

ALog <- vector(mode="list", length=length((years)))
for (i in 1:length(years)){
  A <- as.matrix(allfto[[i]])
  nrows = nrow(A)
  A <- as.matrix(A, nrow=nrows)
  A <- log10(A+1)
  ALog[[i]] = A # list of all logged A matrices
}

# concatenate the A matrices:
allLogConc <- do.call("cbind", ALog) # make one rectangular matrix

# Winsorize the concatenated matrix: 
allLogConcWins <- Winsorize(allLogConc, minval = 0, probs=c(0,0.99))
allLogConcWins <- as.matrix(allLogConcWins, nrow=nrows)

SVD_LogConcWins <- svd(allLogConcWins)


### Compare the 2 ways round of doing things (need to compare SVD_WinsLogConc and SVD_LogConcWins)

# X_hat = V * |Sigma|^(1/2)
# figuring out what to do
V_WC = SVD_WinsLogConc$v
Sigma_WC = diag(SVD_WinsLogConc$d)

V_CW = SVD_LogConcWins$v
Sigma_CW = diag(SVD_LogConcWins$d)

dim(V_WC)
dim(Sigma_WC)

distMeasFunc <- function(SVDlist, d){
  # takes in a list of the 3 components of an SVD: u,d,v and separates these bits out
  V = SVDlist$v
  Sigma_list = SVDlist$d
  U = SVDlist$u
  
  n = dim(U)[1] # this is n from the data input
  T = (dim(V)[1])/(dim(V)[2]) # number of timepoints
  V_d = V[ , 1:d]
  Sigma_list_d = Sigma_list[1:d]
  Sigma_d = diag(Sigma_list_d) # make a diag matrix Sigma
  X_hat = V_d %*% sqrtm(Sigma_d)
  print(dim(X_hat))
  
  # need to then make something that takes the X_hat calculated as the embedding from the SVD and distance measure
  r=0
  for (k in 1:n){
    t=0
    for (j in 1:(T-1)){
      s=0
      for (i in 1:d){
        s = s + (X_hat[(i+j+1),d] - X_hat[(i+j),d])^2
      }
      t = t + s
    }
    r = r + t
  }
  # r=0
  # for (k in 1:7){
  #   t=0
  #   for (j in 1:T){
  #     s=0
  #     for (i in 1:n){
  #       s = s + 1
  #     }
  #     t = t + s
  #   }
  #   r = r + t
  # }
  
  varX = sum(Sigma_list_d)
  val = r/(varX * n * (T-1))
  print(paste0("Without dividing by things at front: ", r))
  print(paste0("With dividing by things at front: ", val))

} # end of function

distMeasFunc(SVD_LogConcWins, 4)
distMeasFunc(SVD_WinsLogConc, 4)
# you do get different singular values depending on the order that you Wins and concat in 
