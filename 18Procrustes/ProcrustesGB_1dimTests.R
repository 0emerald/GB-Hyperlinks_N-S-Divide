library(tidyverse)
library(Matrix)
library(DescTools)
library(sf)
library(ggplot2)
library(grid)
library(plotly)
library(gridExtra)
library(ggrepel)

# first 140 lines are copied in from ProcrustesGB_X_avg_v2.R - some lines from that file are omitted as not needed
# removed the checks as they were for the debugging in the other script

# normalise so values go 0 to 1
normalise=function(x,trunc=0){
  x=as.matrix(x)
  if(all(!is.null(colnames(x)))) colnames(x)=paste0("norm",colnames(x))
  x = apply(x,2,function(x){x + abs(min(x))})
  x = apply(x,2,function(x){x/max(x)})
  x[x < trunc] <- 0
  x
}

# Dan has debugged 
setwd("~/Documents/MiniProject/Network_1_pc/18Procrustes")
# code from 16: matrixSVD script to get in the shp file

## Get the data
namedf=read.csv("LADS_list_GB.csv",header=TRUE) # names of all the LADs in GB, 380 total
colnames(namedf)[1] <- "index"
colnames(namedf)[2] <- "LAD"
namedf$index <- as.numeric(namedf$index)

# read in A matrices
allA= c()
for (year in 2005:2010){
  name = paste0("A_", year)
  filename = paste0("A_LAD_", year, "_GB.mtx")
  A = readMM(filename)
  A = as.matrix(A)
  eval(call("<-", as.name(name), A ))
  allA=cbind(allA, A) # makes the concatenated A
}

# Winsorize the concatenated A matrix
nrows = nrow(allA)
A_wins <- Winsorize(allA, minval = 0, probs=c(0,0.99))
A_wins <- as.matrix(A_wins, nrow=nrows)
# then log
A <- log10(A_wins+1)

getAvg=function(A){
  N=dim(A)[1]
  T=dim(A)[2]/dim(A)[1]
  avgA=matrix(0,N,N)
  for(t in 1:T)avgA=avgA+A[,(1:N)+(t-1)*N]/T
  avgA  
}
avgA=getAvg(A)
svdAvgA=svd(avgA)
svdAvgX=svdAvgA$v


### for the AConcWinsLog data:
SVD_ConcWinsLog <- svd(A)
avgSvdX=t(getAvg(t(SVD_ConcWinsLog$v)))
avgX=avgSvdX[,1:4] %*% (diag(SVD_ConcWinsLog$d)[1:4, 1:4])^0.5


### ADD COLOURS FOR THE LADs DEPENDING ON GEOGRAPHICAL POSITION
LAD11_shp <- read_sf('Local_Authority_Districts_December_2016_GCB_in_the_UK.shp')
LAD11_shp <- LAD11_shp %>%
  mutate(centroids = st_centroid(st_geometry(.))) # find centroids
LAD11_shp <-  cbind(LAD11_shp, st_coordinates(LAD11_shp$centroids)) # add X and Y for centroid
maxX = max(LAD11_shp$X)
maxY = max(LAD11_shp$Y)
minX = min(LAD11_shp$X)
minY = min(LAD11_shp$Y)
LAD11_shp$Xscaled = (LAD11_shp$X - minX)/(maxX-minX)# goes 0 to 1
LAD11_shp$Yscaled = (LAD11_shp$Y - minY)/(maxY-minY)

index = which(LAD11_shp$lad16nm == "Vale of Glamorgan")
# rename 'Vale of Glamorgan' in shp file to 'The Vale of Glamorgan' so that it matches
LAD11_shp[index, 3] <-  'The Vale of Glamorgan' 

# reorder rows of LAD11_shp to match the A matrices
LAD11_shp = merge(LAD11_shp, namedf, by.x="lad16nm", by.y="LAD")
LAD11_shp = LAD11_shp %>%
  arrange(index)


Gfull = data.frame(lad16nm=LAD11_shp$lad16nm,X=LAD11_shp$X,Y=LAD11_shp$Y,
                   Xscaled=LAD11_shp$Xscaled,Yscaled=LAD11_shp$Yscaled)
## scaled for (0,1)
## norm for mean centred with sd 1
Gfull$Xnorm=scale(Gfull$X)
Gfull$Ynorm=scale(Gfull$Y)

G=as.matrix(Gfull[,c("X","Y")])

op = data.frame(lad16nm=LAD11_shp$lad16nm,
                op0=rowMeans(avgA))
op$opRaw=op$op0/max(op$op0)-1e-5
op$op1=0.9*op$opRaw+0.1
op$op2=0.9*op$opRaw^2+0.1
op$op3=0.95*op$opRaw+0.05
op$op4=0.9*order(op$op0)/dim(op)[1]+0.1


### DOING PROCRUSTES TRANSFORMATION

# Try some iterative algorithm - for G, geography
G_tilde_t1 <- cbind(G, rep(0,dim(G)[1]), rep(0,dim(G)[1])) # \tilde{G}_0 = [G G_0 ] # have to set tilde_G_0 outside of for loop
P_tilde_t <- matrix(rep(0, 4^2), nrow =4,ncol=4) # have to decide some P_tilde_0 for it to start from in the distance measure
dist = 10 # so it will start
P_tilde_tm1 <- P_tilde_t
itermax=100
iter=1
while((dist > 1e-10)&(iter<itermax)){
  P_tilde_tm1 <- P_tilde_t #save P_tilde from the last iteration
  Proc <- MCMCpack::procrustes(avgX, G_tilde_t1, translation = TRUE, dilation = TRUE)
  P_tilde_t = Proc$R
  G_told= G_tilde_t1[,3:4]
  G_pred = Proc$X.new[,1:2]
  G_t1 = Proc$X.new[,3:4]
  G_tilde_t1 = cbind(G, G_t1)
  # dist between previous and current P_tilde
  distG <- norm((G-G_pred),type="F")
  dist <- norm((G_t1-G_told),type="F")
  #dist <- norm((P_tilde_t - P_tilde_tm1), type="F")
  print(paste("Iter =",iter,", Distance to G =",distG,", Distance to G_tilde =",dist))
  iter=iter+1
}
Xhat=Proc$s* (avgX %*% Proc$R) + matrix(1,dim(G)[1],1) %*% t(Proc$tt)
all(head(Proc$X.new) == head(Xhat))

X_tilde_full = Xhat[,1:2]
normalisedX_tilde_full = normalise(X_tilde_full,trunc=1/255)



####
# for count of active enterprises 2010
####
#had to look at animations.ipynb in 15 to see how to clean 
ent2010 <- read.csv("count_of_active_enterprises_2010_byLAD_GB_CLEAN.csv")
# remove commas in the numeric fields
ent2010$count.of.active.enterprises <- as.numeric(gsub(",","",ent2010$count.of.active.enterprises))
ent2010 = ent2010[ , 3:4] # keep count and LAD_clean
ent2010 = merge(ent2010, namedf, by.x="LAD_clean", by.y="LAD")
ent2010 = ent2010 %>% 
  arrange(index)


ent2010_tilde_t1 <- cbind(ent2010$count.of.active.enterprises, rep(0,dim(ent2010)[1]), rep(0,dim(ent2010)[1]),rep(0,dim(ent2010)[1])) # \tilde{G}_0 = [G G_0 ] # have to set tilde_G_0 outside of for loop
P_tilde_t <- matrix(rep(0, 4^2), nrow =4,ncol=4) # have to decide some P_tilde_0 for it to start from in the distance measure
dist = 10 # so it will start
P_tilde_tm1 <- P_tilde_t
itermax=100
iter=1
while((dist > 1e-10)&(iter<itermax)){
  P_tilde_tm1 <- P_tilde_t #save P_tilde from the last iteration
  Proc <- MCMCpack::procrustes(avgX, ent2010_tilde_t1, translation = TRUE, dilation = TRUE)
  P_tilde_t = Proc$R
  G_told= ent2010_tilde_t1[,2:4]
  ent2010_pred = Proc$X.new[,1]
  ent2010_t1 = Proc$X.new[,2:4]
  ent2010_tilde_t1 = cbind(ent2010$count.of.active.enterprises, ent2010_t1)
  # dist between previous and current P_tilde
  distG <- norm((as.matrix(ent2010$count.of.active.enterprises) - as.matrix(ent2010_pred)),type="F")
  dist <- norm((ent2010_t1-G_told),type="F")
  #dist <- norm((P_tilde_t - P_tilde_tm1), type="F")
  print(paste("Iter =",iter,", Distance to G =",distG,", Distance to G_tilde =",dist))
  iter=iter+1
}

truth = distG

Xhat=Proc$s* (avgX %*% Proc$R) + matrix(1,dim(ent2010)[1],1) %*% t(Proc$tt)
all(head(Proc$X.new) == head(Xhat))

X_tilde_full = Xhat[,1]
normalisedX_tilde_full = normalise(X_tilde_full,trunc=1/255)

# plot Procusting onto map
ggplot(LAD11_shp) + geom_sf(aes(geometry=centroids), 
                            fill = rgb(red=normalisedX_tilde_full[,1], green=0.5, 
                                       blue=0.5,op$op3), 
                            col=rgb(red=normalisedX_tilde_full[,1], green=0.5, 
                                    blue=0.5,op$op3)) + 
  ggtitle("GB LADs map coloured by \nProcrust PC1-4 into count of active enterprises 2010") 

# plot truth on map
ggplot(LAD11_shp) + geom_sf(aes(geometry=centroids), 
                                     fill = rgb(red=(ent2010$count.of.active.enterprises)/max(ent2010$count.of.active.enterprises), green=0.5, 
                                                blue=0.5, 0.3), 
                                     col=rgb(red=(ent2010$count.of.active.enterprises)/max(ent2010$count.of.active.enterprises), green=0.5, 
                                             blue=0.5, 0.3)) + 
  ggtitle("GB LADs map coloured by count of active enterprises 2010") 

# need to iteratively find the distance between the truth and procrustes, when we change where ent2010$count... is
permDists = c()
for (p in 1:1000){
  
  # sample randomly shuffles the values in the counts vector
  ent2010_tilde_t1 <- cbind(sample(ent2010$count.of.active.enterprises), rep(0,dim(ent2010)[1]), rep(0,dim(ent2010)[1]),rep(0,dim(ent2010)[1])) # \tilde{G}_0 = [G G_0 ] # have to set tilde_G_0 outside of for loop
  P_tilde_t <- matrix(rep(0, 4^2), nrow =4,ncol=4) # have to decide some P_tilde_0 for it to start from in the distance measure
  dist = 10 # so it will start
  P_tilde_tm1 <- P_tilde_t
  itermax=10
  iter=1
  while((dist > 1e-10)&(iter<itermax)){
    P_tilde_tm1 <- P_tilde_t #save P_tilde from the last iteration
    Proc <- MCMCpack::procrustes(avgX, ent2010_tilde_t1, translation = TRUE, dilation = TRUE)
    P_tilde_t = Proc$R
    G_told= ent2010_tilde_t1[,2:4]
    ent2010_pred = Proc$X.new[,1]
    ent2010_t1 = Proc$X.new[,2:4]
    ent2010_tilde_t1 = cbind(ent2010$count.of.active.enterprises, ent2010_t1)
    # dist between previous and current P_tilde
    distG <- norm((as.matrix(ent2010$count.of.active.enterprises) - as.matrix(ent2010_pred)),type="F")
    dist <- norm((ent2010_t1-G_told),type="F")
    #dist <- norm((P_tilde_t - P_tilde_tm1), type="F")
    #print(paste("Iter =",iter,", Distance to G =",distG,", Distance to G_tilde =",dist))
    iter=iter+1
    }
  permDists = c(permDists, distG)
}

length(which(permDists>truth))


####
# test it works
####
distG <- norm((as.matrix(ent2010$count.of.active.enterprises) - as.matrix(rep(0, dim(ent2010)[1]))),type="F")

truth = distG

permDists = c()
for (p in 1:1000){
  
  # sample randomly shuffles the values in the counts vector
  ent2010_tilde_t1 <- cbind(sample(rnorm(n=dim(ent2010)[1])), rep(0,dim(ent2010)[1]), rep(0,dim(ent2010)[1]),rep(0,dim(ent2010)[1])) # \tilde{G}_0 = [G G_0 ] # have to set tilde_G_0 outside of for loop
  P_tilde_t <- matrix(rep(0, 4^2), nrow =4,ncol=4) # have to decide some P_tilde_0 for it to start from in the distance measure
  dist = 10 # so it will start
  P_tilde_tm1 <- P_tilde_t
  itermax=10
  iter=1
  while((dist > 1e-10)&(iter<itermax)){
    P_tilde_tm1 <- P_tilde_t #save P_tilde from the last iteration
    Proc <- MCMCpack::procrustes(avgX, ent2010_tilde_t1, translation = TRUE, dilation = TRUE)
    P_tilde_t = Proc$R
    G_told= ent2010_tilde_t1[,2:4]
    ent2010_pred = Proc$X.new[,1]
    ent2010_t1 = Proc$X.new[,2:4]
    ent2010_tilde_t1 = cbind(ent2010$count.of.active.enterprises, ent2010_t1)
    # dist between previous and current P_tilde
    distG <- norm((as.matrix(ent2010$count.of.active.enterprises) - as.matrix(ent2010_pred)),type="F")
    dist <- norm((ent2010_t1-G_told),type="F")
    #dist <- norm((P_tilde_t - P_tilde_tm1), type="F")
    #print(paste("Iter =",iter,", Distance to G =",distG,", Distance to G_tilde =",dist))
    iter=iter+1
  }
  permDists = c(permDists, distG)
}

length(which(permDists>truth))
# should be around 500 if this works as intended


####
## median pay
####
# IS JUST NOT RENAMED THE VARIABLE, BUT ITS NOT FOR MEDIAN PAY
ent2010 <- read.csv("CleanedMedianMeanPay2010.csv")
ent2010 = ent2010[ , 3:5] # keep Median, mean, and LAD_clean
ent2010 = merge(ent2010, namedf, by.x="LAD_clean", by.y="LAD")
ent2010 = ent2010 %>% 
  arrange(index)
# if there is an x for the value, replaced it with 0, as we don't know the value
ent2010$Median <- as.numeric(gsub("x",0,ent2010$Median))
ent2010$Mean <- as.numeric(gsub("x",0,ent2010$Mean))


ent2010_tilde_t1 <- cbind(ent2010$Median, rep(0,dim(ent2010)[1]), rep(0,dim(ent2010)[1]),rep(0,dim(ent2010)[1])) # \tilde{G}_0 = [G G_0 ] # have to set tilde_G_0 outside of for loop
P_tilde_t <- matrix(rep(0, 4^2), nrow =4,ncol=4) # have to decide some P_tilde_0 for it to start from in the distance measure
dist = 10 # so it will start
P_tilde_tm1 <- P_tilde_t
itermax=100
iter=1
while((dist > 1e-10)&(iter<itermax)){
  P_tilde_tm1 <- P_tilde_t #save P_tilde from the last iteration
  Proc <- MCMCpack::procrustes(avgX, ent2010_tilde_t1, translation = TRUE, dilation = TRUE)
  P_tilde_t = Proc$R
  G_told= ent2010_tilde_t1[,2:4]
  ent2010_pred = Proc$X.new[,1]
  ent2010_t1 = Proc$X.new[,2:4]
  ent2010_tilde_t1 = cbind(ent2010$count.of.active.enterprises, ent2010_t1)
  # dist between previous and current P_tilde
  distG <- norm((as.matrix(ent2010$Median) - as.matrix(ent2010_pred)),type="F")
  dist <- norm((ent2010_t1-G_told),type="F")
  #dist <- norm((P_tilde_t - P_tilde_tm1), type="F")
  print(paste("Iter =",iter,", Distance to G =",distG,", Distance to G_tilde =",dist))
  iter=iter+1
}

truth = distG

Xhat=Proc$s* (avgX %*% Proc$R) + matrix(1,dim(ent2010)[1],1) %*% t(Proc$tt)
all(head(Proc$X.new) == head(Xhat))
X_tilde_full = Xhat[,1]
normalisedX_tilde_full = normalise(X_tilde_full,trunc=1/255)

# need to iteratively find the distance between the truth and procrustes, when we change where ent2010$count... is
permDists = c()
for (p in 1:1000){
  
  # sample randomly shuffles the values in the counts vector
  ent2010_tilde_t1 <- cbind(sample(ent2010$Median), rep(0,dim(ent2010)[1]), rep(0,dim(ent2010)[1]),rep(0,dim(ent2010)[1])) # \tilde{G}_0 = [G G_0 ] # have to set tilde_G_0 outside of for loop
  P_tilde_t <- matrix(rep(0, 4^2), nrow =4,ncol=4) # have to decide some P_tilde_0 for it to start from in the distance measure
  dist = 10 # so it will start
  P_tilde_tm1 <- P_tilde_t
  itermax=10
  iter=1
  while((dist > 1e-10)&(iter<itermax)){
    P_tilde_tm1 <- P_tilde_t #save P_tilde from the last iteration
    Proc <- MCMCpack::procrustes(avgX, ent2010_tilde_t1, translation = TRUE, dilation = TRUE)
    P_tilde_t = Proc$R
    G_told= ent2010_tilde_t1[,2:4]
    ent2010_pred = Proc$X.new[,1]
    ent2010_t1 = Proc$X.new[,2:4]
    ent2010_tilde_t1 = cbind(ent2010$count.of.active.enterprises, ent2010_t1)
    # dist between previous and current P_tilde
    distG <- norm((as.matrix(ent2010$Median) - as.matrix(ent2010_pred)),type="F")
    dist <- norm((ent2010_t1-G_told),type="F")
    #dist <- norm((P_tilde_t - P_tilde_tm1), type="F")
    #print(paste("Iter =",iter,", Distance to G =",distG,", Distance to G_tilde =",dist))
    iter=iter+1
  }
  permDists = c(permDists, distG)
}
length(which(permDists>truth))



####
# geog - lat
###

ent2010_tilde_t1 <- cbind(G[,1], rep(0,dim(ent2010)[1]), rep(0,dim(ent2010)[1]),rep(0,dim(ent2010)[1])) # \tilde{G}_0 = [G G_0 ] # have to set tilde_G_0 outside of for loop
P_tilde_t <- matrix(rep(0, 4^2), nrow =4,ncol=4) # have to decide some P_tilde_0 for it to start from in the distance measure
dist = 10 # so it will start
P_tilde_tm1 <- P_tilde_t
itermax=100
iter=1
while((dist > 1e-10)&(iter<itermax)){
  P_tilde_tm1 <- P_tilde_t #save P_tilde from the last iteration
  Proc <- MCMCpack::procrustes(avgX, ent2010_tilde_t1, translation = TRUE, dilation = TRUE)
  P_tilde_t = Proc$R
  G_told= ent2010_tilde_t1[,2:4]
  ent2010_pred = Proc$X.new[,1]
  ent2010_t1 = Proc$X.new[,2:4]
  ent2010_tilde_t1 = cbind(ent2010$count.of.active.enterprises, ent2010_t1)
  # dist between previous and current P_tilde
  distG <- norm((as.matrix(G[,1]) - as.matrix(ent2010_pred)),type="F")
  dist <- norm((ent2010_t1-G_told),type="F")
  #dist <- norm((P_tilde_t - P_tilde_tm1), type="F")
  print(paste("Iter =",iter,", Distance to G =",distG,", Distance to G_tilde =",dist))
  iter=iter+1
}

truth = distG

Xhat=Proc$s* (avgX %*% Proc$R) + matrix(1,dim(ent2010)[1],1) %*% t(Proc$tt)
all(head(Proc$X.new) == head(Xhat))
X_tilde_full = Xhat[,1]
normalisedX_tilde_full = normalise(X_tilde_full,trunc=1/255)

# need to iteratively find the distance between the truth and procrustes, when we change where ent2010$count... is
permDists = c()
for (p in 1:1000){
  
  # sample randomly shuffles the values in the counts vector
  ent2010_tilde_t1 <- cbind(sample(G[,1]), rep(0,dim(ent2010)[1]), rep(0,dim(ent2010)[1]),rep(0,dim(ent2010)[1])) # \tilde{G}_0 = [G G_0 ] # have to set tilde_G_0 outside of for loop
  P_tilde_t <- matrix(rep(0, 4^2), nrow =4,ncol=4) # have to decide some P_tilde_0 for it to start from in the distance measure
  dist = 10 # so it will start
  P_tilde_tm1 <- P_tilde_t
  itermax=10
  iter=1
  while((dist > 1e-10)&(iter<itermax)){
    P_tilde_tm1 <- P_tilde_t #save P_tilde from the last iteration
    Proc <- MCMCpack::procrustes(avgX, ent2010_tilde_t1, translation = TRUE, dilation = TRUE)
    P_tilde_t = Proc$R
    G_told= ent2010_tilde_t1[,2:4]
    ent2010_pred = Proc$X.new[,1]
    ent2010_t1 = Proc$X.new[,2:4]
    ent2010_tilde_t1 = cbind(ent2010$count.of.active.enterprises, ent2010_t1)
    # dist between previous and current P_tilde
    distG <- norm((as.matrix(G[,1]) - as.matrix(ent2010_pred)),type="F")
    dist <- norm((ent2010_t1-G_told),type="F")
    #dist <- norm((P_tilde_t - P_tilde_tm1), type="F")
    #print(paste("Iter =",iter,", Distance to G =",distG,", Distance to G_tilde =",dist))
    iter=iter+1
  }
  permDists = c(permDists, distG)
}
length(which(permDists>truth))

####
# geog - long
###

ent2010_tilde_t1 <- cbind(G[,2], rep(0,dim(ent2010)[1]), rep(0,dim(ent2010)[1]),rep(0,dim(ent2010)[1])) # \tilde{G}_0 = [G G_0 ] # have to set tilde_G_0 outside of for loop
P_tilde_t <- matrix(rep(0, 4^2), nrow =4,ncol=4) # have to decide some P_tilde_0 for it to start from in the distance measure
dist = 10 # so it will start
P_tilde_tm1 <- P_tilde_t
itermax=100
iter=1
while((dist > 1e-10)&(iter<itermax)){
  P_tilde_tm1 <- P_tilde_t #save P_tilde from the last iteration
  Proc <- MCMCpack::procrustes(avgX, ent2010_tilde_t1, translation = TRUE, dilation = TRUE)
  P_tilde_t = Proc$R
  G_told= ent2010_tilde_t1[,2:4]
  ent2010_pred = Proc$X.new[,1]
  ent2010_t1 = Proc$X.new[,2:4]
  ent2010_tilde_t1 = cbind(ent2010$count.of.active.enterprises, ent2010_t1)
  # dist between previous and current P_tilde
  distG <- norm((as.matrix(G[,2]) - as.matrix(ent2010_pred)),type="F")
  dist <- norm((ent2010_t1-G_told),type="F")
  #dist <- norm((P_tilde_t - P_tilde_tm1), type="F")
  print(paste("Iter =",iter,", Distance to G =",distG,", Distance to G_tilde =",dist))
  iter=iter+1
}

truth = distG

Xhat=Proc$s* (avgX %*% Proc$R) + matrix(1,dim(ent2010)[1],1) %*% t(Proc$tt)
all(head(Proc$X.new) == head(Xhat))
X_tilde_full = Xhat[,1]
normalisedX_tilde_full = normalise(X_tilde_full,trunc=1/255)

# need to iteratively find the distance between the truth and procrustes, when we change where ent2010$count... is
permDists = c()
for (p in 1:1000){
  
  # sample randomly shuffles the values in the counts vector
  ent2010_tilde_t1 <- cbind(sample(G[,2]), rep(0,dim(ent2010)[1]), rep(0,dim(ent2010)[1]),rep(0,dim(ent2010)[1])) # \tilde{G}_0 = [G G_0 ] # have to set tilde_G_0 outside of for loop
  P_tilde_t <- matrix(rep(0, 4^2), nrow =4,ncol=4) # have to decide some P_tilde_0 for it to start from in the distance measure
  dist = 10 # so it will start
  P_tilde_tm1 <- P_tilde_t
  itermax=10
  iter=1
  while((dist > 1e-10)&(iter<itermax)){
    P_tilde_tm1 <- P_tilde_t #save P_tilde from the last iteration
    Proc <- MCMCpack::procrustes(avgX, ent2010_tilde_t1, translation = TRUE, dilation = TRUE)
    P_tilde_t = Proc$R
    G_told= ent2010_tilde_t1[,2:4]
    ent2010_pred = Proc$X.new[,1]
    ent2010_t1 = Proc$X.new[,2:4]
    ent2010_tilde_t1 = cbind(ent2010$count.of.active.enterprises, ent2010_t1)
    # dist between previous and current P_tilde
    distG <- norm((as.matrix(G[,2]) - as.matrix(ent2010_pred)),type="F")
    dist <- norm((ent2010_t1-G_told),type="F")
    #dist <- norm((P_tilde_t - P_tilde_tm1), type="F")
    #print(paste("Iter =",iter,", Distance to G =",distG,", Distance to G_tilde =",dist))
    iter=iter+1
  }
  permDists = c(permDists, distG)
}
length(which(permDists>truth))
