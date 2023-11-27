library(tidyverse)
library(Matrix)
library(DescTools)
library(sf)
library(ggplot2)
library(grid)
library(plotly)
library(gridExtra)
library(ggrepel)

# THE MAJORITY OF THIS CODE COMES FROM ProcrustesGB_X_avg_v2.R

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

op = data.frame(lad16nm=LAD11_shp$lad16nm,
                op0=rowMeans(avgA))
op$opRaw=op$op0/max(op$op0)-1e-5
op$op1=0.9*op$opRaw+0.1
op$op2=0.9*op$opRaw^2+0.1
op$op3=0.95*op$opRaw+0.05
op$op4=0.9*order(op$op0)/dim(op)[1]+0.1
toname=order(op$op0,decreasing=T)[1:20]

# for count of active enterprises 2010
ent2010 <- read.csv("count_of_active_enterprises_2010_byLAD_GB_CLEAN.csv")
# remove commas in the numeric fields
ent2010$count.of.active.enterprises <- as.numeric(gsub(",","",ent2010$count.of.active.enterprises))
ent2010 = ent2010[ , 3:4] # keep count and LAD_clean
ent2010 = merge(ent2010, namedf, by.x="LAD_clean", by.y="LAD")
ent2010 = ent2010 %>% 
  arrange(index)

# G has 3 columns: X, Y, #ents
G=as.matrix(Gfull[,c("X","Y")])
G = as.matrix(cbind(G, ent2010$count.of.active.enterprises))

G_tilde_t1 <- cbind(G, rep(0,dim(G)[1])) # \tilde{G}_0 = [G G_0 ] # have to set tilde_G_0 outside of for loop
P_tilde_t <- matrix(rep(0, 4^2), nrow =4,ncol=4) # have to decide some P_tilde_0 for it to start from in the distance measure
dist = 10 # so it will start
P_tilde_tm1 <- P_tilde_t
itermax=100
iter=1
while((dist > 1e-10)&(iter<itermax)){
  P_tilde_tm1 <- P_tilde_t #save P_tilde from the last iteration
  Proc <- MCMCpack::procrustes(avgX, G_tilde_t1, translation = TRUE, dilation = TRUE)
  P_tilde_t = Proc$R
  G_told= as.matrix(G_tilde_t1[,4])
  G_pred = Proc$X.new[,1:3]
  G_t1 = as.matrix(Proc$X.new[,4])
  G_tilde_t1 = cbind(G, G_t1)
  # dist between previous and current P_tilde
  distG <- norm((G-G_pred),type="F")
  dist <- norm((G_t1-G_told),type="F")
  #dist <- norm((P_tilde_t - P_tilde_tm1), type="F")
  print(paste("Iter =",iter,", Distance to G =",distG,", Distance to G_tilde =",dist))
  iter=iter+1
}
Xhat=Proc$s * (avgX %*% Proc$R) + matrix(1,dim(G)[1],1) %*% t(Proc$tt)
all(head(Proc$X.new) == head(Xhat)) # checks that Xhat and X.new from the algorithm are the same, as you'd want them to


###############################
### look at early/late divide
###############################
# MCMCpack says: X_{new} = s X R + 1 tt'
# if we use the s, R, tt, found using Xavg, we can then replace the 'X' with Xearly and Xlate, to get Xearlyhat, Xlatehat

# make early and late averages of X (2005-2007 and 2008-2010)
N = dim(A_2005)[1]
# early
avgSvdXearly=t(getAvg(t(SVD_ConcWinsLog$v[1:(3*N),]))) # finds average X for the 3 early years
avgXearly=avgSvdXearly[,1:4] %*% (diag(SVD_ConcWinsLog$d)[1:4, 1:4])^0.5

# late
avgSvdXlate=t(getAvg(t(SVD_ConcWinsLog$v[(3*N +1):(6*N),])))
avgXlate=avgSvdXlate[,1:4] %*% (diag(SVD_ConcWinsLog$d)[1:4, 1:4])^0.5


### THIS IS METHOD (a) ABOUT HOW TO FIND THE EARLY AND LATE PROCRUSTED EMBEDDINGS!!
# This is using the P matrix found with the avgX embedding to calculate a Procrustes transform for Xearly and Xlate 
Xhat_early = Proc$s* (avgXearly %*% Proc$R) + matrix(1,dim(G)[1],1) %*% t(Proc$tt)
Xhat_late = Proc$s* (avgXlate %*% Proc$R) + matrix(1,dim(G)[1],1) %*% t(Proc$tt)

# this only keeps the first 2 columns, but we probably want to keep all 4 
X_tilde_full_early = Xhat_early[,1:2]
X_tilde_full_late = Xhat_late[,1:2]
normalisedX_tilde_full_early = normalise(X_tilde_full_early,trunc=1/255)
normalisedX_tilde_full_late = normalise(X_tilde_full_late,trunc=1/255)

normalisedXhat_early = normalise(Xhat_early, trunc=1/255)
normalisedXhat_late = normalise(Xhat_late, trunc=1/255)

write.csv(normalisedXhat_early, "normalisedXhat_early_methodA.csv")
write.csv(normalisedXhat_late, "normalisedXhat_late_methodA.csv")

