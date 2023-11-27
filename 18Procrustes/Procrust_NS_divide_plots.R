library(tidyverse)
library(Matrix)
library(DescTools)
library(sf)
library(ggplot2)
library(grid)
library(plotly)
library(gridExtra)
library(ggrepel)

# this is to make a plot for the paper that is wrt the procrusted dimensions rather than the embedding dimensions of 
# how the N/S divide is changing over time

# normalise so values go 0 to 1
normalise=function(x,trunc=0){
  x=as.matrix(x)
  if(all(!is.null(colnames(x)))) colnames(x)=paste0("norm",colnames(x))
  x = apply(x,2,function(x){x + abs(min(x))})
  x = apply(x,2,function(x){x/max(x)})
  x[x < trunc] <- 0
  x
}

setwd("~/Documents/MiniProject/Network_1_pc/18Procrustes")

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

### for the AConcWinsLog data:
SVD_ConcWinsLog <- svd(A)
avgSvdX=t(getAvg(t(SVD_ConcWinsLog$v)))
avgX=avgSvdX[,1:4] %*% (diag(SVD_ConcWinsLog$d)[1:4, 1:4])^0.5

# embedding of AConcWinsLog
Yhat_all = SVD_ConcWinsLog$v[ , 1:4] %*% (diag(SVD_ConcWinsLog$d)[1:4, 1:4])^0.5

# Find the s, r, tt, from Procrusting into G
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

# this is the whole time averaged embedding procrusted into G
Xhat=Proc$s * (avgX %*% Proc$R) + matrix(1,dim(G)[1],1) %*% t(Proc$tt)

# Procrust each Yhat into G

n = dim(SVD_ConcWinsLog$u)[1]
allYhatProcrusted = c() # save as one thing

for (year in 2005:2010){
  index = year - 2004
  Yhat = Yhat_all[((index-1)*n +1 ):((index)*n) , ]
  Yhat_Procrusted = Proc$s * (Yhat %*% Proc$R) + matrix(1,dim(G)[1],1) %*% t(Proc$tt)
  filename = paste0("ProcrustedEmbedding_", year, ".csv")
  write.csv(Yhat, file=filename, row.names=FALSE)
  allYhatProcrusted=rbind(allYhatProcrusted, Yhat_Procrusted) # makes the concatenated YhatProc
}

# save to a file
write.csv(allYhatProcrusted, "YhatProcrusted_allYears.csv", row.names = FALSE)

