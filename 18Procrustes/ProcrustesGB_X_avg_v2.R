library(tidyverse)
library(Matrix)
library(DescTools)
library(sf)
library(ggplot2)
library(grid)
library(plotly)
library(gridExtra)
library(ggrepel)


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

## Check
plot(svdAvgX[,1],svdAvgX[,2],col=rgb(G[,1],G[,2],1,op$op1),pch=19)
text(svdAvgX[,1],svdAvgX[,2],namedf[,2],cex=0.5)

plot(svdAvgX[,3],svdAvgX[,4],col=rgb(G[,1],G[,2],1,op$op1),pch=19)
text(svdAvgX[,3],svdAvgX[,4],namedf[,2],cex=0.5)

### for the AConcWinsLog data:
SVD_ConcWinsLog <- svd(A)
avgSvdX=t(getAvg(t(SVD_ConcWinsLog$v)))
avgX=avgSvdX[,1:4] %*% (diag(SVD_ConcWinsLog$d)[1:4, 1:4])^0.5

## Check
plot(avgSvdX[,1],avgSvdX[,2],col=rgb(G[,1],G[,2],1,op$op1),pch=19)
text(avgSvdX[,1],avgSvdX[,2],namedf[,2],cex=0.5)

plot(avgSvdX[,3],avgSvdX[,4],col=rgb(G[,1],G[,2],1,op$op1),pch=19)
text(avgSvdX[,3],avgSvdX[,4],namedf[,2],cex=0.5)


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
########
## Check
plot(avgSvdX[,3],avgSvdX[,4],col=rgb(Gfull$Xscaled,Gfull$Yscaled,1,op$op1),pch=19)
text(avgSvdX[,3],avgSvdX[,4],namedf[,2],cex=0.5)

## Simplest rotation: PC3-4
G=as.matrix(Gfull[,c("X","Y")])
Proc <- MCMCpack::procrustes(avgSvdX[,3:4], G, translation = TRUE, dilation = TRUE)
X_tilde_PC34 = Proc$X.new

plot(X_tilde_PC34[,1],X_tilde_PC34[,2],col=rgb(Gfull$Xscaled,Gfull$Yscaled,1,op$op1),pch=19)
text(X_tilde_PC34[,1],X_tilde_PC34[,2],namedf[,2],cex=0.5,main="Procrust PC3-4 onto G")

## Should be worse: PC1-2
Proc <- MCMCpack::procrustes(avgSvdX[,1:2], G, translation = TRUE, dilation = TRUE)
X_tilde_PC12 = Proc$X.new

plot(X_tilde_PC12[,1],X_tilde_PC12[,2],col=rgb(Gfull$Xscaled,Gfull$Yscaled,1,op$op1),pch=19)
text(X_tilde_PC12[,1],X_tilde_PC12[,2],namedf[,2],cex=0.5,main="Procrust PC1-2 onto G")

## Iterative approach using PC1-4
# Try some iterative algorithm
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
Xhat=Proc$s * (avgX %*% Proc$R) + matrix(1,dim(G)[1],1) %*% t(Proc$tt)
all(head(Proc$X.new) == head(Xhat)) # checks that Xhat and X.new from the algorithm are the same, as you'd want them to


 
# ########
# ### ANIMATION STUFF 
# animationData = cbind(SVD_ConcWinsLog$v[ , 1:4] , rep(c(2005:2010), each=n))
# animationData = cbind(animationData, rep(LAD11_shp$Xscaled, 6), rep(LAD11_shp$Yscaled, 6), rep(namedf$LAD, 6))
# animationData = as.data.frame(animationData)
# colnames(animationData) = c("dim1", "dim2", "dim3", "dim4", "year", "Xscaled", "Yscaled", "LAD")
# animationData$col = rgb(red=animationData$Xscaled, blue=animationData$Yscaled, green = 0.2 )
# 
# #### 13/12/2022 - NEW STUFF FOR PROCRUSTES TRANSFORMATION
# animationData[, c("dim1", "dim2", "dim3", "dim4","Xscaled", "Yscaled")] <- sapply(animationData[, c("dim1", "dim2", "dim3", "dim4","Xscaled", "Yscaled")], as.numeric)
G <- as.matrix(Gfull[ , c("Xscaled","Yscaled")])

scores=data.frame(
  # full=norm(G-X_tilde,type="F"),
  full=norm(G-X_tilde_full,type="F"),
  PC34=norm(G-X_tilde_PC34,type="F"),
  PC12=norm(G-X_tilde_PC12,type="F"))
scores

### DOING PROCRUSTES TRANSFORMATION

# Try some iterative algorithm
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
normalisedX_tilde_PC12 = normalise(X_tilde_PC12,trunc=1/255)
normalisedX_tilde_PC34 = normalise(X_tilde_PC34,trunc=1/255)

op = data.frame(lad16nm=LAD11_shp$lad16nm,
                op0=rowMeans(avgA))
op$opRaw=op$op0/max(op$op0)-1e-5
op$op1=0.9*op$opRaw+0.1
op$op2=0.9*op$opRaw^2+0.1
op$op3=0.95*op$opRaw+0.05
op$op4=0.9*order(op$op0)/dim(op)[1]+0.1


pdf("maps.pdf",height=6,width=12)
par(mfrow=c(1,3))
plot(G[,1],G[,2],col=rgb(normalisedX_tilde_PC12[,1],normalisedX_tilde_PC12[,2],
                         1-normalisedX_tilde_PC12[,1]/2-normalisedX_tilde_PC12[,2]/2,op$op3),
     pch=19,main="G coloured by Procrust PC1-2 onto G")
plot(G[,1],G[,2],col=rgb(normalisedX_tilde_PC34[,1],normalisedX_tilde_PC34[,2],
                         1-normalisedX_tilde_PC34[,1]/2-normalisedX_tilde_PC34[,2]/2,op$op3),
     pch=19,main="G coloured by Procrust PC3-4 onto G")
plot(G[,1],G[,2],col=rgb(normalisedX_tilde_full[,1],normalisedX_tilde_full[,2],
                         1-normalisedX_tilde_full[,1]/2-normalisedX_tilde_full[,2]/2,op$op3),
     pch=19,main="G coloured by Procrust PC1-4 into G")
dev.off()


toname=order(op$op0,decreasing=T)[1:20]

pdf("maps2.pdf",height=6,width=12)
par(mfrow=c(1,3))
plot(X_tilde_PC12[,1],X_tilde_PC12[,2],col=rgb(Gfull$Xscaled,Gfull$Yscaled,
                                               1-normalisedX_tilde_PC12[,1]/2-normalisedX_tilde_PC12[,2]/2,op$op3),
     pch=19,main="Procrust PC1-2 onto G coloured by G",sub=paste("Score",scores["PC12"]))
text(X_tilde_PC12[toname,1],X_tilde_PC12[toname,2],namedf[toname,2],cex=0.5)
plot(X_tilde_PC34[,1],X_tilde_PC34[,2],col=rgb(Gfull$Xscaled,Gfull$Yscaled,
                                               1-normalisedX_tilde_PC34[,1]/2-normalisedX_tilde_PC34[,2]/2,op$op3),
     pch=19,main="Procrust PC3-4 onto G coloured by G",sub=paste("Score",scores["PC34"]))
text(X_tilde_PC34[toname,1],X_tilde_PC34[toname,2],namedf[toname,2],cex=0.5)
plot(X_tilde_full[,1],X_tilde_full[,2],col=rgb(Gfull$Xscaled,Gfull$Yscaled,
                                               1-normalisedX_tilde_full[,1]/2-normalisedX_tilde_full[,2]/2,op$op3),
     pch=19,main="Procrust PC1-4 onto G coloured by G",sub=paste("Score",scores["full"]))
text(X_tilde_full[toname,1],X_tilde_full[toname,2],namedf[toname,2],cex=0.5)
dev.off()

# just the third plot from maps2
plot(X_tilde_full[,1],X_tilde_full[,2],col=rgb(Gfull$Xscaled,Gfull$Yscaled,
                                               1-normalisedX_tilde_full[,1]/2-normalisedX_tilde_full[,2]/2,op$op3),
     pch=19,
      main=(bquote(atop("Plot of "~hat(X) ~" obtained by Procrusting", 
      "the time-averaged embedding," ~tilde(X) ~", into geography," ~G))),
     #sub=paste("Score",scores["full"]),
     xlab = "Dimension 1", ylab="Dimension 2"
     )
text(X_tilde_full[toname,1],X_tilde_full[toname,2],namedf[toname,2],cex=0.5)



### TRY AND MAKE THE map.pdf PLOT LOOK LIKE GEOGRAPHY
plot(G[,1],G[,2],col=rgb(normalisedX_tilde_full[,1],normalisedX_tilde_full[,2],
                                    1-normalisedX_tilde_full[,1]/2-normalisedX_tilde_full[,2]/2, op$op3),
                pch=19,main="G coloured by Procrust PC1-4 into G")

plotPrc1_4toG <- ggplot(LAD11_shp) + geom_sf(aes(geometry=centroids), 
                                         fill = rgb(red=normalisedX_tilde_full[,1], green=normalisedX_tilde_full[,2], 
                                                    blue=1-normalisedX_tilde_full[,1]/2-normalisedX_tilde_full[,2]/2,op$op3), 
                                         col=rgb(red=normalisedX_tilde_full[,1], green=normalisedX_tilde_full[,2], 
                                                 blue=1-normalisedX_tilde_full[,1]/2-normalisedX_tilde_full[,2]/2,op$op3)) + 
  ggtitle("GB LADs map coloured by \nProcrust PC1-4 into G") 

plotPrc3_4toG <- ggplot(LAD11_shp) + geom_sf(aes(geometry=centroids), 
                                             fill = rgb(red=normalisedX_tilde_PC34[,1], green=normalisedX_tilde_PC34[,2], 
                                                        blue=1-normalisedX_tilde_PC34[,1]/2-normalisedX_tilde_PC34[,2]/2,op$op3), 
                                             col=rgb(red=normalisedX_tilde_PC34[,1], green=normalisedX_tilde_PC34[,2], 
                                                     blue=1-normalisedX_tilde_PC34[,1]/2-normalisedX_tilde_PC34[,2]/2,op$op3)) + 
  ggtitle("GB LADs map coloured by Procrust PC3-4 into G") 

plotG <- ggplot(LAD11_shp) + geom_sf(aes(geometry=centroids), 
                            fill = rgb(red=G[,1], green=G[,2], 
                                       blue=1-G[,1]/2-G[,2]/2, 0.3), 
                            col=rgb(red=G[,1], green=G[,2], 
                                    blue=1-G[,1]/2-G[,2]/2, 0.3)) + 
  ggtitle("GB LADs map coloured by G") 

grid.arrange(plotPrc3_4toG, plotPrc1_4toG,plotG,ncol=3)


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

# This is using the P matrix found with the avgX embedding to calculate a Procrustes transform for Xearly and Xlate 
Xhat_early = Proc$s* (avgXearly %*% Proc$R) + matrix(1,dim(G)[1],1) %*% t(Proc$tt)
Xhat_late = Proc$s* (avgXlate %*% Proc$R) + matrix(1,dim(G)[1],1) %*% t(Proc$tt)

X_tilde_full_early = Xhat_early[,1:2]
X_tilde_full_late = Xhat_late[,1:2]
normalisedX_tilde_full_early = normalise(X_tilde_full_early,trunc=1/255)
normalisedX_tilde_full_late = normalise(X_tilde_full_late,trunc=1/255)

## 24/04/2023
# saving the normalised Procrustes transformations for X into G = (long,lat) for early and late
write.csv(normalisedX_tilde_full_early, "normalisedX_tilde_full_early.csv")
write.csv(normalisedX_tilde_full_late, "normalisedX_tilde_full_late.csv")

### plots early
par(mfrow=c(1,2))
plot(G[,1],G[,2],col=rgb(normalisedX_tilde_full_early[,1],normalisedX_tilde_full_early[,2],
                         1-normalisedX_tilde_full_early[,1]/2-normalisedX_tilde_full_early[,2]/2,op$op3),
     pch=19,main="G coloured by Procrust\n PC1-4 for X_early into G")

plot(X_tilde_full_early[,1],X_tilde_full_early[,2],col=rgb(Gfull$Xscaled,Gfull$Yscaled,
                                               1-normalisedX_tilde_full_early[,1]/2-normalisedX_tilde_full_early[,2]/2,op$op3),
     pch=19,main="Procrust PC1-4 for X_early onto G coloured by G",sub=paste("Score",scores["full"]))
text(X_tilde_full_early[toname,1],X_tilde_full_early[toname,2],namedf[toname,2],cex=0.5)

plotPrc1_4early_toG <- ggplot(LAD11_shp) + geom_sf(aes(geometry=centroids), 
                                             fill = rgb(red=normalisedX_tilde_full_early[,1], green=normalisedX_tilde_full_early[,2], 
                                                        blue=1-normalisedX_tilde_full_early[,1]/2-normalisedX_tilde_full_early[,2]/2,op$op3), 
                                             col=rgb(red=normalisedX_tilde_full_early[,1], green=normalisedX_tilde_full_early[,2], 
                                                     blue=1-normalisedX_tilde_full_early[,1]/2-normalisedX_tilde_full_early[,2]/2,op$op3)) + 
  ggtitle("GB LADs map coloured by Procrust\n PC1-4 of X early into G") 

grid.arrange(plotPrc1_4early_toG,plotG,ncol=2)


### plots late
par(mfrow=c(1,2))
plot(G[,1],G[,2],col=rgb(normalisedX_tilde_full_late[,1],normalisedX_tilde_full_late[,2],
                         1-normalisedX_tilde_full_late[,1]/2-normalisedX_tilde_full_late[,2]/2,op$op3),
     pch=19,main="G coloured by Procrust \nPC1-4 for X_late into G")

plot(X_tilde_full_late[,1],X_tilde_full_late[,2],col=rgb(Gfull$Xscaled,Gfull$Yscaled,
                                                           1-normalisedX_tilde_full_late[,1]/2-normalisedX_tilde_full_late[,2]/2,op$op3),
     pch=19,main="Procrust PC1-4 for X_late onto G coloured by G",sub=paste("Score",scores["full"]))
text(X_tilde_full_late[toname,1],X_tilde_full_late[toname,2],namedf[toname,2],cex=0.5)

plotPrc1_4late_toG <- ggplot(LAD11_shp) + geom_sf(aes(geometry=centroids), 
                                                   fill = rgb(red=normalisedX_tilde_full_late[,1], green=normalisedX_tilde_full_late[,2], 
                                                              blue=1-normalisedX_tilde_full_late[,1]/2-normalisedX_tilde_full_late[,2]/2,op$op3), 
                                                   col=rgb(red=normalisedX_tilde_full_late[,1], green=normalisedX_tilde_full_late[,2], 
                                                           blue=1-normalisedX_tilde_full_late[,1]/2-normalisedX_tilde_full_late[,2]/2,op$op3)) + 
  ggtitle("GB LADs map coloured by Procrust\n PC1-4 of X late into G") 

grid.arrange(plotPrc1_4late_toG,plotG,ncol=2)

# plot all side by side
grid.arrange(plotPrc1_4early_toG, plotPrc1_4late_toG, plotPrc1_4toG, plotG, ncol=4)

# plot X-avg and reference map side by side, where titles renamed
plotPrc1_4toG_titles <- ggplot(LAD11_shp) + geom_sf(aes(geometry=centroids), 
                                             fill = rgb(red=normalisedX_tilde_full[,1], green=normalisedX_tilde_full[,2], 
                                                        blue=1-normalisedX_tilde_full[,1]/2-normalisedX_tilde_full[,2]/2,op$op3), 
                                             col=rgb(red=normalisedX_tilde_full[,1], green=normalisedX_tilde_full[,2], 
                                                     blue=1-normalisedX_tilde_full[,1]/2-normalisedX_tilde_full[,2]/2,op$op3)) + 
  ggtitle("Local Authorities in Great Britain coloured\n by Procrusting dimensions 1 to 4 of the time\naveraged embedding into geographic location") 
plotG_titles <- ggplot(LAD11_shp) + geom_sf(aes(geometry=centroids), 
                                     fill = rgb(red=G[,1], green=G[,2], 
                                                blue=1-G[,1]/2-G[,2]/2, 0.3), 
                                     col=rgb(red=G[,1], green=G[,2], 
                                             blue=1-G[,1]/2-G[,2]/2, 0.3)) + 
  ggtitle("Local Authorities in Great Britain\ncoloured by geographic location") 


# 05/01/2023 - this little chunk
# plots of the maps2 third plot style, with early and late side by side
# i.e. what X_early looks like when Procrusted into G, and X_late when Procrusted into G
par(mfrow=c(1,2))
plot(X_tilde_full_early[,1],X_tilde_full_early[,2],col=rgb(Gfull$Xscaled,Gfull$Yscaled,
                                                           1-normalisedX_tilde_full_early[,1]/2-normalisedX_tilde_full_early[,2]/2,op$op3),
     pch=19,main="Procrust PC1-4 for X_early onto G coloured by G")#,sub=paste("Score",scores["full"]))
text(X_tilde_full_early[toname,1],X_tilde_full_early[toname,2],namedf[toname,2],cex=0.5)

plot(X_tilde_full_late[,1],X_tilde_full_late[,2],col=rgb(Gfull$Xscaled,Gfull$Yscaled,
                                                           1-normalisedX_tilde_full_late[,1]/2-normalisedX_tilde_full_late[,2]/2,op$op3),
     pch=19,main="Procrust PC1-4 for X_late onto G coloured by G")#,sub=paste("Score",scores["full"]))
text(X_tilde_full_late[toname,1],X_tilde_full_late[toname,2],namedf[toname,2],cex=0.5)
# then colour by just North/South
NorthIndlookup <- read_csv("LAD_NorthIndicator_GB_lookup.csv")
GfullExtra <- merge(Gfull, NorthIndlookup, by.x="lad16nm", by.y="LAD", all.x=TRUE, sort=FALSE) 
GfullExtra <- GfullExtra %>% rename("XscaledG" = "Xscaled")
GfullExtra <- GfullExtra %>% rename("YscaledG" = "Yscaled")
# make an average point for the north and south
nsdf_early = cbind(GfullExtra, X_tilde_full_early) # they are already in matching order 
ndf_early = nsdf_early[nsdf_early$NorthIndicator==1,]
sdf_early = nsdf_early[nsdf_early$NorthIndicator==0,]
nsdf_late = cbind(GfullExtra, X_tilde_full_late) # they are already in matching order
ndf_late = nsdf_late[nsdf_late$NorthIndicator==1,]
sdf_late = nsdf_late[nsdf_late$NorthIndicator==0,]
# find averaage point all 4 scenarios
n_early_avg = c(mean(ndf_early[["Xscaled"]]), mean(ndf_early[["Yscaled"]]))
s_early_avg = c(mean(sdf_early[["Xscaled"]]), mean(sdf_early[["Yscaled"]]))
n_late_avg = c(mean(ndf_late[["Xscaled"]]), mean(ndf_late[["Yscaled"]]))
s_late_avg = c(mean(sdf_late[["Xscaled"]]), mean(sdf_late[["Yscaled"]]))
# plot by north/south colour
par(mfrow=c(1,2))
plot(X_tilde_full_early[,1],X_tilde_full_early[,2],col=rgb(GfullExtra$NorthIndicator,0.3,
                                                           1-GfullExtra$NorthIndicator,op$op3),
     pch=19,main="Procrust PC1-4 for X_early onto G coloured by G",
     xlim = c(-0.05, 1.1), ylim = c(-0.2, 0.8))#,sub=paste("Score",scores["full"]))
text(X_tilde_full_early[toname,1],X_tilde_full_early[toname,2],namedf[toname,2],cex=0.5)
points(n_early_avg[1], n_early_avg[2], col = rgb(1,0.3,0), pch=2, cex=1.5) #early = triangle
points(s_early_avg[1], s_early_avg[2], col = rgb(0,0.3,0), pch=2, cex=1.5)#early
points(n_late_avg[1], n_late_avg[2], col = rgb(1,0.3,0), pch=0,cex=1.5)#late = square
points(s_late_avg[1], s_late_avg[2], col = rgb(0,0.3,0), pch=0, cex=1.5)#late

plot(X_tilde_full_late[,1],X_tilde_full_late[,2],col=rgb(GfullExtra$NorthIndicator,0.3,
                                                         1-GfullExtra$NorthIndicator,op$op3),
     pch=19,main="Procrust PC1-4 for X_late onto G coloured by G",
     xlim = c(-0.05, 1.1), ylim = c(-0.2, 0.8))#,sub=paste("Score",scores["full"]))
text(X_tilde_full_late[toname,1],X_tilde_full_late[toname,2],namedf[toname,2],cex=0.5)
points(n_early_avg[1], n_early_avg[2], col = rgb(1,0.3,0), pch=2, cex=1.5)#early
points(s_early_avg[1], s_early_avg[2], col = rgb(0,0.3,0), pch=2, cex=1.5)#early
points(n_late_avg[1], n_late_avg[2], col = rgb(1,0.3,0), pch=0,cex=1.5)#late
points(s_late_avg[1], s_late_avg[2], col = rgb(0,0.3,0), pch=0, cex=1.5)#late
# end 05/01/2023 chunk


# zoom in on london
LAD_RGN_lookup <- read.csv("LAD_Region_lookup_UK.csv")
London <- merge(LAD11_shp, LAD_RGN_lookup, by.x="lad16nm", by.y="LAD", all.x=TRUE) # perfect join
London <- London %>% arrange(index.x)
london_index <- which(London$RGN11NM== "London")  # have to index while all LADs still present.
London <- London[London$RGN11NM== "London" , ]
London <- London %>%
  mutate(centroids = st_centroid(st_geometry(.))) # find centroids
London <-  cbind(London, st_coordinates(London$centroids)) # add X and Y for centroid
London <- London %>% arrange(index)
london_names = list(London$LAD11NM)

# geom_label_repel from ggrepel package

# early
plotLearly <- ggplot(London) + geom_sf(aes(geometry=geometry),
                         fill = rgb(red=normalisedX_tilde_full[london_index,1], green=normalisedX_tilde_full[london_index,2],
                                    blue=1-normalisedX_tilde_full[london_index,1]/2-normalisedX_tilde_full[london_index,2]/2,op$op3[london_index]),
                         col=rgb(red=normalisedX_tilde_full[london_index,1], green=normalisedX_tilde_full[london_index,2], 
                                 blue=1-normalisedX_tilde_full[london_index,1]/2-normalisedX_tilde_full[london_index,2]/2,op$op3[london_index])) + 
  ggtitle("London LADs map coloured by Procrust\n PC1-4 of X early into G") +
  geom_text_repel(aes(X, Y, label=lad16nm), size=2.5)

# late
plotLlate <- ggplot(London) + geom_sf(aes(geometry=geometry),
                            fill = rgb(red=normalisedX_tilde_full_early[london_index,1], green=normalisedX_tilde_full_early[london_index,2],
                                       blue=1-normalisedX_tilde_full_early[london_index,1]/2-normalisedX_tilde_full_early[london_index,2]/2,op$op3[london_index]),
                            col=rgb(red=normalisedX_tilde_full_early[london_index,1], green=normalisedX_tilde_full_early[london_index,2], 
                                    blue=1-normalisedX_tilde_full_early[london_index,1]/2-normalisedX_tilde_full_early[london_index,2]/2,op$op3[london_index])) + 
  ggtitle("London LADs map coloured by Procrust\n PC1-4 of X late into G") +
  geom_text_repel(aes(X, Y, label=lad16nm), size=2.5)
# aes(geometry = geometry) instead of aes(geometry = centroids) is the solution

# average
plotLavg <- ggplot(London) + geom_sf(aes(geometry=geometry),
                         fill = rgb(red=normalisedX_tilde_full_late[london_index,1], green=normalisedX_tilde_full_late[london_index,2],
                                    blue=1-normalisedX_tilde_full_late[london_index,1]/2-normalisedX_tilde_full_late[london_index,2]/2,op$op3[london_index]),
                         col=rgb(red=normalisedX_tilde_full_late[london_index,1], green=normalisedX_tilde_full_late[london_index,2], 
                                 blue=1-normalisedX_tilde_full_late[london_index,1]/2-normalisedX_tilde_full_late[london_index,2]/2,op$op3[london_index])) + 
  ggtitle("London LADs map coloured by Procrust\n PC1-4 of X avg into G") +
  geom_text_repel(aes(X, Y, label=lad16nm), size=10)

plotLavg

plotL <- ggplot(London) + geom_sf(aes(geometry=geometry),
                                  fill = rgb(red=G[london_index,1], green=G[london_index,2], 
                                             blue=1-G[london_index,1]/2-G[london_index,2]/2, 0.3), 
                                  col=rgb(red=G[london_index,1], green=G[london_index,2], 
                                          blue=1-G[london_index,1]/2-G[london_index,2]/2)) +
  ggtitle("London LADs map coloured by G") +
  geom_text_repel(aes(X, Y, label=lad16nm), size=2.5)

grid.arrange(plotLearly, plotLlate, plotLavg, plotL, ncol=2, nrow=2)

ggplot(London) + geom_sf(aes(geometry=centroids),
                         fill = rgb(red=G[london_index,1], green=G[london_index,2], 
                                    blue=1-G[london_index,1]/2-G[london_index,2]/2, 0.3), 
                         col=rgb(red=G[london_index,1], green=G[london_index,2], 
                                 blue=1-G[london_index,1]/2-G[london_index,2]/2),size=30) +
  ggtitle("London LADs map coloured by G") +
  geom_text_repel(aes(X, Y, label=lad16nm), size=3)


# zoom in on london, east of england, and the south east
GL <- merge(LAD11_shp, LAD_RGN_lookup, by.x="lad16nm", by.y="LAD", all.x=TRUE) # perfect join
GL <- GL %>% arrange(index.x)
GL_index <- which(GL$RGN11NM == "London" | GL$RGN11NM == "South East" | GL$RGN11NM == "East of England")  # have to index while all LADs still present.
GL <- GL[(GL$RGN11NM== "London") | (GL$RGN11NM=="South East") | (GL$RGN11NM == "East of England"), ]
GL <- GL %>%
  mutate(centroids = st_centroid(st_geometry(.))) # find centroids
GL <-  cbind(GL, st_coordinates(GL$centroids)) # add X and Y for centroid
GL <- GL %>% arrange(index)
GL_names = list(GL$LAD11NM)

# early
plotLearly <- ggplot(GL) + geom_sf(aes(geometry=geometry),
                                       fill = rgb(red=normalisedX_tilde_full[GL_index,1], green=normalisedX_tilde_full[GL_index,2],
                                                  blue=1-normalisedX_tilde_full[GL_index,1]/2-normalisedX_tilde_full[GL_index,2]/2,op$op3[london_index]),
                                       col=rgb(red=normalisedX_tilde_full[GL_index,1], green=normalisedX_tilde_full[GL_index,2], 
                                               blue=1-normalisedX_tilde_full[GL_index,1]/2-normalisedX_tilde_full[GL_index,2]/2,op$op3[london_index])) + 
  ggtitle("London, East of England & SE LADs map \ncoloured by Procrust PC1-4 of X early into G") #+
  # geom_sf_text(aes(X, Y, label=lad16nm), size=2.5)

# late
plotLlate <- ggplot(GL) + geom_sf(aes(geometry=geometry),
                                      fill = rgb(red=normalisedX_tilde_full_early[GL_index,1], green=normalisedX_tilde_full_early[GL_index,2],
                                                 blue=1-normalisedX_tilde_full_early[GL_index,1]/2-normalisedX_tilde_full_early[GL_index,2]/2,op$op3[london_index]),
                                      col=rgb(red=normalisedX_tilde_full_early[GL_index,1], green=normalisedX_tilde_full_early[GL_index,2], 
                                              blue=1-normalisedX_tilde_full_early[GL_index,1]/2-normalisedX_tilde_full_early[GL_index,2]/2,op$op3[london_index])) + 
  ggtitle("London, East of England & SE LADs map \ncoloured by Procrust PC1-4 of X late into G")# +
  # geom_sf_text(aes(X, Y, label=lad16nm), size=2.5)
# aes(geometry = geometry) instead of aes(geometry = centroids) is the solution

# average
plotLavg <- ggplot(GL) + geom_sf(aes(geometry=geometry),
                                     fill = rgb(red=normalisedX_tilde_full_late[GL_index,1], green=normalisedX_tilde_full_late[GL_index,2],
                                                blue=1-normalisedX_tilde_full_late[GL_index,1]/2-normalisedX_tilde_full_late[GL_index,2]/2,op$op3[london_index]),
                                     col=rgb(red=normalisedX_tilde_full_late[GL_index,1], green=normalisedX_tilde_full_late[GL_index,2], 
                                             blue=1-normalisedX_tilde_full_late[GL_index,1]/2-normalisedX_tilde_full_late[GL_index,2]/2,op$op3[london_index])) + 
  ggtitle("London, East of England & SE LADs map \ncoloured by Procrust PC1-4 of X avg into G")# +
  # geom_sf_text(aes(X, Y, label=lad16nm), size=2.5)

plotL <- ggplot(GL) + geom_sf(aes(geometry=geometry),
                                  fill = rgb(red=G[GL_index,1], green=G[GL_index,2], 
                                             blue=1-G[GL_index,1]/2-G[GL_index,2]/2, 0.3), 
                                  col=rgb(red=G[GL_index,1], green=G[GL_index,2], 
                                          blue=1-G[GL_index,1]/2-G[GL_index,2]/2)) +
  ggtitle("London, East of England & SE LADs map coloured by G") #+
  # geom_sf_text(aes(X, Y, label=lad16nm), size=1.5)

grid.arrange(plotLearly, plotLlate, plotLavg, plotL, ncol=2, nrow=2)

ggplot(GL) + geom_sf(aes(geometry=centroids),
                         fill = rgb(red=G[GL_index,1], green=G[GL_index,2], 
                                    blue=1-G[GL_index,1]/2-G[GL_index,2]/2, 0.3), 
                         col=rgb(red=G[GL_index,1], green=G[GL_index,2], 
                                 blue=1-G[GL_index,1]/2-G[GL_index,2]/2),size=20) +
  ggtitle("London, East of England & SE LADs map coloured by G")# +
  # geom_sf_text(aes(X, Y, label=lad16nm), size=3)


