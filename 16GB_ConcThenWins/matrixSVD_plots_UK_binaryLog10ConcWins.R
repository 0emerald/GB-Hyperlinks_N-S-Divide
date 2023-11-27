library(tidyverse)
library(Matrix)
library(DescTools)
library(Clarity)
library(sf)
library(ggplot2)
library(grid)
library(plotly)

setwd("~/Documents/MiniProject/Network_1_pc/16GB_ConcThenWins")

####
## THIS SCRIPT MAKES THE GB ANIMATION DATA (via SVD) BUT CONCATENATES THE DATA BEFORE WINSORIZING
# Winsorization is done before concatenation in folder 15
####
## Get the data
namedf=read.csv("LADS_list_UK.csv",header=TRUE) # names of all the LADs in the UK, 380 total
namedf <- as.data.frame(namedf[,-1])
rownames(namedf) <- 0:(nrow(namedf)-1)

# read in A matrices
years=2005:2010
files=paste0("A_LAD_", years, "_UK_v5.mtx")
allf=lapply(files,readMM)

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

# Winsorize the concatenated A matrix
nrows = nrow(alldata)
A_wins <- Winsorize(alldata, minval = 0, probs=c(0,0.99))
A_wins <- as.matrix(A_wins, nrow=nrows)
# then log
A <- log10(A_wins+1)


### for the AConcWinsLog data:
SVD_ConcWinsLog <- svd(A)



### ADD COLOURS FOR THE LADs DEPENDING ON GEOGRAPHICAL POSITION

# 29/11/22 changed the shp file here to try and have it for all UK, inc Scot (not just E&W)
LAD11_shp <- read_sf('Local_Authority_Districts_December_2016_GCB_in_the_UK.shp')
LAD11_shp <- LAD11_shp %>%
  mutate(centroids = st_centroid(st_geometry(.))) # find centroids
LAD11_shp <-  cbind(LAD11_shp, st_coordinates(LAD11_shp$centroids)) # add X and Y for centroid
maxX = max(LAD11_shp$X)
maxY = max(LAD11_shp$Y)
minX = min(LAD11_shp$X)
minY = min(LAD11_shp$Y)

# order the LAD11_shp dataframe so that it is in the same order as the A matrices
LAD_name_order = as.matrix(rownames(plotdata))
LAD_name_order = cbind(LAD_name_order, 1:n)

# rename 'Vale of Glamorgan' in shp file to 'The Vale of Glamorgan' so that it matches
 # manually index cell to change
LAD11_shp[382, 3] <-  'The Vale of Glamorgan' 

LAD11_shp = merge(LAD11_shp, LAD_name_order, by.x="lad16nm", by.y="V1") # by.x="LAd16NM" instead of "lad11nm"
LAD11_shp$V2 = as.numeric(LAD11_shp$V2)
LAD11_shp = LAD11_shp %>%
  arrange(V2)

# test1 = LAD11_shp$lad16nm
# setdiff(LAD_name_order, test1)

# LADs NOT in the shp file: Local_Authority_Districts_December_2016_GCB_in_the_UK.shp
# "The Vale of Glamorgan" 
# BUT IT IS THERE AS "Vle of Glamorgan" !!!

### log10(A1+...+A6) data
dat = as.data.frame(SVD_sumAm$u[,1:5])
dat$X = LAD11_shp$X
dat$Y = LAD11_shp$Y
dat$Xscaled = (dat$X - minX)/(maxX-minX)# goes 0 to 1
dat$Yscaled = (dat$Y - minY)/(maxY-minY)

PC1lim=range(SVD_ConcWinsLog$v[,1])
PC2lim=range(SVD_ConcWinsLog$v[,2])
PC3lim=range(SVD_ConcWinsLog$v[,3])
PC4lim=range(SVD_ConcWinsLog$v[,4])
PC5lim=range(SVD_ConcWinsLog$v[,5])



### ANIMATION STUFF 10/08/2022- -MOVED TO PYTHON
# 30.11.2022 changed the order and added a few lines so this would run correctly
LAD_name_order = as.data.frame(LAD_name_order)
### make dataframes that are suitable to call upon to make an animation for PC3~PC4
animationData = cbind(SVD_ConcWinsLog$v[ , 1:4] , rep(c(2005:2010), each=n))
animationData = cbind(animationData, rep(dat$Xscaled, 6), rep(dat$Yscaled, 6), rep(LAD_name_order$V1, 6))
animationData = as.data.frame(animationData)
colnames(animationData) = c("dim1", "dim2", "dim3", "dim4", "year", "Xscaled", "Yscaled", "LAD")
animationData$col = rgb(red=animationData$Xscaled, blue=animationData$Yscaled, green = 0.2 )

# adds in some extra columns that can be used to label a few nodes
anim = animationData
anim$LAD3_4 = rep("", dim(anim)[1])

for (i in 1:dim(anim)[1]){
  if (anim[i,"LAD"] =="Westminster"){
    anim[i,"LAD3_4"] = "Westminster"
  }
  else if (anim[i,"LAD"] =="Leeds"){
    anim[i,"LAD3_4"] = "Leeds"
  }
  else if (anim[i,"LAD"] =="Cornwall"){
    anim[i,"LAD3_4"] = "Cornwall"
  }
  else if (anim[i,"LAD"] =="Oxford"){
    anim[i,"LAD3_4"] = "Oxford"
  }
  else if (anim[i,"LAD"] =="Highland"){
    anim[i,"LAD3_4"] = "Highland"
  }
  next
}

anim$LAD1_2 = rep("", dim(anim)[1])

for (i in 1:dim(anim)[1]){
  if (anim[i,"LAD"] =="Westminster"){
    anim[i,"LAD1_2"] = "Westminster"
  }
  else if (anim[i,"LAD"] =="Isles of Scilly"){
    anim[i,"LAD1_2"] = "Isles of Scilly"
  }
  else if (anim[i,"LAD"] =="Islington"){
    anim[i,"LAD1_2"] = "Islington"
  }
  else if (anim[i,"LAD"] =="Camden"){
    anim[i,"LAD1_2"] = "Camden"
  }
  
  next
}


## save the csv as the animationDataUK
write.csv(anim, "animationDataGB_ConcWins.csv")

####################################
## PLOT MAKING SECTION 
####################################

pdf("screePlot_A__ConcWins1Log_titledPres_GB.pdf",
    width=8, height=3)
plot(rev(SVD_ConcWinsLog$d), pch=20, ylab="singular value", 
     main="Scree plot of log10(Singular Values) from SVD of Data for GB",
     col=rev(ifelse(SVD_ConcWinsLog$d < SVD_ConcWinsLog$d[4], "black", "red")),
     xlim=c(380,1),xaxt="n") # first 4 dots red
axis(1,at=c(380,280,180,80),label=c(0,100,200,300))
dev.off()

singularValues = SVD_ConcWinsLog$d


## Embedding pLots

pdf("rgb_SVD_A_ConcWins1Log_year_PC3PC4_UK.pdf")
for (t in 1:length(years)){
  plot(SVD_ConcWinsLog$v[n*(t-1)+1:n,3], SVD_ConcWinsLog$v[n*(t-1)+1:n,4],
       xlab="PC3 (+ve to -ve)",ylab="PC4",
       col=rgb(red=dat$Xscaled, blue=dat$Yscaled, green = 0.2 ), 
       pch=20, cex=2, height=20, width=10,
       main=paste0(years[t], " from Concatenating all data then \nWinsorizing c(0, 0.99), log10 PC3 and PC4"),
       xlim=-PC3lim,ylim=PC4lim)
}
dev.off()


### add labels for the outlying places in PC1/PC2
pdf("labelled_rgb_SVD_AWins1Log_year_UK.pdf")
for (t in 1:length(years)){
  plot(SVD_ConcWinsLog$v[n*(t-1)+1:n,1], SVD_ConcWinsLog$v[n*(t-1)+1:n,2],
       xlab="Dimension 1",ylab="Dimension 2",
       col=rgb(red=dat$Xscaled, blue=dat$Yscaled, green = 0.2 ), 
       pch=20, cex=2, height=20, width=10,
       main=paste0(years[t], " from Concatenating all data then \nWinsorizing c(0, 0.99), log10"),
       xlim=PC1lim,ylim=PC2lim)
  outliers = intersect(which(SVD_ConcWinsLog$v[n*(t-1)+1:n,2] > 0.01), which(SVD_ConcWinsLog$v[n*(t-1)+1:n,1] < -0.04))
  text(SVD_ConcWinsLog$v[n*(t-1)+1:n,1][outliers], SVD_ConcWinsLog$v[n*(t-1)+1:n,2][outliers],
       labels=rownames(plotdata)[outliers],cex=0.4, adj = c(-0.2,0))
  
  plot(SVD_ConcWinsLog$v[n*(t-1)+1:n,2], SVD_ConcWinsLog$v[n*(t-1)+1:n,3],
       xlab="Dimension 2",ylab="Dimension 3",
       col=rgb(red=dat$Xscaled, blue=dat$Yscaled, green = 0.2 ), 
       pch=20, cex=2, height=20, width=10,
       main=paste0(years[t], "Concatenate then Winsorized c(0, 0.99), log10"),
       xlim=PC2lim,ylim=PC3lim)
  outliers = which(SVD_ConcWinsLog$v[n*(t-1)+1:n,2] > 0.045)
  text(SVD_ConcWinsLog$v[n*(t-1)+1:n,2][outliers], SVD_ConcWinsLog$v[n*(t-1)+1:n,3][outliers],
       labels=rownames(plotdata)[outliers],cex=0.4, adj = -0.15)
  
  plot(SVD_ConcWinsLog$v[n*(t-1)+1:n,3], SVD_ConcWinsLog$v[n*(t-1)+1:n,4],
       xlab="Dimension 3",ylab="Dimension 4",
       col=rgb(red=dat$Xscaled, blue=dat$Yscaled, green = 0.2 ), 
       pch=20, cex=2, height=20, width=10,
       main=paste0(years[t], " Concatenate then Winsorized c(0, 0.99), log10 \nPC3 and PC4"),
       xlim=-PC3lim,ylim=-PC4lim)
  # add some labels to points - 5 selected, Westminster, Leeds, Cornwall, Oxford, Highland
  numbers = c(1,7,10, 70)
  text(SVD_ConcWinsLog$v[n*(t-1)+1:n,3][numbers], SVD_ConcWinsLog$v[n*(t-1)+1:n,4][numbers],
       labels=rownames(plotdata)[numbers],cex=0.4, adj = c(-0.2,0))
  # minus sign in front of limits to rotate the plot 180 degrees
  
  plot(SVD_ConcWinsLog$v[n*(t-1)+1:n,4], SVD_ConcWinsLog$v[n*(t-1)+1:n,5],
       xlab="Dimension 4",ylab="Dimension 5",
       col=rgb(red=dat$Xscaled, blue=dat$Yscaled, green = 0.2 ), 
       pch=20, cex=2, height=20, width=10,
       main=paste0(years[t], " Concatenate then Winsorized c(0, 0.99), log10"),
       xlim=PC4lim,ylim=PC5lim)
}
dev.off()








