library(tidyverse)
library(Matrix)
library(DescTools)
library(Clarity)
library(sf)
library(ggplot2)
library(grid)
library(plotly)

setwd("~/Documents/MiniProject/Network_1_pc/18Procrustes")
# code from 16: matrixSVD script to get in the shp file

## Get the data
namedf=read.csv("LADS_list_GB.csv",header=TRUE) # names of all the LADs in GB, 380 total
namedf <- as.data.frame(namedf[,-1])
rownames(namedf) <- 0:(nrow(namedf)-1)

# read in A matrices
years=2005:2010
files=paste0("A_LAD_", years, "_GB.mtx")
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
LAD11_shp$Xscaled = (LAD11_shp$X - minX)/(maxX-minX)# goes 0 to 1
LAD11_shp$Yscaled = (LAD11_shp$Y - minY)/(maxY-minY)

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

### ANIMATION STUFF 
# 30.11.2022 changed the order and added a few lines so this would run correctly
LAD_name_order = as.data.frame(LAD_name_order)
### make dataframes that are suitable to call upon to make an animation for PC3~PC4
animationData = cbind(SVD_ConcWinsLog$v[ , 1:4] , rep(c(2005:2010), each=n))
animationData = cbind(animationData, rep(LAD11_shp$Xscaled, 6), rep(LAD11_shp$Yscaled, 6), rep(LAD_name_order$V1, 6))
animationData = as.data.frame(animationData)
colnames(animationData) = c("dim1", "dim2", "dim3", "dim4", "year", "Xscaled", "Yscaled", "LAD")
animationData$col = rgb(red=animationData$Xscaled, blue=animationData$Yscaled, green = 0.2 )

#### 13/12/2022 - NEW STUFF FOR PROCRUSTES TRANSFORMATION
animationData[, c("dim1", "dim2", "dim3", "dim4","Xscaled", "Yscaled")] <- sapply(animationData[, c("dim1", "dim2", "dim3", "dim4","Xscaled", "Yscaled")], as.numeric)
G <- as.matrix(animationData[ , c("Xscaled","Yscaled")])
G <- G[1:380,]
Gtwice <- rbind(G,G)

# early X (2005-2007)
early = animationData[animationData$year %in% c('2005','2006','2007'),]
# late X (2008-2010)
late = animationData[animationData$year %in% c('2008','2009','2010'),]

# average embedding across early timepoints
earlyX <- early %>% group_by(LAD) %>% 
  summarise(across(c("dim1", "dim2", "dim3", "dim4"), mean), 
            .groups = 'drop') %>% 
  as.data.frame()
earlyX <- earlyX[ , -1]
earlyX <- as.matrix(earlyX) #as matrix

# average embedding across late timepoints
lateX <- late %>% group_by(LAD) %>% 
  summarise(across(c("dim1", "dim2", "dim3", "dim4"), mean), 
            .groups = 'drop') %>% 
  as.data.frame()
lateX <- lateX[ , -1]
lateX <- as.matrix(lateX) # as matrix

rowbindX <- rbind(earlyX, lateX)

### DOING PROCRUSTES TRANSFORMATION

#MCMCpack
library(MCMCpack)
G_2Nx4 <- cbind(Gtwice, rep(0,dim(Gtwice)[1]), rep(0,dim(Gtwice)[1]))
Proc <- procrustes(rowbindX, G_2Nx4)
X_tilde <- Proc$X.new

library(bigutilsr)
Proc2 <- procrustes(Gtwice, rowbindX) # reference matrix, matrix to transform - order of inputs to bigutilsr func
Proc2R <- Proc2$R
X_tilde2 <- rowbindX %*% Proc2R

library(IMIFA)
Proc3 <- Procrustes(rowbindX, G_2Nx4)

# MCMCpack and IMIFA wont let G have less cols than X. bigutilsr gives rectangular P.
# smacof wouldn't install package

# 14/12/2022
# Try some iterative algorithm
Gtwice <- rbind(G,G)
X <- rowbindX
G_tilde_t1 <- cbind(Gtwice, rep(0,dim(Gtwice)[1]), rep(0,dim(Gtwice)[1])) # \tilde{G}_0 = [G G_0 ] # have to set tilde_G_0 outside of for loop
P_tilde_t <- as.matrix(rep(0, 4^2), nrow =4) # have to decide some P_tilde_0 for it to start from in the distance measure
dist = 10 # so it will start
P_tilde_tm1 <- P_tilde_t
while(dist > 1e-10){
  P_tilde_tm1 <- P_tilde_t #save P_tilde from the last iteration
  Proc <- MCMCpack::procrustes(X, G_tilde_t1)#, translation = TRUE, dilation = TRUE)
  P_tilde_t = Proc$R
  G_t1 = Proc$X.new[,3:4]
  G_tilde_t1 = cbind(Gtwice, G_t1)
  # dist between previous and current P_tilde
  dist <- norm((P_tilde_t - P_tilde_tm1), type="F")
  print(dist)
}

# X_tilde is the approximation of G 
X_tilde = X %*% P_tilde_t[, 1:2]
norm(Gtwice-X_tilde)

# looking at ratios to see if X_tilde and G look comparable
X_tilde = cbind(X_tilde, X_tilde[,1] / X_tilde[,2])
Gtwice = cbind(Gtwice, Gtwice[,1] / Gtwice[,2])

# normalise after splitting into early and late X
earlynormalisedX_tilde = as.data.frame(X_tilde[1:380, 1:2])
colnames(earlynormalisedX_tilde)[1] <- "normd1"
colnames(earlynormalisedX_tilde)[2] <- "normd2"
maxXtilde_early = max(earlynormalisedX_tilde)
earlynormalisedX_tilde = earlynormalisedX_tilde / maxXtilde_early

latenormalisedX_tilde = as.data.frame(X_tilde[381:760,1:2])
colnames(latenormalisedX_tilde)[1] <- "normd1"
colnames(latenormalisedX_tilde)[1] <- "normd2"
maxXtilde_late = max(latenormalisedX_tilde)
latenormalisedX_tilde = latenormalisedX_tilde / maxXtilde_late

# set values < 1/255 to 0 so that the rgb works
earlynormalisedX_tilde[earlynormalisedX_tilde < 1/255] <- 0
latenormalisedX_tilde[latenormalisedX_tilde < 1/255] <- 0

LAD11_shp_early = cbind(LAD11_shp, earlynormalisedX_tilde)
colnames(LAD11_shp_early)
LAD11_shp_late = cbind(LAD11_shp, latenormalisedX_tilde)

# Read in LAD to opacity (connectivity outward (rowmeans) table)
opacity_lookup_UK = read_csv("opacity_lookup_UK.csv")

# merge to LAD11_shp_early
LAD11_shp_early = merge(LAD11_shp_early, opacity_lookup_UK, by.x = "lad16nm", by.y="LAD", how="left" )
LAD11_shp_early <- LAD11_shp_early %>%
  mutate(centroids = st_centroid(st_geometry(.))) # find centroids

# merge to LAD11_shp_late
LAD11_shp_late = merge(LAD11_shp_late, opacity_lookup_UK, by.x = "lad16nm", by.y="LAD", how="left" )
LAD11_shp_late <- LAD11_shp_late %>%
  mutate(centroids = st_centroid(st_geometry(.))) # find centroids

### map of UK coloured by 2D colourmap - filled
# early
ggplot(LAD11_shp_early) + geom_sf(aes(geometry=centroids), size=1,
                            fill = rgb(red=LAD11_shp_early$normd1, blue=LAD11_shp_early$normd2, green = 0.2), 
                            col=rgb(red=LAD11_shp_early$normd1, blue=LAD11_shp_early$normd2, green = 0.2),lwd=0.1) + 
  ggtitle("GB LADs map filled by 2D Procrustes \nEARLY") 


# late
ggplot(LAD11_shp_late) + geom_sf(aes(geometry=centroids), 
                                  fill = rgb(red=LAD11_shp_late$normd1, blue=LAD11_shp_late$normd2, green = 0.2), 
                                  col=rgb(red=LAD11_shp_late$normd1, blue=LAD11_shp_late$normd2, green = 0.2),lwd=0.1) + 
  ggtitle("GB LADs map filled by 2D Procrustes LATE") 
