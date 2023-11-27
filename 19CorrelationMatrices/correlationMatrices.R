setwd("~/Documents/MiniProject/Network_1_pc/19CorrelationMatrices")

library(tidyverse)
library(Matrix)
library(DescTools)
library(sf)
library(ggplot2)
library(grid)
library(plotly)
library(gridExtra)
library(ggrepel)
library(reshape)


# normalise so values go 0 to 1
normalise=function(x,trunc=0){
  x=as.matrix(x)
  if(all(!is.null(colnames(x)))) colnames(x)=paste0("norm",colnames(x))
  x = apply(x,2,function(x){x + abs(min(x))})
  x = apply(x,2,function(x){x/max(x)})
  x[x < trunc] <- 0
  x
}

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

### Geography data (some is not really needed in this chunk but doesn't matter)

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

### Active enterprise count 2010 data

# for count of active enterprises 2010
#had to look at animations.ipynb in 15 to see how to clean 
ent2010 <- read.csv("count_of_active_enterprises_2010_byLAD_GB_CLEAN.csv")
# remove commas in the numeric fields
ent2010$count.of.active.enterprises <- as.numeric(gsub(",","",ent2010$count.of.active.enterprises))
ent2010 = merge(ent2010, namedf, by.x="LAD_clean", by.y="LAD")
ent2010 = ent2010 %>% 
  arrange(index)

### Pay data

pay2010 <- read.csv("CleanedMedianMeanPay2010.csv")
pay2010 = pay2010[ , 3:5] # keep Median, mean, and LAD_clean
pay2010 = merge(pay2010, namedf, by.x="LAD_clean", by.y="LAD")
pay2010 = pay2010 %>% 
  arrange(index)
# if there is an x for the value, replaced it with 0, as we don't know the value
pay2010$Median <- as.numeric(gsub("x",0,pay2010$Median))
pay2010$Mean <- as.numeric(gsub("x",0,pay2010$Mean))

### Population Density data

popDens2011 = read.csv("popDens_2011census.csv", header= FALSE)
ladwithcode = read.csv("LAD_NorthIndicator_GB_lookup.csv")
ladwithcode = ladwithcode[, c("LAD", "LAD11CD")]
namecodedf = merge(namedf, ladwithcode, by="LAD") 
namecodedf = namecodedf %>% arrange(index)
popDens2011 = merge(namecodedf, popDens2011, by.x ="LAD11CD", by.y = "V1", all.x=TRUE)
# 4 don't match by LAD11CD, but match by name
popDens2011 = popDens2011[, c("LAD", "LAD11CD", "index", "V8")]
colnames(popDens2011)[4] <- "Population Density"
# manually fill 
popDens2011[popDens2011$LAD=="Fife", "Population Density"] = 2.8
popDens2011[popDens2011$LAD=="Perth and Kinross", "Population Density"] = 0.3
popDens2011[popDens2011$LAD=="North Lanarkshire", "Population Density"] = 7.2
popDens2011[popDens2011$LAD=="Glasgow City", "Population Density"] = 34.0
# order
popDens2011 = popDens2011 %>% arrange(index)


### Distance to Westminster (London)

Westminster = Gfull[Gfull$lad16nm =="Westminster",]
EuclidDist_Westminster = function(x,y){
  X_W = Westminster$X
  Y_W = Westminster$Y
  dist = sqrt( (x - X_W)^2 + (y - Y_W)^2 )
}
Gfull$disttoWestminster = EuclidDist_Westminster(Gfull$X, Gfull$Y)

### Qualification level
# census 2011 data: https://www.ons.gov.uk/peoplepopulationandcommunity/educationandchildcare/bulletins/educationenglandandwales/census2021#:~:text=Level%201%3A%20one%20to%20four,A%20Levels%20or%20equivalent%20qualifications
# level 4 and abouve is a HNC, HND, Bachelor's, post-grad qual (highest level)

qual2011 = read.csv("highestLevelofQualification_LAs_GB_census2011.csv")
# Folkestone and Hythe is the same place as Shepway
qual2011[qual2011$geography.code=="E07000112", "geography"]="Shepway"
# taf or taff
qual2011[qual2011$geography.code=="W06000016", "geography"]="Rhondda Cynon Taf"
# no "the"
qual2011[qual2011$geography.code=="W06000014", "geography"]="The Vale of Glamorgan"

qual2011 = merge(qual2011, namecodedf, by.x="geography", by.y="LAD", all=FALSE)
# replace NA with 0 - Scotland doesn't record apprentices or other quals
qual2011[is.na(qual2011)] <- 0
qual2011 = qual2011 %>% arrange(index)

qual2011$propL4orAbove = qual2011[,10]/qual2011[,4]


### Social Grades - census 2011
SG = read.csv("socialGrades_census2011.csv")
# Folkestone and Hythe is the same place as Shepway
SG[SG$geography.code=="E07000112", "geography"]="Shepway"
# taf or taff
SG[SG$geography.code=="W06000016", "geography"]="Rhondda Cynon Taf"
# no "the"
SG[SG$geography.code=="W06000014", "geography"]="The Vale of Glamorgan"
SG = merge(SG, namecodedf,  by.x="geography", by.y="LAD", all=FALSE)
SG[is.na(SG)] <- 0
SG = SG %>% arrange(index)
SG$propAB = SG[ ,5]/SG[ ,4]
# AB = Higher & intermediate managerial, administrative, professional occupations


### Country of Birth - census 2011
cob = read.csv("countryofBirth_census2011.csv")
# Folkestone and Hythe is the same place as Shepway
cob[cob$geography.code=="E07000112", "geography"]="Shepway"
# taf or taff
cob[cob$geography.code=="W06000016", "geography"]="Rhondda Cynon Taf"
# no "the"
cob[cob$geography.code=="W06000014", "geography"]="The Vale of Glamorgan"
cob = merge(cob, namecodedf,  by.x="geography", by.y="LAD", all=FALSE)
cob[is.na(cob)] <- 0
cob = cob %>% arrange(index)
cob$propUK = cob[ ,6]/cob[ ,4]

## CORRELATION MATRICES

# put all the variables to test into one matrix
# can cbind as the rows are all in the same order, according to the index which is defined by LADs_list_GB
variables = cbind(Gfull[, c("X", "Y", "disttoWestminster")] , ent2010$count.of.active.enterprises, 
                  pay2010[, c("Median","Mean")], popDens2011[,c("Population Density")], qual2011$propL4orAbove)
                  #, SG$propAB, cob$propUK)

corMat = cor(avgX, variables)

# via ggplot
x = colnames(variables)
y = paste("dim", c(1:4))
data = expand.grid(X=x, Y=y)
m = matrix(corMat, nrow = length(y))
colnames(m) = c("Longitude", "Latitude", "Euclidean Distance to Westminster", "#Active Enterprises (2010)", "Median Pay (2010)", 
                "Mean Pay (2010)", "Population Density (2011 census)", "Proportion of People with\n a Level 4 or above Qualification\n (2011 census)")#, 
                #"Social Grade A or B", "Proportion of People\n Born in the UK")
rownames(m) = y 
df = melt(m)
colnames(df) <- c("x", "y", "value")
corPlot = ggplot(df , aes(x, y, fill=value)) + geom_tile() +
  coord_fixed() +
  geom_text(aes(label = round(value, 2)), color = "white", size = 4) +
  labs(x="", y="") +
  ggtitle("(a) Correlation Matrix between Embedding Vectors\n and Social/Economic/Geographic Properties") + 
  theme(plot.title = element_text(hjust = 0.7, vjust=0))
  #scale_fill_gradient2(low = "blue", mid="white", high="blue") 

corPlot


###
#### LINEAR MODEL 
###

# G ~ \sum \beta_i X_i
G = Gfull[,c("X", "Y")]
long = G[,1]
lat = G[,2]
lm_long = lm(long ~ avgX)
lm_long
lm_lat = lm(lat ~ avgX)
lm_lat


###
# early/late cor
###
avgSvdX_early=t(getAvg(t(SVD_ConcWinsLog$v[1:(3*380), ])))
avgSvdX_late=t(getAvg(t(SVD_ConcWinsLog$v[(3*380 +1):(6*380), ])))
avgX_early=avgSvdX_early[,1:4] %*% (diag(SVD_ConcWinsLog$d)[1:4, 1:4])^0.5
avgX_late=avgSvdX_late[,1:4] %*% (diag(SVD_ConcWinsLog$d)[1:4, 1:4])^0.5

variables = cbind(Gfull[, c("X", "Y")])

corMat_early = cor(avgX_early, variables)
# via ggplot
x = colnames(variables)
y = paste("dim", c(1:4))
data = expand.grid(X=x, Y=y)
m = matrix(corMat_early, nrow = length(y))
colnames(m) = c("Longitude", "Latitude")
rownames(m) = y 
df = melt(m)
colnames(df) <- c("x", "y", "value")
corPlot_early = ggplot(df , aes(x, y, fill=value)) + geom_tile() +
  coord_fixed() +
  geom_text(aes(label = round(value, 2)), color = "white", size = 5) +
  labs(x="", y="") +
  ggtitle("(b) Correlation Matrix between EARLY Embedding Vectors\n and Social/Economic/Geographic Properties") + 
  theme(plot.title = element_text(hjust = 0.4, vjust=0))

corPlot_early

corMat_late = cor(avgX_late, variables)
# via ggplot
x = colnames(variables)
y = paste("dim", c(1:4))
data = expand.grid(X=x, Y=y)
m = matrix(corMat_late, nrow = length(y))
colnames(m) = c("Longitude", "Latitude")
rownames(m) = y 
df = melt(m)
colnames(df) <- c("x", "y", "value")
corPlot_late = ggplot(df , aes(x, y, fill=value)) + geom_tile() +
  coord_fixed() +
  geom_text(aes(label = round(value, 2)), color = "white", size = 5) +
  labs(x="", y="") +
  ggtitle("(c) Correlation Matrix between LATE Embedding Vectors\n and Social/Economic/Geographic Properties") + 
  theme(plot.title = element_text(hjust = 0.4, vjust=0))

corPlot_late


