library(Matrix)
library(DescTools)

setwd("~/Documents/MiniProject/Network_1_pc/16GB_ConcThenWins")

# make a vector of all the values of all the LAD x LAD GB matrices
# and a matrix
A_vec_all = c()
A_mat_all = c()
for (year in 2005:2010){
  filename = paste0("A_LAD_", year, "_UK_v5.mtx")
  A = readMM(filename)
  A = as.matrix(A)
  A_vec = as.vector(A)
  A_vec_all = c(A_vec_all, A_vec)
  A_mat_all = cbind(A_mat_all, A)
}

n = dim(A_mat_all)[1] 
T=6

## This is from late Feb 2024, trying to see if the data spoke about why to Winsorize.
## From this we see we need to use the embeddings to justify our choice of where to Winsorize

## FIND OUT MEAN AND SPREAD VALUES FOR DATA, WINSORIZED AND UNWINSORIZED
mean(A_vec_all)
sd(A_vec_all)
Wins99raw = Winsorize(A_vec_all, minval = 0, probs=c(0,0.99))
mean(Wins99raw)
sd(Wins99raw)
Wins98raw = Winsorize(A_vec_all, minval = 0, probs=c(0,0.98))
mean(Wins98raw)
sd(Wins98raw)
Wins99_9raw = Winsorize(A_vec_all, minval = 0, probs=c(0,0.999))
mean(Wins99_9raw)
sd(Wins99_9raw)
Wins99_95raw = Winsorize(A_vec_all, minval = 0, probs=c(0,0.9995))
mean(Wins99_95raw)
sd(Wins99_95raw)
Wins99_99raw = Winsorize(A_vec_all, minval = 0, probs=c(0,0.9999))
mean(Wins99_99raw)
sd(Wins99_99raw)


par(mfrow=c(2,3))
# density plot Wins(0, 0.98)
A_wins98 <- Winsorize(log10(A_vec_all+1), minval = 0, probs=c(0,0.98))
plot(density(A_wins98),
     main="Values in log10(A+1) for all years,\n when all data Winsorized at 98th quantile", 
     ylab="''density''", xlab="log10(A+1) value",
     col="purple")

# density plot Wins(0, 0.999)
A_wins999 <- Winsorize(log10(A_vec_all+1), minval = 0, probs=c(0,0.999))
plot(density(A_wins999),
     main="Values in log10(A+1) for all years,\n when all data Winsorized at 99.9th quantile", 
     ylab="''density''", xlab="log10(A+1) value",
     col="purple")

# no wins
plot(density(log10(A_vec_all+1)),
     main="Values in log10(A+1) for all years,\n with no Winsorization", 
     ylab="''density''", xlab="log10(A+1) value",
     col="purple")

# density plot Wins(0, 0.99)
A_wins <- Winsorize(log10(A_vec_all+1), minval = 0, probs=c(0,0.99))
plot(density(A_wins),
     main="Values in log10(A+1) for all years,\n when all data Winsorized at 99th quantile", 
     ylab="''density''", xlab="log10(A+1) value",
     col="purple")

# density plot Wins(0, 0.9995)
A_wins <- Winsorize(log10(A_vec_all+1), minval = 0, probs=c(0,0.9995))
plot(density(A_wins),
     main="Values in log10(A+1) for all years,\n when all data Winsorized at 99.95th quantile", 
     ylab="''density''", xlab="log10(A+1) value",
     col="purple")

# density plot Wins(0, 0.9999)
A_wins <- Winsorize(log10(A_vec_all+1), minval = 0, probs=c(0,0.9999))
plot(density(A_wins),
     main="Values in log10(A+1) for all years,\n when all data Winsorized at 99.99th quantile", 
     ylab="''density''", xlab="log10(A+1) value",
     col="purple")



### Embeddings to decide where to Winsorize
# Need some measure within the embedding

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

d=4

# No Wins - plot for 2010
svd_A <- svd(A_mat_all)
U_A = svd_A$u[1:d, 1:d]
D_A = diag(svd_A$d[1:d])
V_A = svd_A$v[ ,1:d]
# embedding
Y_hat = V_A %*% (D_A)^0.5
# plot
plot(Y_hat[1901:2280,1],Y_hat[1901:2280,2],col=rgb(Gfull$Xscaled,Gfull$Yscaled,1),pch=19,
     xlab="Dimension 1", ylab="Dimension 2", main="NO Winsorization, 2010")

# Wins (0, 0.99) - plot for 2010
A_wins <- Winsorize(log10(A_vec_all+1), minval = 0, probs=c(0,0.99))
A_wins_mat = matrix(A_wins, nrow = n, ncol = n * T)
svd_A <- svd(A_wins_mat)
U_A = svd_A$u[1:d, 1:d]
D_A = diag(svd_A$d[1:d])
V_A = svd_A$v[ ,1:d]
# embedding
Y_hat = V_A %*% (D_A)^0.5
# plot
plot(Y_hat[1901:2280,1],Y_hat[1901:2280,2],col=rgb(Gfull$Xscaled,Gfull$Yscaled,1),pch=19,
     xlab="Dimension 1", ylab="Dimension 2", main="(0, 0.99) Winsorization, 2010")

# Wins (0, 0.98) - plot for 2010
A_wins <- Winsorize(log10(A_vec_all+1), minval = 0, probs=c(0,0.98))
A_wins_mat = matrix(A_wins, nrow = n, ncol = n * T)
svd_A <- svd(A_wins_mat)
U_A = svd_A$u[1:d, 1:d]
D_A = diag(svd_A$d[1:d])
V_A = svd_A$v[ ,1:d]
# embedding
Y_hat = V_A %*% (D_A)^0.5
# plot
plot(Y_hat[1901:2280,1],Y_hat[1901:2280,2],col=rgb(Gfull$Xscaled,Gfull$Yscaled,1),pch=19,
     xlab="Dimension 1", ylab="Dimension 2", main="(0, 0.98) Winsorization, 2010")

# Wins (0, 0.995) - plot for 2010
A_wins <- Winsorize(log10(A_vec_all+1), minval = 0, probs=c(0,0.995))
A_wins_mat = matrix(A_wins, nrow = n, ncol = n * T)
svd_A <- svd(A_wins_mat)
U_A = svd_A$u[1:d, 1:d]
D_A = diag(svd_A$d[1:d])
V_A = svd_A$v[ ,1:d]
# embedding
Y_hat = V_A %*% (D_A)^0.5
# plot
plot(Y_hat[1901:2280,1],Y_hat[1901:2280,2],col=rgb(Gfull$Xscaled,Gfull$Yscaled,1),pch=19,
     xlab="Dimension 1", ylab="Dimension 2", main="(0, 0.995) Winsorization, 2010")


# Wins (0, 0.999) - plot for 2010 - this just scales 144 things in the concat A raw matrix
A_wins <- Winsorize(log10(A_vec_all+1), minval = 0, probs=c(0,0.999))
A_wins_mat = matrix(A_wins, nrow = n, ncol = n * T)
svd_A <- svd(A_wins_mat)
U_A = svd_A$u[1:d, 1:d]
D_A = diag(svd_A$d[1:d])
V_A = svd_A$v[ ,1:d]
# embedding
Y_hat = V_A %*% (D_A)^0.5
# plot
plot(Y_hat[1901:2280,1],Y_hat[1901:2280,2],col=rgb(Gfull$Xscaled,Gfull$Yscaled,1),pch=19,
     xlab="Dimension 1", ylab="Dimension 2", main="(0, 0.999) Winsorization, 2010")


# SOME KIND OF MATHEMATICAL MEASURE
# Calculate L2 norms of each row
l2_norms <- sqrt(rowSums(Y_hat^2))
mean_l2_norms = mean(l2_norms)
max_l2_norms = max(l2_norms)
ratio_mean_max = mean_l2_norms / max_l2_norms
print(ratio_mean_max)
sorted_l2_norms = sort(l2_norms)
sorted_l2_norms[ceil(0.9*n*T)] / max_l2_norms # (0.9 is just a choice to say you don't want 90% 
                        #of the data to be really close to the origin and <10% far away right))
