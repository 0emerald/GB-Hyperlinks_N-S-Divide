library(Matrix)
library(DescTools)
library("gplots")

setwd("~/Documents/MiniProject/Network_1_pc/16GB_ConcThenWins")
### 21/12/2022 - plots for concatenate then Winsorize the data
# also this will use the GB data
# going to make a new script in folder 16
# prev version in 12

# make a vector of all the values of all the LAD x LAD GB matrices
A_vec_all = c()
for (year in 2005:2010){
  filename = paste0("A_LAD_", year, "_UK_v5.mtx")
  A = readMM(filename)
  A = as.matrix(A)
  A_vec = as.vector(A)
  A_vec_all = c(A_vec_all, A_vec)
}


# josh test
library(ggplot2)
df=data.frame(a=log10(A_vec+1))
ggplot(df,aes(x=a))+#+geom_density(bw=0.05) + xlim(c(0,3.5)) +
  geom_histogram(
    #stat="density",
    binwidth = 0.3)


#####
## REAL DATA
plot(density(A_vec_all, bw=0.3), main="Values in A matrices, raw data, across all years", 
     xlab="A value", ylab="''density''")
plot(density(A_vec_all, bw=0.3), ylim=c(0,0.0001), 
     ylab="''density''", xlab = "A value")

quantile99 = quantile(log10(A_vec_all+1), probs = (0.99))

plot(density(log10(A_vec_all +1)), xlim=c(0,4.5), 
     main="Values in log10(A+1) for all years of data", ylab="''density''", xlab="log10(A+1) value",
     col="darkgreen")
abline(v=quantile99, col="purple") # max value of A_wins
text(labels="1%",x=0.9, y=10, col="purple")
text(labels="99%",x=0.3, y=10, col="darkgreen")

A_wins <- Winsorize(log10(A_vec_all+1), minval = 0, probs=c(0,0.99))
plot(density(A_wins),
     main="Values in log10(A+1) for all years,\n when all data Winsorized at 99th quantile", 
     ylab="''density''", xlab="log10(A+1) value",
     col="purple")


## as histograms
hist(A_vec_all, breaks=50, main="Values in A matrices, raw data, across all years", xlab="A value")#, ylab="''density''")
hist(A_vec_all, breaks=50, ylim=c(0,0.0001))#, xlab = "A value")

quantile99 = quantile(log10(A_vec_all+1), probs = (0.99))

hist(log10(A_vec_all +1), xlim=c(0,4.5),#breaks=c(log10(0.5+c(0:10))),
     main="Values in log10(A+1) for all years of data", ylab="Frequency", xlab="log10(A+1) value",
     col="lightskyblue3")
abline(v=quantile99, col="purple") # max value of A_wins
text(labels="1%",x=0.8, y=4e+05, col="purple")
text(labels="99%",x=0.4, y=4e+05, col="lightskyblue4")

A_wins <- Winsorize(log10(A_vec_all+1), minval = 0, probs=c(0,0.99))
hist(A_wins, breaks=50,
     main="Values in log10(A+1) for all years,\n when all data Winsorized at 99th quantile", 
     ylab="''density''", xlab="log10(A+1) value",
     col="lightskyblue3")

hist(A_wins, breaks=log10(0.5+c(0:4)),
     main="Values in log10(A+1) for all years,\n when all data Winsorized at 99th quantile", 
     ylab="Frequency", xlab="log10(A+1) value",
     xlim = c(0, log10(4)),
     col="lightskyblue3")

