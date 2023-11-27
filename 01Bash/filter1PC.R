# Reads all the .tsv files, creates host and suffix (.co.uk),
# keeps .co.uk and creates f (unique pc per host per year) and
# the f.old (observations (i.e. pc with duplicates) per host 
# per year).
# E Tranos code

library(stringr)
library(data.table)
library(plyr)
library(dplyr)
library(stats) #
library(urltools)
library(utils)
library(feather)

setwd("/gpfs/bb/tranose/archive/data")

lookup <-fread("NSPL_NOV_2021_UK.csv")
lookup <- as.data.table(lookup)
setnames(lookup, old=("pcds"), new=("pc"))  #++++++++
lookup <- lookup[,.(pc, lat, long)] #++++++++
setkey(lookup,pc) 

files <- list.files(path="/gpfs/bb/tranose/archive/data", pattern="*.tsv", full.names=T, recursive=FALSE)

a <-data.table()
for(i in files) {
  print(i)
  #dt <- fread(sprintf("bzcat %s | tr -d '\\000'", i), header=F)
  dt <- fread(i, header=F)
  dt <- as.data.table(dt)
  print(dim(dt))
  names(dt)[1]<-"url"   #rename
  names(dt)[2]<-"pc"    #rename
  dt$newurl<-substring(dt$url,16)
  dt$domain.ext <- domain(dt$newurl)
  help <- suffix_extract(dt$domain.ext)
  dt <- cbind(dt, help)
  rm(help)
  dt$domain.ext <- NULL #domain.ext is the same with host
  dt$newurl <-NULL
  dt <- dt[suffix=="co.uk",] #keep suffix = .co.uk
  
  print(dim(dt))
  setkey(dt,pc)
  dt <- dt[lookup, nomatch=0] 
  
  dt$year <- str_sub(dt$url, start= 1, end = 4)
  dt <- dt[,.(url, pc, year, host, domain)] # select columns
  print(dim(dt))

  l = list(a,dt)
  rm(dt)
  a <- rbindlist(l, use.names = T, fill = F, idcol=NULL)
  rm(l)
  print(dim(a))
  print(object.size(a))
}

rm(lookup)
############################### 
# unique pc per host per year #
###############################
f <- a[, uniqueN(pc), by = .(host,year)][order(-V1)]
f.old <- a[, .N, by = .(host,year)][order(-N)] # old f 

# one unique postcode only
f <- f[V1 == 1]

setkey(a, host, year) # table, column
setkey(f, host, year) # 
setkey(f.old, host, year) # 

# perform the join, eliminating not matched rows from Right
dim(a)
dim(f)
dim(f.old)
a<-a[f, nomatch=0]
dim(a)
a<-a[f.old, nomatch=0]
dim(a)

#save(a, file="all.RData", compress = "bzip2") 

all.sample <- sample_n(a, 100000)

#write_feather(a, "all.feather")
#write_feather(all.sample, "all_sample.feather")

fwrite(a, "all.1.csv")
fwrite(all.sample, "all_sample.csv")
