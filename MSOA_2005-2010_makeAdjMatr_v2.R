library(tidyverse)
library(data.table)
library(stats) 
library(urltools)
library(utils)
library(feather)
library(igraph)
library(Matrix)
library(MASS)

## THIS IS VERSION 2:
## ALL HYPERLINKS WHERE HOST1 AND HOST2 ARE THE SAME ARE REMOVED. 

### For all years:
#################################
### HOST LOOKUP TABLE
#################################

# Filter postcode lookup table
lookup_all_hosts <- read.csv("lookup_all.csv")
#some hosts have more than one postcode.
# Just keep the first observation in the table that the host appears for
lookup_all_hosts_1 <- lookup_all_hosts %>%
  distinct(host, .keep_all = TRUE)


#################################
### MSOA TABLE
#################################

# Read in file that links MSOA and postcode
all_MSOA <- read.csv("MSOA_UK_nov2021.csv")

# blank MSOA filled by NA
all_MSOA <- all_MSOA %>%
  mutate(msoa11 = replace(msoa11, msoa11=="", NA )) #%>%
#  drop_na(msoa11)

# For the building of adjacency matrices

# NODES
# vector of just the nodes which are the MSOAs
nodes <- as.matrix(unique(all_MSOA$msoa11))

# turn into a dataframe
nodes <- as.data.frame(nodes)

# remove node for NA for msoa - get rid of any obs with no msoa 
nodes <- nodes %>% drop_na()

# rename column
names(nodes)[1] <- "msoa"

# add index column 
nodes$index <- 1:nrow(nodes)

#################################
### OUTPUT TABLE TO FILL
#################################
output = matrix(rep(NA, 4*(2010-2005+1)) , nrow=(2010-2005+1) , ncol = 4)
# Only for 2005-2010 in version 2
colnames(output) = c("year", "#linkages w missing MSOA", "perc A non-zero", "dim(A)")
output = as.data.frame(output)


##############################
### FUNCTION USED IN FOR LOOP
##############################

# removes ".co.uk" from end of host name, and then takes all the characters after the last first stop
get_host = function(string){
  word(gsub(".co.uk", "", string), -1, sep=fixed("."))
}

### For each year: 

for (year in 2005:2010){

    filepath = paste0(year, "-couk-couk-linkage.tsv")
    couk <- read.delim(filepath, header=FALSE, sep=",")
    
    # 22/4/22 added this in as a safety net
    couk <- couk[ grep( "\\.co\\.uk.*\\|.*\\.co\\.uk", couk$V1) , ]
    
    # remove | delim, order cols as desired
    couk <- cbind(couk, str_split_fixed(couk$V1, "[|]", 3))
    couk <- couk[-1]
    couk <- couk[,c("1", "2", "3", "V2")]
    
    # rename columns
    names(couk)[1] <- "year"
    names(couk)[2] <- "host1"
    names(couk)[3] <- "host2"
    names(couk)[4] <- "hyperlinks"
    
    # added 22/4/22 to format hyperlinks column as numeric not string
    couk$hyperlinks = as.numeric(couk$hyperlinks)
    
    ### VERSION 2 EDITS
    couk <- couk[couk$host1 != couk$host2, ]
    
    # apply the function "get_host" to the dataframe to add two new columns
    # couk <- couk %>%
    #   mutate(host1short = get_host(host1)) %>%
    #   mutate(host2short = get_host(host2))
    
    couk$host1short <- get_host(couk$host1)
    couk$host2short <- get_host(couk$host2)
    na.omit(couk)
    
    couk$host1short = tolower(couk$host1short)
    couk$host2short = tolower(couk$host2short)
    
    couk <- couk[couk$host1short != couk$host2short , ]
      
    couk <- couk[, -c(5,6)]
    
    # JOIN LINKAGE DATA AND HOST LOOKUP (POSTCODE) TABLE
    # If host1 has a postcode it is kept
    couk_host1 <- merge(x = couk , y = lookup_all_hosts_1, by.x = "host1", by.y="host", all = FALSE)
    # all = FALSE does an inner join
    # If host2 has a  postcode, and the host1 has already been identified as having a postcode
    all_couk_hostsboth <- merge(x = couk_host1, y = lookup_all_hosts_1, by.x = "host2", by.y="host", all = FALSE)
    
    # rename columns
    names(all_couk_hostsboth)[6] <- "pc2"
    names(all_couk_hostsboth)[5] <- "pc1"
    
    # ADD MSOA INFORMATION TO THE DATA
    
    # FIRST HOST
    df <- merge(x = all_couk_hostsboth , y = all_MSOA, by.x = "pc1", by.y="pcds", all.x = TRUE)
    # all.x = TRUE does a left outer join, so we add in MSOA to the host1
    # rename column added of msoa for host1
    names(df)[7] <- "msoa1"
    
    ### SECOND HOST
    df <- merge(x = df , y = all_MSOA, by.x = "pc2", by.y="pcds", all.x = TRUE)
    # rename column added of msoa for host1
    names(df)[8] <- "msoa2"
    
    # REMOVE OBSERVATIONS THAT HAVE A HOST PC WITH NO MSOA ASSOCIATED (UNLIKELY TO HAPPEN)
    df <- df %>% drop_na(msoa1, msoa2)
    
    # BUILD ADJACENCY MATRIX
    
    # EDGES
    # want to get the number of edges between nodes, equal to the number of hyperlinks between nodes
    V = dim(nodes)[1]
    LINKS = dim(df)[1]
    
    # adjacency Matrix to fill with values
    A = Matrix(0, nrow=V, ncol = V)  
    
    counter = 0 
    
    for (i in 1:LINKS){
      
      # nodes the ith link connects (ith entry in df)
      n1 <- df[i,"msoa1"]
      n2 <- df[i,"msoa2"]
      
      #hyperlinks
      h <- as.integer(df[i,"hyperlinks"])
      
      # index values for the two hosts in the ith edge
      i1 <- as.integer(nodes[nodes$msoa == n1, ]["index"])#linker
      i2 <- as.integer(nodes[nodes$msoa == n2, ]["index"])#linked to
      
      if(is.na(i1) || is.na(i2)){
        next
      }
      else{
        A[i1,i2] = A[i1,i2] + h
        counter = counter + 1
        
        #if (i1 != i2){
          #A[i2,i1] = A[i2, i1] + h # symmetric A, and not on diag
        #}
      }
    }
    
    # percentage of non-zero entries
    perc_non_zero_A <- 100*length(A@x)/V^2
    
    # fill up output table
    output[(year-1995), "year"] = year
    output[(year-1995), "#linkages w missing MSOA"] = LINKS - counter
    output[(year-1995), "perc A non-zero"] = perc_non_zero_A
    output[(year-1995), "dim(A)"] = dim(A)[1]
    
    # save matrix A
    A_filepath_name = paste0("A_for_", year, "_v2.mtx")
    writeMM(obj = A, file=A_filepath_name)
    
    rm(list=c("couk", "df","LINKS", "A", "counter", "couk_host1", "all_couk_hostsboth")) # clear some memory


} # end of year for loop

print(output)
write.csv(output, file="2005-2010_makeAdjMatr_v2_info.csv")
