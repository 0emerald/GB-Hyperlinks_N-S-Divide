# read in each year's linkage tsv file
couk1996 <- read.delim("1996-co-uk-linkage.tsv", header=FALSE)

# add this to grep now
# this filters that the two hosts are both ".co.uk"
couk1996 <- couk1996[ grep( "\\.co\\.uk.*\\|.*\\.co\\.uk", couk1996$V1) , ]

write.csv(couk1996,"1996-couk-couk-linkage.tsv", row.names = FALSE)

# remove data from memory now
rm(list=c("couk1996"))

# 1997
couk1997 <- read.delim("1997-co-uk-linkage.tsv", header=FALSE)
couk1997 <- couk1997[ grep( "\\.co\\.uk.*\\|.*\\.co\\.uk", couk1997$V1) , ]
write.csv(couk1997,"1997-couk-couk-linkage.tsv", row.names = FALSE)
rm(list=c("couk1997"))

# 1998
couk1998 <- read.delim("1998-co-uk-linkage.tsv", header=FALSE)
couk1998 <- couk1998[ grep( "\\.co\\.uk.*\\|.*\\.co\\.uk", couk1998$V1) , ]
write.csv(couk1998,"1998-couk-couk-linkage.tsv", row.names = FALSE)
rm(list=c("couk1998"))

# 1999
couk1999 <- read.delim("1999-co-uk-linkage.tsv", header=FALSE)
couk1999 <- couk1999[ grep( "\\.co\\.uk.*\\|.*\\.co\\.uk", couk1999$V1) , ]
write.csv(couk1999,"1999-couk-couk-linkage.tsv", row.names = FALSE)
rm(list=c("couk1999"))

# 2000
couk2000 <- read.delim("2000-co-uk-linkage.tsv", header=FALSE)
couk2000 <- couk2000[ grep( "\\.co\\.uk.*\\|.*\\.co\\.uk", couk2000$V1) , ]
write.csv(couk2000,"2000-couk-couk-linkage.tsv", row.names = FALSE)
rm(list=c("couk2000"))

# 2001
couk2001 <- read.delim("2001-co-uk-linkage.tsv", header=FALSE)
couk2001 <- couk2001[ grep( "\\.co\\.uk.*\\|.*\\.co\\.uk", couk2001$V1) , ]
write.csv(couk2001,"2001-couk-couk-linkage.tsv", row.names = FALSE)
rm(list=c("couk2001"))

# 2002
couk2002 <- read.delim("2002-co-uk-linkage.tsv", header=FALSE)
couk2002 <- couk2002[ grep( "\\.co\\.uk.*\\|.*\\.co\\.uk", couk2002$V1) , ]
write.csv(couk2002,"2002-couk-couk-linkage.tsv", row.names = FALSE)
rm(list=c("couk2002"))

# 2003
couk2003 <- read.delim("2003-co-uk-linkage.tsv", header=FALSE)
couk2004 <- couk2003[ grep( "\\.co\\.uk.*\\|.*\\.co\\.uk", couk2003$V1) , ]
write.csv(couk2003,"2003-couk-couk-linkage.tsv", row.names = FALSE)
rm(list=c("couk2003"))

# 2005
couk2005 <- read.delim("2005-co-uk-linkage.tsv", header=FALSE)
couk2005 <- couk2005[ grep( "\\.co\\.uk.*\\|.*\\.co\\.uk", couk2005$V1) , ]
write.csv(couk2005,"2005-couk-couk-linkage.tsv", row.names = FALSE)
rm(list=c("couk2005"))

# 2006
couk2006 <- read.delim("2006-co-uk-linkage.tsv", header=FALSE)
couk2006 <- couk2006[ grep( "\\.co\\.uk.*\\|.*\\.co\\.uk", couk2006$V1) , ]
write.csv(couk2006,"2006-couk-couk-linkage.tsv", row.names = FALSE)
rm(list=c("couk2006"))

# 2007
couk2007 <- read.delim("2007-co-uk-linkage.tsv", header=FALSE)
couk2007 <- couk2007[ grep( "\\.co\\.uk.*\\|.*\\.co\\.uk", couk2007$V1) , ]
write.csv(couk2007,"2007-couk-couk-linkage.tsv", row.names = FALSE)
rm(list=c("couk2007"))

# 2008
couk2008 <- read.delim("2008-co-uk-linkage.tsv", header=FALSE)
couk2008 <- couk2008[ grep( "\\.co\\.uk.*\\|.*\\.co\\.uk", couk2008$V1) , ]
write.csv(couk2008,"2008-couk-couk-linkage.tsv", row.names = FALSE)
rm(list=c("couk2008"))

# 2009
couk2009 <- read.delim("2009-co-uk-linkage.tsv", header=FALSE)
couk2009 <- couk2009[ grep( "\\.co\\.uk.*\\|.*\\.co\\.uk", couk2009$V1) , ]
write.csv(couk2009,"2009-couk-couk-linkage.tsv", row.names = FALSE)
rm(list=c("couk2009"))

# 2010
couk2010 <- read.delim("2010-co-uk-linkage.tsv", header=FALSE)
couk2010 <- couk2010[ grep( "\\.co\\.uk.*\\|.*\\.co\\.uk", couk2010$V1) , ]
write.csv(couk2010,"2010-couk-couk-linkage.tsv", row.names = FALSE)
rm(list=c("couk2010"))









