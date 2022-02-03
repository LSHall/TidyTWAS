## Read in TWAS results
# requires an empty R workspace as line 39 combines all the dataframes in the workspace into a list. 
# Read in all FUSION output files
library(tidyverse)
library(data.table)
library(xlsx)
library(purrr)
library(plyr)
library(stringr)
library(sqldf)
# specific TWAS phenotype of interest (change to MDD, SCZ, or PD as appropriate)
phenotype<-"MDD"

# Set working directory for TWAS phenotype
for(phen in phenotype){
setwd(paste0("~/Bristol/Results/TWAS/",phen))
}

# Read in all TWAS files
filenames <- list.files()
names<-gsub("Brain_|.allChr.txt", "", filenames)

# obtain number of characters required for substring(tissue, N) ;  MDD = 4, cognition = 10, anhedonia = 10, PD = 3, SCZ = 4)
pheno.nchars<-as.numeric(nchar(str_extract(getwd(), "(?<=/)[^/]*$")))+1

# read the TWAS results for each phenotype ~ tissue combination into the workspace
for(i in names){
  tissue<-paste(i,".allChr",sep="")
  filepath <- file.path(paste(dir(pattern = substring(tissue, pheno.nchars))))
  assign(i, read.delim(filepath,
                       colClasses=c(rep("character",3),rep("numeric",4),"character","numeric","character",rep("numeric",5),"character",rep("numeric",4)),
                       sep = "\t"))
  }

# combine dataframes into a list for ease of data management
l.df <- mget(ls()[sapply(ls(), function(i) class(get(i))) == "data.frame"])

# order all dataframes by Chromosome and P0 (gene start position)
l.df.sorted <- lapply(l.df, function(df){
    df[with(df, order(CHR,P0)),]
})

# create bonferroni P-value column by correcting for the number of non-missing TWAS P values
l.df.sorted.bonf <- lapply(l.df.sorted, function(df){
  cbind(df, P.bonf = p.adjust(p = df$TWAS.P,method = "bonferroni",n = length(df$TWAS.P[which(!is.na(df$TWAS.P))])))
})

# extract rows with bonferroni significant P-values to a single dataframe
l.df.bonf.only <- lapply(l.df.sorted.bonf, function(df){
  df[(df$P.bonf<=0.05),]
})

# includes rows with NA, so remove those rows
l.df.bonf.only.noNA <- lapply(l.df.sorted.bonf, function(df){
    df[(!is.na(df$TWAS.P) & df$P.bonf<=0.05),]
})

# convert to data.frame
TWAS <- ldply(l.df.bonf.only.noNA, data.frame)

# Convert data.frame to data.table
TWAS<-setDT(TWAS)

# I have already checked that all the genes are on the same strand (so P1 is bigger than P0) by 
# checking in the dataframe that has all of the genes in all of the tissues in it
# if all the genes are on the same strand, the table of the ifelse statement should be all 1, or all 0. 
# That means either P1 is always bigger than P0, or P0 is always bigger than P1 
# If there was a mix of 0s and 1s, that would suggest a mix of forward and reverse strands
	# test <- ldply(l.df, data.frame)
	# test$strand<-ifelse(test = (test$P1-test$P0)>0,yes = 1,no = 0)
	# table(test$strand)

# Create 4 columns: KBdistP0P0 KBdistP0P1 KBdistP1P0 and KBdistP1P1
# to work out the distance in base pairs between the two genes start (P0) and stop points (P1) within each tissue
TWAS <- TWAS[, `:=`(KBdistP0P0=abs(P0 - shift(P0, 1L, type="lag")),
                  KBdistP0P1=abs(P1 - shift(P0, 1L, type="lag")),
                  KBdistP1P0=abs(P0 - shift(P1, 1L, type="lag")),
                  KBdistP1P1=abs(P1 - shift(P1, 1L, type="lag"))), keyby=PANEL]
					 
####### adding GWAS loci info to SMR gene data #######

# Create an empty column called locus
TWAS$locus<-NA

# Make row one of locus for each tissue equal to 1
add_locus<-TWAS[, .(ID=ID[1]), by=PANEL]
add_locus$add_locus<-1
TWAS<-left_join(TWAS,add_locus,by=c("PANEL","ID"))
TWAS$locus[is.na(TWAS$locus)] <- TWAS$add_locus[is.na(TWAS$locus)]
TWAS$add_locus<-NULL

# Need 4 logical vectors so that I can ask if ANY of them are true
# This requires the table to be ordered by PANEL, CHR, P0, to make sure that it is by tissue, and not across tissues by genomic positioning
# Create a logical vector stating whether gene (row n+1) and gene (row n) are on the same chromosome and within 500KB
indP0P0 <- (diff(TWAS$CHR) == 0 & TWAS$KBdistP0P0[-1] < 5e+5)
indP0P1 <- (diff(TWAS$CHR) == 0 & TWAS$KBdistP0P1[-1] < 5e+5)
indP1P0 <- (diff(TWAS$CHR) == 0 & TWAS$KBdistP1P0[-1] < 5e+5)
indP1P1 <- (diff(TWAS$CHR) == 0 & TWAS$KBdistP1P1[-1] < 5e+5)

# create vector that sums the indices, and returns a 1 if the sum of the indices is over 1
# as this indicates that the start or end of one of the gene1 is within 500kb or the start or end of gene2
indsum<-data.frame(indP0P0=indP0P0,indP0P1=indP0P1,indP1P0=indP1P0,indP1P1=indP1P1,total=numeric(length(TWAS$KBdistP0P0)-1))
indsum<-data.frame(lapply(indsum, function(x) as.numeric(as.logical(x))))
indsum$total<-rowSums(indsum[,1:4], na.rm=T)
indsum$sameLocus<-ifelse(indsum$total>=1,TRUE,FALSE)
# Populate the locus column so that if ind = T, make locus number the same, if ind = F, add 1
vapply(2:nrow(TWAS), function (k) TWAS$locus[k] <<- TWAS$locus[k-1] + 1 - indsum$sameLocus[k-1], numeric(1))

# Drop any added columns (.id, and KBdistP0P0, KBdistP0P1, KBdistP1P0, KBdistP1P1
TWAS$tissue<-TWAS$PANEL
TWAS<-TWAS[,c("locus","tissue","PANEL","FILE","ID","CHR","P0","P1","HSQ","BEST.GWAS.ID","BEST.GWAS.Z","EQTL.ID","EQTL.R2","EQTL.Z","EQTL.GWAS.Z","NSNP","NWGT","MODEL","MODELCV.R2","MODELCV.PV","TWAS.Z","TWAS.P")]

# Write out loci with more than one gene to file for --input flag in TWAS conditional analysis
cond.loci <- TWAS %>%
     group_by(locus,tissue) %>%
     filter(n() > 1) %>%
	 nest_by()

mapply(function(x, nm) write.table(x, nm,sep = "\t",quote = FALSE,col.names = TRUE,row.names = FALSE),
        cond.loci$data, paste0(phenotype,".",cond.loci$tissue,".loc", cond.loci$locus, ".txt"))