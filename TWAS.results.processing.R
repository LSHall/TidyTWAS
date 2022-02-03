# TWAS results
# Read in all FUSION output files
library(tidyverse)
library(data.table)
library(xlsx)
library(purrr)
library(plyr)
library(stringr)
library(biomaRt)
library(sqldf)
# import mart from biomaRt
ensembl.GRCh37<-useMart(biomart = "ENSEMBL_MART_ENSEMBL", host = "grch37.ensembl.org", path = "/biomart/martservice", dataset = "hsapiens_gene_ensembl")

# specific TWAS phenotype of interest
phenotype<-"PD"

# Set working directory for TWAS phenotype
for(phen in phenotype){
setwd(paste0("~/Bristol/Results/TWAS/",phen))
}

# Read in all TWAS files
filenames <- list.files()
names<-gsub("Brain_|.allChr.txt", "", filenames)

# obtain number of characters required for substring(tissue, N) ;  MDD = 4, cog = 10, anh = 10, PD = 3)
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


# order all dataframes by their TWAS P value column
l.df.sorted <- lapply(l.df, function(df){
  df[order(df$TWAS.P),]
})

# create bonferroni P-value column by correcting for the number of non-missing TWAS P values
l.df.sorted.bonf <- lapply(l.df.sorted, function(df){
  cbind(df, P.bonf = p.adjust(p = df$TWAS.P,method = "bonferroni",n = length(df$TWAS.P[which(!is.na(df$TWAS.P))])))
})

# extract rows with bonferroni significant P-values to a single dataframe
l.df.bonf.only <- lapply(l.df.sorted.bonf, function(df){
  df[(df$P.bonf<=0.05),]
})

TWAS <- ldply(l.df.bonf.only, data.frame)

# includes rows with NA, so remove those rows
TWAS<-TWAS[!is.na(TWAS$TWAS.P),]

# for today, remove rows pertaining to blood tissue
#TWAS<- TWAS[!grepl("BLOOD", TWAS$FILE),]

# and fill in PANEL column for CMC and PEC data manually
TWAS$PANEL <- ifelse(grepl("CMC", TWAS$.id), "CMC", TWAS$PANEL)
TWAS$PANEL <- ifelse(grepl("PEC", TWAS$.id), "PEC", TWAS$PANEL)
TWAS$PANEL <- ifelse(grepl("NTR", TWAS$.id), "NTR", TWAS$PANEL)
TWAS$PANEL <- ifelse(grepl("YFS", TWAS$.id), "YFS", TWAS$PANEL)

# add in ensembl gene IDs for tissues with hgnc IDs
ensembl.lookup<-unique(TWAS$ID[!grepl(pattern="ENSG",x = TWAS$ID)])
link2ensembl<-getBM(attributes = c('ensembl_gene_id','hgnc_symbol','entrezgene_description'), 
                    filters = 'hgnc_symbol', 
                    values = ensembl.lookup, 
                    mart = ensembl.GRCh37)

# collapse instances of multiple ensembl IDs (or gene descriptions) so only one row per hgnc gene ID
link2ensembl<-sqldf("select ensembl_gene_id, hgnc_symbol, group_concat(entrezgene_description) entrezgene_description from link2ensembl group by ensembl_gene_id, hgnc_symbol", method = "raw")


# create ensembl gene ID column
TWAS<-add_column(.data = TWAS, ensembl_gene_id = if_else(grepl("ENSG",TWAS$ID),TWAS$ID,link2ensembl$ensembl_gene_id[match(TWAS$ID,link2ensembl$hgnc_symbol)]),.before = "CHR")
  
# add in hgnc gene IDs for tissues with ensembl IDs
hgnc.lookup<-unique(TWAS$ID[grepl(pattern="ENSG",x = TWAS$ID)])
link2hgnc<-getBM(attributes = c('ensembl_gene_id','hgnc_symbol','entrezgene_description'), 
                 filters = 'ensembl_gene_id', 
                 values = hgnc.lookup, 
                 mart = ensembl.GRCh37)

link2hgnc<-sqldf("select ensembl_gene_id, hgnc_symbol, group_concat(entrezgene_description) entrezgene_description from link2hgnc group by ensembl_gene_id, hgnc_symbol", method = "raw")

# combine gene info
geneinfolink<-rbind(link2ensembl,link2hgnc)

# collapse instances of multiple hgnc IDs so only one row per ensembl gene ID

# update hgnc symbol for PEC genes which had ensembl gene ids in that column
TWAS$ID = if_else(grepl("ENSG",TWAS$ID),link2hgnc$hgnc_symbol[match(TWAS$ID,link2hgnc$ensembl_gene_id)],TWAS$ID)

# add gene description to results
TWAS<-add_column(.data = TWAS, entrezgene_description = geneinfolink$entrezgene_description[match(TWAS$ID,geneinfolink$hgnc_symbol)], .after = "TWAS.P")



