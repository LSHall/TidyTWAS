#!/usr/bin/Rscript
# This script was written by Lynsey Hall
start.time <- Sys.time()
suppressMessages(library("optparse"))

option_list = list(
make_option("--twas", action="store", default=NA, type='character',
		help="Path to genome-wide TWAS results [required]"),
make_option("--phen1", action="store", default=NA, type='character',
		help="Name of GWAS phenotype used in TWAS [required]"),
make_option("--phen2", action="store", default=NA, type='character',
		help="Name of GWAS phenotype used in TWAS [optional]"),	
make_option("--phen3", action="store", default=NA, type='character',
		help="Name of GWAS phenotype used in TWAS [optional]"),
make_option("--phen4", action="store", default=NA, type='character',
		help="Name of GWAS phenotype used in TWAS [optional]"),
make_option("--phen5", action="store", default=NA, type='character',
		help="Name of GWAS phenotype used in TWAS [optional]"),		
make_option("--build", action="store", default=37, type='numeric',
		help="Name of genome build used in TWAS, 37 or 38 [required]"),
)		

opt = parse_args(OptionParser(option_list=option_list))

# Load libraries
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(xlsx))
suppressMessages(library(purrr))
suppressMessages(library(plyr))
suppressMessages(library(stringr))
suppressMessages(library(biomaRt))
suppressMessages(library(sqldf))


# import mart from biomaRt
# if genome build GRCh37 import GRCH37 mart from biomaRt
if(opt$build == 37){
	ensembl<-useMart(biomart = "ENSEMBL_MART_ENSEMBL", host = "grch37.ensembl.org", path = "/biomart/martservice", dataset = "hsapiens_gene_ensembl")   
	} else {
    # otherwise load GRCh38
    ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
		}

# set working directory to folder with TWAS results
setwd(paste0((opt$twas)))

# establish how many phenotypes have been assigned
nphens<-c(opt$phen1,opt$phen2,opt$phen3,opt$phen4,opt$phen5)
nphens<-nphens[!is.na(nphens)]

# For each phenotype listed
# for(i in nphens){
# Read in all TWAS files
filenames.phen1 <- list.files(paste0((opt$phen1)))
names.phen1<-gsub("Brain_|.allChr.txt", "", filenames.phen1)

# obtain number of characters required for substring(tissue, N) ;  MDD = 4, cog = 10, anh = 10, PD = 3)
phen1.nchars<-nchar(paste0(opt$phen1))+1

# read the TWAS results for each phenotype ~ tissue combination into the workspace
for(i in names.phen1){
  tissue<-paste(i,".allChr",sep="")
filepath <- file.path(paste0(opt$phen1,'/',dir(paste0(opt$phen1),pattern = substring(tissue, phen1.nchars))))
  assign(i, read.delim(filepath,
                       colClasses=c(rep("character",3),rep("numeric",4),"character","numeric","character",rep("numeric",5),"character",rep("numeric",4)),
                       sep = "\t"))
  }

# Add in phenotype column to each dataframe before merging into list.


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



