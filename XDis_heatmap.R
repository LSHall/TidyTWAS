## Create heatmaps from XDis_multiphen_heatmap_setup.txt
# load requisite packages
if (!require("gplots")) {
     install.packages("gplots", dependencies = TRUE)
     library(gplots)
 }
if (!require("RColorBrewer")) {
     install.packages("RColorBrewer", dependencies = TRUE)
     library(RColorBrewer)
 }

# Read in gene/transcript ID linkers
genes<-read.xlsx("Manuscripts/FOETAL/TWAS_sig_multiphen.xlsx",sheetIndex = 1,colIndex = 1:2,colClasses = "character",header = T,stringsAsFactors=F)
trans<-read.xlsx("Manuscripts/FOETAL/TWAS_sig_multiphen.xlsx",sheetIndex = 2,colIndex = 1:3,colClasses = "character",header = T,stringsAsFactors=F)

# read in sig genes 
data<-read.table("Manuscripts/FOETAL/multiphen_sig_gene_Z.txt",header=T,stringsAsFactors=F)
# make gene IDs rnames
rnames <- data[,1]                           
# retain Z scores, round to 2 d.p and convert to matrix format
mat_data <- round(data.matrix(data[,2:ncol(data)]),2)
# assign gene IDs as row names
rownames(mat_data) <- rnames  
# create colour paletter  
my_palette <- colorRampPalette(c("dodgerblue2", "white", "red"))(n = 299)
# plot heatmap
heatmap.2(mat_data, cellnote = mat_data, main = "",  notecol="black", density.info="none", trace="none", col=my_palette, dendrogram="none", Colv="NA", Rowv="NA", xlab = "Trait",ylab = "Gene ID",key.xlab = "TWAS Z-score",margins = c(6,13),key.title = "",notecex=1.4,keysize = 1,cexCol = 1.5) 
# add title
title("Gene-level TWAS Z-scores across traits", line = -1, cex.main=2)

# repeat for transcript level data
data<-read.table("Manuscripts/FOETAL/multiphen_sig_transcript_Z.txt",header=T,stringsAsFactors=F)
rnames <- data[,1]
mat_data <- round(data.matrix(data[,2:ncol(data)]),2)
rownames(mat_data) <- rnames    
heatmap.2(mat_data, cellnote = mat_data, main = "",  notecol="black", density.info="none", trace="none", col=my_palette, dendrogram="none", Colv="NA", Rowv="NA", xlab = "Trait",ylab = "Transcript ID",key.xlab = "TWAS Z-score",margins = c(6,13),key.title = "",notecex=1.4,keysize = 1,cexCol = 1.5) 
title("Transcript-level TWAS Z-scores across traits", line = -1, cex.main=2)