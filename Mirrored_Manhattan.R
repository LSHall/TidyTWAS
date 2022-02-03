manhattan_mirror = function(dataframe, colors=c("gray10", "gray50"), ymax="max", ymin=ymin, cex.x.axis=1, limitchromosomes=1:23, suggestiveline=-log10(1e-5), genomewideline=-log10(5e-8), annotate=NULL, highlight=NULL, ...) {
d=dataframe
if (!("CHR" %in% names(d) & "BP" %in% names(d) & "Z" %in% names(d))) stop("Make sure your data frame contains columns CHR, BP, and Z")
if (any(limitchromosomes)) d=d[d$CHR %in% limitchromosomes, ]
d=na.omit(d[order(d$CHR, d$BP), ]) # remove na's and sort
d$pos=NA
ticks=NULL
lastbase=0
colors <- rep(colors,max(d$CHR))[1:max(d$CHR)]
if (ymax=="max") ymax<-ceiling(max(abs(d$Z)))
ymin=-ymax
#if (ymax<8) ymax<-8
numchroms=length(unique(d$CHR))
if (numchroms==1) {
d$pos=d$BP
ticks=floor(length(d$pos))/2+1
} else {
for (i in unique(d$CHR)) {
if (i==1) {
d[d$CHR==i, ]$pos=d[d$CHR==i, ]$BP
} else {
lastbase=lastbase+tail(subset(d,CHR==i-1)$BP, 1)
d[d$CHR==i, ]$pos=d[d$CHR==i, ]$BP+lastbase
}
ticks=c(ticks, d[d$CHR==i, ]$pos[floor(length(d[d$CHR==i, ]$pos)/2)+1])
}
}
if (numchroms==1) {
with(d, plot(pos, Z, ylim=c(ymin,ymax), ylab="Z-score", xlab=paste("Physical position",unique(d$CHR),"position"), ...))
} else {
with(d, plot(pos, Z, ylim=c(ymin,ymax), ylab="Z-score", xlab="Physical position", xaxt="n", type="n", ...))
axis(1, at=ticks, lab=unique(d$CHR), cex.axis=cex.x.axis)
icol=1
for (i in unique(d$CHR)) {
with(d[d$CHR==i, ],points(pos, Z, col=colors[icol], ...))
icol=icol+1
}
}
if (!is.null(annotate)) {
d.annotate=d[which(d$SNP %in% annotate), ]
with(d.annotate, points(pos, Z, col="red", ...)) 
}
if (suggestiveline) abline(h=suggestiveline, col="blue")
if (suggestiveline) abline(h=-suggestiveline, col="blue")
if (genomewideline) abline(h=genomewideline, col="red")
if (genomewideline) abline(h=-genomewideline, col="red")
abline(h=0, col="black")
 if (!is.null(highlight)) {
        d.highlight = d[which(d$SNP %in% highlight), ]
        with(d.highlight, textxy(pos, Z, 
                offset = 0.625, labs = d.highlight$SNP, cex = 0.45), 
                ...)
        }
        }