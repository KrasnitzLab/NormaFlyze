#### Use default R = 3.3 for DNAcopy_1.50.1

library("DNAcopy", lib.loc="/mnt/wigclust5/data/safe/kendall/DNAcopy_1.50.1")

lowess.gc <- function(jtkx, jtky) {
        jtklow <- lowess(jtkx, log(jtky), f=0.05)
        jtkz <- approx(jtklow$x, jtklow$y, jtkx)
        return(exp(log(jtky) - jtkz$y))
}

cbs.segment.uber01.hg19dm6 <- function(outdir, indir, varbin.gc, thisUber, abspos.col, alpha, nperm, undo.SD, min.width) {

	gc <- read.table(varbin.gc, header=T)

	chrom.numeric <- substring(gc$bin.chrom, 4)
	chrom.numeric[which(gc$bin.chrom == "chrX")] <- "23"
	chrom.numeric[which(gc$bin.chrom == "chrY")] <- "24"
	chrom.numeric[which(gc$bin.chrom == "dm6chr2L")] <- "25"
	chrom.numeric[which(gc$bin.chrom == "dm6chr2R")] <- "26"
	chrom.numeric[which(gc$bin.chrom == "dm6chr3L")] <- "27"
	chrom.numeric[which(gc$bin.chrom == "dm6chr3R")] <- "28"
	chrom.numeric[which(gc$bin.chrom == "dm6chr4")] <- "29"
	chrom.numeric[which(gc$bin.chrom == "dm6chrX")] <- "30"
	chrom.numeric[which(gc$bin.chrom == "dm6chrY")] <- "31"
	chrom.numeric <- as.numeric(chrom.numeric)

	thisUberRatio <- thisUber
	thisUberRatioQuantal <- thisUber
	thisUberSeg <- thisUber
	thisUberSegQuantal <- thisUber

	for (i in 1:ncol(thisUber)) {
		sample.name <- dimnames(thisUber)[[2]][i]
		cat(i, dimnames(thisUber)[[2]][i], "\n")
		
		thisBincount <- thisUber[, i] + 1
		thisRatio <- thisBincount / mean(thisBincount)
		thisLowratio <- lowess.gc(gc$gc.content, thisRatio)
		
		thisUberRatio[, i] <- thisLowratio
		
		
		cat(min(thisLowratio), max(thisLowratio), mean(thisLowratio), "\n")
		cat(min(chrom.numeric), max(chrom.numeric), mean(chrom.numeric), "\n")
		cat(min(gc$bin.start.chrompos), max(gc$bin.start.chrompos), mean(gc$bin.start.chrompos), "\n")
		
		
		set.seed(25) 
		CNA.object <- CNA(log(thisLowratio), chrom.numeric, gc$bin.start.chrompos, data.type="logratio", sampleid=dimnames(thisUber)[[2]][i]) 
		smoothed.CNA.object <- smooth.CNA(CNA.object) 
		segment.smoothed.CNA.object <- segment(smoothed.CNA.object, alpha=alpha, nperm=nperm, undo.splits="sdundo", undo.SD=undo.SD, min.width=2) 
		thisShort <- segment.smoothed.CNA.object[[2]]

		m <- matrix(data=0, nrow=nrow(thisUber), ncol=1)	
		prevEnd <- 0
		for (j in 1:nrow(thisShort)) {
			thisStart <- prevEnd + 1
			thisEnd <- prevEnd + thisShort$num.mark[j]
			m[thisStart:thisEnd, 1] <- exp(thisShort$seg.mean[j])
			prevEnd = thisEnd
		}

		thisUberSeg[, i] <- m[, 1]
		write.table(thisShort, sep="\t", file=paste(outdir, "/", sample.name, ".hg19dm6.5k.k100.varbin.short.txt", sep=""), quote=F, row.names=F) 

		chr <- chrom.numeric
		chr.shift <- c(chr[-1], chr[length(chr)])
		vlines <- c(1, abspos.col[which(chr != chr.shift) + 1], abspos.col[nrow(thisUber)])
		hlines <- c(0.5, 1.0, 1.5, 2.0)
		chr.text <- c(1:22, "X", "Y", "2L", "2R", "3L", "3R", "4", "X", "Y")
		vlines.shift <- c(vlines[-1], 4*10^9)
		chr.at <- vlines + (vlines.shift - vlines) / 2
		x.at <- c(0, 0.5, 1, 1.5, 2, 2.5, 3) * 10^9
		x.labels <- c("0", "0.5", "1.0", "1.5", "2.0", "2.5", "3.0")

		png(paste(outdir, "/", sample.name, ".5k.wg.png", sep=""), height=800, width=1200)
		plot(x=abspos.col, y=thisLowratio, log="y", main=paste(sample.name, ""), xaxt="n", xlab="Genome Position Gb", ylab="Ratio", col="#CCCCCC")
		axis(1, at=x.at, labels=x.labels)
		lines(x=abspos.col, y=thisLowratio, col="#CCCCCC")
		points(x=abspos.col, y=thisUberSeg[, i], col="#0000AA")
		lines(x=abspos.col, y=thisUberSeg[, i], col="#0000AA")
		abline(h=hlines)
		abline(v=vlines)
		mtext(chr.text, at = chr.at)
		dev.off()
		
		thisGrid <- seq(1.5, 5.5, by=0.05)
		thisOuter <- thisUberSeg[, i] %o% thisGrid
		thisOuterRound <- round(thisOuter)
		thisOuterDiff <- (thisOuter - thisOuterRound) ^ 2
		thisOuterColsums <- colSums(thisOuterDiff, na.rm = FALSE, dims = 1)
		thisMultiplier <- thisGrid[which.min(thisOuterColsums)]
		thisError <- min(thisOuterColsums)
#####		thisShredded <- length(which(thisRatio$seg.mean.LOWESS[which(chrom.numeric < 23)] < 0.1)) / length(which(chrom.numeric < 23))

		thisLowratioQuantal <- thisLowratio * thisMultiplier
		thisSegQuantal <- thisUberSeg[, i] * thisMultiplier

		thisUberRatioQuantal[, i] <- thisLowratioQuantal
		thisUberSegQuantal[, i] <- thisSegQuantal
	
		hlines <- c(1, 2, 3, 4, 5, 6)
		
		png(paste(outdir, "/", sample.name, ".5k.wg.quantal.png", sep=""), height=800, width=1200)
		plot(x=abspos.col, y=thisLowratioQuantal, log="y", main=paste(sample.name, ""), xaxt="n", xlab="Genome Position Gb", ylab="Ratio", col="#CCCCCC")
		axis(1, at=x.at, labels=x.labels)
		lines(x=abspos.col, y=thisLowratioQuantal, col="#CCCCCC")
		points(x=abspos.col, y=thisSegQuantal, col="#0000AA")
		lines(x=abspos.col, y=thisSegQuantal, col="#0000AA")
		abline(h=hlines)
		abline(v=vlines)
		mtext(chr.text, at = chr.at)
		dev.off()

		thisRatioOut <- data.frame(gc[, 1:3], "bincount"=thisBincount, "ratio" = thisRatio, "gc.content" = gc$gc.content, "lowratio" = thisLowratio, "seg.mean.LOWESS" = thisUberSeg[, i], "ratio.quantal" = thisLowratioQuantal, "seg.quantal" = thisSegQuantal)
		write.table(thisRatioOut, sep="\t", file=paste(outdir, "/", sample.name, ".hg19dm6.5k.k100.varbin.data.txt", sep=""), quote=F, row.names=F)
		
	}
	
	return(list(thisUberRatio, thisUberRatioQuantal, thisUberSeg, thisUberSegQuantal))
}


#####  ALEX hg19dm6
a.varbin <- read.table("/mnt/wigclust5/data/safe/kendall/nla3_88/b2/BEI014298.bt1.vt.r2.varbin.5k.txt", sep="\t", header=F, as.is=T, stringsAsFactors=F)
uber.ratio.nla3.88 <- read.table("/mnt/wigclust5/data/safe/kendall/nla3_88/uber.bt1.vt.r1.varbin.5k.txt", sep="\t", header=T, as.is=T, stringsAsFactors=F)
nla3.88.count.mat <- as.matrix(uber.ratio.nla3.88)
nla3.88.results <- cbs.segment.uber01.hg19dm6(outdir="/mnt/wigclust1/data/safe/kostic/SNS_data", indir="", varbin.gc="/mnt/wigclust5/data/safe/kendall/sequences/hg19dm6.varbin.gc.content.5k.bowtie.k100.txt", thisUber=nla3.88.count.mat, abspos.col=a.varbin$V3, alpha=0.02, nperm=1000, undo.SD=0.5, min.width=3)
nla3.88.ratio.mat <- nla3.88.results[[1]]
nla3.88.ratio.quantal.mat <- nla3.88.results[[2]]
nla3.88.seg.mat <- nla3.88.results[[3]]
nla3.88.seg.quantal.mat <- nla3.88.results[[4]]
write.table(nla3.88.ratio.mat, "/mnt/wigclust1/data/safe/kostic/SNS_data/uber.bt1.vt.r1.varbin.5k.ratio.R.txt", quote=F, row.names=F, sep="\t")
write.table(nla3.88.ratio.quantal.mat, "/mnt/wigclust1/data/safe/kostic/SNS_data/uber.bt1.vt.r1.varbin.5k.ratio.quantal.R.txt", quote=F, row.names=F, sep="\t")
write.table(nla3.88.seg.mat, "/mnt/wigclust1/data/safe/kostic/SNS_data/uber.bt1.vt.r1.varbin.5k.seg.R.txt", quote=F, row.names=F, sep="\t")
write.table(nla3.88.seg.quantal.mat, "/mnt/wigclust1/data/safe/kostic/SNS_data/uber.bt1.vt.r1.varbin.5k.seg.quantal.R.txt", quote=F, row.names=F, sep="\t")


x <- nla3.88.ratio.mat
x <- nla3.88.seg.mat
mean(x[4790:4959, 1])
mean(x[4790:4959, 2])
mean(x[4790:4959, 3])
mean(x[4790:4959, 4])
mean(x[4790:4959, 5])
mean(x[4790:4959, 6])
mean(x[4790:4959, 7])
mean(x[4790:4959, 8])
mean(x[4790:4959, 9])

source("/mnt/wigstore3/user/kendall/rscripts/heatmap.jk02.R")

z <- nla3.88.seg.mat
tree.colors <- c(rep("white", ncol(z)))

#centromeres.5k <- c(851, 2213, 3917, 5026, 6352, 7701, 8898, 9866, 10868, 11648, 12651, 13467, 14152, 14850, 15478, 16245, 16701, 17186, 17784, 18181, 18425, 18666, 19291)
cell.colors <- tree.colors
x1 <- a.varbin$V1
x2 <- c(x1[-1], x1[length(x1)])
chrombreaks <- c(0, which(x1 != x2), length(x1))
chromlabs <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")
chromlabspos <- chrombreaks[-length(chrombreaks)] + c((chrombreaks[-1] - chrombreaks[-length(chrombreaks)]) / 2)
lr <- rep("", nrow(z))
lr[chromlabspos] <- chromlabs

png("/mnt/wigclust5/data/safe/kendall/nla3_88/nla3.88.seg.mat.heatmap01.png", height=1400, width=900)

heatmap.jk02(z, Rowv=NA, main="nla3.88.seg.mat Heatmap", scale="none", labRow=lr, cexRow=1.5, cexCol=1.5, 
breaks=c(0, 0.9, 1.10, 9999), col=c("blue3", "white", "red3"),
ColSideColors=cell.colors,
legend.text=c(""), 
legend.colors=c("white"),
distance.function="manhattan",
clustering.method="ward",
add.expr = abline(h=sort(c(chrombreaks)), lty=c(rep(c(1), 23), 1, 1)) )
dev.off()

