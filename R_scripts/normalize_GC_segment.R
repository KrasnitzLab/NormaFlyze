#### Use default R = 3.3 for DNAcopy_1.50.1

library("DNAcopy", lib.loc="/mnt/wigclust5/data/safe/kendall/DNAcopy_1.50.1")

lowess.gc <- function(jtkx, jtky) {
        jtklow <- lowess(jtkx, log(jtky), f=0.05)
        jtkz <- approx(jtklow$x, jtklow$y, jtkx)
        return(exp(log(jtky) - jtkz$y))
}
lowess.norm <- function(jtkx, jtky, partition) {
        loweights<-c(rep(1,partition),rep(0,5000-partition))
		pseudo<-0.001
		logbad <-log(jtky+pseudo)
		lfit<-loess(logbad~jtkx,span=0.05,degree=1)
		norm_ratio<-exp(lfit$residuals)
		return(norm_ratio)
}

cbs.segment.uber01.hg19dm6 <- function(outdir, indir, varbin.gc, varbin.len, thisUber, abspos.col, alpha, nperm, undo.SD, min.width) {

	gc <- read.table(varbin.gc, header=T)
	lengths <- read.table(varbin.len, header = T)
	chrom.numeric <- substring(gc$bin.chrom, 4)
	chrom.numeric[which(gc$bin.chrom == "chrX")] <- "23"
	chrom.numeric[which(gc$bin.chrom == "chrY")] <- "24"
	chrom.numeric[which(gc$bin.chrom == "chr2L_dm")] <- "25"
	chrom.numeric[which(gc$bin.chrom == "chr2R_dm")] <- "26"
	chrom.numeric[which(gc$bin.chrom == "chr3L_dm")] <- "27"
	chrom.numeric[which(gc$bin.chrom == "chr3R_dm")] <- "28"
	chrom.numeric[which(gc$bin.chrom == "chr4_dm")] <- "29"
	chrom.numeric[which(gc$bin.chrom == "chrX_dm")] <- "30"
	chrom.numeric[which(gc$bin.chrom == "chrY_dm")] <- "31"
	chrom.numeric <- as.numeric(chrom.numeric)

########################

	flyStart <- min(which(chrom.numeric=="25"))

	thisUberRatio <- thisUber
	thisUberRatioQuantal <- thisUber
	thisUberSeg <- thisUber
	thisUberSegQuantal <- thisUber

#######################

	for (i in 1:ncol(thisUber)) {
		sample.name <- dimnames(thisUber)[[2]][i]
		cat(i, dimnames(thisUber)[[2]][i], "\n")
		
		#get the column at i
		thisBincount <- thisUber[, i] + 1
		thisRatio <- thisBincount / mean(thisBincount)
		####normalize for GC content and median bin fragment length
		thisLowratio_1 <- lowess.gc(gc$gc.content, thisRatio)
		thisLowratio <- lowess.norm(lengths$median.len, thisLowratio_1, flyStart)
		thisUberRatio[, i] <- thisLowratio

		cat("thisLowratio stats: ", min(thisLowratio), max(thisLowratio), mean(thisLowratio), "\n")

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

		cat("thisUberseg stats: ", min(thisUberSeg[, i]), max(thisUberSeg[, i]), mean(thisUberSeg[, i]), "\n")

		chr <- chrom.numeric
		chr.shift <- c(chr[-1], chr[length(chr)])
		vlines <- c(1, abspos.col[which(chr != chr.shift) + 1], abspos.col[nrow(thisUber)])
		hlines <- c(0.25, 0.5, 1.0, 1.5, 2.0)
		chr.text <- c(1:22, "X", "Y", "2L\n\n", "2R\n", "3L", "3R\n\n", "4\n", "X", "Y\n\n")
		vlines.shift <- c(vlines[-1], 4*10^9)
		chr.at <- vlines + (vlines.shift - vlines) / 2
		x.at <- c(0, 0.5, 1, 1.5, 2, 2.5, 3) * 10^9
		x.labels <- c("0", "0.5", "1.0", "1.5", "2.0", "2.5", "3.0")

		# makes the graph of the non-scaled data
		png(paste(outdir, "/", sample.name, ".5k.len.norm.png", sep=""), height=800, width=1200)
		plot(x=abspos.col, y=thisLowratio, ylim=c(0,6), main=paste(sample.name, ""), xaxt="n", xlab="Genome Position Gb", ylab="Ratio", col="#CCCCCC")
		axis(1, at=x.at, labels=x.labels)
		lines(x=abspos.col, y=thisLowratio, col="#CCCCCC")
		points(x=abspos.col, y=thisUberSeg[, i], col="#0000AA")
		lines(x=abspos.col, y=thisUberSeg[, i], col="#0000AA")
		abline(h=hlines)
		abline(v=vlines)
		mtext(chr.text, at = chr.at)
		dev.off()
		
		#get multiplier by minimizing the residuals
		thisGrid <- seq(1.5, 5.5, by=0.05)
		thisOuter <- thisUberSeg[, i] %o% thisGrid
		thisOuterRound <- round(thisOuter)
		thisOuterDiff <- (thisOuter - thisOuterRound) ^ 2
		thisOuterColsums <- colSums(thisOuterDiff, na.rm = FALSE, dims = 1)
		thisMultiplier <- thisGrid[which.min(thisOuterColsums)]
		cat("multiplier", thisMultiplier, "\n")
		thisError <- min(thisOuterColsums)

		thisLowratioQuantal <- thisLowratio * thisMultiplier
		thisSegQuantal <- thisUberSeg[, i] * thisMultiplier

		thisUberRatioQuantal[, i] <- thisLowratioQuantal
		thisUberSegQuantal[, i] <- thisSegQuantal
	
		hlines <- c(0, 1, 2, 3, 4, 5, 6)
		
		png(paste(outdir, "/", sample.name, ".5k.GC_length_loess.png", sep=""), height=800, width=1200)
		plot(x=abspos.col, y=thisLowratioQuantal, ylim=c(0,10), main=paste(sample.name, ""), xaxt="n", xlab="Genome Position Gb", ylab="Ratio", col="#CCCCCC")
		axis(1, at=x.at, labels=x.labels)
		lines(x=abspos.col, y=thisLowratioQuantal, col="#CCCCCC")
		points(x=abspos.col, y=thisSegQuantal, col="#0000AA")
		lines(x=abspos.col, y=thisSegQuantal, col="#0000AA")
		abline(h=hlines)
		abline(v=vlines)
		mtext(chr.text, at = chr.at)
		dev.off()

		thisRatioOut <- data.frame(gc[, 1:3], "bincount"=thisBincount, "ratio" = thisRatio, "gc.content" = gc$gc.content, "lowratio" = thisLowratio, "seg.mean.LOWESS" = thisUberSeg[, i], "ratio.quantal" = thisLowratioQuantal, "seg.quantal" = thisSegQuantal)
		write.table(thisRatioOut, sep="\t", file=paste(outdir, "/", sample.name, ".hg19dm6.varbin.data.txt", sep=""), quote=F, row.names=F)

	}

	return(list(thisUberRatio, thisUberRatioQuantal, thisUberSeg, thisUberSegQuantal))
}

#use this just to get the column containing absolute start positions of bins
a.varbin <- read.table("/mnt/wigclust1/data/safe/kostic/SNS_data_2/b15_5_varbin_count.txt", sep="\t", header=F, as.is=T, stringsAsFactors=F)
#####change file to one with correct date
uber.ratio.nla3.88 <- read.table("/mnt/wigclust1/data/safe/kostic/SNS_data_2/range125_600_uber_varbin_count_data.txt", sep="\t", header=T, as.is=T, stringsAsFactors=F)
nla3.88.count.mat <- as.matrix(uber.ratio.nla3.88)
nla3.88.results <- cbs.segment.uber01.hg19dm6(outdir="/mnt/wigclust1/data/safe/kostic/SNS_data_2", indir="", varbin.gc="/mnt/wigclust1/data/safe/kostic/bin_mapping/range125_600_GC.txt", varbin.len = "/mnt/wigclust1/data/safe/kostic/bin_mapping/range125_600_LEN.txt", thisUber=nla3.88.count.mat, abspos.col=a.varbin$V3, alpha=0.02, nperm=1000, undo.SD=0.5, min.width=3)
