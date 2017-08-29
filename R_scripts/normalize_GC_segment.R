#### Use default R = 3.3 for DNAcopy_1.50.1

library("DNAcopy", lib.loc="/mnt/wigclust5/data/safe/kendall/DNAcopy_1.50.1")

lowess.gc <- function(jtkx, jtky) {
        jtklow <- lowess(jtkx, log(jtky), f=0.05)
        jtkz <- approx(jtklow$x, jtklow$y, jtkx)
        return(exp(log(jtky) - jtkz$y))
}
lowess.norm <- function(jtkx, jtky) {
        # jtklow <- lowess(jtkx, log(jtky), f=0.05)
        # jtkz <- approx(jtklow$x, jtklow$y, jtkx)
        # return(exp(log(jtky) - jtkz$y))
        loweights<-c(rep(1,4850),rep(0,150))
		pseudo<-0.001
		logbad <-log(jtky+pseudo)
		lfit<-loess(logbad~jtkx,weights=loweights,span=0.05,degree=1)
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

	thisUberRatio <- thisUber
	thisUberRatioQuantal <- thisUber
	thisUberSeg <- thisUber
	thisUberSegQuantal <- thisUber

#######################

	chroms <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31")
	hg_means <- matrix(nrow=ncol(thisUber), ncol=length(chroms))
	colnames(hg_means) <- chroms
	rownames(hg_means) <- colnames(thisUber)
	hg_meds <- hg_means
	seg_means <- hg_means
	seg_meds <- hg_means

########################

	for (i in 1:ncol(thisUber)) {
		sample.name <- dimnames(thisUber)[[2]][i]
		cat(i, dimnames(thisUber)[[2]][i], "\n")
		
		#get the column at i
		thisBincount <- thisUber[, i] + 1
		thisRatio <- thisBincount / mean(thisBincount)
		####normalize for GC content and median bin fragment length
		thisLowratio_1 <- lowess.norm(lengths$median.len.ratio, thisRatio)
		thisLowratio <- lowess.gc(gc$gc.content, thisLowratio_1)
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
		# png(paste(outdir, "/", sample.name, ".5k.len.norm.png", sep=""), height=800, width=1200)
		# plot(x=abspos.col, y=thisLowratio, ylim=c(0,10), main=paste(sample.name, ""), xaxt="n", xlab="Genome Position Gb", ylab="Ratio", col="#CCCCCC")
		# axis(1, at=x.at, labels=x.labels)
		# lines(x=abspos.col, y=thisLowratio, col="#CCCCCC")
		# points(x=abspos.col, y=thisUberSeg[, i], col="#0000AA")
		# lines(x=abspos.col, y=thisUberSeg[, i], col="#0000AA")
		# abline(h=hlines)
		# abline(v=vlines)
		# mtext(chr.text, at = chr.at)
		# dev.off()

		############
		startL <- min(which(chrom.numeric=="25"))
		flySum <- sum(thisUber[startL:5000,i])
		humanSum <- sum(thisUber[1:startL-1, i])
		cat("fly/human number of reads: ", flySum/humanSum, "\n")

		flySum <- sum(thisLowratio_1[startL:5000])
		humanSum <- sum(thisLowratio_1[1:startL-1])
		cat("GC norm fly/human number of reads: ", flySum/humanSum, "\n")

		flySum <- sum(thisLowratio[startL:5000])
		humanSum <- sum(thisLowratio[1:startL-1])
		cat("len norm fly/human number of reads: ", flySum/humanSum, "\n")
		
		humanFragLen <- mean(lengths[1:startL-1,5])
		flyFragLen <- mean(lengths[startL:5000,5])
		cat("mean len human: ", humanFragLen, " mean len fly: ", flyFragLen,"\n")
		##########

		
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

		########################
		#get the means		
		startLoc = 1
		for (k in 1:length(chroms)){
			end <- max(which(chrom.numeric==chroms[k]))
			#mean and median ratio for each chromosme using unsegmented data
			hg_means[i,k] <- mean(thisLowratioQuantal[startLoc:end]) 
			hg_meds[i,k] <- median(thisLowratioQuantal[startLoc:end]) 
			#mean and median ratio for each chrom using segmented data
			seg_means[i,k] <- mean(thisSegQuantal[startLoc:end]) 
			seg_meds[i,k] <- median(thisSegQuantal[startLoc:end]) 
			startLoc = end + 1
		}
		########################

		thisUberRatioQuantal[, i] <- thisLowratioQuantal
		thisUberSegQuantal[, i] <- thisSegQuantal
	
		hlines <- c(0, 1, 2, 3, 4, 5, 6)
		
		png(paste(outdir, "/", sample.name, ".5k.len.quantal.png", sep=""), height=800, width=1200)
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
##############
	
	# write.table(hg_means, file = paste(outdir, "/", "means_allchroms_125_600.txt", sep=""), quote=F, sep ="\t", row.names=T, col.names=T)
	# write.table(hg_meds, file = paste(outdir, "/", "medians_allchroms_125_600.txt", sep=""), quote=F, sep ="\t", row.names=T, col.names=T)
	# write.table(seg_means, file = paste(outdir, "/", "means_seg_allchroms_125_600.txt", sep=""), quote=F, sep ="\t", row.names=T, col.names=T)
	# write.table(seg_meds, file = paste(outdir, "/", "medians_seg_allchroms_125_600.txt", sep=""), quote=F, sep ="\t", row.names=T, col.names=T)

##############	
	return(list(thisUberRatio, thisUberRatioQuantal, thisUberSeg, thisUberSegQuantal))
}

#use this just to get the column containing absolute start positions of bins
a.varbin <- read.table("/mnt/wigclust1/data/safe/kostic/SNS_data_2/b15_5_varbin_count.txt", sep="\t", header=F, as.is=T, stringsAsFactors=F)
#####change file to one with correct date
uber.ratio.nla3.88 <- read.table("/mnt/wigclust1/data/safe/kostic/SNS_data_2/range125_600_uber_varbin_count_data.txt", sep="\t", header=T, as.is=T, stringsAsFactors=F)
nla3.88.count.mat <- as.matrix(uber.ratio.nla3.88)
#######change to correct dated files
nla3.88.results <- cbs.segment.uber01.hg19dm6(outdir="/mnt/wigclust1/data/safe/kostic/SNS_data_2", indir="", varbin.gc="/mnt/wigclust1/data/safe/kostic/bin_mapping/range125_600_GC.txt", varbin.len = "/mnt/wigclust1/data/safe/kostic/bin_mapping/range125_600_LEN.txt", thisUber=nla3.88.count.mat, abspos.col=a.varbin$V3, alpha=0.02, nperm=1000, undo.SD=0.5, min.width=3)
nla3.88.ratio.mat <- nla3.88.results[[1]]
nla3.88.ratio.quantal.mat <- nla3.88.results[[2]]
nla3.88.seg.mat <- nla3.88.results[[3]]
nla3.88.seg.quantal.mat <- nla3.88.results[[4]]
# write.table(nla3.88.ratio.mat, "/mnt/wigclust1/data/safe/kostic/SNS_data/uber.bt1.vt.r1.varbin.5k.ratio.R.txt", quote=F, row.names=F, sep="\t")
# write.table(nla3.88.ratio.quantal.mat, "/mnt/wigclust1/data/safe/kostic/SNS_data/uber.bt1.vt.r1.varbin.5k.ratio.quantal.R.txt", quote=F, row.names=F, sep="\t")
# write.table(nla3.88.seg.mat, "/mnt/wigclust1/data/safe/kostic/SNS_data/uber.bt1.vt.r1.varbin.5k.seg.R.txt", quote=F, row.names=F, sep="\t")
dt <- data.frame(a.varbin[, 1], nla3.88.seg.quantal.mat)
write.table(dt, "/mnt/wigclust1/data/safe/kostic/SNS_data_2/uber_varbin_seg_quantal.txt", quote=F, row.names=F, sep="\t")



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

cell.colors <- tree.colors
x1 <- a.varbin$V1
x2 <- c(x1[-1], x1[length(x1)])
chrombreaks <- c(0, which(x1 != x2), length(x1))
chromlabs <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")
chromlabspos <- chrombreaks[-length(chrombreaks)] + c((chrombreaks[-1] - chrombreaks[-length(chrombreaks)]) / 2)
lr <- rep("", nrow(z))
lr[chromlabspos] <- chromlabs

png("/mnt/wigclust1/data/safe/kostic/R_scripts/nla3.88.seg.mat.heatmap01.png", height=1400, width=900)
heatmap.jk02(z, Rowv=NA, main="nla3.88.seg.mat Heatmap", scale="none", labRow=lr, cexRow=1.5, cexCol=1.5, 
breaks=c(0, 0.9, 1.10, 9999), col=c("blue3", "white", "red3"),
ColSideColors=cell.colors,
legend.text=c(""), 
legend.colors=c("white"),
distance.function="manhattan",
clustering.method="ward",
add.expr = abline(h=sort(c(chrombreaks)), lty=c(rep(c(1), 23), 1, 1)) )
dev.off()