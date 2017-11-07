
library("DNAcopy", lib.loc="/mnt/wigclust5/data/safe/kendall/DNAcopy_1.50.1")

segment_uber_hg19dm6_quantile <- function(outdir, indir, varbin, thisUber, alpha, nperm, undo.SD, min.width){

	bins <- read.table(varbin, header=T)

	genomic_bins <- nrow(bins)
	quantiles <- nrow(thisUber) / genomic_bins

	thisUberSeg <- data.frame(matrix(0, nrow = genomic_bins, ncol = ncol(thisUber)))

	chrom.numeric <- substring(bins[,1], 4)
	chrom.numeric[which(bins[,1] == "chrX")] <- "23"
	chrom.numeric[which(bins[,1] == "chrY")] <- "24"
	chrom.numeric[which(bins[,1] == "chr2L_dm")] <- "25"
	chrom.numeric[which(bins[,1] == "chr2R_dm")] <- "26"
	chrom.numeric[which(bins[,1] == "chr3L_dm")] <- "27"
	chrom.numeric[which(bins[,1] == "chr3R_dm")] <- "28"
	chrom.numeric[which(bins[,1] == "chr4_dm")] <- "29"
	chrom.numeric[which(bins[,1] == "chrX_dm")] <- "30"
	chrom.numeric[which(bins[,1] == "chrY_dm")] <- "31"
	chrom.numeric <- as.numeric(chrom.numeric)

	abspos.col <- bins$bin.start.abspos
	bin.start.pos <- bins$bin.start.chrompos

	end_hg_auto <- quantiles * max(which(chrom.numeric=="22"))

	for (i in 1:ncol(thisUber)){

		sample.name <- dimnames(thisUber)[[2]][i]
		cat(i, dimnames(thisUber)[[2]][i], "\n")
		P_j <- 1/quantiles

		N_counts <- numeric(genomic_bins)
		pj_list <- numeric(quantiles)

		for (j in 1:quantiles){
			#find pj using only human data
			bins_j <- thisUber[seq(j, end_hg_auto, quantiles), i]
			nj <- sum(bins_j)
			if(nj >= 200){
				pj <- (nj + 1)/ (sum(thisUber[1: end_hg_auto,i]) + 2)
				pj_list[j] <- pj
			}
		}

		nj_zero <- which(pj_list == 0)
		end <- nrow(thisUber)

		remove <- c()
		for(j in 1:length(nj_zero)){
			remove <- c(remove, seq(nj_zero[j], end, quantiles))
		}
		
		thisUberEdit <- thisUber[-remove,i]
		pj_list <- pj_list[-nj_zero]
		new_quantiles <- length(pj_list)

		#for each genomic bin (so each set of 100 consecutive rows)
		for (b in 1:genomic_bins){
			#1-100, 101-200, 201-300
			start <- (b-1)*new_quantiles + 1
			end <- b * new_quantiles
			n_jb <- thisUberEdit[start:end]
			#Nb is a vector of 5000 raw bin counts
			Nb <- sum(P_j * n_jb / pj_list)
			if(is.na(Nb)){
				cat("na in nratio\n")
			}
			N_counts[b] <- Nb
		}

		N_ratio <- N_counts / mean(N_counts)
		
		new_end_hg_auto <- max(which(chrom.numeric=="22"))
		start_dm <- min(which(chrom.numeric== "25"))
		end_dm_auto <- max(which(chrom.numeric=="29"))

		set.seed(25) 
		CNA.object <- CNA(N_ratio, chrom.numeric, bin.start.pos, data.type="logratio", sampleid=dimnames(thisUber)[[2]][i]) 
		smoothed.CNA.object <- smooth.CNA(CNA.object) 
		segment.smoothed.CNA.object <- segment(smoothed.CNA.object, alpha=alpha, nperm=nperm, undo.splits="sdundo", undo.SD=undo.SD, min.width=2) 
		thisShort <- segment.smoothed.CNA.object[[2]]

		m <- matrix(data=0, nrow=genomic_bins, ncol=1)	
		prevEnd <- 0
		for (j in 1:nrow(thisShort)) {
			thisStart <- prevEnd + 1
			thisEnd <- prevEnd + thisShort$num.mark[j]
			m[thisStart:thisEnd, 1] <- thisShort$seg.mean[j]
			prevEnd = thisEnd
		}

		thisUberSeg[, i] <- m[, 1]
		diff <- abs(thisUberSeg[, i] - N_ratio)

		chr <- chrom.numeric
		chr.shift <- c(chr[-1], chr[length(chr)])
		vlines <- c(1, abspos.col[which(chr != chr.shift) + 1], abspos.col[genomic_bins])
		hlines <- c(0.25, 0.5, 1.0, 1.5, 2.0)
		chr.text <- c(1:22, "X", "Y", "2L\n\n", "2R\n", "3L", "3R\n\n", "4\n", "X", "Y\n\n")
		vlines.shift <- c(vlines[-1], 4*10^9)
		chr.at <- vlines + (vlines.shift - vlines) / 2
		x.at <- c(0, 0.5, 1, 1.5, 2, 2.5, 3) * 10^9
		x.labels <- c("0", "0.5", "1.0", "1.5", "2.0", "2.5", "3.0")


##################
		# chr <- chrom.numeric
		# chr.shift <- c(chr[-1], chr[length(chr)])
		# vlines <- c(1, abspos.col[which(chr != chr.shift) + 1], abspos.col[nrow(thisUber)])
		# hlines <- c(0.25, 0.5, 1.0, 1.5, 2.0)
		# chr.text <- c(1:22, "X", "Y", "2L\n\n", "2R\n", "3L", "3R\n\n", "4\n", "X", "Y\n\n")
		# vlines.shift <- c(vlines[-1], 4*10^9)
		# chr.at <- vlines + (vlines.shift - vlines) / 2
		# print(chr.at)
		# x.at <- c(0, 0.5, 1, 1.5, 2, 2.5, 3) * 10^9
		# x.labels <- c("0", "0.5", "1.0", "1.5", "2.0", "2.5", "3.0")
###################


		thisGrid <- seq(0.5, 5.5, by=0.05)
		thisOuter <- thisUberSeg[1:new_end_hg_auto, i] %o% thisGrid
		thisOuterRound <- round(thisOuter)
		thisOuterDiff <- (thisOuter - thisOuterRound) ^ 2
		thisOuterColsums <- colSums(thisOuterDiff, na.rm = FALSE, dims = 1)
		thisMultiplier <- thisGrid[which.min(thisOuterColsums)]
		#cat("multiplier", thisMultiplier, "\n")
		thisError <- min(thisOuterColsums)

		thisLowratioQuantal <- t(N_ratio)# * thisMultiplier)
		thisSegQuantal <- thisUberSeg[, i] #* thisMultiplier

		hlines <- c(0, 1, 2, 3, 4, 5, 6)

		# png(paste(outdir, "/", sample.name, ".10by5.quantile.png", sep=""), height=800, width=1200)
		# #.5k.quant.nomult.png
		# plot(x=abspos.col, y=N_ratio, ylim=c(0,10), main=paste(sample.name, ""), xaxt="n", xlab="Genome Position Gb", ylab="Ratio", col="#CCCCCC")
		# axis(1, at=x.at, labels=x.labels)
		# lines(x=abspos.col, y=N_ratio, col="#CCCCCC")
		# points(x=abspos.col, y=thisUberSeg[, i], col="#0000AA")
		# lines(x=abspos.col, y=thisUberSeg[, i], col="#0000AA")
		# abline(h=hlines)
		# abline(v=vlines)
		# mtext(chr.text, at = chr.at)
		# dev.off()

		png(paste(outdir, "/", sample.name, ".10by10.quantile.png", sep=""), height=800, width=1500)
		plot(x=abspos.col, y=thisLowratioQuantal, ylim=c(0,10), main=paste(sample.name, ""), xaxt="n", xlab="Genome Position Gb", ylab="Ratio", col="#CCCCCC")
		axis(1, at=x.at, labels=x.labels)
		lines(x=abspos.col, y=thisLowratioQuantal, col="#CCCCCC")
		points(x=abspos.col, y=thisSegQuantal, col="#0000AA")
		lines(x=abspos.col, y=thisSegQuantal, col="#0000AA")
		abline(h=hlines)
		abline(v=vlines)
		mtext(chr.text, at = chr.at)
		dev.off()
	}
}

df <- read.table("/mnt/wigclust1/data/safe/kostic/SNS_data_2/GClen_uber_10by10_count_data.txt", sep="\t", header=T, as.is=T, stringsAsFactors=F)
uber_counts_nla3_quantile <- df[,-c(12,13,16,17)]
uber_counts_mat <- as.matrix(uber_counts_nla3_quantile)
nla3.88.results <- segment_uber_hg19dm6_quantile(outdir="/mnt/wigclust1/data/safe/kostic/SNS_data_2", indir="", varbin="/mnt/wigclust1/data/safe/kostic/bin_mapping/range125_600_GC.txt", thisUber=uber_counts_mat, alpha=0.02, nperm=1000, undo.SD=0.5, min.width=3)