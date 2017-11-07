library("DNAcopy", lib.loc="/mnt/wigclust5/data/safe/kendall/DNAcopy_1.50.1")
# for each gc/len bin i, sum up bin counts across all the genomic bins b --> ni
# for each gc/len bin i, find pi
# p_i = (n_i + 1)/(Sum_j n_j + 100), sumj nj is the sum of counts for all gc/len bins

segment_uber_hg19dm6_quantile <- function(outdir, indir, varbin, thisUber, alpha, nperm, undo.SD, min.width){

	bins <- read.table(varbin, header=F)

	genomic_bins <- nrow(bins)
	quantiles <- nrow(thisUber) / genomic_bins

	#cat("num genomic and quantile bins: ", genomic_bins, quantiles, "\n")

	thisUberSeg <- data.frame(matrix(0, nrow = genomic_bins, ncol = ncol(thisUber)))

	#rep(seq_len(nrow(bins)), each=quantiles),]
	#bins.expanded <- bins[rep(seq_len(nrow(bins)), each=quantiles),]
	#rep(row.names(bins), quantiles),]

	chrom.numeric <- substring(bins[,1], 4)
	chrom.numeric[which(bins[,1] == "chrX")] <- "23"
	chrom.numeric[which(bins[,1] == "chrY")] <- "24"
	chrom.numeric <- as.numeric(chrom.numeric)

	abspos.col <- bins[,3]
	bin.start.pos <- bins[,2]

	end_hg_auto <- quantiles * max(which(chrom.numeric=="22"))

	#print(bins.expanded[1:200,])

	#median difference in seg_ratio - ratio for each well
	mn_diff <- thisUber[1:2,]

	for (i in 1:ncol(thisUber)){
		sample.name <- dimnames(thisUber)[[2]][i]
		cat(i, dimnames(thisUber)[[2]][i], "\n")
		P_j <- 1/quantiles

		N_counts <- numeric(genomic_bins)
		pj_list <- numeric(quantiles)

		for (j in 1:quantiles){
			#find pj using only human data
			bins_j <- thisUber[seq(j, nrow(thisUber), quantiles), i]
			#bins_j <- thisUber[seq(j, nrow(thisUber), quantiles), i]
			nj <- sum(bins_j)
			if(nj >= 200){
				pj <- (nj + 1)/ (sum(thisUber[1: nrow(thisUber),i]) + 100)
				pj_list[j] <- pj	
			}
		}

		################### removal stuff
		nj_zero <- which(pj_list == 0)
		end <- nrow(thisUber)

		print(nj_zero)
		remove <- c()
		for(j in 1:length(nj_zero)){
			remove <- c(remove, seq(nj_zero[j], end, quantiles))
		}

		
		thisUberEdit <- thisUber[-remove,i]
		pj_list <- pj_list[-nj_zero]
		new_quantiles <- length(pj_list)
		#thisUberEdit <- thisUber[,i]
		###############################################

		#for each genomic bin (so each set of 100 consecutive rows)
		for (b in 1:genomic_bins){
			#1-100, 101-200, 201-300
			start <- (b-1)*new_quantiles + 1
			end <- b * new_quantiles
			n_jb <- thisUberEdit[start:end]
			#Nb is a vector of 5000 raw bin counts
			Nb <- sum(P_j * n_jb / pj_list)
			if(Nb == 0){
				cat("count is 0 ", start, end, "\n")
			}
			N_counts[b] <- Nb
		}

		new_end_hg_auto <- max(which(chrom.numeric=="22"))
	
		N_ratio <- N_counts/mean(N_counts)
		print(mean(N_ratio[1:new_end_hg_auto]))
		print(mean(N_ratio[new_end_hg_auto:length(N_counts)]))

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
		print(sum(thisShort$num.mark))
		thisUberSeg[, i] <- m[, 1]
		print(mean(m[1:new_end_hg_auto,1]))
		print(mean(m[new_end_hg_auto:nrow(m),1]))

		diff <- abs(thisUberSeg[, i] - N_ratio)
		mn_diff[1,i] <- median(diff)

		chr <- chrom.numeric
		chr.shift <- c(chr[-1], chr[length(chr)])
		vlines <- c(1, abspos.col[which(chr != chr.shift) + 1], abspos.col[genomic_bins])
		hlines <- c(0.25, 0.5, 1.0, 1.5, 2.0)
		chr.text <- c(1:22, "X", "Y")
		vlines.shift <- c(vlines[-1], 4*10^9)
		chr.at <- vlines + (vlines.shift - vlines) / 2
		x.at <- c(0, 0.5, 1, 1.5, 2, 2.5, 3) * 10^9
		x.labels <- c("0", "0.5", "1.0", "1.5", "2.0", "2.5", "3.0")


		thisGrid <- seq(0.5, 5.5, by=0.05)
		thisOuter <- thisUberSeg[, i] %o% thisGrid
		thisOuterRound <- round(thisOuter)
		thisOuterDiff <- (thisOuter - thisOuterRound) ^ 2
		thisOuterColsums <- colSums(thisOuterDiff, na.rm = FALSE, dims = 1)
		thisMultiplier <- thisGrid[which.min(thisOuterColsums)]
		cat("multiplier", thisMultiplier, "\n")
		thisError <- min(thisOuterColsums)

		thisLowratioQuantal <- t(N_ratio) * thisMultiplier
		thisSegQuantal <- thisUberSeg[, i] * thisMultiplier

		hlines <- c(0, 1, 2, 3, 4, 5, 6)

		png(paste(outdir, "/", sample.name, ".5k.quantile.norm.png", sep=""), height=800, width=1200)
		plot(x=abspos.col, y=N_ratio, ylim=c(0,10), main=paste(sample.name, ""), xaxt="n", xlab="Genome Position Gb", ylab="Ratio", col="#CCCCCC")
		axis(1, at=x.at, labels=x.labels)
		lines(x=abspos.col, y=N_ratio, col="#CCCCCC")
		points(x=abspos.col, y=thisUberSeg[, i], col="#0000AA")
		lines(x=abspos.col, y=thisUberSeg[, i], col="#0000AA")
		abline(h=hlines)
		abline(v=vlines)
		mtext(chr.text, at = chr.at)
		dev.off()

		# png(paste(outdir, "/", sample.name, ".5k.quantile.mult.png", sep=""), height=800, width=1200)
		# plot(x=abspos.col, y=thisLowratioQuantal, ylim=c(0,10), main=paste(sample.name, ""), xaxt="n", xlab="Genome Position Gb", ylab="Ratio", col="#CCCCCC")
		# axis(1, at=x.at, labels=x.labels)
		# lines(x=abspos.col, y=thisLowratioQuantal, col="#CCCCCC")
		# points(x=abspos.col, y=thisSegQuantal, col="#0000AA")
		# lines(x=abspos.col, y=thisSegQuantal, col="#0000AA")
		# abline(h=hlines)
		# abline(v=vlines)
		# mtext(chr.text, at = chr.at)
		# dev.off()
	}

	print(mn_diff[1,])
	print(mean(mn_diff[1,]))
	#write.table(mn_diff,file = "/mnt/wigclust1/data/safe/kostic/SNS_data_2/hg19_only/seg_ratio_diff.txt", quote=F, sep = "\t", row.names=F, col.names=T)

}

uber_counts_nla3_quantile <- read.table("/mnt/wigclust1/data/safe/kostic/SNS_data_2/hg19_only/GClen_uber_hg19_count_data.txt", sep="\t", header=T, as.is=T, stringsAsFactors=F)
#GClen_uber_hg19_count_data.txt
uber_counts_mat <- as.matrix(uber_counts_nla3_quantile)
nla3.88.results <- segment_uber_hg19dm6_quantile(outdir="/mnt/wigclust1/data/safe/kostic/SNS_data_2/hg19_only", indir="", varbin="/mnt/wigclust1/data/safe/kostic/bin_mapping/hg19_bin_boundaries_sorted_125_600.txt", thisUber=uber_counts_mat, alpha=0.02, nperm=1000, undo.SD=0.5, min.width=3)
