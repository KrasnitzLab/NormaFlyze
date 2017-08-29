d1 <- read.table("/mnt/wigclust1/data/safe/kostic/bin_mapping/hybrid_bin_boundaries_sorted_125_600.txt", sep="\t", header=F, as.is=T, stringsAsFactors=F)

humanEnd <- max(which(d1[,1]=="chrY"))
print(d1[humanEnd,])
humanSum <- sum(d1[1:humanEnd,6])
print(d1[5000,6])
flySum <- sum(humanEnd+1:5000,6)
print(flySum/humanSum)