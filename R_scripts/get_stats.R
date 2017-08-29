#d1 <- read.table("/mnt/wigclust1/data/safe/kostic/SNS_data/means_allchroms_125_600.txt", sep="\t", header=T, as.is=T, stringsAsFactors=F)
#d2 <- read.table("/mnt/wigclust1/data/safe/kostic/SNS_data_2/means_allchroms_125_600.txt", sep="\t", header=T, as.is=T, stringsAsFactors=F)

d1 <- read.table("/mnt/wigclust1/data/safe/kostic/SNS_data/medians_allchroms_125_600.txt", sep="\t", header=T, as.is=T, stringsAsFactors=F)
d2 <- read.table("/mnt/wigclust1/data/safe/kostic/SNS_data/medians_allchroms_125_600.txt", sep="\t", header=T, as.is=T, stringsAsFactors=F)

s1 <- read.table("/mnt/wigclust1/data/safe/kostic/SNS_data/medians_seg_allchroms_125_600.txt", sep="\t", header=T, as.is=T, stringsAsFactors=F)
s2 <- read.table("/mnt/wigclust1/data/safe/kostic/SNS_data/medians_seg_allchroms_125_600.txt", sep="\t", header=T, as.is=T, stringsAsFactors=F)

#"/mnt/wigclust1/data/safe/kostic/SNS_data/range125_600_uber_varbin_count_data.txt"

data <- rbind(d1, d2)
seg_data <- rbind(s1, s2)

toDrop <- c("BEI014304", "BEI014307", "BEI014308", "BEI014311", "BEI014312", "BEI014315", "BEI014316", "BEI014317", "Fly_1cell_G7", "Fly_1cell_G8", "SKBR3_Fly_F2", "SKBR3_Fly_F3")

datawithoutSingles <- data[!row.names(data)%in%toDrop,]
segNoSingles <- seg_data[!row.names(data)%in%toDrop,]

means <- integer(ncol(datawithoutSingles))
medians <- integer(ncol(datawithoutSingles))

s_means <- integer(ncol(segNoSingles))
s_meds <- integer(ncol(segNoSingles))


for(i in 1:ncol(datawithoutSingles)){
	set <- datawithoutSingles[,i]
	seg_set <- segNoSingles[,i]
	medians[i] <- median(set)
	means[i] <- mean(set)

	s_meds[i] <- median(seg_set)
	s_means[i] <- mean(seg_set)
}

c2L <- as.vector(segNoSingles$X25)
c2R <- as.vector(segNoSingles$X26)
chrom2LR <- c(c2L, c2R)
med2 <- median(chrom2LR)
c3L <- as.vector(segNoSingles$X27)
c3R <- as.vector(segNoSingles$X28)
chrom3LR <- c(c3L, c3R)
med2 <- median(chrom3LR)
cat("chrom 3LR ", chrom3LR)

# print(medians)
# print(means)
df1 <- as.data.frame(t(means))
df2 <- as.data.frame(t(medians))
sm <- as.data.frame(t(s_means))
smd <- as.data.frame(t(s_meds))
df <- rbind(df1, df2, sm, smd)
row.names(df) <- c("CN means", "CN medians", "seg means", "seg medians")
colnames(df) <- colnames(datawithoutSingles)

write.table(df, file = paste("/mnt/wigclust1/data/safe/kostic/chrom_CNs.txt", sep=""), quote=F, sep ="\t", row.names=T, col.names=T)


hg <- df[4,1:22]
dm <- df[4,25:31]
hg_cn_mean <- mean(as.numeric(hg))

dm_ratios <- integer(length(dm))

for(i in 1:ncol(dm)){
	dm_ratios[i] <- dm[i]/hg_cn_mean
}

dms <- as.numeric(dm_ratios)
dm_r <- data.frame(t(dms))
colnames(dm_r) <- colnames(dm)

write.table(dm_r, file = paste("/mnt/wigclust1/data/safe/kostic/dm_hg_seg_med_ratios.txt"), quote=F, sep ="\t", row.names=F, col.names=T)
