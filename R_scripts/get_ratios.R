# files1 <- list.files(path="/mnt/wigclust1/data/safe/kostic/SNS_data", pattern="*hg19dm6.varbin.data.txt", full.names=T, recursive=FALSE)
# files2 <- list.files(path="/mnt/wigclust1/data/safe/kostic/SNS_data_2", pattern="*hg19dm6.varbin.data.txt", full.names=T, recursive=FALSE)

d1 <- read.table("/mnt/wigclust1/data/safe/kostic/SNS_data/uber_varbin_seg_quantal.txt", sep="\t", header=T, as.is=T, stringsAsFactors=F)
d2 <- read.table("/mnt/wigclust1/data/safe/kostic/SNS_data_2/uber_varbin_seg_quantal.txt", sep="\t", header=T, as.is=T, stringsAsFactors=F)


dt <- merge(d1, d2)

toDrop <- c("BEI014304", "BEI014307", "BEI014308", "BEI014311", "BEI014312", "BEI014315", "BEI014316", "BEI014317", "Fly_1cell_G7", "Fly_1cell_G8", "SKBR3_Fly_F2", "SKBR3_Fly_F3")
t <- dt[,!colnames(dt)%in%toDrop]

# files <- c(files1, files2)
cat("numcols ", ncol(t))

end_hg_auto <- max(which(t[,1]=="chr22"))
all_hg_auto <- t[0:end_hg_auto,]
vec <- as.vector(as.matrix(all_hg_auto[,2:ncol(t)]))
med_hg <- median(vec)

start2L <- min(which(t[,1]=="chr2L_dm"))
start2R <- min(which(t[,1]=="chr2R_dm"))
end2R <- max(which(t[,1]=="chr2R_dm"))
start3L <- min(which(t[,1]=="chr3L_dm"))
start3R <- min(which(t[,1]=="chr3R_dm"))
end3R <- max(which(t[,1]=="chr3R_dm"))

dm_2 <- t[start2L:end2R,]
dm_2L <- t[start2L:start2R-1,]
dm_2R <- t[start2R:end2R,]

dm_3 <- t[start3L:end3R,]
dm_3L <- t[start3L:start3R-1,]
dm_3R <- t[start3R:end3R,]

all_dm_auto <- rbind(dm_2, dm_3)

vec2 <- as.vector(as.matrix(dm_2[,2:ncol(t)]))
med2 <- median(vec2)
print(mean(vec2))

vec3 <- as.vector(as.matrix(dm_3[,2:ncol(t)]))
med3 <- median(vec3)
print(mean(vec3))

vec <- as.vector(as.matrix(all_dm_auto[,2:ncol(t)]))
med23 <- median(vec)


out <- data.frame(median(as.vector(as.matrix(dm_2L[,2:ncol(t)]))), median(as.vector(as.matrix(dm_2R[,2:ncol(t)]))), med2, median(as.vector(as.matrix(dm_3L[,2:ncol(t)]))), median(as.vector(as.matrix(dm_3R[,2:ncol(t)]))), med3, med_hg)
colnames(out) <- c("dm_2L", "dm_2R", "dm_2", "dm_3L", "dm_3R", "dm_3", "hg_auto")

row2 <- data.frame(out$dm_2L / out$hg_auto, out$dm_2R / out$hg_auto, out$dm_2 / out$hg_auto, out$dm_3L / out$hg_auto, out$dm_3R / out$hg_auto, out$dm_3 / out$hg_auto, out$hg_auto / out$hg_auto)
colnames(row2) <- c("dm_2L", "dm_2R", "dm_2", "dm_3L", "dm_3R", "dm_3", "hg_auto")

final <- rbind(out, row2)
rownames(final) <- c("medians", "ratios")

write.table(final, "/mnt/wigclust1/data/safe/kostic/dm_hg_data.txt", sep="\t", quote=F, row.names=T, col.names=T)


