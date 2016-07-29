library(pegas)
library(seqinr)
library(RColorBrewer)
### script to estimate the haplotypes in each population
Vspecies <- "Mammuthus_primigenius"
for(sp in seq_along(Vspecies)){
        temp_hap_DB <- DATABASE[which(DATABASE$Species == Vspecies[sp] &  nchar(DATABASE$Sequence) > 1),]
        temp_points_hap_DB <- as.data.frame(matrix(nrow=nrow(temp_hap_DB), ncol=11))
        colnames(temp_points_hap_DB) <- c("Vspecies[sp]","Longitude", "Latitude", "Time_sample", "Time_bin", "Layer", "Sequence","Temp","Prec","Biome", "Haplotype")
        Vlayer <- c(seq(0, 21000,  by=1000), seq(22000, 50000, by=2000))
        for (row in seq_along(temp_hap_DB[,1])){
                temp_points_hap_DB[row,c(1,2,3,4,5,6,7)] <- c(temp_hap_DB$Species[row], 
                                                              temp_hap_DB$Longitude[row],
                                                              temp_hap_DB$Latitude[row],
                                                              temp_hap_DB$Mean_Age[row],
                                                              min(Vlayer[which(Vlayer > temp_hap_DB$Mean_Age[row])]),
                                                              min(which(Vlayer > temp_hap_DB$Mean_Age[row])),
                                                              temp_hap_DB$Sequence[row])
                temp_points_hap_DB[row,8]<- extract(Tmp, layer=as.numeric(temp_points_hap_DB[row,6]), nl=1, y=matrix(as.numeric( temp_points_hap_DB[row,c(2,3)]), nrow=1, ncol=2))
                temp_points_hap_DB[row,9]<- extract(Prc, layer=as.numeric(temp_points_hap_DB[row,6]), nl=1, y=matrix(as.numeric( temp_points_hap_DB[row,c(2,3)]), nrow=1, ncol=2))
                temp_points_hap_DB[row,10]<- extract(Bio, layer=as.numeric(temp_points_hap_DB[row,6]), nl=1, y=matrix(as.numeric( temp_points_hap_DB[row,c(2,3)]), nrow=1, ncol=2))
                
        }
}
### Clean the data: remove NAs and -1000 values
setwd("/Users/afr/Desktop/PhD/Thesis/Chapter3/Ch3/Ch3_results/")
write.table(temp_points_hap_DB, file=paste(Vspecies, "points.txt", sep="_"), sep = "\t")
### IMPORTANT SAFETY POINT backup_hap <- temp_points_hap_DB
if(sum(is.na(temp_points_hap_DB$Longitude)) > 0){
temp_points_hap_DB <- temp_points_hap_DB[-which(is.na(temp_points_hap_DB$Longitude)),]
}
if(sum(is.na(temp_points_hap_DB$Biome)) > 0){
temp_points_hap_DB <- temp_points_hap_DB[-which(is.na(temp_points_hap_DB$Biome)),]
}
if(sum(temp_points_hap_DB$Biome == -1000) > 0){
temp_points_hap_DB <- temp_points_hap_DB[-which(temp_points_hap_DB$Biome == -1000),]
}
if(sum(temp_points_hap_DB$Time_bin == "Inf") > 0){
temp_points_hap_DB <- temp_points_hap_DB[-which(temp_points_hap_DB$Time_bin == "Inf"),]
}

### Create a DNA object
write.fasta(as.list(temp_points_hap_DB$Sequence), names = paste(seq(1,length.out = length(temp_points_hap_DB$Longitude)), temp_points_hap_DB$Biome, sep="_"), file.out = "Mp_hap_biome")
setwd("/Users/afr/Desktop/PhD/Thesis/Chapter3/Ch3/Ch3_data/Fasta_files/")
temp_align_biome <- read.dna("Mp_align_biome.fasta", format = "fasta")
temp_hap_biome <- haplotype(temp_align_biome)
### remove sequences with lots of Ns (gaps)
strings <- strsplit(dimnames(temp_align_biome)[1][[1]], split="_")
rows_good_seqs  <- as.numeric(unlist(strings)[2*(1:length(dimnames(temp_align_biome)[1][[1]]))-1])
for (seq in 1:dim(temp_align_biome)[1]){
        temp_points_hap_DB$Sequence[rows_good_seqs[seq]] <- paste(as.character(temp_align_biome[1,]), collapse = "")
}
temp_points_hap_DB <- temp_points_hap_DB[rows_good_seqs,]
### assign haplotypes to each sample
hap_index <- 1
for (hap in seq_along(attr(temp_hap_biome, which = "index"))){
        samples <- attr(temp_hap_biome, which = "index")[[hap]]
        temp_points_hap_DB$Haplotype[samples] <- hap_index
        hap_index <- hap_index + 1
}

### Subset the data in time bins
temp_hap_preLGM <- temp_points_hap_DB[temp_points_hap_DB$Time_bin >= 25000,]
temp_hap_LGM <- temp_points_hap_DB[temp_points_hap_DB$Time_bin < 25000 & temp_points_hap_DB$Time_bin >= 15000,]
temp_hap_posLGM <- temp_points_hap_DB[temp_points_hap_DB$Time_bin < 15000,]
### Histograms for haplotypes
Bio_hap_preLGM <- temp_hap_preLGM$Biome
hBio_hap_preLGM <- hist(Bio_hap_preLGM, breaks=seq(3,27,by=1), plot=F)
Bio_hap_LGM <- temp_hap_LGM$Biome
hBio_hap_LGM <- hist(Bio_hap_LGM, breaks=seq(3,27,by=1), plot=F)
Bio_hap_posLGM <- temp_hap_posLGM$Biome
hBio_hap_posLGM <- hist(Bio_hap_posLGM, breaks=seq(3,27,by=1), plot=F)

### Plotting haplotypes
        ### Create color palette
unique(temp_hap_preLGM$Biome)
unique(temp_hap_LGM$Biome)
unique(temp_hap_posLGM$Biome)

plot.new()
par(mar=c(12,0,17,0))
haplo_col <- rainbow(n = 87)
Layout <- layout(matrix(c(0,3,0,2,0,1,0),ncol=7, nrow=1),widths =c(4,10,1,10,1,10,4), heights = c(1,1,1))
plot(NULL, type = "n", xlim = c(0, max(hBio_hap_preLGM$counts)), ylim = c(min(hBio_hap_preLGM$breaks), max(hBio_hap_preLGM$breaks)), axes=FALSE, frame=T, xlab="Pre-LGM", xaxs="i", yaxs="i", cex.lab=2)
for (biome in seq_along(unique(temp_hap_preLGM$Biome))){
        temp_hap_bio_rect <- temp_hap_preLGM[which(temp_hap_preLGM$Biome == unique(temp_hap_preLGM$Biome)[biome]),]
        base <- 0
        for(hap_col in seq_along(temp_hap_bio_rect$Haplotype)){
                rect(xleft = base,xright =  base+1, ytop = temp_hap_bio_rect$Biome[hap_col], ybottom = temp_hap_bio_rect$Biome[hap_col]+1, col=haplo_col[temp_hap_bio_rect$Haplotype[hap_col]])
                base <- base+1
        }
}
axis(side=1)

plot(NULL, type = "n", xlim = c(0, max(hBio_hap_LGM$counts)), ylim = c(min(hBio_hap_LGM$breaks), max(hBio_hap_LGM$breaks)), axes=FALSE, frame=T, xlab="LGM", xaxs="i", yaxs="i", cex.lab=2)
for (biome in seq_along(unique(temp_hap_LGM$Biome))){
        temp_hap_bio_rect <- temp_hap_LGM[which(temp_hap_LGM$Biome == unique(temp_hap_LGM$Biome)[biome]),]
        base <- 0
        for(hap_col in seq_along(temp_hap_bio_rect$Haplotype)){
                rect(xleft = base,xright =  base+1, ytop = temp_hap_bio_rect$Biome[hap_col], ybottom = temp_hap_bio_rect$Biome[hap_col]+1, col=haplo_col[temp_hap_bio_rect$Haplotype[hap_col]])
                base <- base+1
        }
}
axis(side=1)
mtext(Vspecies, side = 3, line=5)

plot(NULL, type = "n", xlim = c(0, max(hBio_hap_posLGM$counts)), ylim = c(min(hBio_hap_posLGM$breaks), max(hBio_hap_posLGM$breaks)), axes=FALSE, frame=T, xlab="post-LGM", xaxs="i", yaxs="i", cex.lab=2)
for (biome in seq_along(unique(temp_hap_posLGM$Biome))){
        temp_hap_bio_rect <- temp_hap_posLGM[which(temp_hap_posLGM$Biome == unique(temp_hap_posLGM$Biome)[biome]),]
        base <- 0
        for(hap_col in seq_along(temp_hap_bio_rect$Haplotype)){
                rect(xleft = base,xright =  base+1, ytop = temp_hap_bio_rect$Biome[hap_col], ybottom = temp_hap_bio_rect$Biome[hap_col]+1, col=haplo_col[temp_hap_bio_rect$Haplotype[hap_col]])
                base <- base+1
        }
}
axis(side=1)
axis(side=2, at=seq(4,27,by=1), labels = seq(4,27,by=1),las=2)
