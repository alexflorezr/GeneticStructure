rm(list=ls())
setwd()
Dataset <- read.delim(file="~/Desktop/PhD/Thesis/Raw_data/Clean_database/DATABASE_16-06-16.txt", sep = "\t", stringsAsFactors = F)
Species <- unique(Dataset$Species)
Bins <- c(15000, 25000, 50000)
seqs_per_bin <- as.data.frame(matrix(data=as.numeric(), nrow = length(Species), ncol = length(Bins)))
colnames(seqs_per_bin) <- c("post-LGM", "LGM", "pre-LGM")
row.names(seqs_per_bin) <- Species
for(sp in seq_along(Species)){
        Dataset_sp <- Dataset[which(Dataset$Species == Species[sp] & nchar(Dataset$Sequence) > 0),]
        if (dim(Dataset_sp)[1] > 0){
                lower_bound <- 0
                for (bin in seq_along(Bins)){
                        seqs_per_bin[sp,bin]  <- sum(Dataset_sp$Mean_Age > lower_bound & Dataset_sp$Mean_Age <= Bins[bin])
                        lower_bound <- Bins[bin]       
                }
        }else{
                seqs_per_bin[sp,] <- NA
        }      
}
seqs_per_bin <- na.omit(seqs_per_bin)
which(sum(is.na()))
### create a z matrix from the x,y matrix
seqs_per_bin_category <- seqs_per_bin[]
seqs_per_bin_category[seqs_per_bin_category  == 0] <- -1
seqs_per_bin_category[seqs_per_bin_category  >0 & seqs_per_bin_category <= 9] <- -2
seqs_per_bin_category[seqs_per_bin_category  >=10 & seqs_per_bin_category <= 20] <- -3
seqs_per_bin_category[seqs_per_bin_category  >21 ] <- -4

### plot
par(mar=c(6,10,1,5))
image(t(as.matrix(seqs_per_bin_category)), col=rev(c("#FFFFFF", "#A8A8A8", "#4A4A4A", "#050505")), xaxs = "i", axes=F, xlim=c(-0.2,1.5))
#abline(v=1.4/Bins, col="#FFFFFF")
#unit_x <- 1/length(Bins)
#min_x <- unit_x/2
#axis(1, at = seq(0-min_x,(1+min_x + (unit_x)), by= 1/length(Bins[-c(1,2)])), labels = Bins, las=2)
mtext(colnames(seqs_per_bin), side=1, at=c(0,0.5,1), cex=0.7)
mtext(row.names(seqs_per_bin_category), side=2, at=seq(0,1, by=1/length(row.names(seqs_per_bin_category)[-1])), las=2, cex=0.7)
abline(h=seq(0,1, by=1/length(row.names(seqs_per_bin_category)[-1]))+(1/length(row.names(seqs_per_bin_category)[-1])/2), col="#FFFFFF")
legend("topright", legend = c("0", "0< n <10", "10< n <21", "<21"))
