ch3_seq_sampling(Dataset2)
### Sequence sampling per time period for all the species in Database
ch3_seq_sampling <- function(Dataset){
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
mtext(colnames(seqs_per_bin), side=1, at=c(0,0.5,1), line=1.5)
#strsplit(row.names(seqs_per_bin_category), split="_")
mtext(row.names(seqs_per_bin_category), side=2, at=seq(0,1, by=1/length(row.names(seqs_per_bin_category)[-1])), las=2, cex=0.7)
abline(h=seq(0,1, by=1/length(row.names(seqs_per_bin_category)[-1]))+(1/length(row.names(seqs_per_bin_category)[-1])/2), col="#FFFFFF", lwd=2)
legend(x=0.5, y=1, legend = c(" 0", " 0< n <10", "10< n <21", "<21"), bty = "n", fill =c("#FFFFFF", "#A8A8A8", "#4A4A4A", "#050505"), col =c("#FFFFFF", "#A8A8A8", "#4A4A4A", "#050505"))
}
### Create and save a dataframe for further Ch3 analyses
        ### Import paleobiomes
        ### Biomes using Bio4_CO2
        setwd("/Users/afr/Desktop/PhD/Thesis/Chapter3/Ch3/Ch3_data/Bio50k/")
        Bio <- stack(dir())
        ### Biomes using Koppen
        setwd("/Users/afr/Desktop/PhD/Thesis/Chapter3/Ch3/Ch3_data/Kop50k/")
        Kop <- stack(dir())
ch3_data.frame <- function(Dataset, Bio, Kop){
require(raster)
require(rgdal)
temp_hap_DB <- DATABASE[nchar(DATABASE$Sequence) > 1,]
temp_points_hap_DB <- as.data.frame(matrix(nrow=nrow(temp_hap_DB), ncol=10))
colnames(temp_points_hap_DB) <- c("Vspecies[sp]","Longitude", "Latitude", "Time_sample", "Time_bin", "Layer", "Sequence","Bio4","Kopp", "Haplotype")
Vlayer <- c(seq(0, 21000,  by=1000), seq(22000, 50000, by=2000))
for (row in seq_along(temp_hap_DB[,1])){
        temp_points_hap_DB[row,c(1,2,3,4,5,6,7)] <- c(temp_hap_DB$Species[row], 
                                                temp_hap_DB$Longitude[row],
                                                temp_hap_DB$Latitude[row],
                                                temp_hap_DB$Mean_Age[row],
                                                min(Vlayer[which(Vlayer > temp_hap_DB$Mean_Age[row])]),
                                                min(which(Vlayer > temp_hap_DB$Mean_Age[row])),
                                                temp_hap_DB$Sequence[row])
                temp_points_hap_DB[row,8]<- extract(Bio, layer=as.numeric(temp_points_hap_DB[row,6]), nl=1, y=matrix(as.numeric( temp_points_hap_DB[row,c(2,3)]), nrow=1, ncol=2))
                temp_points_hap_DB[row,9]<- extract(Kop, layer=as.numeric(temp_points_hap_DB[row,6]), nl=1, y=matrix(as.numeric( temp_points_hap_DB[row,c(2,3)]), nrow=1, ncol=2))
                
}
back_up <- temp_points_hap_DB
}
ch3_data.frame(ch3_database, Bio = Bio, Kop = Kop)
### save it raw
setwd("/Users/afr/Desktop/PhD/Thesis/Chapter3/Ch3/Ch3_results/")
write.table(temp_points_hap_DB, file="all_sp_points.txt", sep = "\t")

### Clean the data: remove NAs and -1000 values
### IMPORTANT SAFETY POINT backup_hap <- temp_points_hap_DB
### remove rows with no coordinates
if(sum(is.na(temp_points_hap_DB$Longitude)) > 0){
        temp_points_hap_DB <- temp_points_hap_DB[-which(is.na(temp_points_hap_DB$Longitude)),]
}
### remove rows with Biomes equal to NA
if(sum(temp_points_hap_DB$Bio4 == -1000) > 0 + sum(temp_points_hap_DB$kopp == 0) > 0){
        temp_points_hap_DB <- temp_points_hap_DB[-which(temp_points_hap_DB$Biome == -1000),]
        temp_points_hap_DB <- temp_points_hap_DB[-which(temp_points_hap_DB$kopp == 0),]
}
if(sum(temp_points_hap_DB$Time_bin == "Inf") > 0){
        temp_points_hap_DB <- temp_points_hap_DB[-which(temp_points_hap_DB$Time_bin == "Inf"),]
}
}



### Haplotype barplots per time period for all the species using either Bio4 or Koppen biomes