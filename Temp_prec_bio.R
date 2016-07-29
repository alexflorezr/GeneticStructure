rm(list=ls())
library(raster)
library(rgdal)
library(VennDiagram)
### Import the paleoclimatic data
### Temperature
setwd("/Users/afr/Desktop/PhD/Thesis/Chapter3/Ch3/Ch3_data/Temp50k/")
Tmp <- stack(dir())
### Precipitation
setwd("/Users/afr/Desktop/PhD/Thesis/Chapter3/Ch3/Ch3_data/Prec50k/")
Prc <- stack(dir())
### Biomes using Bio4_CO2
setwd("/Users/afr/Desktop/PhD/Thesis/Chapter3/Ch3/Ch3_data/Bio50k/")
Bio <- stack(dir())
### Import the database
DATABASE <- read.delim(file.choose(), header = T, sep = "\t", stringsAsFactors = F)
### Extract the values of climate (temp and prec) and biomes for every species
Vspecies <- unique(DATABASE$Species)
Vspecies <- "Mammuthus_primigenius" 
for(sp in seq_along(Vspecies)){
        temp_DB <- DATABASE[DATABASE$Species == Vspecies[sp], ]
        temp_points_DB <- as.data.frame(matrix(nrow=nrow(temp_DB), ncol=10))
        colnames(temp_points_DB) <- c("Vspecies[sp]","Longitude", "Latitude", "Time_sample", "Time_bin", "Layer", "Type","Temp","Prec","Biome")
        Vlayer <- c(seq(0, 21000,  by=1000), seq(22000, 50000, by=2000))
        for (row in seq_along(temp_DB[,1])){
                temp_points_DB[row,c(1,2,3,4,5,6,7)] <- c(temp_DB$Species[row], 
                        temp_DB$Longitude[row],
                        temp_DB$Latitude[row],
                        temp_DB$Mean_Age[row],
                        min(Vlayer[which(Vlayer > temp_DB$Mean_Age[row])]),
                        min(which(Vlayer > temp_DB$Mean_Age[row])),
                        ifelse(nchar(temp_DB$Sequence[row]) > 1, "Seq", "Fossil")
                )
        }
        for(row2 in seq_along(temp_points_DB$Longitude)){
                temp_points_DB[row2,8]<- extract(Tmp, layer=as.numeric(temp_points_DB[row2,6]), nl=1, y=matrix(as.numeric(temp_points_DB[row2,c(2,3)]), nrow=1, ncol=2))
                temp_points_DB[row2,9]<- extract(Prc, layer=as.numeric(temp_points_DB[row2,6]), nl=1, y=matrix(as.numeric(temp_points_DB[row2,c(2,3)]), nrow=1, ncol=2))
                temp_points_DB[row2,10]<- extract(Bio, layer=as.numeric(temp_points_DB[row2,6]), nl=1, y=matrix(as.numeric(temp_points_DB[row2,c(2,3)]), nrow=1, ncol=2))
        }
}
### Remove bad data in temp_points_DB
backup <- temp_points_DB
temp_points_DB <- temp_points_DB[-1,]
temp_points_DB <- temp_points_DB[-seq(1855,1868, by=1),]
temp_points_DB <- temp_points_DB[-which(is.na(temp_points_DB$Longitude)),]
### Subsetting the data in time periods (preLGM, LGM, posLGM)
temp_preLGM <- temp_points_DB[temp_points_DB$Time_bin >= 25000,]
temp_LGM <- temp_points_DB[temp_points_DB$Time_bin < 25000 & temp_points_DB$Time_bin >= 15000,]
temp_posLGM <- temp_points_DB[temp_points_DB$Time_bin < 15000,]
### Histograms for temperature
hist(temp_preLGM$Temp, col="black", freq=T)
hist(temp_LGM$Temp, col="red", add=T)
hist(temp_posLGM$Temp, col="blue", add=T)
### Histograms for precipitation
hist(temp_preLGM$Prec, freq=T )
hist(temp_LGM$Prec, col="red", add=T)
hist(temp_posLGM$Prec, col="blue", add=T)

### Histograms for the biomes
Bio_all_preLGM <- temp_preLGM$Biome
Bio_all_preLGM <- Bio_all_preLGM[-which(is.na(Bio_all_preLGM))]
Bio_all_preLGM <- Bio_all_preLGM[-which(Bio_all_preLGM == -1000)]
hBio_preLGM <- hist(Bio_all_preLGM, breaks=seq(3,27,by=1), plot=F)

Bio_all_LGM <- temp_LGM$Biome
Bio_all_LGM <- Bio_all_LGM[-which(is.na(Bio_all_LGM))]
hBio_LGM <- hist(Bio_all_LGM, breaks=seq(3,27,by=1), plot=F)

Bio_all_posLGM <- temp_posLGM$Biome
Bio_all_posLGM <- Bio_all_posLGM[-which(is.na(Bio_all_posLGM))]
Bio_all_posLGM <- Bio_all_posLGM[-which(Bio_all_posLGM == -1000)]
hBio_posLGM <- hist(Bio_all_posLGM, breaks=seq(3,27,by=1), plot=F)

### Histograms for the sequences
Bio_seq_preLGM <- temp_preLGM$Biome[temp_preLGM$Type == "Seq"]
Bio_seq_preLGM <- Bio_seq_preLGM[-which(Bio_seq_preLGM == -1000)]
hBio_seq_preLGM <- hist(Bio_seq_preLGM, breaks=seq(3,27,by=1), plot=F)

Bio_seq_LGM <- temp_LGM$Biome[temp_LGM$Type == "Seq"]
Bio_seq_LGM <- Bio_seq_LGM[-which(is.na(Bio_seq_LGM))]
hBio_seq_LGM <- hist(Bio_seq_LGM, breaks=seq(3,27,by=1), plot=F)

Bio_seq_posLGM <- temp_posLGM$Biome[temp_posLGM$Type == "Seq"]
Bio_seq_posLGM <- Bio_seq_posLGM[-which(Bio_seq_posLGM == -1000)]
hBio_seq_posLGM <- hist(Bio_seq_posLGM, breaks=seq(3,27,by=1), plot=F)

### Histograms (density) for the region for each period (preLGM, LGM, posLGM)
e <- extent(-180,180,30, 90)
Bio_Holartic <- crop(Bio, e)
Bio_Holartic_preLGM <- Bio_Holartic[[seq(37,25, by=-1)]]
VBio_Holartic_preLGM <- values(Bio_Holartic_preLGM)
VBio_Holartic_preLGM <- VBio_Holartic_preLGM[-which(VBio_Holartic_preLGM == -1000)]
density_preLGM <- density(VBio_Holartic_preLGM, from=3, to=27)
plot(density_preLGM)

Bio_Holartic_LGM <- Bio_Holartic[[seq(24,16, by=-1)]]
VBio_Holartic_LGM <- values(Bio_Holartic_LGM)
VBio_Holartic_LGM <- VBio_Holartic_LGM[-which(VBio_Holartic_LGM == -1000)]
density_LGM <- density(VBio_Holartic_LGM, from=3, to=27)
plot(density_LGM)
Bio_Holartic_posLGM <- Bio_Holartic[[seq(15,1, by=-1)]]
VBio_Holartic_posLGM <- values(Bio_Holartic_posLGM)
VBio_Holartic_posLGM <- VBio_Holartic_posLGM[-which(VBio_Holartic_posLGM == -1000)]
density_posLGM <- density(VBio_Holartic_posLGM, from=3, to=27)
plot(density_posLGM)

### Plotting 
par(mar=c(7,0,4,0))
Layout <- layout(matrix(c(0,3,0,2,0,1,0),ncol=7, nrow=1),widths =c(4,10,1,10,1,10,4), heights = c(1,1,1))
plot(NULL, type = "n", xlim = c(0, max(hBio_preLGM$counts)), ylim = c(min(hBio_preLGM$breaks), max(hBio_preLGM$breaks)), axes=FALSE, frame=T, xlab="Pre-LGM", xaxs="i", yaxs="i")
rect(0, hBio_preLGM$breaks[1:(length(hBio_preLGM$breaks) - 1)], hBio_preLGM$counts, hBio_preLGM$breaks[2:length(hBio_preLGM$breaks)], col="#CDAA7D")
rect(0, hBio_seq_preLGM$breaks[1:(length(hBio_seq_preLGM$breaks) - 1)], hBio_seq_preLGM$counts, hBio_seq_preLGM$breaks[2:length(hBio_seq_preLGM$breaks)], col="#8B7355")
axis(side=1)
par(new=T)
plot(y=density_preLGM$x, x=density_preLGM$y, type='l', axes=F, xlab=NA, xaxs="i", yaxs="i", lwd=2, ylim=c(3,27))
axis(side = 3)

plot(NULL, type = "n", xlim = c(0, max(hBio_LGM$counts)), ylim = c(min(hBio_LGM$breaks), max(hBio_LGM$breaks)), axes=FALSE, frame=T, xlab="LGM", xaxs="i", yaxs="i")
rect(0, hBio_LGM$breaks[1:(length(hBio_LGM$breaks))], hBio_LGM$counts, hBio_LGM$breaks[2:length(hBio_LGM$breaks)], col="#BFEFFF")
rect(0, hBio_seq_LGM$breaks[1:(length(hBio_seq_LGM$breaks) - 1)], hBio_seq_LGM$counts, hBio_seq_LGM$breaks[2:length(hBio_seq_LGM$breaks)], col="#87CEFA")
axis(side=1)
par(new=T)
plot(y=density_LGM$x, x=density_LGM$y, type='l', axes=F, xlab=NA, xaxs="i", yaxs="i", lwd=2, ylim=c(3,27))
axis(side = 3)

plot(NULL, type = "n", xlim = c(0, max(hBio_posLGM$counts)), ylim = c(min(hBio_posLGM$breaks), max(hBio_posLGM$breaks)), axes=FALSE, frame=T, xlab="post-LGM", xaxs="i", yaxs="i")
rect(0, hBio_posLGM$breaks[1:(length(hBio_posLGM$breaks) - 1)], hBio_posLGM$counts, hBio_posLGM$breaks[2:length(hBio_posLGM$breaks)], col="#9ACD32")
rect(0, hBio_seq_posLGM$breaks[1:(length(hBio_seq_posLGM$breaks) - 1)], hBio_seq_posLGM$counts, hBio_seq_posLGM$breaks[2:length(hBio_seq_posLGM$breaks)], col="#698B22")
axis(side=1)
axis(side=2, at=seq(4,27,by=1), labels = seq(4,27,by=1),las=2)
par(new=T)
plot(y=density_posLGM$x, x=density_posLGM$y, type='l', axes=F, xlab=NA, xaxs="i", yaxs="i", lwd=2, ylim=c(3,27))
axis(side = 3)
axis(side = 2)

### Venn diagrams for the biomes used
### Using the unique values
Bio_preLGM <- unique(temp_preLGM$Biome)
Bio_preLGM <- Bio_preLGM[-which(is.na(Bio_preLGM))]
Bio_preLGM <- Bio_preLGM[-which(Bio_preLGM == -1000)]
Bio_LGM <- unique(temp_LGM$Biome)
Bio_LGM <- Bio_LGM[-which(is.na(Bio_LGM))]
Bio_posLGM <- unique(temp_posLGM$Biome)
Bio_posLGM <- Bio_posLGM[-which(is.na(Bio_posLGM))]
Bio_posLGM <- Bio_posLGM[-which(Bio_posLGM == -1000)]
overlap <- calculate.overlap(list(Bio_preLGM,Bio_LGM,Bio_posLGM))
venn.diagram <- draw.triple.venn(length(c(overlap$a1, overlap$a2, overlap$a4, overlap$a5)),
                                 length(c(overlap$a2, overlap$a3, overlap$a5, overlap$a6)),
                                 length(c(overlap$a4, overlap$a5, overlap$a6, overlap$a7)),
                                 n12 = length(c(overlap$a2,overlap$a5)),
                                 n23 = length(c(overlap$a5,overlap$a6)), 
                                 n13 = length(c(overlap$a4,overlap$a5)), 
                                 n123 = length(overlap$a5), 
                                 c("PreLGM", "LGM", "PosLGM"), fill = c("#8B7355", "#98F5FF", "#EE7600"))

