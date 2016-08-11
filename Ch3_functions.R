### Import paleobiomes
### Biomes using Bio4_CO2 and Koppen
setwd("/Users/afr/Desktop/3rd_chapter/Ch3_data/Bio50k/")
Bio <- stack(dir())
### Biomes using Koppen
setwd("/Users/afr/Desktop/3rd_chapter/Ch3_data/Kop50k/")
Kop <- stack(dir())

#### The function: Start
ch3_data.frame <- function(Dataset, Bio, Kop){
        require(raster)
        require(rgdal)
        temp_hap_DB <- Dataset
        temp_points_hap_DB <- as.data.frame(matrix(nrow=nrow(temp_hap_DB), ncol=10))
        colnames(temp_points_hap_DB) <- c("Species","Longitude", "Latitude", "Time_sample", "Time_bin", "Layer", "Sequence","Bio4","Kopp", "Haplotype")
        Vlayer <- c(seq(0, 21000,  by=1000), seq(22000, 50000, by=2000))
        for (row in seq_along(temp_hap_DB[,1])){
                temp_points_hap_DB[row,c(1,2,3,4,5,6,7)] <- c(temp_hap_DB$Species[row], 
                                                              temp_hap_DB$Longitude[row],
                                                              temp_hap_DB$Latitude[row],
                                                              temp_hap_DB$Mean_Age[row],
                                                              min(Vlayer[which(Vlayer > temp_hap_DB$Mean_Age[row])]),
                                                              min(which(Vlayer > temp_hap_DB$Mean_Age[row])),
                                                              temp_hap_DB$Sequence[row])
        }
        for (layers in seq_along(unique(temp_points_hap_DB$Layer))){
                Vlayers <- sort(as.numeric(unique(temp_points_hap_DB$Layer)))
                Vnrows <- sum(as.numeric(temp_points_hap_DB$Layer) == Vlayers[layers])
                y  <- matrix(as.numeric(), ncol=2, nrow = Vnrows)
                y[,1] <- as.numeric(temp_points_hap_DB[which(as.numeric(temp_points_hap_DB$Layer) == Vlayers[layers]),2])
                y[,2] <- as.numeric(temp_points_hap_DB[which(as.numeric(temp_points_hap_DB$Layer) == Vlayers[layers]),3])
                temp_points_hap_DB[which(as.numeric(temp_points_hap_DB$Layer) == Vlayers[layers]),8]<- extract(Bio, layer=Vlayers[layers], nl=1, y)
                temp_points_hap_DB[which(as.numeric(temp_points_hap_DB$Layer) == Vlayers[layers]),9]<- extract(Kop, layer=Vlayers[layers], nl=1, y)
        }
        return(temp_points_hap_DB)
}
#### The function: End

ch3_seq_data_frame <- ch3_data.frame(ch3_database, Bio = Bio, Kop = Kop)
### SAFETY POINT >>> save the data frame to not run the script again and make a copy in the R environment
ch3_seq_data_frame <- ch3_seq_data_frame[-which(is.na(ch3_seq_data_frame$Kopp)),]
setwd("/Users/afr/Desktop/3rd_chapter/Ch3_results/")
write.table(ch3_data_frame, file="all_sp_points.txt", sep = "\t")
### remove rows with NA values in Bio4 or Kopp
ch3_seq_data_frame <- ch3_seq_data_frame[-which(is.na(ch3_seq_data_frame$Bio4)),]
### replace -1000 values in Bio4 by 0
ch3_seq_data_frame$Bio4[which(ch3_seq_data_frame$Bio4 == -1000)] <- 0
### Subsetting the data in time periods (preLGM, LGM, posLGM) and in seq anf fossil

#### The function: Start
ch3_seq_biome.barplot <- function(ch3_data_frame, Kop){
        tropical <- colorRampPalette(c("pink","#8B2323"))(4)
        arid <- colorRampPalette(c("orange", "orange4"))(4)
        warm <- colorRampPalette(c("#90EE90", "darkgreen"))(9)
        snow <- colorRampPalette(c("midnightblue","#1C86EE", "skyblue1", "#AEEEEE", "#6C7B8B"))(12)
        polar <- colorRampPalette(c("thistle", "thistle4"))(2) 
        bio_colors <- c("white", tropical, arid, warm, snow, polar)
        bio_names <- c("None","Af" ,"Am" ,"As" ,"Aw" ,"BWk" ,"BWh" ,"BSk" ,"BSh" ,"Cfa" ,"Cfb" ,"Cfc" ,"Csa" ,"Csb" ,"Csc" ,"Cwa" ,"Cwb" ,"Cwc" ,"Dfa" ,"Dfb" ,"Dfc" ,"Dfd" ,"Dsa" ,"Dsb" ,"Dsc" ,"Dsd" ,"Dwa" ,"Dwb" ,"Dwc" ,"Dwd" ,"EF" ,"ET")
        
        for (species in seq_along(unique(ch3_data_frame$Species))){
                Vspecies <- unique(ch3_data_frame$Species)
                temp_all_sp <- ch3_data_frame[ch3_data_frame$Species == Vspecies[species],]
                ### subset sequences and fossils
                temp_seq_sp <- temp_all_sp[which(nchar(temp_all_sp$Sequence) > 1), ]
                temp_fossil_sp <- temp_all_sp[which(nchar(temp_all_sp$Sequence) < 1), ]
                ### subset per time bin
                temp_seq_sp_preLGM <- temp_seq_sp[temp_seq_sp$Time_bin >= 25000,]
                temp_seq_sp_LGM <- temp_seq_sp[temp_seq_sp$Time_bin < 25000 & temp_seq_sp$Time_bin >= 15000,]
                temp_seq_sp_posLGM <- temp_seq_sp[temp_seq_sp$Time_bin < 15000,]
                temp_fossil_sp_preLGM <- temp_fossil_sp[temp_fossil_sp$Time_bin >= 25000,]
                temp_fossil_sp_LGM <- temp_fossil_sp[temp_fossil_sp$Time_bin < 25000 & temp_fossil_sp$Time_bin >= 15000,]
                temp_fossil_sp_posLGM <- temp_fossil_sp[temp_fossil_sp$Time_bin < 15000,] 
                ### Histograms for the biomes
                #Bio_min_limit <- min(temp_all_sp$Bio4)
                #Bio_max_limit <- max(temp_all_sp$Bio4)
                Kop_min_limit <- 0
                Kop_max_limit <- 31
                #Bio_seq_preLGM <- temp_seq_sp_preLGM$Bio4
                Kop_seq_preLGM <- temp_seq_sp_preLGM$Kopp
                #hBio_seq_preLGM <- hist(Bio_seq_preLGM, breaks=seq(Bio_min_limit, Bio_max_limit, by=1), plot=F)
                hKop_seq_preLGM <- hist(Kop_seq_preLGM, breaks=seq(Kop_min_limit, Kop_max_limit, by=1), plot=F)
                
                #Bio_fossil_preLGM <- temp_fossil_sp_preLGM$Bio4
                Kop_fossil_preLGM <- temp_fossil_sp_preLGM$Kopp
                #hBio_fossil_preLGM <- hist(Bio_fossil_preLGM, breaks=seq(Bio_min_limit, Bio_max_limit, by=1), plot=F)
                hKop_fossil_preLGM <- hist(Kop_fossil_preLGM, breaks=seq(Kop_min_limit, Kop_max_limit, by=1), plot=F)
                
                #Bio_seq_LGM <- temp_seq_sp_LGM$Bio4
                Kop_seq_LGM <- temp_seq_sp_LGM$Kopp
                #hBio_seq_LGM <- hist(Bio_seq_LGM, breaks=seq(Bio_min_limit, Bio_max_limit, by=1), plot=F)
                hKop_seq_LGM <- hist(Kop_seq_LGM, breaks=seq(Kop_min_limit, Kop_max_limit, by=1), plot=F)
                
                #Bio_fossil_LGM <- temp_fossil_sp_LGM$Bio4
                Kop_fossil_LGM <- temp_fossil_sp_LGM$Kopp
                #hBio_fossil_LGM <- hist(Bio_fossil_LGM, breaks=seq(Bio_min_limit, Bio_max_limit, by=1), plot=F)
                hKop_fossil_LGM <- hist(Kop_fossil_LGM, breaks=seq(Kop_min_limit, Kop_max_limit, by=1), plot=F)
                
                #Bio_seq_posLGM <- temp_seq_sp_posLGM$Bio4
                Kop_seq_posLGM <- temp_seq_sp_posLGM$Kopp
                #hBio_seq_posLGM <- hist(Bio_seq_posLGM, breaks=seq(Bio_min_limit, Bio_max_limit, by=1), plot=F)
                hKop_seq_posLGM <- hist(Kop_seq_posLGM, breaks=seq(Kop_min_limit, Kop_max_limit, by=1), plot=F)
                
                #Bio_seq_posLGM <- temp_seq_sp_posLGM$Bio4
                Kop_seq_posLGM <- temp_seq_sp_posLGM$Kopp
                #hBio_seq_posLGM <- hist(Bio_seq_posLGM, breaks=seq(Bio_min_limit, Bio_max_limit, by=1), plot=F)
                hKop_seq_posLGM <- hist(Kop_seq_posLGM, breaks=seq(Kop_min_limit, Kop_max_limit, by=1), plot=F)
                
                #Bio_fossil_posLGM <- temp_fossil_sp_posLGM$Bio4
                Kop_fossil_posLGM <- temp_fossil_sp_posLGM$Kopp
                #hBio_fossil_posLGM <- hist(Bio_fossil_posLGM, breaks=seq(Bio_min_limit, Bio_max_limit, by=1), plot=F)
                hKop_fossil_posLGM <- hist(Kop_fossil_posLGM, breaks=seq(Kop_min_limit, Kop_max_limit, by=1), plot=F)
                
                ### Histograms (density) for the region for each period (preLGM, LGM, posLGM)
                e <- extent(-180,180,30, 90)
                #Bio_Holartic <- crop(Bio, e)
                #Bio_Holartic_preLGM <- Bio_Holartic[[seq(37,25, by=-1)]]
                #VBio_Holartic_preLGM <- values(Bio_Holartic_preLGM)
                #VBio_Holartic_preLGM <- VBio_Holartic_preLGM[-which(VBio_Holartic_preLGM == -1000)]
                #Bio_density_preLGM <- density(VBio_Holartic_preLGM, from=3, to=Bio_max_limit+1)
                
                #Bio_Holartic_LGM <- Bio_Holartic[[seq(24,16, by=-1)]]
                #VBio_Holartic_LGM <- values(Bio_Holartic_LGM)
                #VBio_Holartic_LGM <- VBio_Holartic_LGM[-which(VBio_Holartic_LGM == -1000)]
                #Bio_density_LGM <- density(VBio_Holartic_LGM, from=1, to=Bio_max_limit+1)
                
                #Bio_Holartic_posLGM <- Bio_Holartic[[seq(15,1, by=-1)]]
                #VBio_Holartic_posLGM <- values(Bio_Holartic_posLGM)
                #VBio_Holartic_posLGM <- VBio_Holartic_posLGM[-which(VBio_Holartic_posLGM == -1000)]
                #Bio_density_posLGM <- density(VBio_Holartic_posLGM, from=1, to=Bio_max_limit+1)
                
                Kop_Holartic <- crop(Kop, e)
                Kop_Holartic_preLGM <- Kop_Holartic[[seq(37,25, by=-1)]]
                VKop_Holartic_preLGM <- values(Kop_Holartic_preLGM)
                VKop_Holartic_preLGM <- VKop_Holartic_preLGM[-which(VKop_Holartic_preLGM == 0)]
                Kop_density_preLGM <- density(VKop_Holartic_preLGM, from=1, to=32)
                
                Kop_Holartic_LGM <- Kop_Holartic[[seq(24,16, by=-1)]]
                VKop_Holartic_LGM <- values(Kop_Holartic_LGM)
                VKop_Holartic_LGM <- VKop_Holartic_LGM[-which(VKop_Holartic_LGM == 0)]
                Kop_density_LGM <- density(VKop_Holartic_LGM, from=1, to=32)
                
                Kop_Holartic_posLGM <- Kop_Holartic[[seq(15,1, by=-1)]]
                VKop_Holartic_posLGM <- values(Kop_Holartic_posLGM)
                VKop_Holartic_posLGM <- VKop_Holartic_posLGM[-which(VKop_Holartic_posLGM == 0)]
                Kop_density_posLGM <- density(VKop_Holartic_posLGM, from=1, to=32)
                
                ### Plotting
                setwd("/Users/afr/Desktop/")
                sp <- paste(substr(unlist(strsplit(Vspecies[species], split="_")), start = 1, stop=1), collapse = "")  
                pdf(paste(sp, "_seq_fossil_biomes.pdf"), paper = "a4")
                par(mar=c(4,0,5,0))
                Layout <- layout(matrix(c(5,4,0,3,0,2,0,1,0),ncol=9, nrow=1),widths =c(3,2,0.3,10,1,10,1,10,3), heights = c(1,1,1))
                #max_bio_pre <- max(c(max(hBio_fossil_preLGM$counts), max(hBio_seq_preLGM$counts)))
                #plot(NULL, type = "n", xlim = c(0, max_bio_pre+(max_bio_pre/30)), ylim = c(min(hBio_fossil_preLGM$breaks), max(hBio_fossil_preLGM$breaks)+1), axes=FALSE, frame=T, xlab=NA, xaxs="i", yaxs="i")
                #rect(0, hBio_fossil_preLGM$breaks[1:(length(hBio_fossil_preLGM$breaks) - 1)], hBio_fossil_preLGM$counts, hBio_fossil_preLGM$breaks[2:length(hBio_fossil_preLGM$breaks)], col="#CDAA7D", border="white", lwd=0.5)
                #rect(0, hBio_seq_preLGM$breaks[1:(length(hBio_seq_preLGM$breaks) - 1)], hBio_seq_preLGM$counts, hBio_seq_preLGM$breaks[2:length(hBio_seq_preLGM$breaks)], col="#8B7355", border = "white", lwd=0.5)
                #axis(side=1)
                #par(new=T)
                #plot(y=Bio_density_preLGM$x, x=Bio_density_preLGM$y, type='l', axes=F, xlab=NA, xaxs="i", yaxs="i", lwd=1, ylim=c(0,27))
                #axis(side = 3)
                
                max_kop_pre <- max(c(max(hKop_fossil_preLGM$counts), max(hKop_seq_preLGM$counts)))
                plot(NULL, type = "n", xlim = c(0, max_kop_pre+(max_kop_pre/30)), ylim = c(min(hKop_fossil_preLGM$breaks), max(hKop_fossil_preLGM$breaks)+1), axes=FALSE, frame=T, xaxs="i", yaxs="i", ylab=NA, xlab=NA)
                rect(0, hKop_fossil_preLGM$breaks[1:(length(hKop_fossil_preLGM$breaks) - 1)], hKop_fossil_preLGM$counts, hKop_fossil_preLGM$breaks[2:length(hKop_fossil_preLGM$breaks)], col="#228B22", border="white" , lwd=0.5)
                rect(0, hKop_seq_preLGM$breaks[1:(length(hKop_seq_preLGM$breaks) - 1)], hKop_seq_preLGM$counts, hKop_seq_preLGM$breaks[2:length(hKop_seq_preLGM$breaks)], col="#90EE90", border ="white", lwd=0.5 )
                axis(side=1)
                mtext("Pre-LGM", side=1, line=3, cex=1.25)
                par(new=T)
                plot(y=Kop_density_preLGM$x, x=Kop_density_preLGM$y, type='l', axes=F, xlab=NA, xaxs="i", yaxs="i", lwd=1, ylim=c(0,32), ylab=NA)
                axis(side = 3)
                
                ### LGM
                #max_bio_lgm <- max(c(max(hBio_fossil_LGM$counts), max(hBio_seq_LGM$counts)))
                #plot(NULL, type = "n", xlim = c(0, max_bio_lgm+(max_bio_lgm/30)), ylim = c(min(hBio_fossil_LGM$breaks), max(hBio_fossil_LGM$breaks)+1), axes=FALSE, frame=T, xlab=NA, xaxs="i", yaxs="i")
                #rect(0, hBio_fossil_LGM$breaks[1:(length(hBio_fossil_LGM$breaks) - 1)], hBio_fossil_LGM$counts, hBio_fossil_LGM$breaks[2:length(hBio_fossil_LGM$breaks)], col="#BFEFFF", border = "white", lwd=0.5)
                #rect(0, hBio_seq_LGM$breaks[1:(length(hBio_seq_LGM$breaks) - 1)], hBio_seq_LGM$counts, hBio_seq_LGM$breaks[2:length(hBio_seq_LGM$breaks)], col="#87CEFA", border = "white", lwd=0.5)
                #axis(side=1)
                #par(new=T)
                #plot(y=Bio_density_LGM$x, x=Bio_density_LGM$y, type='l', axes=F, xlab=NA, xaxs="i", yaxs="i", lwd=1, ylim=c(0,27))
                #axis(side = 3)
                
                
                max_kop_lgm <- max(c(max(hKop_fossil_LGM$counts), max(hKop_seq_LGM$counts)))
                plot(NULL, type = "n", xlim = c(0, max_kop_lgm+(max_kop_lgm/30)), ylim = c(min(hKop_fossil_LGM$breaks), max(hKop_fossil_LGM$breaks)+1), axes=FALSE, frame=T, xaxs="i", yaxs="i", ylab=NA, xlab=NA)
                rect(0, hKop_fossil_LGM$breaks[1:(length(hKop_fossil_LGM$breaks) - 1)], hKop_fossil_LGM$counts, hKop_fossil_LGM$breaks[2:length(hKop_fossil_LGM$breaks)], col="#009ACD", border = "white" , lwd=0.5)
                rect(0, hKop_seq_LGM$breaks[1:(length(hKop_seq_LGM$breaks) - 1)], hKop_seq_LGM$counts, hKop_seq_LGM$breaks[2:length(hKop_seq_LGM$breaks)], col="#87CEFA", border="white", lwd=0.5)
                axis(side=1)
                mtext("LGM", side=1, line=3, cex=1.25)
                mtext(str_replace(Vspecies[species], "_", " "), side=3, line=3, cex=1.25)
                par(new=T)
                plot(y=Kop_density_LGM$x, x=Kop_density_LGM$y, type='l', axes=F, xlab=NA, xaxs="i", yaxs="i", lwd=1, ylim=c(0,32))
                axis(side = 3)
                
                ### posLGM
                #max_bio_poslgm <- max(c(max(hBio_fossil_posLGM$counts), max(hBio_seq_posLGM$counts)))
                #plot(NULL, type = "n", xlim = c(0, max_bio_poslgm+(max_bio_poslgm/30)), ylim = c(min(hBio_fossil_posLGM$breaks), max(hBio_fossil_posLGM$breaks)+1), axes=FALSE, frame=T, xlab=NA, xaxs="i", yaxs="i")
                #rect(0, hBio_fossil_posLGM$breaks[1:(length(hBio_fossil_posLGM$breaks) - 1)], hBio_fossil_posLGM$counts, hBio_fossil_posLGM$breaks[2:length(hBio_fossil_posLGM$breaks)], col="#9ACD32", border="white", lwd=0.5)
                #rect(0, hBio_seq_posLGM$breaks[1:(length(hBio_seq_posLGM$breaks) - 1)], hBio_seq_posLGM$counts, hBio_seq_posLGM$breaks[2:length(hBio_seq_posLGM$breaks)], col="#698B22", border = "white", lwd=0.5)
                #axis(side=1)
                #mtext("Biomes Bio4 CO2", side=2, line=3, cex=1.25)
                #par(new=T)
                #plot(y=Bio_density_posLGM$x, x=Bio_density_posLGM$y, type='l', axes=F, xlab=NA, xaxs="i", yaxs="i", lwd=1, ylim=c(0,27))
                #axis(side = 3)
                #axis(side = 2, at = seq(min(hBio_fossil_posLGM$breaks), max(hBio_fossil_posLGM$breaks)+1, by=1), las=2)
                
                max_kop_poslgm <- max(c(max(hKop_fossil_posLGM$counts), max(hKop_seq_posLGM$counts)))
                plot(NULL, type = "n", xlim = c(0, max_kop_poslgm+(max_kop_poslgm/30)), ylim = c(min(hKop_fossil_posLGM$breaks), max(hKop_fossil_posLGM$breaks)+1), axes=FALSE, frame=T, xaxs="i", yaxs="i", ylab=NA, xlab=NA)
                rect(0, hKop_fossil_posLGM$breaks[1:(length(hKop_fossil_posLGM$breaks) - 1)], hKop_fossil_posLGM$counts, hKop_fossil_posLGM$breaks[2:length(hKop_fossil_posLGM$breaks)], col="#CD3700", border ="white", lwd=0.5 )
                rect(0, hKop_seq_posLGM$breaks[1:(length(hKop_seq_posLGM$breaks) - 1)], hKop_seq_posLGM$counts, hKop_seq_posLGM$breaks[2:length(hKop_seq_posLGM$breaks)], col="#FF7F24", border="white", lwd=0.5)
                axis(side=1)
                mtext("Pos-LGM", side=1, line=3, cex=1.25)
                mtext("Biomes Köppen-Geiger", side=2.5, line=5, cex=1.25)
                par(new=T)
                plot(y=Kop_density_posLGM$x, x=Kop_density_posLGM$y, type='l', axes=F, xlab=NA, xaxs="i", yaxs="i", lwd=1, ylim=c(0,32), ylab=NA)
                axis(side = 3)

                ### Biomes
                plot(NULL, type = "n", xlim = c(0, max_posLGM), ylim = c(min(breaks), max(breaks)+1), axes=FALSE, frame=T, xlab=NA,ylab=NA, xaxs="i", yaxs="i", cex.lab=2)
                for (bio_col in seq_along(breaks)){
                        rect(xleft = 0,xright =  base+max_posLGM, ytop = breaks[bio_col], ybottom = breaks[bio_col]+1, col=bio_colors[bio_col], border = "white")
                }
                axis(side=2, at = breaks+0.5, labels = bio_names, las=2, hadj=0, line = -1.5, cex=0.7, tick = F)
                
                dev.off()
        }
}
#### The function: End
ch3_seq_fossil_bio_all <- read.table(file.choose(), header = T, stringsAsFactors = F, sep = "\t")
ch3_data_frame <- ch3_seq_fossil_bio_all
ch3_seq_biome.barplot(ch3_data_frame =ch3_seq_fossil_bio_all, Kop = Kop)

### create fasta files from the ch3_data_frame
setwd("/Users/afr/Desktop/3rd_chapter/Ch3_results/Ch3_fasta/")

#### The function: Start
ch3_fasta.files <- function(ch3_data_frame){
        require(seqinr)
        for (species in seq_along(unique(ch3_data_frame$Species))){
                Vspecies <- unique(ch3_data_frame$Species)
                temp_all_sp <- ch3_data_frame[ch3_data_frame$Species == Vspecies[species],]
                temp_seq_sp <- temp_all_sp[which(nchar(temp_all_sp$Sequence) > 1), ]
                sp <- paste(substr(unlist(strsplit(Vspecies[species], split="_")), start = 1, stop=1), collapse = "")  
                write.fasta(as.list(temp_seq_sp$Sequence), names = paste(sp,seq(1,length.out = length(temp_points_hap_DB$Longitude)), as.numeric(temp_seq_sp$Time_bin)/1000, temp_seq_sp$Bio4, temp_seq_sp$Kopp, sep="_"), file.out = paste(sp, "_hap_bio_kop.fasta", sep = ""))
                write.table(temp_seq_sp, file=paste(sp, "bio_kop_df.txt", sep=""), sep="\t", row.names = F)
        }
}
#### The function: End
ch3_fasta.files(ch3_data_frame =ch3_data_frame)

### Function to assign haplotypes to Biomes
"/Users/afr/Desktop/3rd_chapter/Ch3_results/Ch3_fasta_clean/"
#### The function: Start
ch3_assing.hap <- function(directory){
        require(ape)
        setwd(directory)
        fasta_files <- dir(pattern = ".fasta")
        table_files <- dir(pattern = "_df.txt")
        for (file in seq_along(fasta_files)){
                temp_fasta_sp <- read.dna(fasta_files[file], format = "fasta")
                temp_ch3_df_sp <- read.table(table_files[file], stringsAsFactors = F, header = T)
                temp_split <- strsplit(dimnames(temp_fasta_sp)[[1]], split = "_")
                clean <- as.numeric(unlist(temp_split)[(seq(1, length(dimnames(temp_fasta_sp)[[1]]), by=1)*5)-3])
                temp_hap <- haplotype(temp_fasta_sp)
                temp_ch3_df_sp <- temp_ch3_df_sp[clean,]
                hap_index <- 1
                sp <- substring(fasta_files[file], first = 1, last = 2)
                for (hap in seq_along(attr(temp_hap, which = "index"))){
                        samples <- attr(temp_hap, which = "index")[[hap]]
                        temp_ch3_df_sp$Haplotype[samples] <- hap_index
                        hap_index <- hap_index + 1
                }
                #compare the fasta names and the table names
                dimnames(temp_fasta_sp)[[1]]
                names_fasta <- paste(as.numeric(unlist(temp_split)[c((seq(1, length(dimnames(temp_fasta_sp)[[1]]), by=1)*5)-2)]),
                                     as.numeric(unlist(temp_split)[c((seq(1, length(dimnames(temp_fasta_sp)[[1]]), by=1)*5)-1)]), 
                                     as.numeric(unlist(temp_split)[c((seq(1, length(dimnames(temp_fasta_sp)[[1]]), by=1)*5)-0)]),
                                     sep="_")
                names_table <- paste(temp_ch3_df_sp$Time_bin/1000,temp_ch3_df_sp$Bio4, temp_ch3_df_sp$Kopp, sep="_")
                if(sum(names_table == names_fasta) == length(names_fasta)){
                        "Proceed, df and fasta match"
                        write.table(temp_ch3_df_sp, file = paste(sp, "_hap_bio_kop.txt", sep=""), sep="\t", row.names = F)
                }else{
                        "Wrong df and fasta matching"
                }
        }
} 
#### The function: End
ch3_assing.hap("/Users/afr/Desktop/3rd_chapter/Ch3_results/Ch3_fasta_clean/")


### Function for the map for all the species
#### The function: Start
ch3_maps.plot <- function(ch3_hap_bio_kop){
        require(stringr)
        require(maps)
        require(mapproj)
        map_data <- read.table(ch3_hap_bio_kop, header=T, stringsAsFactors = F)
        map_data$Period <- NA
        map_data$Col_period <- NA
        map_data$Period[which(map_data$Time_bin > 25000)] <- "pre-LGM"
        map_data$Col_period[which(map_data$Time_bin > 25000)] <- "#8B7355"
        map_data$Period[which(map_data$Time_bin <= 25000 & map_data$Time_bin > 15000)] <- "LGM"
        map_data$Col_period[which(map_data$Time_bin <= 25000 & map_data$Time_bin > 15000)] <- "#87CEFA"
        map_data$Period[which(map_data$Time_bin <= 15000)] <- "pos-LGM"
        map_data$Col_period[which(map_data$Time_bin <= 15000)] <- "#698B22"
        table_periods <- table(map_data$Period)
        #table(map_data$Period)
        species <- str_replace(unique(map_data$Species),"_", " ")
        sp <- paste(substring(unlist(strsplit(unique(map_data$Species), split="_")), first = 1, last = 1), collapse = "")
        pdf(paste(sp, "_ch3_maps.pdf", sep=""), paper = "a4")
        map("world",proj="azequidist", resolution = 0.2, ylim = c(33,90), lwd=1.5, mar=c(0,2,2,2), interior = F) 
        points(mapproject(as.numeric(map_data$Longitude), as.numeric(map_data$Latitude)), 
               cex=1.4, bg=paste(map_data$Col_period, "99", sep=""), col=map_data$Col_period, pch=21)
        #draw.circle(0,0,1.17, lwd=2)
        mtext(species, side=3, line=1)
        mtext(paste("n = ", sum(!is.na(map_data$Longitude))), line=0)
        legend( -1.2,-0.55, pch=16, c(paste("pre-LGM", " (", table_periods[3], ")", sep = ""),
                                      paste("LGM", " (", table_periods[1], ")", sep = ""),
                                      paste("pos-LGM", " (", table_periods[2], ")", sep = "")),
                col=c("#8B735599","#87CEFA99", "#698B2299"), bty = "n", pt.cex = 2, pt.lwd = 3, title="Time period", title.adj = c(0,0), y.intersp = 1, x.intersp = 0.5)
        dev.off()
}
#### The function: End
for(file in dir(pattern = "kop.txt")){
        ch3_maps.plot(file)
}

### Function to plot the haplotypes per biome per time period (using Koppen)
#### The function: Start
ch3_hap_kopp.barplot <- function(ch3_hap_bio_kop, seq_fossil_bio){
        require(stringr)
        ## Colors for biomes
        tropical <- colorRampPalette(c("pink","#8B2323"))(4)
        arid <- colorRampPalette(c("orange", "orange4"))(4)
        warm <- colorRampPalette(c("#90EE90", "darkgreen"))(9)
        snow <- colorRampPalette(c("midnightblue","#1C86EE", "skyblue1", "#AEEEEE", "#6C7B8B"))(12)
        polar <- colorRampPalette(c("thistle", "thistle4"))(2) 
        bio_colors <- c("white", tropical, arid, warm, snow, polar)
        bio_names <- c("None","Af" ,"Am" ,"As" ,"Aw" ,"BWk" ,"BWh" ,"BSk" ,"BSh" ,"Cfa" ,"Cfb" ,"Cfc" ,"Csa" ,"Csb" ,"Csc" ,"Cwa" ,"Cwb" ,"Cwc" ,"Dfa" ,"Dfb" ,"Dfc" ,"Dfd" ,"Dsa" ,"Dsb" ,"Dsc" ,"Dsd" ,"Dwa" ,"Dwb" ,"Dwc" ,"Dwd" ,"EF" ,"ET")
        ### time periods
        hap_bio_kop_data <- read.table(ch3_hap_bio_kop, header=T, stringsAsFactors = F)
        seq_fossil_bio_data <- read.table(seq_fossil_bio, header=T, stringsAsFactors = F)
        species <- str_replace(unique(hap_bio_kop_data$Species), pattern = "_", replacement = " ")
        sp <- paste(substring(unlist(strsplit(unique(hap_bio_kop_data$Species), split="_")), first = 1, last = 1), collapse = "")
        seq_fossil_bio_data_sp <- seq_fossil_bio_data[seq_fossil_bio_data$Species == unique(hap_bio_kop_data$Species),]
        seq_fossil_bio_data_sp <- seq_fossil_bio_data_sp[which(nchar(seq_fossil_bio_data_sp$Sequence) < 1),]
        pdf(paste(sp, "_ch3_hap_kopp.pdf", sep=""), paper = "a4")
        par(mar=c(12,0,17,0), lwd=0.7)
        Haplo_col_period <- as.data.frame(matrix(nrow = dim(hap_bio_kop_data)[1], ncol=4))
        colnames(Haplo_col_period) <- c("Haplotype", "Time_bin", "Haplo_color", "color_name")
        Haplo_col_period$Haplotype  <-  hap_bio_kop_data$Haplotype
        Haplo_col_period$Time_bin <- hap_bio_kop_data$Time_bin
        Haplo_col_period <- Haplo_col_period[order(Haplo_col_period$Time_bin, decreasing = T),]
        duplicated_haplo <- duplicated(Haplo_col_period$Haplotype)
        Haplo_col_period <- Haplo_col_period[-which(duplicated_haplo),]
        Haplo_col_period_preLGM <- Haplo_col_period[which(Haplo_col_period$Time_bin > 25000),]
        Haplo_col_period_LGM <- Haplo_col_period[which(Haplo_col_period$Time_bin <= 25000 & Haplo_col_period$Time_bin > 15000),]
        Haplo_col_period_posLGM <- Haplo_col_period[which(Haplo_col_period$Time_bin <= 15000),]
        
        #col_pre <- rainbow(length(Haplo_col_period_preLGM$Haplotype), start = 0.3, end = 0.45)
        col_pre <- colorRampPalette(c("green4", "lightgreen"))(length(Haplo_col_period_preLGM$Haplotype))
        Haplo_col_period_preLGM$Haplo_color <- col_pre
        Haplo_col_period_preLGM$color_name <- "greens"
        
        #col_LGM <- rainbow(length(Haplo_col_period_LGM$Haplotype), start = 0.5, end = 0.6)
        col_LGM <- colorRampPalette(c("dodgerblue3", "lightskyblue"))(length(Haplo_col_period_LGM$Haplotype))
        Haplo_col_period_LGM$Haplo_color <- col_LGM
        Haplo_col_period_LGM$color_name <- "blues"
        
        #col_posLGM <- rainbow(length(Haplo_col_period_posLGM$Haplotype), start = 0, end = 0.17)
        col_posLGM <- colorRampPalette(c("darkorange4", "darkorange"))(length(Haplo_col_period_posLGM$Haplotype))
        Haplo_col_period_posLGM$Haplo_color <- col_posLGM
        Haplo_col_period_posLGM$color_name <- "oranges"
        
        to_add_colors <- rbind(Haplo_col_period_posLGM, Haplo_col_period_LGM, Haplo_col_period_preLGM)
        hap_bio_kop_data$hap_colors <- to_add_colors$Haplo_color[match(hap_bio_kop_data$Haplotype, to_add_colors$Haplotype)]
        hap_bio_kop_data$name_colors  <- to_add_colors$color_name[match(hap_bio_kop_data$Haplotype, to_add_colors$Haplotype)]
        Layout <- layout(matrix(c(5,4,0,3,0,2,0,1,0),ncol=9, nrow=1),widths =c(3,2,0.3,10,1,10,1,10,3), heights = c(1,1,1))
        #layout.show(Layout)
        breaks <- seq(0, 31, by=1)
        ### preLGM
        hap_bio_kop_data_preLGM <- hap_bio_kop_data[which(hap_bio_kop_data$Time_bin > 25000),]
        seq_fossil_bio_data_sp_preLGM <- seq_fossil_bio_data_sp[which(seq_fossil_bio_data_sp$Time_bin > 25000),]
        if(dim(hap_bio_kop_data_preLGM)[1] < 1){
                plot(NULL, type = "n", xlim = c(0, 1), ylim = c(min(breaks), max(breaks)+1), axes=FALSE, frame=T, xlab=NA,ylab=NA, xaxs="i", yaxs="i", cex.lab=2)
        }else{
                max_pre <- max(table(hap_bio_kop_data_preLGM$Kopp))
                plot(NULL, type = "n", xlim = c(0, max_pre), ylim = c(min(breaks), max(breaks)+1), axes=FALSE, frame=T, xlab=NA, xaxs="i", yaxs="i", cex.lab=2)
                
                for (biome in seq_along(unique(hap_bio_kop_data_preLGM$Kopp))){
                        temp_hap_bio_rect <- hap_bio_kop_data_preLGM[which(hap_bio_kop_data_preLGM$Kopp == unique(hap_bio_kop_data_preLGM$Kopp)[biome]),]
                        temp_hap_bio_rect <- temp_hap_bio_rect[order(temp_hap_bio_rect$hap_colors),]
                        base <- 0
                        for(hap_col in seq_along(temp_hap_bio_rect$Haplotype)){
                                rect(xleft = base,xright =  base+1, ytop = temp_hap_bio_rect$Kopp[hap_col], ybottom = temp_hap_bio_rect$Kopp[hap_col]+1, col=temp_hap_bio_rect$hap_colors[hap_col])
                                base <- base+1
                        }
                }
                #only fossil 
                matched <- which(is.na(match(sort(unique(seq_fossil_bio_data_sp_preLGM$Kopp)),sort(unique(hap_bio_kop_data_preLGM$Kopp)))))
                only_fossil <- sort(unique(seq_fossil_bio_data_sp_preLGM$Kopp))[matched]
                for (biome_fossil in seq_along(only_fossil)){
                        base <- 0
                        rect(xleft = base,xright =  base+max_pre, ytop = only_fossil[biome_fossil], ybottom = only_fossil[biome_fossil]+1, col="#69696950", border = NA)
                        }
                }
        abline(v = 4, lwd=0.5, lty=2)
        axis(side=1)
        n_pre <- dim(hap_bio_kop_data_preLGM)[1]
        mtext(paste("(", n_pre, ")"), side=1, line=4.5, cex=0.7)
        
        ### LGM
        hap_bio_kop_data_LGM <- hap_bio_kop_data[which(hap_bio_kop_data$Time_bin <= 25000 & hap_bio_kop_data$Time_bin > 15000),]
        seq_fossil_bio_data_sp_LGM <- seq_fossil_bio_data_sp[which(seq_fossil_bio_data_sp$Time_bin <= 25000 & seq_fossil_bio_data_sp$Time_bin > 15000),]
        if(dim(hap_bio_kop_data_LGM)[1] < 1){
                plot(NULL, type = "n", xlim = c(0, 1), ylim = c(min(breaks), max(breaks)+1), axes=FALSE, frame=T, xlab=NA, xaxs="i", yaxs="i", cex.lab=2)
        }else{
                max_LGM <- max(table(hap_bio_kop_data_LGM$Kopp))
                plot(NULL, type = "n", xlim = c(0, max_LGM), ylim = c(min(breaks), max(breaks)+1), axes=FALSE, frame=T, xlab=NA,ylab=NA, xaxs="i", yaxs="i", cex.lab=2)
                for (biome in seq_along(unique(seq_fossil_bio_data$Kopp))){
                        temp_hap_bio_rect <- hap_bio_kop_data_LGM[which(hap_bio_kop_data_LGM$Kopp == unique(hap_bio_kop_data_LGM$Kopp)[biome]),]
                        temp_hap_bio_rect <- temp_hap_bio_rect[order(temp_hap_bio_rect$hap_colors),]
                        base <- 0
                        for(hap_col in seq_along(temp_hap_bio_rect$Haplotype)){
                                rect(xleft = base,xright =  base+1, ytop = temp_hap_bio_rect$Kopp[hap_col], ybottom = temp_hap_bio_rect$Kopp[hap_col]+1,col=temp_hap_bio_rect$hap_colors[hap_col])
                                base <- base+1
                        }
                }
                matched <- which(is.na(match(sort(unique(seq_fossil_bio_data_sp_LGM$Kopp)),sort(unique(hap_bio_kop_data_preLGM$Kopp)))))
                only_fossil <- sort(unique(seq_fossil_bio_data_sp_LGM$Kopp))[matched]
                for (biome_fossil in seq_along(only_fossil)){
                        base <- 0
                        rect(xleft = base,xright =  base+max_pre, ytop = only_fossil[biome_fossil], ybottom = only_fossil[biome_fossil]+1, col="#69696950", border = NA)
                }
                }
        abline(v = 4, lwd=0.5, lty=2)
        axis(side=1)
        n_lgm <- dim(hap_bio_kop_data_LGM)[1]
        mtext(paste("(", n_lgm, ")"), side=1, line=4.5, cex=0.7)
        mtext(species, side=3, line=4.5, cex=1.25)
        ### Pos LGM
        hap_bio_kop_data_posLGM <- hap_bio_kop_data[which(hap_bio_kop_data$Time_bin <= 15000),]
        seq_fossil_bio_data_sp_posLGM <- seq_fossil_bio_data_sp[which(seq_fossil_bio_data_sp$Time_bin <= 15000),]
        if(dim(hap_bio_kop_data_posLGM)[1] < 1){
                plot(NULL, type = "n", xlim = c(0, 1), ylim = c(min(breaks), max(breaks)+1), axes=FALSE, frame=T, xlab=NA, xaxs="i", yaxs="i", cex.lab=2)
        }else{
                max_posLGM <- max(table(hap_bio_kop_data_posLGM$Kopp))
                plot(NULL, type = "n", xlim = c(0, max_posLGM), ylim = c(min(breaks), max(breaks)+1), axes=FALSE, frame=T, xlab=NA,ylab=NA, xaxs="i", yaxs="i", cex.lab=2)
                for (biome in seq_along(unique(hap_bio_kop_data_posLGM$Kopp))){
                        temp_hap_bio_rect <- hap_bio_kop_data_posLGM[which(hap_bio_kop_data_posLGM$Kopp == unique(hap_bio_kop_data_posLGM$Kopp)[biome]),]
                        temp_hap_bio_rect <- temp_hap_bio_rect[order(temp_hap_bio_rect$hap_colors),]
                        base <- 0
                        for(hap_col in seq_along(temp_hap_bio_rect$Haplotype)){
                                rect(xleft = base,xright =  base+1, ytop = temp_hap_bio_rect$Kopp[hap_col], ybottom = temp_hap_bio_rect$Kopp[hap_col]+1, col=temp_hap_bio_rect$hap_colors[hap_col])
                                base <- base+1
                        }
                }
                matched <- which(is.na(match(sort(unique(seq_fossil_bio_data_sp_posLGM$Kopp)),sort(unique(hap_bio_kop_data_preLGM$Kopp)))))
                only_fossil <- sort(unique(seq_fossil_bio_data_sp_posLGM$Kopp))[matched]
                for (biome_fossil in seq_along(only_fossil)){
                        base <- 0
                        rect(xleft = base,xright =  base+max_pre, ytop = only_fossil[biome_fossil], ybottom = only_fossil[biome_fossil]+1, col="#69696950", border = NA)
                }
        }
        abline(v = 4, lwd=0.5, lty=2)
        axis(side=1)
        n_pos <- dim(hap_bio_kop_data_posLGM)[1]
        mtext(paste("(", n_pos, ")"), side=1, line=4.5, cex=0.7)
        ### Biomes
        plot(NULL, type = "n", xlim = c(0, max_posLGM), ylim = c(min(breaks), max(breaks)+1), axes=FALSE, frame=T, xlab=NA,ylab=NA, xaxs="i", yaxs="i", cex.lab=2)
        for (bio_col in seq_along(breaks)){
                rect(xleft = 0,xright =  base+max_posLGM, ytop = breaks[bio_col], ybottom = breaks[bio_col]+1, col=bio_colors[bio_col], border = "white")
        }
        axis(side=2, at = breaks+0.5, labels = bio_names, las=2)
        dev.off()
        }
#### The function: End

setwd("~/Desktop/PhD/Thesis/Chapter3/Ch3/Ch3_figures/Ch3_figures_data/Ch3_hap_barplots_data/")
setwd("/Users/afr/Desktop/Ch3_fasta_haplo_plot/")
for(file in dir(pattern = "kop.txt")){
        ch3_hap_kopp.barplot(file, "ch3_seq_fossil_biome_df.txt")
}
### function to estimate the fst among time periods
ch3_periods.fst <- function(directory){
        data_files <- dir(pattern = "df.txt")
        fasta_files <- dir(pattern = "clean.fasta")
        for(files in seq_along(data_files)){
                sp_data <- strsplit(data_files[files], split = "_")[[1]][1]
                sp_fasta <- strsplit(fasta_files[files], split = "_")[[1]][1]
                if(sp_data == sp_fasta){
                ch3_periods_data_fst <- data_files[files]
                ch3_periods_fasta_fst <- fasta_files[files]
        ### output data frame for periods overall
        temp_out_diff_overall <- as.data.frame(matrix(ncol=4, nrow = 0))
        colnames(temp_out_diff_overall) <- c("sp","diff_metric","estimate","p.val") 
        ### output data frame for periods pairwise
        temp_out_diff_pairwise <- as.data.frame(matrix(ncol=12, nrow = 0))
        colnames(temp_out_diff_pairwise) <-  c("sp","pair.label","strata.1","strata.2","n.1","n.2","Fst","Fst.p.val","PHIst","PHIst.p.val","Chi2","Chi2.p.val")
        ### output data frame for nucleotide divergence within
        temp_out_div_within <- as.data.frame(matrix(ncol=8, nrow = 0))
        colnames(temp_out_div_within) <-  c("sp","period","mean","pct.0","pct.0.025","pct.0.5","pct.0.975","pct.1")
        ### output data frame for nucleotide divergence between
        temp_out_div_between <- as.data.frame(matrix(ncol=10, nrow = 0))
        colnames(temp_out_div_between) <-  c("sp" ,"strata.1", "strata.2","dA","mean","pct.0","pct.0.025","pct.0.5","pct.0.975","pct.1")
        Periods<- c("pre_LGM", "LGM", "pos_LGM")
        temp_sp_kop_data <- read.table(ch3_periods_data_fst, header=T, stringsAsFactors = F)
        temp_sp_kop_data$Period[which(temp_sp_kop_data$Time_bin > 25000)] <- Periods[1]
        temp_sp_kop_data$Period[which(temp_sp_kop_data$Time_bin <= 25000 & temp_sp_kop_data$Time_bin > 15000)] <- Periods[2]
        temp_sp_kop_data$Period[which(temp_sp_kop_data$Time_bin <= 15000)] <- Periods[3]
        ### temp_align_clean <- read.fasta("Bs_align_clean.fasta")
        temp_align_clean  <- read.dna(ch3_periods_fasta_fst, format = "fasta", as.character = FALSE, as.matrix=NULL)
        temp_sp_kop_data <- temp_sp_kop_data[as.numeric(unlist(temp_split)[c((seq(1, length(labels(temp_align_clean)), by=1) *5)-3)]),]
        temp_sp_kop_data$UniqueInd <- labels(temp_align_clean)
        sp <- strsplit(ch3_periods_data_fst, split = "_")[[1]][1]
        haps_period <- FindHaplo(align = temp_align_clean)
        haps_period <- as.data.frame(haps_period)
        colnames(haps_period) <-c("UniqueInd", "haplotype")
        temp_data_stats <- temp_sp_kop_data
        matched <- match(haps_period$UniqueInd, temp_data_stats$UniqueInd)
        temp_data_stats <- temp_data_stats[matched,]
        temp_haps_time_Kopp <- cbind(temp_data_stats, haps_period)
        ### add a column to specify the biome group
        sp_haps_period <-GetHaplo(align = temp_align_clean, saveFile=T, outname=paste(sp,"_periods_fst",".fasta", sep=""), format="fasta", seqsNames="Inf.Hap") #haps are now a DNAbin
        sp_ghaps_period <-read.fasta(paste(sp,"_periods_fst",".fasta", sep="")) #imports haps as gtypes
        sp_gtype_period <- gtypes(gen.data=data.frame(temp_haps_time_Kopp$UniqueInd,temp_haps_time_Kopp$Period,temp_haps_time_Kopp$haplotype),id.col=1,strata.col=2,locus.col=3,dna.seq=sp_ghaps_period)
        ### Estimate the genetic parameters       
        temp_periods_diff <- pop.diff.test(sp_gtype_period)
                ### output for diff overall
                temp_sp_periods_overall <- as.data.frame(cbind(sp, rownames(temp_periods_diff$overall$result), temp_periods_diff$overall$result))
                colnames(temp_sp_periods_overall)[2] <- "diff_metric"
                temp_out_diff_overall <- rbind(temp_out_diff_pairwise,temp_sp_periods_overall)
                ### output for diff pairwise
                temp_sp_periods_pairwise <- cbind(sp,temp_periods_diff$pairwise$result)                                   
                temp_out_diff_pairwise <- rbind(temp_out_diff_pairwise,temp_sp_periods_pairwise)
        #Nucleotide divergence
                temp_nuc_divergence_all <- nucleotide.divergence(sp_gtype_period)
                ### within
                temp_nuc_within <- as.data.frame(cbind(row.names(temp_nuc_divergence_all$within), temp_nuc_divergence_all$within), row.names = F)
                colnames(temp_nuc_within)[1] <- "period"
                temp_out_div_within <- rbind(temp_out_div_within,cbind(sp,temp_nuc_within))
                ### between
                temp_nuc_between <- as.data.frame(cbind(sp, temp_nuc_divergence_all$between))
                temp_out_div_between <- rbind(temp_out_div_between,cbind(sp,temp_nuc_between))
                }else{
                print("The data files do not match the fasta files")
                }
        }
        setwd("/Users/afr/Desktop/")
        write.table(temp_out_diff_overall, file = "ch3_diff_overall.txt", sep="\t", row.names = F)
        write.table(temp_out_diff_pairwise, file = "ch3_diff_pairwise.txt", sep="\t", row.names = F)
        write.table(temp_out_div_within, file = "ch3_div_within.txt", sep="\t", row.names = F)
        write.table(temp_out_div_between, file = "ch3_div_between.txt", sep="\t", row.names = F)
}
ch3_periods.fst("/Users/afr/Desktop/Ch3_fst_period_data/")

### function to plot phist triangles among the time periods
ch3_periods_fst.triangles <- function(temp_out_diff_biome_pairwise){
        require(igraph)
        require(plotrix)
        temp_table_diff_network <- read.table(temp_out_diff_biome_pairwise, header = T, stringsAsFactors = F )
        species <- unique(temp_table_diff_network$sp)
        for(sp in seq_along(species)){
                pdf(paste(species[sp], "ch3_fst_triangles.pdf", sep=""), paper = "a4")
                temp_table_diff_network_sp <- temp_table_diff_network[which(temp_table_diff_network$sp == species[sp]),]
                time_periods <- unique(c(temp_table_diff_network_sp$strata.1, temp_table_diff_network_sp$strata.2))
                ### vertices
                Vertices <- as.data.frame(matrix(nrow = 3, ncol=3))
                colnames(Vertices) <- c("Time_period","Size", "Color")
                Vertices$Time_period <- time_periods
                Vertices$Size <- 2
                #Vertices[1,] <- temp_table_diff_network_sp[1,c(3,5)]
                #Vertices[c(2,3),] <- temp_table_diff_network_sp[c(1,2),c(4,6)]
                Vertices$Color <- c("#87CEFA", "#FF7F24", "#228B22")
                ### edges
                Edges <- temp_table_diff_network_sp[,c(3,4,9,10)]
                Edges$significance <- NA
                Edges$significance[which(Edges$PHIst.p.val <= 0.05)] <- "1"
                Edges$significance[which(Edges$PHIst.p.val > 0.05)] <- "3"
                temp_net <- graph.data.frame(directed = F,d = Edges, vertices = Vertices)
                V(temp_net)$size <- V(temp_net)$Size*25
                V(temp_net)$color <- V(temp_net)$Color
                E(temp_net)$width <- E(temp_net)$PHIst*50
                plot(temp_net, vertex.frame.color=Vertices$Color, edge.color="black", edge.lty=as.numeric(Edges$significance),
                     layout=layout_on_grid, vertex.label=Vertices$Time_period, vertex.label.family="Helvetica", vertex.label.color="black")
                mtext(species[sp], side=2, las=2, line=-3)
                dev.off()
                }
}
setwd("/Users/afr/Desktop/PhD/Thesis/Chapter3/Ch3/Ch3_figures/Ch3_figures_data/Ch3_phist_triangles_data")
ch3_periods_fst.triangles(temp_out_diff_biome_pairwise = "ch3_diff_pairwise.txt")
temp_out_diff_biome_pairwise <- "ch3_diff_pairwise.txt"
temp_out_div_biome_between <- "ch3_div_between.txt"
