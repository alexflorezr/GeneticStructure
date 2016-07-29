# Upload all the packages needed for the script
require(ape)
require(strataG)
require(sidier)
require(pegas)
require(seqinr)
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
### Biomes using Koppen
setwd("/Users/afr/Desktop/PhD/Thesis/Chapter3/Ch3/Ch3_data/Kop50k/")
Kop <- stack(dir())
### Import the database
DATABASE <- read.delim(file.choose(), header = T, sep = "\t", stringsAsFactors = F)
### Extract the values of climate (temp and prec) and biomes for every species
Periods<- c("pre_LGM", "LGM", "pos_LGM")
Vspecies <- "Mammuthus_primigenius"
for(sp in seq_along(Vspecies)){
        temp_hap_DB <- DATABASE[which(DATABASE$Species == Vspecies[sp] &  nchar(DATABASE$Sequence) > 1),]
        temp_points_hap_DB <- as.data.frame(matrix(nrow=nrow(temp_hap_DB), ncol=12))
        colnames(temp_points_hap_DB) <- c("Vspecies[sp]","Longitude", "Latitude", "Time_sample", "Time_bin", "Layer", "Sequence","Temp","Prec","Biome_Bio","Biome_Kop", "Haplotype")
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
                temp_points_hap_DB[row,11]<- extract(Kop, layer=as.numeric(temp_points_hap_DB[row,6]), nl=1, y=matrix(as.numeric( temp_points_hap_DB[row,c(2,3)]), nrow=1, ncol=2))
                
        }
}
setwd("/Users/afr/Desktop/PhD/Thesis/Chapter3/Ch3/Ch3_results/")
write.table(temp_points_hap_DB, file=paste(Vspecies, "points.txt", sep="_"), sep = "\t")
### IMPORTANT SAFETY POINT backup_hap <- temp_points_hap_DB
### Create a DNA object including all the sequences in the alignment
        write.fasta(as.list(temp_points_hap_DB$Sequence), names = paste(seq(1,length.out = length(temp_points_hap_DB$Longitude)), as.numeric(temp_points_hap_DB$Time_bin)/1000, temp_points_hap_DB$Biome_Bio, temp_points_hap_DB$Biome_Kop, sep="_"), file.out = "Mp_hap_biome")
        setwd("/Users/afr/Desktop/PhD/Thesis/Chapter3/Ch3/Ch3_data/Fasta_files/")
        temp_align_biome <- read.dna( "Mp_alig_no_N.fasta", format = "fasta", as.character = FALSE, as.matrix=NULL)
        temp_haps <- FindHaplo(align=temp_align_biome)
        colnames(temp_haps) <-c("UniqueInd", "haplotype")
        temp_data_stats <- as.data.frame(matrix(nrow=dim(temp_align_biome)[1], ncol=5))
        colnames(temp_data_stats) <- c("Index", "Time", "Period", "Bio", "Kop")
        temp_data_stats$Index <- dimnames(temp_align_biome)[[1]]
        time_strings <- strsplit(dimnames(temp_align_biome)[[1]], split="_")
        temp_data_stats$Index <- dimnames(temp_align_biome)[1][[1]]
        temp_data_stats$Time <- unlist(time_strings)[2*seq(1,length(dimnames(temp_align_biome)[1][[1]])*2, by=2)]
        temp_data_stats$Period[which(as.numeric(temp_data_stats$Time) > 25)] <- Periods[1]
        temp_data_stats$Period[which(as.numeric(temp_data_stats$Time) < 25 & as.numeric(temp_data_stats$Time) >= 15)] <- Periods[2]
        temp_data_stats$Period[which(as.numeric(temp_data_stats$Time) < 15)] <- Periods[3]
        temp_data_stats$Bio <- unlist(time_strings)[2*seq(1,length(dimnames(temp_align_biome)[1][[1]])*2, by=2)+1]
        temp_data_stats$Kop <- unlist(time_strings)[2*seq(1,length(dimnames(temp_align_biome)[1][[1]])*2, by=2)+2]
        temp_haps_time_bio <- cbind(temp_data_stats, temp_haps)
        sp_haps  <- GetHaplo(align=temp_align_biome, saveFile=T, outname=paste("Mp","_Pops.fasta", sep=""), format="fasta", seqsNames="Inf.Hap") #haps are now a DNAbin
        sp_ghaps <- read.fasta(paste("Mp","_Pops.fasta", sep="")) #imports haps as gtypes
        sp_gtype <- gtypes(gen.data=data.frame(temp_haps_time_bio$UniqueInd,temp_haps_time_bio$Period,temp_haps_time_bio$haplotype),id.col=1,strata.col=2,locus.col=3,dna.seq=sp_ghaps)
### Create a DNA object for each time period
Periods<- c("pre_LGM", "LGM", "pos_LGM")
for(period in seq_along(Periods)){
                temp_fasta_stats_period <- temp_align_biome[match(temp_haps_time_bio$UniqueInd[temp_haps_time_bio$Period == Periods[period]], labels(temp_align_biome)),]
                haps_period <-FindHaplo(align=temp_fasta_stats_period, saveFile=F)
                haps_period <- as.data.frame(haps_period)
                colnames(haps_period)<-c("UniqueInd", "haplotype")
                haps_period$UniqueInd <- levels(haps_period$UniqueInd)
                temp_data_stats_period <- temp_data_stats[temp_data_stats$Period == Periods[period],]
                matched <- match(haps_period$UniqueInd, temp_data_stats_period[,1])
                temp_data_stats_period_m <- temp_data_stats_period[matched,]
                temp_haps_time_bio_period <- cbind(temp_data_stats_period_m, haps_period)
                        ### add a column to specify the biome group 
                        temp_haps_time_bio_period$Bio_clust[temp_haps_time_bio_period$Bio <= 15] <- "Tens"
                        temp_haps_time_bio_period$Bio_clust[temp_haps_time_bio_period$Bio > 15] <- "Twenties"
                sp_haps_period <-GetHaplo(align=temp_fasta_stats_period, saveFile=T, outname=paste("Mp","_Pops_",  Periods[period],".fasta", sep=""), format="fasta", seqsNames="Inf.Hap") #haps are now a DNAbin
                sp_ghaps_period <-read.fasta(paste("Mp","_Pops_",  Periods[period],".fasta", sep="")) #imports haps as gtypes
                sp_gtype_period <- gtypes(gen.data=data.frame(temp_haps_time_bio_period$UniqueInd,temp_haps_time_bio_period$Bio_clust,temp_haps_time_bio_period$haplotype),id.col=1,strata.col=2,locus.col=3,dna.seq=sp_ghaps_period)
}
if (pop_diff == TRUE){
        temp_pop_all <- pop.diff.test(sp_gtype)
        temp_pop_diff <-pop.diff.test(sp_gtype_period)
        #FST
        temp_FST <- as.vector(temp_pop_diff$overall$result[1,])
        temp_FST_all <- as.vector(temp_pop_all$overall$result[1,])
        #Fi-st
        temp_phist <- as.vector(temp_pop_diff$overall$result[2,])
        #chi squared
        Temp_chi2 <-  as.vector(temp_pop_diff$overall$result[3,])
        #Fixed differences
        temp_fixed <- fixed.differences(sp_gtype, count.indels = F,bases = c("a", "c", "g", "t"))$num.fixed[,3]
        temp_fixed_period <- fixed.differences(sp_gtype_period, count.indels = F,bases = c("a", "c", "g", "t"))$num.fixed[,3]
        #Nucleotide divergence
        temp_nuc_divergence_all <- nucleotide.divergence(sp_gtype)
        temp_nuc_divergence_period <- nucleotide.divergence(sp_gtype_period)
        temp_nuc_divergence <- temp_nuc_divergence_all$between$mean
        #Nei's DA
        temp_nei_DA <- temp_nuc_divergence_all$between$dA
        #Shared Haplotypes
        temp_shared <- shared.haps(sp_gtype)$shared.haps
}