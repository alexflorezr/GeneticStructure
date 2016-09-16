### Biomes using Bio4_CO2
setwd("/Users/afr/Desktop/PhD/Thesis/Chapter3/Ch3/Ch3_data/Kop50k/")
Kop_Holartic <- stack(dir())
e <- extent(-180,180,30, 90)
kop_Hol_croped <- crop(Kop_Holartic, e)
### Cell area
cell_area <- raster("/Users/afr/Desktop/PhD/Collaborations/Rangel/Area_per_cell/image_area.txt")
### 
#time_raster <- rev(substr(names(kop_Hol_croped), start = 2, stop =4 ))
nrows <- prod(dim(kop_Hol_croped))
time_kop_area <- as.data.frame(matrix(nrow = nrows, ncol = 4))
colnames(time_kop_area) <- c("Time", "Koppen", "Area", "Cell_number")
row_start <- 1
row_end <- prod(dim(kop_Hol_croped)[c(1,2)])
for (raster in 1:dim(kop_Hol_croped)[3]){
        temp_time <- as.numeric(substr(names(kop_Hol_croped[[raster]]), start = 2, stop = 4))
        total_cells <- prod(dim(kop_Hol_croped[[raster]]))
        time_kop_area[row_start:row_end,c(4,2)] <- extract(kop_Hol_croped[[raster]], e, cellnumbers=T)
        time_kop_area[row_start:row_end,1] <- temp_time*1000
        time_kop_area[row_start:row_end,3] <- extract(cell_area, time_kop_area[row_start:row_end,4])
        row_start <- row_start+total_cells
        row_end <- row_end+total_cells
}
write.table(time_kop_area, file = "time_kop_area.txt",row.names = F, sep = "\t")

time <- c(seq(50000, 22000, -2000), seq(21000, 0, -1000))
time_total_areaXbiome <- data.frame(matrix(nrow=31*37, ncol=3))
colnames(time_total_areaXbiome) <- c("Time", "Koppen", "Total_area")
row_start <- 1
row_end <- 31
for (bin in seq_along(time)){
        temp_time <- time_kop_area[which(time_kop_area$Time == time[bin]),]
        time_total_areaXbiome[row_start:row_end,1] <- time[bin]
        time_total_areaXbiome[row_start:row_end,c(2,3)] <- aggregate(temp_time$Area, by=list(biome=temp_time$Koppen), FUN=sum)[2:32,]
        row_start <- row_end+1    
        row_end <- row_end+31
        }
 
time_total_area_ordered <- time_total_areaXbiome[order(time_total_areaXbiome$Koppen, decreasing = T),]
my_fun=function(vec){ as.numeric(vec[3]) / sum(time_total_area_ordered$Total_area[which(time_total_area_ordered$Koppen == as.numeric(vec[2]))]) *100 }
time_total_area_ordered$prop=apply(time_total_area_ordered, 1 , my_fun)
### area plot ##
pdf(file="biome_area_time.pdf", paper = "a4")
ggplot(time_total_area_ordered, aes(x=Time, y=Total_area, fill=as.factor(time_total_area_ordered$Koppen))) + 
        geom_area(colour="white", size=.2, alpha=0.99) +
        scale_fill_manual(values=bio_colors) +
        theme_bw()+
        theme(panel.grid.major = element_blank(),
              panel.background = element_blank()) 
dev.off()        
getwd()

### barplots ###
for_barplots <- matrix(nrow=31, ncol=37)
rownames(for_barplots) <-  seq(31,1,-1)
colnames(for_barplots) <- c(seq(50000, 22000, -2000), seq(21000, 0, -1000))
for (col in 1:dim(for_barplots)[2]){
        koppen_raw <- time_total_area_ordered$Koppen[which(time_total_area_ordered$Time == as.numeric(colnames(for_barplots)[col]))]
        if(sum(koppen_raw == rownames(for_barplots)) == 31){
                for_barplots[,col] <- time_total_area_ordered$Total_area[which(time_total_area_ordered$Time == as.numeric(colnames(for_barplots)[col]))]
        }    
}
par(lwd=0.4)
tropical <- colorRampPalette(c("pink","#8B2323"))(4)
arid <- colorRampPalette(c("orange", "orange4"))(4)
warm <- colorRampPalette(c("#90EE90", "darkgreen"))(9)
snow <- colorRampPalette(c("midnightblue","#1C86EE", "skyblue1", "#AEEEEE", "#6C7B8B"))(12)
polar <- colorRampPalette(c("thistle", "thistle4"))(2) 
bio_colors <- c(tropical, arid, warm, snow, polar)
barplot(for_barplots,col=rev(bio_colors))
