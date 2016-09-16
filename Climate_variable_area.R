setwd("/Users/afr/Desktop/PhD/Thesis/Chapter3/Ch3/Ch3_data/Temp50k/")
temp <- stack(dir())
setwd("/Users/afr/Desktop/PhD/Thesis/Chapter3/Ch3/Ch3_data/Prec50k/")
Climate_variable <- stack(dir())
cell_area <- raster("/Users/afr/Desktop/PhD/Collaborations/Rangel/Area_per_cell/image_area.txt")
Temp_colors <- colorRampPalette(c( "#CD0000", "#FF4040","#FFC1C1", "#87CEFA","#00B2EE", "#1874CD"))(10)
Prec_colors <- colorRampPalette(c( "#228B22", "#698B69","#C1FFC1", "#DEB887","#8B7355", "#8B4513"))(10)
Temp_colors <- Prec_colors
Climate_Variable_area <- function(Climate_variable, Colors, cell_area, file_name){
require(ggplot2)
Tmp <- Climate_variable
e <- extent(-180,180,30, 90)
Tmp_hol <- crop(Tmp, e)
### create temperature categories ###
Tmp_hol_categories <- Tmp_hol
Tmp_quantile <- quantile(as.vector(Tmp_hol_categories), probs=seq(0.1, 1, 0.1))
Tmp_categories <- 1:10
for(bin in 1:dim(Tmp_hol_categories)[3]){
        for (quant in seq_along(Tmp_quantile)){
                if(quant == 1){
                        Tmp_hol_categories[[bin]][which(values(Tmp_hol[[bin]]) <= Tmp_quantile[quant])] <- Tmp_categories[quant]
                }
                if(is.element(quant, 2:10)){
                        Tmp_hol_categories[[bin]][which(values(Tmp_hol[[bin]]) > Tmp_quantile[quant-1] & values(Tmp_hol[[bin]]) <= Tmp_quantile[quant])] <- Tmp_categories[quant]  
                }
        }
}
### create a table with the temperature categories and area ###
nrows <- prod(dim(Tmp_hol_categories))
temp_cat_area <- as.data.frame(matrix(nrow = nrows, ncol = 4))
colnames(temp_cat_area) <- c("Time", "Temp_cat", "Area", "Cell_number")
row_start <- 1
row_end <- prod(dim(Tmp_hol_categories)[c(1,2)])
for (raster in 1:dim(Tmp_hol_categories)[3]){
        temp_time <- as.numeric(substr(names(Tmp_hol_categories[[raster]]), start = 2, stop = 4))
        total_cells <- prod(dim(Tmp_hol_categories[[raster]]))
        temp_cat_area[row_start:row_end,c(4,2)] <- extract(Tmp_hol_categories[[raster]], e, cellnumbers=T)
        temp_cat_area[row_start:row_end,1] <- temp_time*1000
        temp_cat_area[row_start:row_end,3] <- extract(cell_area, temp_cat_area[row_start:row_end,4])
        row_start <- row_start+total_cells
        row_end <- row_end+total_cells
}

### estimate the total area per temperature category ###
time <- c(seq(50000, 22000, -2000), seq(21000, 0, -1000))
time_area_tempCategory <- data.frame(matrix(nrow=10*37, ncol=3))
colnames(time_area_tempCategory) <- c("Time", "TempCategory", "Total_area")
row_start <- 1
row_end <- 10
for (bin in seq_along(time)){
        temp_time <- temp_cat_area[which(temp_cat_area$Time == time[bin]),]
        time_area_tempCategory[row_start:row_end,1] <- time[bin]
        temp_aggregate <- aggregate(temp_time$Area, by=list(TempCategory=temp_time$Temp_cat), FUN=sum)
        if(dim(temp_aggregate)[1] < length(Tmp_categories)){
                time_area_tempCategory[row_start:row_end,c(2,3)] <- rbind(c(1,0), temp_aggregate)
        }else{
                time_area_tempCategory[row_start:row_end,c(2,3)] <- temp_aggregate
        }
        row_start <- row_end+1    
        row_end <- row_end+10
}
time_area_tempCategory_ordered <- time_area_tempCategory[order(time_area_tempCategory$TempCategory, decreasing = F),]
### area plot ##
#pdf(file=file_name, paper = "a4")
Climate_var_label <- vector()
for (quant in seq_along(Tmp_quantile)){
        Climate_var_label[1] <-  c(paste("<", round(Tmp_quantile[1], digits = 2)))
        if(is.element(quant, 2:10)){
                Climate_var_label[quant]   <- paste(round(Tmp_quantile[quant-1], digits = 2), " - ",round(Tmp_quantile[quant], digits = 2))   
        }
        }

ggplot(time_area_tempCategory_ordered, aes(x=Time, y=Total_area, fill=as.factor(time_area_tempCategory_ordered$TempCategory))) + 
        geom_area(colour="white", size=.2, alpha=0.8) +
        scale_fill_manual(values=rev(Temp_colors), breaks=1:10,
                          labels=Climate_var_label) +
        scale_alpha_discrete(range=seq(0,1,length.out=10)) +
        theme_bw()+
        theme(panel.grid.major = element_blank(),
              panel.background = element_blank())+
        guides(fill = guide_legend(reverse=TRUE, title="Quantiles"))

### for polygons ##
id <- paste(rep(1:10, each=36), rep(1:36, 10), sep = ".")
ids <- rep(ids, each=4)
for_polygons <- data.frame(matrix(nrow=length(ids), ncol=3))
colnames(for_polygons) <- c("ids", "x", "y")
for_polygons$ids <- ids
time_bins <- unique(time_area_tempCategory_ordered$Time)
for_x <- vector()
for(time in seq_along(time_bins)[-37]){
        for_x <- c(for_x, c(time_bins[time], time_bins[time+1], time_bins[time],time_bins[time+1]))
}
for_polygons$x <- rep(for_x, 10)
for_y <- vector()
categories <- unique(time_area_tempCategory_ordered$TempCategory)
for(cat in seq_along(categories)){
        temp_cat <- time_area_tempCategory_ordered[which(time_area_tempCategory_ordered$TempCategory == categories[cat]),]
        temp_cat1 <- time_area_tempCategory_ordered[which(time_area_tempCategory_ordered$TempCategory == categories[cat+1]),]
        for(time in seq_along(time_bins)[-37]){
                if(cat == 1){
                        for_y <- c(for_y, c(0,0, temp_cat$Total_area[time], temp_cat$Total_area[time]))
                }else{
                        for_y <- c(for_y, c(temp_cat$Total_area[time], temp_cat$Total_area[time], temp_cat1$Total_area[time+1], temp_cat1$Total_area[time+1]))
                }
        }
}
for_polygons$y <- for_y
### values(to use transparency for the fossil record)
values <- data.frame(
        ids = ids,
        value = seq(1, 1440, 1)
)
for_polygons_merged <- merge(values, for_polygons, by=c("ids"))
p <- ggplot(for_polygons_merged, aes(x=x, y=y)) + geom_polygon(aes(fill=value, group=ids))
p
ggplot(for_polygons_merged, aes(x=x, y=y)) + 
        geom_polygon() 
        scale_fill_manual(values=rev(Temp_colors), breaks=1:10,
                          labels=Climate_var_label) +
        scale_alpha_discrete(range=seq(0,1,length.out=10)) +
        theme_bw()+
        theme(panel.grid.major = element_blank(),
              panel.background = element_blank())+
        guides(fill = guide_legend(reverse=TRUE, title="Quantiles"))

#dev.off()        
}

Climate_Variable_area(Climate_variable = Climate_variable, Colors = Temp_colors, cell_area = cell_area, file_name = "temp_cat_area.pdf")
