library(maps)
library(mapproj)
library(plotrix)
library(stringr)
pdf("ch3_figure_1_map.pdf", paper = "a4")
map("world",proj="azequidist", resolution = 0.2, ylim = c(33,90), lwd=1.5, mar=c(0,2,2,2), interior = F) 
points(mapproject(as.numeric(ch3_data$Longitude), as.numeric(ch3_data$Latitude)), 
       bg="#8B735590", 
       cex=1.1, pch=21, 
       col="black", lwd=0.5)
draw.circle(0,0,1.17, lwd=2)
mtext(paste("n = ", dim(ch3_data)[1], sep = ""), line=1)
dev.off()
getwd()
setwd("/Users/afr/Desktop/")
kk <- read.delim(file.choose(), header=T, stringsAsFactors = T, sep="\t")
kk_mp <- kk[which(kk$Species == "Mammuthus_primigenius"),]
max(kk_mp$Cli_vel_tmp)
for(time in seq_along(unique(kk_mp$Time_bin))){
        kk_mp_time <- kk_mp[which(kk_mp$Time_bin == unique(kk_mp$Time_bin)[time]),]
        kk_median <- median(kk_mp_time$Cli_vel_tmp, na.rm = T)
        print(paste(unique(kk_mp$Time_bin)[time], round(kk_median, digits = 3), sep="    "))
        
}