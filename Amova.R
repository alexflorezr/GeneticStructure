library(ade4)
data(humDNAm)
amovahum <- amova(humDNAm$samples, sqrt(humDNAm$distances), humDNAm$structures)
### create sample object
nchar(temp_points_hap_DB$Sequence)
temp_samples <- table(temp_points_hap_DB$Haplotype,temp_points_hap_DB$Biome)
temp_samples_df <- as.data.frame(matrix(nrow=81, ncol=15))
for(row in seq_along(temp_samples_df[,1])){
        temp_samples_df[row,] <- temp_samples[row,]
}
colnames(temp_samples_df) <- colnames(temp_samples)
rownames(temp_samples_df) <- rownames(temp_samples)
### create distance object
temp_dist <- dist.dna(temp_hap_biome, model = "raw")
attr(temp_dist, "Labels") <- as.character(seq(1, dim(temp_hap_biome)[1]))

### create structure object
temp_structure <- as.data.frame(colnames(temp_samples), stringsAsFactors=F)
colnames(temp_structure) <- "Biomes"
temp_structure$Biomes <- ifelse(as.numeric(temp_structure$Biomes) <= 16, "bio_tens", "bio_twenties")
### amova
temp_amova <- amova(temp_samples_df, sqrt(temp_dist), temp_structure)
