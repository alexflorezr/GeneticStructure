### create haplotype networks
temp_align_biome <- read.dna("Mp_align_biome.fasta", format = "fasta")
temp_haplo <- haplotype(temp_align_biome)
temp_net <- haploNet(temp_haplo)
plot(temp_net, labels=T, size=attr(temp_net, "freq"), cex=0.5,show.mutation=F, bg=haplo_col, col="white", scale.ratio=3, asp=1)

