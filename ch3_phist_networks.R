ch3_biomes_fst.network <- function(temp_out_diff_biome_pairwise){
        require(igraph)
        require(plotrix)
        #temp_out_diff_biome_pairwise <- "ch3_diff_biome_pairwise_network.txt"
        #temp_out_diff_biome_overall <- "ch3_diff_biome_overall_network.txt"
        #seq_fossil_bio <- "ch3_seq_fossil_biome_df.txt"
        #par(mfrow=c(4,3), mar=c(2,2,0,0))
        temp_diff_pairwise_network <- read.table(temp_out_diff_biome_pairwise, header = T, stringsAsFactors = F, sep="\t" )
        #temp_seq_fossil_bio <- read.table(seq_fossil_bio, header = T, stringsAsFactors = F )
        #sp_all <- strsplit(temp_seq_fossil_bio$Species, split = "_")
        #sp_all_first <- substring(unlist(sp_all), first = 1, last = 1)
        #temp_seq_fossil_bio$Sp <- paste(sp_all_first[seq(1,length(sp_all_first), by=2)], sp_all_first[seq(2,length(sp_all_first), by=2)], sep = "")
        species <- unique(temp_diff_pairwise_network$sp)
        ### biome colors
        bio_colors <- as.data.frame(matrix(nrow=32, ncol=2))
        colnames(bio_colors) <- c("biome", "bio_color")
        bio_colors$biome <- biomes_all_sp
        tropical <- colorRampPalette(c("pink","#8B2323"))(4)
        arid <- colorRampPalette(c("orange", "orange4"))(4)
        warm <- colorRampPalette(c("#90EE90", "darkgreen"))(9)
        snow <- colorRampPalette(c("midnightblue","#1C86EE", "skyblue1", "#AEEEEE", "#6C7B8B"))(12)
        polar <- colorRampPalette(c("thistle", "thistle4"))(2) 
        bio_colors$bio_color <- c("white", tropical, arid, warm, snow, polar)
        plot(bio_colors$biome,y=rep(0.90,length(bio_colors$biome)) ,ylim=c(0.9, 1.1), bg=bio_colors$bio_color, pch=21, cex=3)
        for(sp in seq_along(species)){
                temp_diff_pairwise_network_sp <- temp_diff_pairwise_network[which(temp_diff_pairwise_network$sp == species[sp]),]
                #temp_seq_fossil_bio_sp <- temp_seq_fossil_bio[which(temp_seq_fossil_bio$Sp == species[sp]),]
                time_periods_sp <- unique(temp_diff_pairwise_network_sp$Periods)
                temp_biomes_all_periods_sp <- sort(unique(c(temp_diff_pairwise_network_sp$strata.1, temp_diff_pairwise_network_sp$strata.2)))
                comb_biomes_all_periods_sp <- combn(temp_biomes_all_periods_sp, m=2)
                for(bin in seq_along(time_periods)){
                        pdf(paste(species[sp], time_periods[bin], "_ch3_fst_network.pdf", sep=""), paper = "a4")
                        temp_diff_pairwise_network_sp_biomes <- temp_diff_pairwise_network_sp[which(temp_diff_pairwise_network_sp$Periods == time_periods[bin]),]
                        temp_diff_pairwise_network_sp_biomes$paste <- paste(temp_diff_pairwise_network_sp_biomes$strata.1, temp_diff_pairwise_network_sp_biomes$strata.2, sep = "_")
                        ### vertices
                        Vertices <- as.data.frame(matrix(nrow =length(temp_biomes_all_periods_sp) , ncol=3))
                        colnames(Vertices) <- c("Biome", "Color_in", "Color_out")
                        Vertices[,1] <- temp_biomes_all_periods_sp
                        Vertices$Color_out <- bio_colors$bio_color[match(Vertices$Biome, bio_colors$biome)]
                        biome_sp_period <- unique(c(temp_diff_pairwise_network_sp_biomes$strata.1, temp_diff_pairwise_network_sp_biomes$strata.2))
                        matched <- match(biome_sp_period,Vertices$Biome)
                        Vertices$Color_in[matched] <- Vertices$Color_out[matched]
                        Vertices$Color_in[-matched] <- "white"
                        Vertices$Size <- 2
                        ### edges
                        Edges <- as.data.frame(matrix(nrow=length(comb_biomes_all_periods_sp[1,]), ncol = 4))
                        colnames(Edges) <- c("node_one", "node_two", "phist", "p")
                        Edges$node_one <- comb_biomes_all_periods_sp[1,]
                        Edges$node_two <- comb_biomes_all_periods_sp[2,]
                        Edges$paste <- paste(Edges$node_one, Edges$node_two, sep = "_")
                        matched   <- match(temp_diff_pairwise_network_sp_biomes$paste, Edges$paste)
                        Edges[matched,c(3,4)] <- temp_diff_pairwise_network_sp_biomes[, c(10, 11)]
                        Edges[-matched,c(3,4)] <- NA
                        Edges$significance <- NA
                        Edges$significance <- ifelse(Edges$phist <= 0.05, "black", "#D2B48C80")
                        #Edges$line_color <- ifelse(Edges$phist == 0, "grey", "black")
                        #max_edge <- Edges()
                        temp_net <- graph.data.frame(directed = F,d = Edges, vertices = Vertices)
                        V(temp_net)$size <- V(temp_net)$Size*15
                        V(temp_net)$color <- V(temp_net)$Color_in
                        E(temp_net)$width <- E(temp_net)$phist*35
                        plot(temp_net, vertex.frame.color=Vertices$Color_out, edge.color=Edges$significance, 
                              main=paste(species[sp],time_periods[bin], sep=" "),
                             layout=layout_in_circle, vertex.label=NA)
                        dev.off()
                        }
        }
}
ch3_biomes_fst.network("ch3_diff_biome_pairwise_network.txt")
        edge.lty=as.numeric(Edges$significance)
        ch3_biomes.fst("/Users/afr/Desktop/Ch3_fst_period_data/")
        ch3_biomes.fst("/Users/afr/Desktop/Ch3_fst_period_data/")
        ### vertices
        Vertices <- temp_diff_pairwise_network_sp_biomes[,c(4,6)]
        Vertices <- Vertices[-which(duplicated(Vertices)),]
        Vertices <- rbind(Vertices, c(NA,NA))
        Vertices[dim(Vertices)[1],] <- as.vector(temp_diff_pairwise_network_sp_biomes[dim(temp_diff_pairwise_network_sp_biomes)[1],c(5,7)])
        colrs <- terrain.colors(n=length(unique(nams)))
        Vertices$colors <- colrs
        ### edges
        Edges <- cbind(temp_diff_pairwise_network_sp_biomes[,c(4,5,10,11)], temp_table_div_network_sp_biomes[,c(5,6)])
        Edges$edge_colors <- color.scale(Edges$dA, extremes = c("lightgrey", "black"))
        Edges$significance <- NA
        Edges$significance[which(Edges$PHIst.p.val <= 0.05)] <- "1"
        Edges$significance[which(Edges$PHIst.p.val > 0.05)] <- "2"
        temp_net <- graph.data.frame(directed = F,d = Edges, vertices = Vertices)
        V(temp_net)$size <- V(temp_net)$n.1*2
        V(temp_net)$color <- V(temp_net)$colors
        V(temp_net)$vertex.frame.color <- V(temp_net)$colors
        E(temp_net)$width <- E(temp_net)$PHIst*100
        E(temp_net)$arrow.size <- 2
        
        
        temp_net <- simplify(temp_net, remove.multiple = F, remove.loops = T) 
        
        plot(temp_net, vertex.frame.color=colrs, edge.color="black", edge.lty=as.numeric(Edges$significance))
        
        
        # Generate colors base on media type:
        colrs <- rainbow(n=length(unique(nams)))
        V(temp_net)$color <- colrs[V(temp_net)$media.type]
        
        # Compute node degrees (#links) and use that to set node size:
        deg <- degree(temp_net, mode="all")
        V(temp_net)$size <- deg*3
        # We could also use the audience size value:
        V(temp_net)$size <- V(temp_net)$n.1
        
        # The labels are currently node IDs.
        # Setting them to NA will render no labels:
        V(net)$label <- NA
        
        # Set edge width based on weight:
        E(temp_net)$width <- E(temp_net)$Fst*210
        
        #change arrow size and edge color:
        E(net)$arrow.size <- .2
        E(net)$edge.color <- "gray80"
        E(net)$width <- 1+E(net)$weight/12
        plot(temp_net) 
        
        temp_table_diff_network_sp_biomes$Fst, temp_table_diff_network_sp_biomes$PHIst
        
        
        nodes2 <- read.csv("Dataset2-Media-User-Example-NODES.csv", header=T, as.is=T)
        links2 <- read.csv("Dataset2-Media-User-Example-EDGES.csv", header=T, row.names=1)
        