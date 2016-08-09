ch3_biomes_fst.network <- function(temp_out_diff_biome_pairwise){
        require(igraph)
        require(plotrix)
        #par(mfrow=c(4,3), mar=c(2,2,0,0))
        temp_table_diff_network <- read.table(temp_out_diff_biome_pairwise, header = T, stringsAsFactors = F )
        species <- unique(temp_table_network$sp)
        unique_biomes <- sort(unique(c(temp_table_diff_network$strata.1, temp_table_diff_network$strata.2)))
        bio_colors <- as.data.frame(matrix(nrow=32, ncol=2))
        colnames(bio_colors) <- c("biome", "bio_color")
        bio_colors$biome <- seq(0,31, 1)
        #bio_colors$bio_color <- terrain.colors(length(bio_colors$biome))
        #bio_colors$bio_color <- colorRampPalette(c("pink", "red"))(length(bio_colors$biome))
        bio_colors$bio_color <- rev(rich.colors(32))
        for(sp in seq_along(species)){
                temp_table_diff_network_sp <- temp_table_diff_network[which(temp_table_diff_network$sp == species[sp]),]
                time_periods <- unique(temp_table_diff_network_sp$Periods.period.)
                unique_biomes_sp <- sort(unique(c(temp_table_diff_network_sp$strata.1, temp_table_diff_network_sp$strata.2)))
                comb_biomes <- combn(unique_biomes_sp, m=2)
                for(bin in seq_along(time_periods)){
                        pdf(paste(species[sp], time_periods[bin], "_ch3_fst_network.pdf", sep=""), paper = "a4")
                        temp_table_diff_network_sp_biomes <- temp_table_diff_network_sp[which(temp_table_diff_network_sp$Periods.period. == time_periods[bin]),]
                        temp_table_diff_network_sp_biomes$paste <- paste(temp_table_diff_network_sp_biomes$strata.1, temp_table_diff_network_sp_biomes$strata.2, sep = "_")
                        ### vertices
                        Vertices <- as.data.frame(matrix(nrow =length(unique_biomes_sp) , ncol=3))
                        colnames(Vertices) <- c("Biome", "Color_in", "Color_out")
                        Vertices[,1] <- unique_biomes_sp
                        Vertices$Color_out <- bio_colors$bio_color[match(Vertices$Biome, bio_colors$biome)]
                        biome_sp_period <- unique(c(temp_table_diff_network_sp_biomes$strata.1, temp_table_diff_network_sp_biomes$strata.2))
                        matched <- match(biome_sp_period,Vertices$Biome)
                        Vertices$Color_in[matched] <- Vertices$Color_out[matched]
                        Vertices$Color_in[-matched] <- NA
                        Vertices$Size <- 2
                        ### edges
                        Edges <- as.data.frame(matrix(nrow=length(comb_biomes[1,]), ncol = 4))
                        colnames(Edges) <- c("node_one", "node_two", "phist", "p")
                        Edges$node_one <- comb_biomes[1,]
                        Edges$node_two <- comb_biomes[2,]
                        Edges$paste <- paste(Edges$node_one, Edges$node_two, sep = "_")
                        matched   <- match(temp_table_diff_network_sp_biomes$paste, Edges$paste)
                        Edges[matched,c(3,4)] <- temp_table_diff_network_sp_biomes[, c(10, 11)]
                        Edges[-matched,c(3,4)] <- 0
                        Edges$significance <- NA
                        Edges$significance <- ifelse(Edges$phist <= 0.05, "1", "2")
                        Edges$line_color <- ifelse(Edges$phist == 0, "white", "black")
                        #max_edge <- Edges()
                        temp_net <- graph.data.frame(directed = F,d = Edges, vertices = Vertices)
                        V(temp_net)$size <- V(temp_net)$Size*20
                        V(temp_net)$color <- V(temp_net)$Color_in
                        E(temp_net)$width <- E(temp_net)$phist*100
                        plot(temp_net, vertex.frame.color=Vertices$Color_out, edge.color=Edges$line_color, 
                             edge.lty=as.numeric(Edges$significance), main=paste(species[sp],
                                time_periods[bin], sep=" "), layout=layout_in_circle)
                        dev.off()
                        }
        }
        
        ch3_biomes.fst("/Users/afr/Desktop/Ch3_fst_period_data/")
        ch3_biomes.fst("/Users/afr/Desktop/Ch3_fst_period_data/")
        ### vertices
        Vertices <- temp_table_diff_network_sp_biomes[,c(4,6)]
        Vertices <- Vertices[-which(duplicated(Vertices)),]
        Vertices <- rbind(Vertices, c(NA,NA))
        Vertices[dim(Vertices)[1],] <- as.vector(temp_table_diff_network_sp_biomes[dim(temp_table_diff_network_sp_biomes)[1],c(5,7)])
        colrs <- terrain.colors(n=length(unique(nams)))
        Vertices$colors <- colrs
        ### edges
        Edges <- cbind(temp_table_diff_network_sp_biomes[,c(4,5,10,11)], temp_table_div_network_sp_biomes[,c(5,6)])
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
        