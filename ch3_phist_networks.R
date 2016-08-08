ch3_biomes_fst.network <- function(temp_out_diff_biome_pairwise, temp_out_div_biome_between){
        require(igraph)
        require(plotrix)
        #par(mfrow=c(4,3), mar=c(2,2,0,0))
        temp_table_diff_network <- read.table(temp_out_diff_biome_pairwise, header = T, stringsAsFactors = F )
        temp_table_div_network <- read.table(temp_out_div_biome_between, header = T, stringsAsFactors = F )
        species <- unique(temp_table_network$sp)
        for(sp in seq_along(species)){
                temp_table_diff_network_sp <- temp_table_diff_network[which(temp_table_diff_network$sp == species[sp]),]
                temp_table_div_network_sp <- temp_table_div_network[which(temp_table_div_network$sp == species[sp]),]
                time_periods <- unique(temp_table_network_sp$Periods.period.)
                for(bin in seq_along(time_periods)){
                        temp_table_diff_network_sp_biomes <- temp_table_diff_network_sp[which(temp_table_diff_network_sp$Periods.period. == time_periods[bin]),]
                        temp_table_div_network_sp_biomes <- temp_table_div_network_sp[which(temp_table_div_network_sp$Periods.period. == time_periods[bin]),]
                        cbind(temp_table_diff_network_sp_biomes[,c(4,5)], temp_table_div_network_sp_biomes[,c(3,4)])
                        #nams <- with(temp_table_diff_network_sp_biomes, unique(c(as.character(strata.1), (as.character(strata.2)))))
                        #dist_network <- with(temp_table_network_sp_biomes, structure(Fst, Size=length(nams), Labels = nams, Diag =FALSE, UPPER=FALSE, method="user", class="dist"))
                        #temp_graph <- qgraph(dist_network, layout="spring", vsize=20)
                        
                }
        }
        
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
        