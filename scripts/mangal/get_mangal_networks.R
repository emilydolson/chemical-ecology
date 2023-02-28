#export PKG_CONFIG_PATH=$(pkg-config --variable pc_path pkg-config)${PKG_CONFIG_PATH:+:}${PKG_CONFIG_PATH}
#remotes::install_github("ropensci/rmangal")
#install.packages('tidygraph')

library("rmangal")
library(tidygraph)

#diet, Diet
networks <- search_networks(query = 'diet') %>% get_collection
for (i in 1:length(networks)){
    start_nodes <- c()
	for (k in 1:5){
		net_tbl <- networks[[i]]$interactions

		if ((net_tbl$attribute.name[1] == "presence/absence") | (net_tbl$network_id[1] > 4000)){
			next
		}

		print(net_tbl$network_id[1])
		print(unique(net_tbl$type))

		normalize <- FALSE
		if (max(net_tbl$value) > 1){
		    normalize <- TRUE
		}

		all_nodes_temp <- unique(c(net_tbl$node_from, net_tbl$node_to))
		len_all_nodes <- length(all_nodes_temp)
		all_nodes <- c(1:len_all_nodes)
		names(all_nodes) <- sort(all_nodes_temp)
		nodes_from <- all_nodes[as.character(net_tbl$node_from)]
		nodes_to <- all_nodes[as.character(net_tbl$node_to)]
		adj_matrix <- matrix(replicate(len_all_nodes, 0), nrow=len_all_nodes, ncol=len_all_nodes)

		for (j in 1:nrow(net_tbl)){
		    value <- net_tbl$value[j]
			if (normalize){
				adj_matrix[nodes_to[j], nodes_from[j]] <- (value - min(net_tbl$value)) / (max(net_tbl$value) - min(net_tbl$value))
			}
			else{
		    	adj_matrix[nodes_to[j], nodes_from[j]] <- value
			}
			if (((net_tbl$type[j] == 'predation') | (net_tbl$type[j] == 'herbivory') | (net_tbl$type[j] == 'parasitism')) & (identical(nodes_to[j], nodes_from[j]) == 0)){
				neg_value <- rnorm(1, adj_matrix[nodes_to[j], nodes_from[j]], 0.1)
				if (neg_value > 0){
					neg_value <- -neg_value
				}
				adj_matrix[nodes_from[j], nodes_to[j]] <- neg_value
			}
		}

		attribute_name <- net_tbl$attribute.name[1]
		attribute_name <- gsub("/", "-", attribute_name)
		attribute_name <- gsub(" ", "-", attribute_name)
		network_id <- net_tbl$network_id[1]
		write.table(adj_matrix, paste(network_id, "_", attribute_name, "_", k, ".csv", sep=""), sep=",", row.names=FALSE, col.names=FALSE)
		
		#start bipartite conversion
		set_A <- c()
		set_B <- c()
		start_node <- sample(all_nodes[! all_nodes %in% start_nodes], 1)
		start_nodes <- c(start_nodes, start_node)
		set_A <- c(set_A, start_node)
		uni_nodes <- all_nodes[all_nodes != start_node]
		
		for (n in sample(uni_nodes)) {
            if (adj_matrix[start_node, n] == 0 & adj_matrix[n, start_node] == 0) {
                set_B <- c(set_B, n)
                uni_nodes <- uni_nodes[uni_nodes != n]
                break
            }
		}
		if (length(set_B) == 0) {
            n <- sample(uni_nodes, 1)
            set_B <- c(set_B, n)
            uni_nodes <- uni_nodes[uni_nodes != n]
		}
		
        while (length(uni_nodes) > 0) {
            if (length(uni_nodes) > 1) {
                sample_node <- sample(uni_nodes, 1)
            } else {
                sample_node <- uni_nodes
            }
            
            sum_A <- 0
            for (n in set_A) {
                if (adj_matrix[sample_node, n] != 0 | adj_matrix[n, sample_node] != 0) {
                    sum_A <- sum_A + 1
                }
            }
            sum_B <- 0
            for (n in set_B) {
                if (adj_matrix[sample_node, n] != 0 | adj_matrix[n, sample_node] != 0) {
                    sum_B <- sum_B + 1
                }
            }
            
            if (sum_A > sum_B) {
                set_A <- c(set_A, sample_node)
            } else if (sum_B > sum_A) {
                set_B <- c(set_B, sample_node)
            } else {
                if (runif(1, 0, 1) > 0.5) {
                    set_A <- c(set_A, sample_node)
                } else {
                    set_B <- c(set_B, sample_node)
                }
            }
            
            uni_nodes <- uni_nodes[uni_nodes != sample_node]
        }
		
		bi_adj <- adj_matrix
		for (a1 in set_A) {
            for (a2 in set_A) {
                if (a1 != a2) {
                    bi_adj[a1,a2] <- 0
                    bi_adj[a2,a1] <- 0
                }
            }
		}
        for (b1 in set_B) {
            for (b2 in set_B) {
                if (b1 != b2) {
                    bi_adj[b1,b2] <- 0
                    bi_adj[b2,b1] <- 0
                }
            }
		}
		
		write.table(bi_adj, paste(network_id, "_", attribute_name, "_", k, "_bi", ".csv", sep=""), sep=",", row.names=FALSE, col.names=FALSE)
	}
}
