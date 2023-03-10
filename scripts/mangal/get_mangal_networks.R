#export PKG_CONFIG_PATH=$(pkg-config --variable pc_path pkg-config)${PKG_CONFIG_PATH:+:}${PKG_CONFIG_PATH}
#remotes::install_github("ropensci/rmangal")
#install.packages('tidygraph')

library("rmangal")
library(tidygraph)

#diet, Diet
networks <- search_networks(query = 'diet') %>% get_collection
for (i in 1:length(networks)){
	net_tbl <- networks[[i]]$interactions

	if ((net_tbl$attribute.name[1] == "presence/absence") | (net_tbl$network_id[1] > 4000) | length(unique(c(net_tbl$node_from, net_tbl$node_to))) > 20){
		print('next')
		next
	}

	print(net_tbl$network_id[1])
	print(unique(net_tbl$type))

	for (k in 1:5){
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
	}
}
