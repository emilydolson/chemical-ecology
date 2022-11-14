#export PKG_CONFIG_PATH=$(pkg-config --variable pc_path pkg-config)${PKG_CONFIG_PATH:+:}${PKG_CONFIG_PATH}
#remotes::install_github("ropensci/rmangal")
#install.packages('tidygraph')

library("rmangal")
library(tidygraph)

networks <- search_networks(query = 'diet') %>% get_collection
for (i in 1:length(networks)){
    net_tbl <- networks[[i]]$interactions

    normalize <- FALSE
    if (max(net_tbl$value) > 1){
        normalize <- TRUE
        sum_values <- sum(net_tbl$value)
    }

    nodes_from <- as.numeric(factor(net_tbl$node_from))
    nodes_to <- as.numeric(factor(net_tbl$node_to))
    len_all_nodes <- length(unique(c(net_tbl$node_from, net_tbl$node_to)))
    adj_matrix <- matrix(replicate(len_all_nodes, 0), nrow=len_all_nodes, ncol=len_all_nodes)

    for (j in 1:nrow(net_tbl)){
        value <- net_tbl$value[j]
        if (net_tbl$type[j] == "predation"){
            value <- -value
        }

        if (normalize){
            adj_matrix[nodes_from[j], nodes_to[j]] <- value / sum_values
        }
        else {
            adj_matrix[nodes_from[j], nodes_to[j]] <- value
        }
    }

    attribute_name <- net_tbl$attribute.name[1]
    attribute_name <- gsub("/", "-", attribute_name)
    attribute_name <- gsub(" ", "-", attribute_name)
    network_id <- net_tbl$network_id[1]
    write.table(adj_matrix, paste(network_id, "_", attribute_name, ".csv", sep=""), sep=",", row.names=FALSE, col.names=FALSE)
}
