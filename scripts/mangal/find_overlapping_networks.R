#find networks with both "positive" and "negative" interaction types
#results: no networks with both mutualism and parasitism, nor mutualism and predation, nor mutualism and herbivory, nor mutualism and competition

library("rmangal")
library(tidygraph)
library(geosphere)

networks <- search_networks(query = 'diet') %>% get_collection
for (i in 1:length(networks)){
    for (j in 1:length(networks)){
        if (i != j){
            network1 <- networks[[i]]$network
            network2 <- networks[[j]]$network
            if (!is.na(network1$geom_type) && !is.na(network2$geom_type)){
                if (network1$geom_type == "Point" && network2$geom_type == "Point"){
                    distance <- distm(c(network1$geom_lon[[1]], network1$geom_lat[[1]]), c(network2$geom_lon[[1]], network2$geom_lat[[1]]), fun=distHaversine)[[1]]
                    if (distance < 100){
                        print(c(network1$network_id, network2$network_id, distance))
                    }
                }
            }
        }
    }
}

if (FALSE){
interactions_pos <- search_interactions(type = "mutualism")
interactions_neg <- search_interactions(type = "competition")
overlapping_ids <- Reduce(intersect, list(interactions_pos$network_id, interactions_neg$network_id))
print(overlapping_ids)
}