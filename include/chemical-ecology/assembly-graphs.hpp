#pragma once

#include <iostream>
#include "emp/Evolve/World.hpp"
#include "emp/math/distances.hpp"
#include "emp/math/random_utils.hpp"
#include "emp/datastructs/Graph.hpp"
#include "emp/bits/BitArray.hpp"
#include "emp/data/DataNode.hpp"
#include "chemical-ecology/config_setup.hpp"
#include "emp/datastructs/map_utils.hpp"
#include <string>
#include <math.h>
#include "pagerank.h"

/**************************************
This file contains all code pertaining to the calculation of assembly graphs, 
fitness landscapes, and pageranks.

After the switch to stochastic models, this code is no longer in use, 
but it remains a possible area of exploration in the future.

Functions were called with the follwing 3 lines in the Run() fuction of a-eco:

//emp::Graph assemblyGraph = CalculateCommunityAssemblyGraph();
//emp::WeightedGraph wAssembly = calculateWeightedAssembly(assemblyGraph, config->PROB_CLEAR(), config->SEEDING_PROB());
//std::map<std::string, float> assembly_pr_map = CalculateWeightedPageRank(wAssembly);

***************************************/

emp::Graph CalculateCommunityAssemblyGraph() {
    // long long int is 64 bits, so this can handle
    // communities of up to 64 species
    long long unsigned int n_communities = pow(2, N_TYPES);
    struct community_info {
        bool stable;
        long long unsigned int transitions_to;
        community_info(bool s=false, long long int t=0) {
        stable = s;
        transitions_to = t;
        }
    };

    std::map<long long unsigned int, community_info> communities;

    // community_info t = community_info(true, 0);
    communities.emplace(0, community_info(true, 0));
    long long unsigned int n_stable = 1;

    for (long long unsigned int community_id = 1; community_id < pow(2,N_TYPES); community_id++) {
        std::bitset<64> temp(community_id);
        emp::BitArray<64> comm(temp);

        // Number of time steps to run test world for
        int time_limit = 25;
        world_t test_world;
        // world is 1x1
        test_world.resize(1);
        // Initialize empty test world
        test_world[0].resize(config->N_TYPES(), 0);
        for (int pos = comm.FindOne(); pos >= 0; pos = comm.FindOne(pos+1)) {
        // std::cout << comm.ToBinaryString() << " " << pos
        test_world[0][pos] = 1;
        }

        // Make a next_world for the test world to
        // handle simulated updates in our test environment
        world_t next_world;
        next_world.resize(1);
        next_world[0].resize(config->N_TYPES(), 0);

        // Run DoGrowth for 20 time steps to hopefully
        // get cell to equilibrium
        // (if cell is an organism, this is analogous to 
        // development)
        for (int i = 0; i < time_limit; i++) {
        DoGrowth(0, test_world, next_world);
        std::swap(test_world, next_world);
        std::fill(next_world[0].begin(), next_world[0].end(), 0);
        }

        emp::BitArray<64> new_comm;
        for (int i = 0; i < config->N_TYPES(); i++) { 
        if (test_world[0][i]){
            new_comm.Toggle(i);
        }
        }
        communities.emplace(community_id, community_info(comm==new_comm, new_comm.GetUInt64(0)));
        if (communities[community_id].stable) {
        n_stable++;
        }
    }

    for (auto c : communities) {
        //std::cout << c.first << " " << c.second.stable << " " << c.second.transitions_to << std::endl;

        if (!communities[c.second.transitions_to].stable) {
        n_stable++;
        communities[c.second.transitions_to].stable = true;
        }
    }

    emp::Graph g(n_stable);
    std::map<int, int> node_map;
    int curr_node = 0;
    for (long long unsigned int i = 0; i < pow(2,N_TYPES); i++) {
        // Would be great to use n_communities instead of
        // 64, but we don't know it at compile time
        // Will end up with a bunch of leading 0s.
        if (!communities[i].stable) {
        continue;
        }
        node_map[i] = curr_node;
        std::bitset<64> temp(i);
        std::string temp_str = temp.to_string();
        temp_str.erase(0, 64-N_TYPES);
        g.SetLabel(curr_node, temp_str);
        curr_node++;
    }

    for (long long unsigned int i = 0; i < pow(2,N_TYPES); i++) {
        if (!communities[i].stable) {
        continue;
        }
        std::bitset<64> temp(i);
        emp::BitArray<64> comm(temp);
        curr_node = node_map[i];

        emp::BitArray<64> missing_species = comm;
        missing_species.Toggle();
        for (int pos = missing_species.FindOne(); pos >= 0 && pos < config->N_TYPES(); pos = missing_species.FindOne(pos+1)) { 
        comm.Toggle(pos);
        long long unsigned int id = comm.GetUInt64(0);
        if (communities[id].stable) {
            if (curr_node != node_map[id]) {
            g.AddEdge(curr_node, node_map[id]);
            //std::cout << g.GetLabel(curr_node) << ", " << g.GetLabel(node_map[id]) << ", " << pos << std::endl;
            }
        } else {
            emp_assert(emp::Has(node_map, communities[id].transitions_to));
            if (curr_node != node_map[communities[id].transitions_to]) {
            g.AddEdge(curr_node, node_map[communities[id].transitions_to]);
            //std::cout << g.GetLabel(curr_node)  << ", " << g.GetLabel(node_map[communities[id].transitions_to])  << ", " << pos << std::endl;         
            } 
        }
        comm.Toggle(pos);
        }
    }

    return g;

}

emp::Graph CalculateCommunityLevelFitnessLandscape(std::string fitness_measure) {
    emp::Graph g = CalculateCommunityAssemblyGraph();
    struct fitnesses{
        double growth_rate;
        double biomass; 
        double heredity;
        double invasion_ability;
        double resiliance;
    };

    std::map<std::string, fitnesses> found_fitnesses;
    emp::vector<emp::Graph::Node> all_nodes = g.GetNodes();

    for(size_t curr_pos = 0; curr_pos < g.GetSize(); curr_pos++){
        emp::Graph::Node n = g.GetNode(curr_pos);
        std::string label = n.GetLabel();
        bool found = false;
        fitnesses curr_node_fitness;
        fitnesses adjacent_node_fitness;

        if(emp::Has(found_fitnesses, label)){
        found = true;
        curr_node_fitness = found_fitnesses[label];
        }
        if(found == false){
        emp::vector<double> community;
        for(char& c : label){
            community.push_back((double)c - 48);
        }
        //Need to reverse the vector, since the bit strings have reversed order from the world vectors
        std::reverse(community.begin(), community.end());
        CellData data = doGetFitness(community);
        //TODO Use equillib growth rate?
        curr_node_fitness.growth_rate = data.equilib_growth_rate;
        curr_node_fitness.biomass = data.biomass;
        curr_node_fitness.heredity = data.heredity;
        curr_node_fitness.invasion_ability = data.invasion_ability;
        curr_node_fitness.resiliance = n.GetDegree();
        found_fitnesses[label] = curr_node_fitness;
        }
        emp::BitVector out_nodes = n.GetEdgeSet();
        for(int pos = out_nodes.FindOne(); pos >= 0 && pos < out_nodes.size(); pos = out_nodes.FindOne(pos+1)) {
        bool adj_found = false;
        std::string adj_label = g.GetLabel(pos);
        if(emp::Has(found_fitnesses, adj_label)){
            adj_found = true;
            adjacent_node_fitness = found_fitnesses[adj_label];
        }
        if(adj_found == false){
            emp::vector<double> adj_community;
            for(char& c : adj_label){
            adj_community.push_back((double)c - 48);
            }
            //reverse community to get world vector
            std::reverse(adj_community.begin(), adj_community.end());
            CellData adj_data = doGetFitness(adj_community);
            adjacent_node_fitness.growth_rate = adj_data.equilib_growth_rate;
            adjacent_node_fitness.biomass = adj_data.biomass;
            adjacent_node_fitness.heredity = adj_data.heredity;
            adjacent_node_fitness.invasion_ability = adj_data.invasion_ability;
            adjacent_node_fitness.resiliance = g.GetDegree(pos);
            found_fitnesses[adj_label] = adjacent_node_fitness;
        }
        if(fitness_measure.compare("Biomass") == 0){
            if(curr_node_fitness.biomass > adjacent_node_fitness.biomass){
            g.SetEdge(curr_pos, pos, false);
            }
        }
        if(fitness_measure.compare("Growth_Rate") == 0){
            if(curr_node_fitness.growth_rate > adjacent_node_fitness.growth_rate){
            g.SetEdge(curr_pos, pos, false);
            }
        }
        if(fitness_measure.compare("Heredity") == 0){
            if(curr_node_fitness.heredity > adjacent_node_fitness.heredity){
            g.SetEdge(curr_pos, pos, false);
            }
        }
        if(fitness_measure.compare("Invasion_Ability") == 0){
            if(curr_node_fitness.invasion_ability > adjacent_node_fitness.invasion_ability){
            g.SetEdge(curr_pos, pos, false);
            }
        }
        if(fitness_measure.compare("Resiliance") == 0){
            //Need to do < here -- resiliance is the out degree. Lower out degree is better
            if(curr_node_fitness.resiliance < adjacent_node_fitness.resiliance){
            g.SetEdge(curr_pos, pos, false);
            }
        }
        }
    }
    return g;
}

double calcAdaptabilityScore(std::map<std::string, double> finalCommunities, std::string fitness_measure) {
    emp::Graph fitnessGraph = CalculateCommunityLevelFitnessLandscape(fitness_measure);
    std::map<std::string, float> fitness_pr_map = calculatePageRank(fitnessGraph);

    double fitness_score = 0;
    for(auto& [key, val] : finalCommunities){
        std::string node = key;
        float proportion = val;
        if (fitness_pr_map.find(node) != fitness_pr_map.end()) {
        fitness_score += (proportion * fitness_pr_map[node]);
        }
        else {
        fitness_score *= pow(10.0, -10.0);
        }
    }

    return fitness_score;
}

emp::WeightedGraph calculateWeightedAssembly(emp::Graph g, double prob_clear, double seed_prob){
    size_t num_communities = pow(2, N_TYPES);
    emp::WeightedGraph wAssembly(num_communities);

    emp::vector<emp::Graph::Node> all_nodes = g.GetNodes();
    for(emp::Graph::Node n: all_nodes){
        std::string label = n.GetLabel();
        std::size_t num = std::stoi(label, 0, 2);
        emp::BitVector out_nodes = n.GetEdgeSet();
        size_t out_degree = n.GetDegree();
        //weight is the probability of a seeding event of the needed type for the transition 
        double out_weight = seed_prob/(double)N_TYPES;
        for (int pos = out_nodes.FindOne(); pos >= 0 && pos < g.GetSize(); pos = out_nodes.FindOne(pos+1)){
        int adj_num = std::stoi(g.GetLabel(pos), 0, 2);
        wAssembly.AddEdge(num, adj_num, out_weight);
        }
        //cell must be not seeded and cleared
        double clear_weight = prob_clear*(1-out_weight);
        //original assembly node points to 0
        if (n.HasEdge(0)) {
        out_degree -= 1;
        }
        //Add an edge for clear probability
        wAssembly.AddEdge(num, 0, clear_weight);
        //Add self edge, to account for chance the cell is neither cleared nor seeded 
        double self_weight = 1 - clear_weight - (out_weight*out_degree);
        if (num == 0) {
        self_weight += clear_weight;
        }
        wAssembly.AddEdge(num, num, self_weight);
    }

    return wAssembly;
}

std::map<std::string, float> calculatePageRank(emp::Graph g) {
    Table t;

    t.set_trace(false);
    t.set_numeric(false);
    t.set_delim(" ");
    //t.read_graph(g);
    t.pagerank();

    std::map<std::string, float> map = t.get_pr_map();
    return map;
}

std::map<std::string, float> calculatePageRank(emp::WeightedGraph edge_weights) {
    Table t;

    t.set_trace(false);
    t.set_numeric(false);
    t.set_delim(" ");
    t.read_graph(edge_weights);
    t.set_edge_weights(edge_weights);
    t.weighted_pagerank();

    std::map<std::string, float> map = t.get_pr_map();
    return map;
}

std::map<std::string, float> CalculateWeightedPageRank(emp::WeightedGraph edge_weights, double alpha = 1) {
    //initalize pagerank table
    std::map<int, float> pagerank;
    int N = edge_weights.GetSize();
    emp::vector<int> not_dangling;
    emp::vector<emp::vector<double> > better_edge_weights = edge_weights.GetWeights();
    for (int i = 0; i < N; i++) {
        if (std::accumulate(better_edge_weights[i].begin(), better_edge_weights[i].end(), 0.0) != 0)
        not_dangling.push_back(i);
        else
        pagerank[i] = 0;
    }
    for (int i : not_dangling) {
        pagerank[i] = 1.0/(not_dangling.size());
    }

    //iterate
    double tol = 1e-6;
    double delta = 100;
    int iterations = 0;
    while ((iterations < 500) && (delta > tol)) {
        std::map<int, float> next_pagerank;
        for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            next_pagerank[j] += alpha * pagerank[i] * edge_weights.GetWeight(i,j);
        }
        }

        delta = 0.0;
        for (int i = 0; i < N; i++) {
            delta += abs(next_pagerank[i] - pagerank[i]);
        }

        pagerank = next_pagerank;
        iterations += 1;
    }

    //convert nodes into bitstrings
    std::map<std::string, float> final_pagerank;
    for (int i = 0; i < N; i++) {
        std::string label = std::bitset<64>(i).to_string();
        final_pagerank[label.erase(0, 64-N_TYPES)] = pagerank[i];
    }

    return final_pagerank;
}