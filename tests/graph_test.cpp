#define CATCH_CONFIG_MAIN

#include "Catch/single_include/catch2/catch.hpp"

#include <iostream>

#include "emp/base/vector.hpp"
#include "emp/datastructs/Graph.hpp"

#include "chemical-ecology/a-eco.hpp"
#include "chemical-ecology/ExampleConfig.hpp"

chemical_ecology::Config cfg;

TEST_CASE("Test Assembly Self Equality")
{   
    AEcoWorld world;
    world.Setup(cfg);
    emp::Graph g = world.CalculateCommunityAssemblyGraph();
    emp::Graph g2 = world.CalculateCommunityAssemblyGraph();
    assert(g.GetSize() == g2.GetSize());
    assert(g.GetEdgeCount() == g2.GetEdgeCount());
    for(size_t i = 0; i < g.GetSize(); i++){
        emp::Graph::Node n = g.GetNode(i);
        emp::Graph::Node n2 = g2.GetNode(i);
        assert(n.GetEdgeSet() == n2.GetEdgeSet());
    }
}

TEST_CASE("Test Assembly and Fitness Node Equality")
{   
    AEcoWorld world;
    world.Setup(cfg);
    emp::vector<std::string> measures;
    measures.push_back("Growth_Rate");
    measures.push_back("Biomass");
    measures.push_back("Heredity");
    measures.push_back("Invasion_Ability");
    measures.push_back("Resiliance");
    emp::Graph g = world.CalculateCommunityAssemblyGraph();
    for(int i = 0; i < 5; i++){
        emp::Graph g2 = world.CalculateCommunityLevelFitnessLandscape(measures[i]);
        //Fitness graphs should have equal number of nodes to the assembly graph
        //Even though we dont print nodes with 0 in/out degree, they should still show up here
        assert(g.GetSize() >= g2.GetSize());
    }
}

TEST_CASE("Resiliance edge test")
{   
    AEcoWorld world;
    world.Setup(cfg);
    emp::Graph g = world.CalculateCommunityLevelFitnessLandscape("Resiliance");
    double curr_out;
    double adj_out;
    for(size_t curr_pos = 0; curr_pos < g.GetSize(); curr_pos++){
    emp::Graph::Node n = g.GetNode(curr_pos);
    curr_out = n.GetDegree();
    emp::BitVector out_nodes = n.GetEdgeSet();
    for(int pos = out_nodes.FindOne(); pos >= 0 && pos < out_nodes.size(); pos = out_nodes.FindOne(pos+1)){
        adj_out = g.GetDegree(pos);
        assert(curr_out > adj_out);
        }
    }
}