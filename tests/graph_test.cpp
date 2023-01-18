#define CATCH_CONFIG_MAIN

#include "Catch/single_include/catch2/catch.hpp"

#include <iostream>

#include "emp/base/vector.hpp"
#include "emp/datastructs/Graph.hpp"

#include "chemical-ecology/a-eco.hpp"
#include "chemical-ecology/ExampleConfig.hpp"

chemical_ecology::Config cfg;

TEST_CASE("Test Assembly Size")
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
    std::cout << g.GetSize();
    std::cout << "\n";
}

TEST_CASE("Test Assembly and Fitness Equality")
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
        assert(g.GetSize() == g2.GetSize());
    }
    std::cout << g.GetSize();
    std::cout << "\n";
    std::cout << g.GetEdgeCount();
}

TEST_CASE("Dominant Community")
{   
    AEcoWorld world;
    world.Setup(cfg);
    emp::vector<int> community {0, 0, 0, 1, 1, 1, 1, 1, 1};
    emp::vector<int> community2 {1, 1, 1, 0, 0, 0, 0, 0, 0};
    CellData data = world.doGetFitness(community);
    CellData data2 = world.doGetFitness(community2);
    std::cout << data.species << "\n" << data.equilib_growth_rate << "\n" << data.biomass<<"\n"<<data.heredity<<"\n"<<data.invasion_ability<<"\n\n";
    std::cout << data2.species << "\n" << data2.equilib_growth_rate << "\n" << data2.biomass<<"\n"<<data2.heredity<<"\n"<<data2.invasion_ability<<"\n";
    assert(data.equilib_growth_rate > data2.equilib_growth_rate);
}

TEST_CASE("Fittest edge"){
    bool found = false;
    AEcoWorld world;
    world.Setup(cfg);
    emp::Graph g = world.CalculateCommunityLevelFitnessLandscape("Growth_Rate");
    for(size_t i = 0; i < g.GetSize(); i++){
        emp::Graph::Node n = g.GetNode(i);
        if(n.GetLabel() == "111111000"){
            found = true;
            assert(found);
            std::cout << g.GetInDegree(i);
            emp::BitVector out_nodes = n.GetEdgeSet();
            for(int pos = out_nodes.FindOne(); pos >= 0 && pos < g.GetSize(); pos = out_nodes.FindOne(pos+1)) {
                std::cout << g.GetLabel(pos) << "\n";
            }   
        }
    }
    std::cout << "done";
}

TEST_CASE("Assembly edge"){
    bool found = false;
    AEcoWorld world;
    world.Setup(cfg);
    emp::Graph g = world.CalculateCommunityAssemblyGraph();
    for(size_t i = 0; i < g.GetSize(); i++){
        emp::Graph::Node n = g.GetNode(i);
        if(n.GetLabel() == "111111000"){
            found = true;
            assert(found);
            std::cout << "here" << g.GetInDegree(i);
            emp::BitVector out_nodes = n.GetEdgeSet();
            for(int pos = out_nodes.FindOne(); pos >= 0 && pos < g.GetSize(); pos = out_nodes.FindOne(pos+1)) {
                std::cout << g.GetLabel(pos) << "\n";
            }   
        }
    }
    std::cout << "done";
}

TEST_CASE("Test Edges"){
    emp::vector<int> community {0, 0, 0, 1, 1, 1, 1, 1, 1};
    emp::vector<int> community2 {1, 0, 0, 0, 0, 0, 1, 1, 1};
    emp::vector<int> community3 {0, 1, 0, 0, 0, 0, 1, 1, 1};
    emp::vector<int> community4 {0, 0, 1, 0, 0, 0, 1, 1, 1};
    AEcoWorld world;
    world.Setup(cfg);
    CellData data = world.doGetFitness(community);
    CellData data2 = world.doGetFitness(community2);
    CellData data3 = world.doGetFitness(community3);
    CellData data4 = world.doGetFitness(community4);
    std::cout << data.species << "\nGrowth rate: " << data.equilib_growth_rate << "\n";
    std::cout << data2.species << "\nGrowth rate: " << data2.equilib_growth_rate << "\n";
    std::cout << data3.species << "\nGrowth rate: " << data3.equilib_growth_rate << "\n";
    std::cout << data4.species << "\nGrowth rate: " << data4.equilib_growth_rate << "\n";
}