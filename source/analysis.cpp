//  This file is part of Artificial Ecology for Chemical Ecology Project
//  Copyright (C) Emily Dolson, 2022.
//  Released under MIT license; see LICENSE

#include <iostream>
#include <map>
#include <fstream>

#include "emp/base/vector.hpp"
#include "emp/datastructs/Graph.hpp"

#include "chemical-ecology/a-eco.hpp"
#include "chemical-ecology/ExampleConfig.hpp"

chemical_ecology::Config cfg;

void printSoupWorld(AEcoWorld world, ofstream& File, chemical_ecology::Config * config)
{
  std::map<std::string, double> soupResults = world.getSoupWorlds(100, 1000, config->PROB_CLEAR(), config->SEEDING_PROB());
  for(auto& [key, val] : soupResults){
    File << key << " " << val << std::endl;
  }
}

void printPageRank(AEcoWorld world, emp::Graph g, ofstream& File, emp::WeightedGraph wAssembly)
{
  std::map<std::string, float> map = world.CalculateWeightedPageRank(wAssembly);
  std::map<std::string, float>::iterator it = map.begin();
  while (it != map.end())
  {
      File << it->first << " " << it->second << std::endl;
      ++it;
  }
}

void printGraph(emp::Graph g, ofstream& File, emp::WeightedGraph wAssembly, chemical_ecology::Config * config)
{
  for (size_t from = 0; from < wAssembly.GetSize(); from++) {
    for (size_t to = 0; to < wAssembly.GetSize(); to++) {
      if (wAssembly.HasEdge(from, to) == false) continue;
      File << std::bitset<64>(from).to_string().erase(0, 64-config->N_TYPES()) << " " << std::bitset<64>(to).to_string().erase(0, 64-config->N_TYPES()) << " " << wAssembly.GetWeight(from, to) << std::endl;
    }
  }
}

int main(int argc, char* argv[])
{ 
  // Set up a configuration panel for native application
  setup_config_native(cfg, argc, argv);
  cfg.Write(std::cout);
  AEcoWorld world;
  world.Setup(cfg);
  chemical_ecology::Config * config = &cfg;

  // Initiate files to write to
  ofstream GraphFile("community_fitness.txt");
  ofstream PageRankFile("page_rank.txt");
  ofstream FinalCommunitiesFile("final_communities.txt");
  ofstream SoupWorldFile("soup_world.txt");

  // Calculate community assembly and fitness graphs
  std::map<std::string, emp::Graph> graphs;
  emp::Graph g1 = world.CalculateCommunityAssemblyGraph();
  emp::WeightedGraph wAssembly = world.calculateWeightedAssembly(g1, config->PROB_CLEAR(), config->SEEDING_PROB());
  graphs["Community Assembly"] = g1;
  std::vector<std::string> fitnessTypes{"Biomass", "Growth_Rate", "Heredity", "Invasion_Ability", "Resiliance"};
  /*for(int i = 0; i < 5; i++)
  {
    emp::Graph g2 = world.CalculateCommunityLevelFitnessLandscape(fitnessTypes[i]);
    graphs[fitnessTypes[i]] = g2;
  }*/

  // Write graphs and the node pageranks
  for(auto& [graphName, graph] : graphs)
  {
    GraphFile << "***" << graphName << "***" << std::endl;
    printGraph(graph, GraphFile, wAssembly, config);
    PageRankFile << "***" << graphName << "***" << std::endl;
    printPageRank(world, graph, PageRankFile, wAssembly);
  }

  // Write soup world
  printSoupWorld(world, SoupWorldFile, config);

  // Write final communities
  world.Run();
  emp::vector<emp::vector<double>> stable_world = world.stableUpdate(1000);
  std::map<std::string, double> finalCommunities = world.getFinalCommunities(stable_world);
  for(auto& [node, proportion] : finalCommunities)
  {
    FinalCommunitiesFile << node << " " << proportion << std::endl;
  }

  // End
  GraphFile.close();
  PageRankFile.close();
  FinalCommunitiesFile.close();
  return 0;
}