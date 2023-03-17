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

void printPageRank(AEcoWorld world, emp::Graph g, ofstream& File)
{
  std::map<std::string, float> map = world.calculatePageRank(g);
  std::map<std::string, float>::iterator it = map.begin();
  while (it != map.end())
  {
      File << it->first << " " << it->second << std::endl;
      ++it;
  }
}

void printGraph(emp::Graph g, ofstream& File)
{
  emp::vector<emp::Graph::Node> all_nodes = g.GetNodes();
  for(emp::Graph::Node n: all_nodes)
  {
    emp::BitVector out_nodes = n.GetEdgeSet();
    for (int pos = out_nodes.FindOne(); pos >= 0 && pos < g.GetSize(); pos = out_nodes.FindOne(pos+1))
    {
      File << n.GetLabel() << " " << g.GetLabel(pos) << std::endl;
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

  // Initiate files to write to
  ofstream GraphFile("community_fitness.txt");
  ofstream PageRankFile("page_rank.txt");
  ofstream FinalCommunitiesFile("final_communities.txt");

  // Calculate community assembly and fitness graphs
  std::map<std::string, emp::Graph> graphs;
  emp::Graph g1 = world.CalculateCommunityAssemblyGraph();
  graphs["Community Assembly"] = g1;
  std::vector<std::string> fitnessTypes{"Biomass", "Growth_Rate", "Heredity", "Invasion_Ability", "Resiliance"};
  for(int i = 0; i < 5; i++)
  {
    emp::Graph g2 = world.CalculateCommunityLevelFitnessLandscape(fitnessTypes[i]);
    graphs[fitnessTypes[i]] = g2;
  }

  // Write graphs and the node pageranks
  for(auto& [graphName, graph] : graphs)
  {
    GraphFile << "***" << graphName << "***" << std::endl;
    printGraph(graph, GraphFile);
    PageRankFile << "***" << graphName << "***" << std::endl;
    printPageRank(world, graph, PageRankFile);
  }

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