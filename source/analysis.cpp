//  This file is part of Artificial Ecology for Chemical Ecology Project
//  Copyright (C) Emily Dolson, 2022.
//  Released under MIT license; see LICENSE

#include <iostream>

#include "emp/base/vector.hpp"
#include "emp/datastructs/Graph.hpp"

#include "chemical-ecology/a-eco.hpp"
#include "chemical-ecology/ExampleConfig.hpp"

chemical_ecology::Config cfg;

void printAdaptiveness(AEcoWorld world, std::set<std::string> finalCommunities, emp::Graph g)
{
  std::map<std::string, float> map = world.calculatePageRank(g);
  std::map<std::string, float>::iterator it = map.begin();
  while (it != map.end())
  {
    if (finalCommunities.find(it->first) != finalCommunities.end())
      std::cout << it->first << " " << it->second << std::endl;
    else
      std::cout << it->first << " " << "0" << std::endl;
    ++it;
  }
}

void printPageRank(AEcoWorld world, emp::Graph g)
{
  std::map<std::string, float> map = world.calculatePageRank(g);
  std::map<std::string, float>::iterator it = map.begin();
  while (it != map.end())
  {
      std::cout << it->first << " " << it->second << std::endl;
      ++it;
  }
}

void printGraph(emp::Graph g)
{
  emp::vector<emp::Graph::Node> all_nodes = g.GetNodes();
  for(emp::Graph::Node n: all_nodes)
  {
    emp::BitVector out_nodes = n.GetEdgeSet();
    for (int pos = out_nodes.FindOne(); pos >= 0 && pos < g.GetSize(); pos = out_nodes.FindOne(pos+1))
    {
      std::cout << n.GetLabel() << " " << g.GetLabel(pos) << std::endl;
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
  //world.Run();

  //emp::vector<emp::vector<int>> stable_world = world.stableUpdate(20);
  //std::set<std::string> finalCommunities = world.getFinalCommunities(stable_world);

  std::cout << "***Community Assembly***" << std::endl;
  emp::Graph g = world.CalculateCommunityAssemblyGraph();
  //printAdaptiveness(world, finalCommunities, g);
  printGraph(g);

  std::vector<std::string> fitnessTypes{"Biomass", "Growth_Rate", "Heredity", "Invasion_Ability", "Resiliance"};
  for(int i = 0; i < 5; i++)
  {
    std::cout << "***" << fitnessTypes[i] << "***" << std::endl;
    emp::Graph g2 = world.CalculateCommunityLevelFitnessLandscape(fitnessTypes[i]);
    printGraph(g2);
    //printAdaptiveness(world, finalCommunities, g2);
  }

  world.WriteInteractionMatrix("interaction.dat");

  return 0;
}