//  This file is part of Artificial Ecology for Chemical Ecology Project
//  Copyright (C) Emily Dolson, 2022.
//  Released under MIT license; see LICENSE

#include <iostream>

#include "emp/base/vector.hpp"
#include "emp/datastructs/Graph.hpp"

#include "chemical-ecology/a-eco.hpp"
#include "chemical-ecology/ExampleConfig.hpp"

#include "pagerank.h"

chemical_ecology::Config cfg;

void printPageRank(emp::Graph g)
{
    Table t;

    t.set_trace(false);
    t.set_numeric(false);
    t.set_delim(" ");
    t.read_graph(g);
    t.pagerank();

    std::map<std::string, float> map = t.get_pr_map();
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

  std::cout << "***Community Assembly***" << std::endl;

  emp::Graph g = world.CalculateCommunityAssemblyGraph();
  //printPageRank(g);
  printGraph(g);

  std::vector<std::string> fitnessTypes{"Biomass", "Growth_Rate", "Heredity", "Invasion_Ability", "Resiliance"};
  for(int i = 0; i < 5; i++)
  {
    std::cout << "***" << fitnessTypes[i] << "***" << std::endl;
    emp::Graph g2 = world.CalculateCommunityLevelFitnessLandscape(fitnessTypes[i]);
    printGraph(g2);
  }

  world.WriteInteractionMatrix("interaction.dat");

  return 0;
}