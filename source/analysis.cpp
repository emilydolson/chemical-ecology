//  This file is part of Artificial Ecology for Chemical Ecology Project
//  Copyright (C) Emily Dolson, 2022.
//  Released under MIT license; see LICENSE

#include <iostream>

#include "emp/base/vector.hpp"

#include "chemical-ecology/a-eco.hpp"
#include "chemical-ecology/ExampleConfig.hpp"

chemical_ecology::Config cfg;

int main(int argc, char* argv[])
{ 
  // Set up a configuration panel for native application
  setup_config_native(cfg, argc, argv);
  cfg.Write(std::cout);

  AEcoWorld world;
  world.Setup(cfg);
  emp::Graph g = world.CalculateCommunityAssemblyGraph();
  g.PrintDirected();
  world.WriteInteractionMatrix("interaction.dat");

  return 0;
}
