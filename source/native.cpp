//  This file is part of Artificial Ecology for Chemical Ecology Project
//  Copyright (C) Emily Dolson, 2021.
//  Released under MIT license; see LICENSE

#include <iostream>

#include "emp/base/vector.hpp"

#include "chemical-ecology/AEcoWorld.hpp"
#include "chemical-ecology/Config.hpp"
#include "chemical-ecology/utils/config_setup.hpp"

// This is the main function for the NATIVE version of Artificial Ecology for Chemical Ecology Project.

chemical_ecology::Config cfg;

int main(int argc, char* argv[])
{
  // Set up a configuration panel for native application
  chemical_ecology::utils::setup_config_native(cfg, argc, argv);
  if(cfg.V()){
    cfg.Write(std::cout);
  }

  chemical_ecology::AEcoWorld world;
  world.Setup(cfg);
  world.Run();

  return 0;
}
