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

// A struct to hold data about individual cells 
struct CellData {
  double invasion_ability;  // Ability of this cell to move across an empty world
  double heredity; // Extent to which offspring cells look same
  double growth_rate;  // How fast this cell produces biomass
  double equilib_growth_rate;  // How fast cell produces biomass
                              // after growing to equilibrium
  double biomass;
  int count = 1;  // For counting the number of cells occupied by the same community
  emp::vector<double> species; // The counts of each type in this cell
};
// Old datafile code to track fitness proxies
// std::map<emp::vector<double>, CellData> dom_map;
/*
// emp::DataNode<double, emp::data::Stats> growth_rate_node;
// emp::DataNode<double, emp::data::Stats> heredity_node;
// emp::DataNode<double, emp::data::Stats> invasion_node;
//
// emp::DataFile data_file;
// Add columns to data file
// Time column will print contents of curr_update
data_file.AddVar(curr_update, "Time", "Time");
// Add columns for stats (mean, variance, etc.) from the 
// data nodes. Clear data from nodes.
data_file.AddStats(biomass_node, "Biomass", "Biomass", true);
// Add columns calculated by running functions that look up values
// in dom_map
data_file.AddFun((std::function<std::string()>)[this](){return emp::to_string(fittest);}, "fittestCommunity", "fittest community");
// Add function that runs before each row is printed in the datafile
// to ensure fitnesses have been calculated (it's expensive so we don't
// want to do it every update if we aren't going to use that data)
data_file.AddPreFun([this](){CalcAllFitness();});
// Set data_file to print a new row every 10 time steps
data_file.SetTimingRepeat(10);
// Print column names to data_file
data_file.PrintHeaderKeys();
*/

// Calculate the fitnesses of all cells
// (useful for preparing to print data)
void CalcAllFitness() {
  for (size_t pos = 0; pos < world.size(); pos++) {
    // Calculate fitness
    CellData res = GetFitness(pos);
    // Keep track of counts of each community
    if (emp::Has(dom_map, res.species)) {
      dom_map[res.species].count++;
    } else {
      dom_map[res.species] = res;
    }
    // Keep track of fittest cell we've seen so far
    if (fittest.size() == 0 || dom_map[res.species].equilib_growth_rate > dom_map[fittest].equilib_growth_rate) {
      fittest = res.species;
    }
  }

  // Track most populus community
  dominant = (*std::max_element(dom_map.begin(), dom_map.end(), 
              [](const std::pair<emp::vector<double>, CellData> & v1, const std::pair<emp::vector<double>, CellData> & v2 ){return v1.second.count < v2.second.count;})).first;
}

//This function allows us to calculate the fitness with just take a community composition.
CellData doGetFitness(emp::vector<double> community) {
  CellData data;
  // Size of test world 
  // Test world is a num_cells x 1 grid
  // (i.e. a strip of cells)
  int num_cells = 5;
  // Number of time steps to run test world for
  int time_limit = 250;
  world_t test_world;
  test_world.resize(num_cells);
  // Initialize empty test world
  for (emp::vector<double> & v : test_world) {
    v.resize(config->N_TYPES(), 0);
  }

  // Put contents of cell being evaluated
  // in left-most spot in test world
  test_world[0] = community;

  // Make a next_world for the test world to
  // handle simulated updates in our test environment
  world_t next_world;
  next_world.resize(num_cells);
  for (emp::vector<double> & v : next_world) {
    v.resize(config->N_TYPES(), 0);
  }

  // Run DoGrowth for 20 time steps to hopefully
  // get cell to equilibrium
  // (if cell is an organism, this is analogous to 
  // development)
  for (int i = 0; i < 20; i++) {
    DoGrowth(0, test_world, next_world);
    std::swap(test_world, next_world);
    std::fill(next_world[0].begin(), next_world[0].end(), 0);
  }

  // Record equilibrium community
  emp::vector<double> starting_point = test_world[0];
  // Calculate pure-math growth rate
  data.growth_rate = doCalcGrowthRate(community);
  // Calculate equilibrium growth rate
  data.equilib_growth_rate = doCalcGrowthRate(starting_point);
  // Record community composition
  data.species = starting_point;

  // Add equilibrium growth rate to growth_rate data node
  growth_rate_node.Add(data.equilib_growth_rate);
  //Add biomass to biomass data node
  data.biomass = accumulate(data.species.begin(), data.species.end(), 0.0);
  biomass_node.Add(accumulate(data.species.begin(), data.species.end(), 0.0));

  // Now we're going to simulate up to time_limit number
  // of time steps to measure heredity and 
  // propogation-based fitness
  int time = 1;
  // Stop if hit the time limit or if the right-most cell is
  // able to do group-level reproduction (meaning the community
  // has spread across the entire test world and could keep going)

  while (time < time_limit) {

    // Group level repro will never occur with the values we're using
    // Instead check if the final cell is similar enough to the beginning cell
    // using heredity
    double curr_heredity = 1.0 - emp::EuclideanDistance(test_world[num_cells - 1], starting_point) / ((double)MAX_POP*N_TYPES);
    if(curr_heredity > .999){
      break;
    }

    // Handle within cell growth for each cell in test world
    for (int pos = 0; pos < num_cells; pos++) {
      DoGrowth(pos, test_world, next_world);
    }

    // Handle movement of biomass between cells
    for (int pos = 0; pos < num_cells; pos++) {
      // Because world is a non-toroidal rectangle
      // calculating adjacent positions is easier
      emp::vector<int> adj;
      if (pos > 0) {
        adj.push_back(pos - 1);
      }
      if (pos < num_cells - 1) {
        adj.push_back(pos + 1);
      }

      // DoRepro works as normal except that we don't
      // allow seeding or clearing because the goal is
      // to evaluate this cell's community in a clean
      // environment
      DoRepro(pos, adj, test_world, next_world, 0, 0, config->DIFFUSION(), false);
    }
    // We've finished calculating next_world so we can swap it
    // into test_world
    std::swap(test_world, next_world);
    std::fill(next_world[0].begin(), next_world[0].end(), 0);
    // Track time
    time++;
  }
  // Set invasion ability equal to number of time steps taken to traverse world, inverted, so that higher is better
  data.invasion_ability = 1/(double)time;
  invasion_node.Add(data.invasion_ability);
  // Heredity is how similar the community in the final cell is to 
  // the equilibrium community
  data.heredity = 1.0 - emp::EuclideanDistance(test_world[num_cells - 1], starting_point) / ((double)MAX_POP*N_TYPES);

  // Keep track of the heredity values we've seen
  heredity_node.Add(data.heredity);

  return data;
}


//Get fitness of one cell
CellData GetFitness(size_t orig_pos){
  return doGetFitness(world[orig_pos]);
}