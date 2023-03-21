#pragma once

#include <iostream>
#include "emp/Evolve/World.hpp"
#include "emp/math/distances.hpp"
#include "emp/datastructs/Graph.hpp"
#include "emp/bits/BitArray.hpp"
#include "emp/data/DataNode.hpp"
#include "chemical-ecology/config_setup.hpp"
#include "emp/datastructs/map_utils.hpp"
#include <string>
#include <math.h>
#include "pagerank.h"

// TERMINOLOGY NOTES:
//
// - Cell = location in the world grid
// - Type = a type of thing that can live in a cell
//          (e.g. a chemical or a species)
// - ??? = a set of types that commonly co-occur
//         (don't currently have word for this but need one)
//          I'm going to use the word "community" for now,
//          but I don't love that one

// A struct to hold data about individual cells 
struct CellData {
  double invasion_ability;  // Ability of this cell to move across an empty world
  double heredity; // Extent to which offspring cells look same
  double growth_rate;  // How fast this cell produces biomass
  double synergy;  // Attempt at measuring facilitation in cell
  double equilib_growth_rate;  // How fast cell produces biomass
                              // after growing to equilibrium
  double biomass;
  int count = 1;  // For counting the number of cells occupied by the same community
  emp::vector<double> species; // The counts of each type in this cell
};

// A class to handle running the simple ecology
class AEcoWorld {
  private:
    // The world is a vector of cells (where each cell is
    // represented as a vector ints representing the count of
    // each type in each cell). 
    
    // Although the world is stored as a flat vector of cells
    // it represents a grid of cells
    using world_t = emp::vector<emp::vector<double> >;

    // The matrix of interactions between types
    // The diagonal of this matrix represents the 
    // intrinsic growth rate (r) of each type
    emp::vector<emp::vector<double> > interactions;

    // A random number generator for all our random number
    // generating needs
    emp::Random rnd;

    // All configuration information is stored in config
    chemical_ecology::Config * config = nullptr;

    // These values are set in the config but we have local
    // copies for efficiency in accessing their values
    int N_TYPES;
    double MAX_POP;
    int curr_update;

    // Set up data tracking
    emp::DataFile data_file;
    emp::DataFile score_file;
    emp::DataNode<double, emp::data::Stats> growth_rate_node;
    emp::DataNode<double, emp::data::Stats> synergy_node;
    emp::DataNode<double, emp::data::Stats> heredity_node;
    emp::DataNode<double, emp::data::Stats> invasion_node;
    emp::DataNode<double, emp::data::Stats> biomass_node;
    std::map<emp::vector<double>, CellData> dom_map;
    emp::vector<double> fittest;
    emp::vector<double> dominant;

  public:
  // Default constructor just has to set name of output file
  AEcoWorld() : data_file("a-eco_data.csv"), score_file("scores.csv") {;}

  // Initialize vector that keeps track of grid
  world_t world;

  // Getter for interaction matrix
  emp::vector<emp::vector<double> > GetInteractions() {
    return interactions;
  }

  // Set interaction value between types x and y to w
  void SetInteraction(int x, int y, double w) {
    interactions[x][y] = w;
  }

  // Set the interaction matrix to contain random values
  // (with probabilities determined by configs)
  void SetupRandomInteractions() {
    // interaction matrix will ultimately have size
    // N_TYPES x N_TYPES
    interactions.resize(N_TYPES);

    for (int i = 0; i < N_TYPES; i++) {
      interactions[i].resize(N_TYPES);
      for (int j = 0; j < N_TYPES; j++) {
        // PROB_INTERACTION determines probability of there being
        // an interaction between a given pair of types.
        // Controls sparsity of interaction matrix
        if (rnd.P(config->PROB_INTERACTION())){
          // If there's an interaction, it is a random double
          // between -INTERACTION_MAGNITUDE and INTERACTION_MAGNITUDE
          interactions[i][j] = rnd.GetDouble(config->INTERACTION_MAGNITUDE() * -1, config->INTERACTION_MAGNITUDE());
        }
      }
    }
  }

  // Load an interaction matrix from the specified file
  void LoadInteractionMatrix(std::string filename) {    
    emp::File infile(filename);
    emp::vector<emp::vector<double>> interaction_data = infile.ToData<double>();

    interactions.resize(N_TYPES);
    for (int i = 0; i < N_TYPES; i++) {
      interactions[i].resize(N_TYPES);
      for (int j = 0; j < N_TYPES; j++) {
        interactions[i][j] = interaction_data[i][j];
      }
    }
  }

  // Store the current interaction matrix in a file
  void WriteInteractionMatrix(std::string filename) {
    emp::File outfile;
    
    for (int i = 0; i < N_TYPES; i++) {
      outfile += emp::join(interactions[i], ",");
    }
    
    outfile.Write(filename);
  }

  // Setup the world according to the specified configuration
  void Setup(chemical_ecology::Config & cfg) {
    // Store cfg for future reference
    config = &cfg;

    // Set local config variables based on given configuration
    N_TYPES = config->N_TYPES();
    MAX_POP = double(config->MAX_POP());

    // Set seed to configured value for reproducibility
    rnd.ResetSeed(config->SEED());

    // world vector needs a spot for each cell in the grid
    world.resize(config->WORLD_X() * config->WORLD_Y());
    // Initialize world vector
    for (emp::vector<double> & v : world) {
      v.resize(N_TYPES);
      for (double & count : v) {
        // The quantity of each type in each cell is either 0 or 1
        // The probability of it being 1 is controlled by SEEDING_PROB
        count = rnd.P(config->SEEDING_PROB());
      }
    }

    // Setup interaction matrix based on the method
    // specified in the configuration file
    if (config->INTERACTION_SOURCE() == "") {
      SetupRandomInteractions();
    } else {
      LoadInteractionMatrix(config->INTERACTION_SOURCE());
    }

    // Add columns to data file
    // Time column will print contents of curr_update
    data_file.AddVar(curr_update, "Time", "Time");
    // Add columns for stats (mean, variance, etc.) from the 
    // three data nodes. Clear data from nodes.
    data_file.AddStats(growth_rate_node, "Growth_Rate", "Growth_Rate", true);
    data_file.AddStats(synergy_node, "Synergy", "Synergy", true);
    data_file.AddStats(heredity_node, "Heredity", "Heredity", true);
    data_file.AddStats(invasion_node, "Invasion_Ability", "Invasion_Ability", true);
    data_file.AddStats(biomass_node, "Biomass", "Biomass", true);
    // Add columns calculated by running functions that look up values
    // in dom_map
    data_file.AddFun((std::function<int()>)[this](){return dom_map[fittest].count;}, "fittest_count", "fittest_count");
    data_file.AddFun((std::function<int()>)[this](){return dom_map[fittest].heredity;}, "fittest_heredity", "fittest_heredity");
    data_file.AddFun((std::function<std::string()>)[this](){return emp::to_string(fittest);}, "fittestCommunity", "fittest community");
    data_file.AddFun((std::function<int()>)[this](){return dom_map[dominant].count;}, "dominant_count", "dominant_count");
    data_file.AddFun((std::function<int()>)[this](){return dom_map[dominant].heredity;}, "dominant_heredity", "dominant_heredity");
    data_file.AddFun((std::function<std::string()>)[this](){return emp::to_string(dominant);}, "dominantCommunity", "dominant community");
    // Add function that runs before each row is printed in the datafile
    // to ensure fitnesses have been calculated (it's expensive so we don't
    // want to do it every update if we aren't going to use that data)
    data_file.AddPreFun([this](){CalcAllFitness();});
    // Set data_file to print a new row every 10 time steps
    data_file.SetTimingRepeat(10);
    // Print column names to data_file
    data_file.PrintHeaderKeys();
  }

  // Handle an individual time step
  // ud = which time step we're on
  void Update(int ud) {

    // Create a new world object to store the values 
    // for the next time step
    world_t next_world;
    next_world.resize(config->WORLD_X() * config->WORLD_Y());
    for (emp::vector<double> & v : next_world) {
      v.resize(N_TYPES, 0);
    }

    // Handle population growth for each cell
    for (size_t pos = 0; pos < world.size(); pos++) {
      DoGrowth(pos, world, next_world);
    }

    // Handle everything that allows biomass to move from
    // one cell to another. Do so for each cell
    for (int pos = 0; pos < (int)world.size(); pos++) {

      // Figure out which cells are above, below, left
      // and right of the focal cell.
      // This process is a little arduous because it
      // needs to handle toroidal wraparound
      int x = pos % config->WORLD_X(); // x coordinate
      int y = pos / config->WORLD_Y(); // y coordinate
      int left = pos - 1;  // Assuming no wraparound    
      int right = pos + 1; // Assuming no wraparound

      // Check whether we need to adjust left and right
      // cell ids for wraparound
      if (x == 0) {
        left += config->WORLD_X(); 
      } else if (x == config->WORLD_X() - 1) {
        right -= config->WORLD_X();
      }
 
      // Calculate up and down assuming no wraparound
      int up = pos - config->WORLD_X();
      int down = pos + config->WORLD_X();

      // Adjust for wraparound if necessary
      if (y == 0) {
        up += config->WORLD_X() * (config->WORLD_Y());
      } else if (y == config->WORLD_Y() - 1) {
        down -= config->WORLD_X() * (config->WORLD_Y());
      }

      // Actually call function that handles between-cell
      // movement
      emp::vector<int> adj = {up, down, left, right};
      DoRepro(pos, adj, world, next_world, config->SEEDING_PROB(), config->PROB_CLEAR());
    }

    // Update our time tracker
    curr_update = ud;
    
    // Give data_file the opportunity to write to the file
    data_file.Update(curr_update);

    // We're done calculating the type counts for the next
    // time step. We can now swap our counts for the next
    // time step into the main world variable 
    std::swap(world, next_world);

    // Clean-up data trackers
    dom_map.clear();
    fittest.clear();
    dominant.clear();
  }


  //This function should be called to create a stable copy of the world for graph calculations
  world_t stableUpdate(int num_updates){
    world_t stable_world;
    world_t next_stable_world;

    stable_world.resize(config->WORLD_X() * config->WORLD_Y());
    for (emp::vector<double> & v : next_stable_world) {
      v.resize(N_TYPES, 0);
    }

    next_stable_world.resize(config->WORLD_X() * config->WORLD_Y());
    for (emp::vector<double> & v : next_stable_world) {
      v.resize(N_TYPES, 0);
    }

    //copy current world into stable world
    for(size_t i = 0; i < world.size(); i++){
      stable_world[i].assign(world[i].begin(), world[i].end());
    }

    for(int i = 0; i < num_updates; i++){
      // Handle population growth for each cell
      for (size_t pos = 0; pos < world.size(); pos++) {
        DoGrowth(pos, stable_world, next_stable_world);
      }

      std::swap(stable_world, next_stable_world);
      dom_map.clear();
      fittest.clear();
      dominant.clear();
    }
    return stable_world;
  }

  // Handle the process of running the program through
  // all time steps
  void Run() {
    // Call update the specified number of times
    for (int i = 0; i < config->UPDATES(); i++) {
      Update(i);
    }
    
    //temp fix for very small communities
    world_t stable_world_temp = stableUpdate(1000);
    world_t stable_world;
    for(size_t i = 0; i < stable_world_temp.size(); i++){
      float sum_of_elems = 0;
      for (auto& n : world[i])
        sum_of_elems += n;
      if (sum_of_elems > 500)
        stable_world.push_back(stable_world_temp[i]);
    }
    
    std::map<std::string, double> finalCommunities = getFinalCommunities(stable_world);

    emp::Graph assemblyGraph = CalculateCommunityAssemblyGraph();
    std::map<std::string, float> assembly_pr_map = calculatePageRank(assemblyGraph);
    double assembly_score = 0;
    for(auto& [key, val] : finalCommunities){
      std::string node = key;
      float proportion = val;
      if (assembly_pr_map.find(node) != assembly_pr_map.end()) {
        assembly_score += (proportion * assembly_pr_map[node]);
      }
      else {
        assembly_score *= pow(10.0, -10.0);
      }
    }

    double biomass_score = calcAdaptabilityScore(finalCommunities, "Biomass");
    double growth_rate_score = calcAdaptabilityScore(finalCommunities, "Growth_Rate");
    double heredity_score = calcAdaptabilityScore(finalCommunities, "Heredity");
    double invasion_score = calcAdaptabilityScore(finalCommunities, "Invasion_Ability");
    double resiliance_score = calcAdaptabilityScore(finalCommunities, "Resiliance");

    score_file.AddVar(biomass_score, "Biomass_Score", "Biomass_Score");
    score_file.AddVar(growth_rate_score, "Growth_Rate_Score", "Growth_Rate_Score");
    score_file.AddVar(heredity_score, "Heredity_Score", "Heredity_Score");
    score_file.AddVar(invasion_score, "Invasion_Ability_Score", "Invasion_Ability_Score");
    score_file.AddVar(resiliance_score, "Resiliance_Score", "Resiliance_Score");
    score_file.AddVar(assembly_score, "Assembly_Score", "Assembly_Score");

    score_file.PrintHeaderKeys();

    score_file.Update();
    
    // Print out final state
    /*std::cout << "World Vectors:" << std::endl;
    for (auto & v : world) {
      std::cout << emp::to_string(v) << std::endl;
    }
    std::cout << "Stable World Vectors:" << std::endl;
    for (auto & v : stable_world) {
      std::cout << emp::to_string(v) << std::endl;
    }*/

    // Store interaction matrix in a file in case we
    // want to do stuff with it later
    //WriteInteractionMatrix("interaction_matrix.dat");
  }

  double doCalcGrowthRate(emp::vector<double> community){
    double growth_rate = 0;
    for (int i = 0; i < N_TYPES; i++) {
      double modifier = 0;
      for (int j = 0; j < N_TYPES; j++) {
        // Each type contributes to the growth rate modifier
        // of each other type based on the product of its
        // interaction value with that type and its population size
        modifier += interactions[i][j] * community[j];
      }

      // Add this type's overall growth rate to the cell-level
      // growth-rate
      growth_rate += (modifier*(community[i]/MAX_POP)); // * ((double)(MAX_POP - pos[i])/MAX_POP));
    }
    return growth_rate;
  }

  // Calculate the growth rate (one measurement of fitness)
  // for a given cell
  double CalcGrowthRate(size_t pos, world_t & curr_world) {
    return doCalcGrowthRate(curr_world[pos]);
  }

  // Calculate the amount of cell-level growth that is due to 
  // facilitation (i.e. beneficial interactions)
  double CalcSynergy(size_t pos, world_t & curr_world) {
    double growth_rate = 0;
    for (int i = 0; i < N_TYPES; i++) {
      double modifier = 0;
      for (int j = 0; j < N_TYPES; j++) {
        if (i != j) {
          // Calculate is same as growth rate except we exclude
          // the type's interaction with itself (i.e. intrinsic
          // growth rate)
          modifier += interactions[i][j] * curr_world[pos][j];
        }
      }

      growth_rate += modifier;
      // growth_rate += ceil(modifier*((double)curr_world[pos][i]/(double)MAX_POP)); // * ((double)(MAX_POP - pos[i])/MAX_POP));
    }
    return growth_rate;
  }

  // Handles population growth of each type within a cell
  void DoGrowth(size_t pos, world_t & curr_world, world_t & next_world) {
    for (int i = 0; i < N_TYPES; i++) {
      double modifier = 0;
      for (int j = 0; j < N_TYPES; j++) {
        // Sum up growth rate modifier for current type 
        modifier += interactions[i][j] * curr_world[pos][j];
      }

      // Grow linearly until we hit the cap
      double new_pop = modifier*curr_world[pos][i]; // * ((double)(MAX_POP - pos[i])/MAX_POP));
      // Population size cannot be negative
      next_world[pos][i] = std::max(curr_world[pos][i] + new_pop, 0.0);
      // Population size capped at MAX_POP
      next_world[pos][i] = std::min(next_world[pos][i], MAX_POP);
    }
  }

  // Allows the assembly graph to handle growth
  // void GraphDoGrowth(size_t pos, world_t & curr_world, world_t & next_world) {
  //   for (int i = 0; i < N_TYPES; i++) {
  //     double modifier = 0;
  //     for (int j = 0; j < N_TYPES; j++) {
  //       // Sum up growth rate modifier for current type 
  //       modifier += interactions[i][j] * curr_world[pos][j];
  //     }

  //     // Grow linearly until we hit the cap
  //     int new_pop = floor(modifier*curr_world[pos][i]); // * ((double)(MAX_POP - pos[i])/MAX_POP));
  //     // Population size cannot be negative
  //     next_world[pos][i] = std::max(curr_world[pos][i] + new_pop, 0.0);
  //     // Population size capped at MAX_POP
  //     next_world[pos][i] = std::min(next_world[pos][i], MAX_POP);
  //   }
  // }

  // Check whether the total population of a cell is large enough
  // for group level replication
  bool IsReproReady(size_t pos, world_t & w) {
    return emp::Sum(w[pos]) > config->MAX_POP()*(config->REPRO_THRESHOLD());
  }

  // Handle movement of biomass between cells
  void DoRepro(size_t pos, emp::vector<int> & adj, world_t & curr_world, world_t & next_world, double seed_prob, double prob_clear) {

    // Check whether conditions for group-level replication are met
    if (IsReproReady(pos, curr_world)) {
      // Choose a random cell out of those that are adjacent
      // to replicate into
      int direction = adj[rnd.GetInt(adj.size())];
      for (int i = 0; i < N_TYPES; i++) {
        // Add a portion (configured by REPRO_DILUTION) of the quantity of the type
        // in the focal cell to the cell we're replicating into
        //
        // NOTE: An important decision here is whether to clear the cell first.
        // We have chosen not to, but can revisit that choice
        next_world[direction][i] += curr_world[pos][i] * config->REPRO_DILUTION();
      }
    }

    // Each cell has a chance of being cleared on every time step
    if (rnd.P(prob_clear)) {
      for (int i = 0; i < N_TYPES; i++) {
        next_world[pos][i] = 0;
      }
    }

    // Handle diffusion
    for (int direction : adj) {
      for (int i = 0; i < N_TYPES; i++) {
        // Calculate amount diffusing
        double avail = curr_world[pos][i] * config->DIFFUSION();
        
        // Send 1/4 of it in the direction we're currently processing
        next_world[direction][i] +=  avail / 4;
        // if (avail/4 > 0 && adj.size() == 4) {
        //   std::cout << pos <<  " " << i << " Diffusing: " << next_world[direction][i] << " " << avail/4 << std::endl;
        // }

        // make sure final value is legal number
        next_world[direction][i] = std::min(next_world[direction][i], MAX_POP);
        next_world[direction][i] = std::max(next_world[direction][i], 0.0);
      }
    }


    //subtract diffusion
    for (int i = 0; i < N_TYPES; i++) {
      next_world[pos][i] -= curr_world[pos][i] * config->DIFFUSION();
      // We can't have negative population sizes
      next_world[pos][i] = std::max(next_world[pos][i], 0.0);
    }
    if (rnd.P(seed_prob)) {
      int cell = rnd.GetInt(N_TYPES);
      next_world[pos][cell]++;
    }

  }

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
    // Calculate synergy
    data.synergy = CalcSynergy(0, test_world);
    // Calculate equilibrium growth rate
    data.equilib_growth_rate = doCalcGrowthRate(starting_point);
    // Record community composition
    data.species = starting_point;

    // Add equilibrium growth rate to growth_rate data node
    growth_rate_node.Add(data.equilib_growth_rate);
    // Add synergy to synergy data node
    synergy_node.Add(data.synergy);
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

    while (!IsReproReady(num_cells - 1, test_world) && time < time_limit) {

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
        DoRepro(pos, adj, test_world, next_world, 0, 0);
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

  std::map<std::string, float> calculatePageRank(emp::Graph g) {
    Table t;

    t.set_trace(false);
    t.set_numeric(false);
    t.set_delim(" ");
    t.read_graph(g);
    t.pagerank();

    std::map<std::string, float> map = t.get_pr_map();
    return map;
  }
  

std::map<std::string, double> getFinalCommunities(world_t stable_world) {
    double size = stable_world.size();
    std::map<std::string, double> finalCommunities;
    for(emp::vector<double> cell: stable_world){
      std::string temp = "";
      for(double species: cell){
        if(species > 0.001){
          temp.append("1");
        }
        else{
          temp.append("0");
        }
      }
      std::reverse(temp.begin(), temp.end());
      if(finalCommunities.find(temp) == finalCommunities.end()){
          finalCommunities.insert({temp, 1});
      }
      else{
        finalCommunities[temp] += 1;
      }
    }
    for(auto& [key, val] : finalCommunities){
      val = val/size;
    }
    return finalCommunities;
  }

  // Getter for current update/time step
  int GetTime() {
    return curr_update;
  }
};