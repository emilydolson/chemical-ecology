#pragma once

#include <iostream>
#include "emp/Evolve/World.hpp"
#include "emp/math/distances.hpp"
#include "emp/data/DataNode.hpp"
#include "chemical-ecology/config_setup.hpp"

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
  double fitness;  // Fitness of this cell
  double heredity; // Extent to which offspring cells look same
  double growth_rate;  // How fast this cell produces biomass
  double synergy;  // Attempt at measuring facilitation in cell
  double equilib_growth_rate;  // How fast cell produces biomass
                              // after growing to equilibrium
  int count = 1;  // For counting the number of cells occupied by the same community
  emp::vector<int> species; // The counts of each type in this cell
};

// A class to handle running the simple ecology
class AEcoWorld {
  private:
    // The world is a vector of cells (where each cell is
    // represented as a vector ints representing the count of
    // each type in each cell). 
    
    // Although the world is stored as a flat vector of cells
    // it represents a grid of cells
    using world_t = emp::vector<emp::vector<int> >;

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
    int MAX_POP;
    int curr_update;

    // Set up data tracking
    emp::DataFile data_file;
    emp::DataNode<double, emp::data::Stats> fitness_node;
    emp::DataNode<double, emp::data::Stats> synergy_node;
    emp::DataNode<double, emp::data::Stats> heredity_node;
    std::map<emp::vector<int>, CellData> dom_map;
    emp::vector<int> fittest;
    emp::vector<int> dominant;

  public:
  // Default constructor just has to set name of output file
  AEcoWorld() : data_file("a-eco_data.csv") {;}

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
      // std::cout << "Vector in interaction matrix at pos " << i << " in SetupRandomInteractions():" << std::endl;
      // std::cout << emp::to_string(interactions[i]) << std::endl;
    }
    // std::cout << std::endl;
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
      // std::cout << "Vector in interaction matrix at pos " << i << " in LoadInteractionMatrix():" << std::endl;
      // std::cout << emp::to_string(interactions[i]) << std::endl;
    }
    // std::cout << std::endl;
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
    MAX_POP = config->MAX_POP();

    // Set seed to configured value for reproducibility
    rnd.ResetSeed(config->SEED());

    // world vector needs a spot for each cell in the grid
    world.resize(config->WORLD_X() * config->WORLD_Y());
    // Initialize world vector
    for (emp::vector<int> & v : world) {
      v.resize(N_TYPES);
      for (int & count : v) {
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
    data_file.AddStats(fitness_node, "Fitness", "Fitness", true);
    data_file.AddStats(synergy_node, "Synergy", "Synergy", true);
    data_file.AddStats(heredity_node, "Heredity", "Heredity", true);
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
    for (emp::vector<int> & v : next_world) {
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

  // Handle the process of running the program through
  // all time steps
  void Run() {
    // Call update the specified number of times
    for (int i = 0; i < config->UPDATES(); i++) {
      Update(i);
      //std::cout << i << std::endl;
    }

    // Print out final state
    std::cout << "World Vectors:" << std::endl;
    for (auto & v : world) {
      std::cout << emp::to_string(v) << std::endl;
    }

    // Store interaction matrix in a file in case we
    // want to do stuff with it later
    WriteInteractionMatrix("interaction_matrix.dat");
  }

  // Calculate the growth rate (one measurement of fitness)
  // for a given cell
  double CalcGrowthRate(size_t pos, world_t & curr_world) {
    double growth_rate = 0;
    for (int i = 0; i < N_TYPES; i++) {
      double modifier = 0;
      for (int j = 0; j < N_TYPES; j++) {
        // Each type contributes to the growth rate modifier
        // of each other type based on the product of its
        // interaction value with that type and its population size
        modifier += interactions[i][j] * curr_world[pos][j];
      }

      // Add this type's overall growth rate to the cell-level
      // growth-rate
      growth_rate += ceil(modifier*((double)curr_world[pos][i]/(double)MAX_POP)); // * ((double)(MAX_POP - pos[i])/MAX_POP));
    }
    return growth_rate;
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
      int new_pop = ceil(modifier*curr_world[pos][i]); // * ((double)(MAX_POP - pos[i])/MAX_POP));
      // std::cout << pos[i] << " " << new_pop << " " << modifier << " " << ((double)(MAX_POP - pos[i])/MAX_POP) <<  std::endl;
      // Population size cannot be negative
      next_world[pos][i] = std::max(curr_world[pos][i] + new_pop, 0);
      // Population size capped at MAX_POP
      next_world[pos][i] = std::min(next_world[pos][i], MAX_POP);
    }
    // std::cout << emp::to_string(pos) << std::endl;
  }

  // Check whether the total population of a cell is large enough
  // for group level replication
  bool IsReproReady(size_t pos, world_t & w) {
    return emp::Sum(w[pos]) > config->MAX_POP()*(config->REPRO_THRESHOLD());
  }

  // Handle movement of biomass between cells
  void DoRepro(size_t pos, emp::vector<int> & adj, world_t & curr_world, world_t & next_world, double seed_prob, double prob_clear) {

    // Check whether conditions for group-level replication are met
    if (IsReproReady(pos, curr_world)) {
      //std::cout << "group repro" << std::endl;
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
        next_world[direction][i] = std::max(next_world[direction][i], 0);
      }
    }

    // Subtract material that diffused away from the current cell
    for (int i = 0; i < N_TYPES; i++) {
      next_world[pos][i] -= curr_world[pos][i] * config->DIFFUSION();
      // We can't have negative population sizes
      next_world[pos][i] = std::max(next_world[pos][i], 0);
      // Check whether we should randomly add one individual
      // of this type
      // From the community level perspective this is basically
      // a mutation
      // Note that this means we don't have mutations that remove
      // a member from a community, which is why PROB_CLEAR
      // is such an important parameter
      if (rnd.P(seed_prob)) {
        next_world[pos][i]++;
      }
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
                [](const std::pair<emp::vector<int>, CellData> & v1, const std::pair<emp::vector<int>, CellData> & v2 ){return v1.second.count < v2.second.count;})).first;
  }

  // Calculate fitness of individual cell
  CellData GetFitness(size_t orig_pos) {
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
    for (emp::vector<int> & v : test_world) {
      v.resize(config->N_TYPES(), 0);
    }

    // Put contents of cell being evaluated
    // in left-most spot in test world
    test_world[0] = world[orig_pos];

    // Make a next_world for the test world to
    // handle simulated updates in our test environment
    world_t next_world;
    next_world.resize(num_cells);
    for (emp::vector<int> & v : next_world) {
      v.resize(config->N_TYPES(), 0);
    }

    // Run DoGrowth for 20 time steps to hopefully
    // get cell to equilibrium
    // (if cell is an organism, this is analogous to 
    // development)
    for (int i = 0; i < 20; i++) {
      DoGrowth(0, test_world, next_world);
      std::swap(test_world, next_world);
    }

    // Record equilibrium community
    emp::vector<int> starting_point = test_world[0];
    // Calculate pure-math growth rate
    data.growth_rate = CalcGrowthRate(orig_pos, world);
    // Calculate synergy
    data.synergy = CalcSynergy(0, test_world);
    // Calculate equilibrium growth rate
    data.equilib_growth_rate = CalcGrowthRate(0, test_world);
    // Record community composition
    data.species = starting_point;

    // Add equilibrium growth rate to fitness data node
    fitness_node.Add(data.equilib_growth_rate);
    // Add synergy to synergy data node
    synergy_node.Add(data.synergy);

    // Now we're going to simulate up to time_limit number
    // of time steps to measure heredity and 
    // propogation-based fitness
    int time = 0;
    // Stop if hit the time limit or if the right-most cell is
    // able to do group-level reproduction (meaning the community
    // has spread across the entire test world and could keep going)
    while (!IsReproReady(num_cells - 1, test_world) && time < time_limit) {

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
      // Track time
      time++;
    }
    // Lower fitness is better by this metric - indicates
    // ability to spread from one end of the world to the other
    // more quickly
    // Note that even though this is called "fitness" it is currently
    // not the main fitness metric that we are using
    data.fitness = time;
    // Heredity is how similar the community in the final cell is to 
    // the equilibrium community
    data.heredity = 1.0 - emp::EuclideanDistance(test_world[num_cells - 1], starting_point) / ((double)MAX_POP*N_TYPES);

    // Keep track of the heredity values we've seen
    heredity_node.Add(data.heredity);

    return data;
  }

  // Getter for current update/time step
  int GetTime() {
    return curr_update;
  }

};