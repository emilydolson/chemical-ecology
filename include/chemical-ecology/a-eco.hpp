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
    int curr_update2;

    // Set up data tracking
    emp::DataFile data_file;
    emp::DataFile stochastic_data_file;
    emp::DataNode<double, emp::data::Stats> growth_rate_node;
    emp::DataNode<double, emp::data::Stats> heredity_node;
    emp::DataNode<double, emp::data::Stats> invasion_node;
    emp::DataNode<double, emp::data::Stats> biomass_node;
    std::map<emp::vector<double>, CellData> dom_map;
    emp::vector<double> fittest;
    emp::vector<double> dominant;
    world_t worldState;
    world_t stochasticWorldState;
    std::string worldType;

  public:
  // Default constructor just has to set name of output file
  AEcoWorld() : data_file("a-eco_data.csv"), stochastic_data_file("a-eco_model_data.csv") {;}

  // Initialize vector that keeps track of grid
  world_t world;
  
  // List of any isolated communities
  emp::vector<emp::vector<int>> subCommunities;

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
    /*
    data_file.AddStats(growth_rate_node, "Growth_Rate", "Growth_Rate", true);
    data_file.AddStats(heredity_node, "Heredity", "Heredity", true);
    data_file.AddStats(invasion_node, "Invasion_Ability", "Invasion_Ability", true);
    */
    data_file.AddStats(biomass_node, "Biomass", "Biomass", true);
    // Add columns calculated by running functions that look up values
    // in dom_map
    /*
    data_file.AddFun((std::function<int()>)[this](){return dom_map[fittest].count;}, "fittest_count", "fittest_count");
    data_file.AddFun((std::function<int()>)[this](){return dom_map[fittest].heredity;}, "fittest_heredity", "fittest_heredity");
    */
    data_file.AddFun((std::function<std::string()>)[this](){return emp::to_string(fittest);}, "fittestCommunity", "fittest community");
    /*
    data_file.AddFun((std::function<int()>)[this](){return dom_map[dominant].count;}, "dominant_count", "dominant_count");
    data_file.AddFun((std::function<int()>)[this](){return dom_map[dominant].heredity;}, "dominant_heredity", "dominant_heredity");
    data_file.AddFun((std::function<std::string()>)[this](){return emp::to_string(dominant);}, "dominantCommunity", "dominant community");
    data_file.AddFun((std::function<std::string()>)[this](){return emp::to_string(worldState);}, "worldState", "world state");
    */
    // Add function that runs before each row is printed in the datafile
    // to ensure fitnesses have been calculated (it's expensive so we don't
    // want to do it every update if we aren't going to use that data)
    data_file.AddPreFun([this](){CalcAllFitness();});
    // Set data_file to print a new row every 10 time steps
    data_file.SetTimingRepeat(10);
    // Print column names to data_file
    data_file.PrintHeaderKeys();

    stochastic_data_file.AddVar(curr_update2, "Time", "Time");
    //Takes all communities in the current world, and prints them to a file along with their stochastic world proportions
    stochastic_data_file.AddFun((std::function<std::string()>)[this](){return emp::to_string(stochasticWorldState);}, "stochasticWorldState", "stochastic world state");
    stochastic_data_file.AddVar(worldType, "worldType", "world type");
    stochastic_data_file.PrintHeaderKeys();

    // Make sure you get sub communities after setting up the matrix
    // otherwise the matrix will be empty when this is called
    subCommunities = findSubCommunities();
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

    // We need to handle cell updates in a random order
    // so that cells at the end of the world do not have an
    // advantage when group repro triggers
    int unordered[(int)world.size()];
    for (int i = 0; i < (int)world.size(); i++) {
      unordered[i] = i;
    }

    // Handle everything that allows biomass to move from
    // one cell to another. Do so for each cell

    for (int i = 0; i < (int)world.size(); i++) {

      int pos = unordered[i];
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
      DoRepro(pos, adj, world, next_world, config->SEEDING_PROB(), config->PROB_CLEAR(), config->DIFFUSION(), false);
    }

    // Update our time tracker
    curr_update = ud;
    
    // Give data_file the opportunity to write to the file
    worldState = next_world;
    data_file.Update(curr_update);

    // We're done calculating the type counts for the next
    // time step. We can now swap our counts for the next
    // time step into the main world variable 
    std::swap(world, next_world);

    // Clean-up data trackers
    dom_map.clear();
    fittest.clear();
    dominant.clear();
    worldState.clear();
  }

  //Depth first search allows us to recursively find sub communities
  void dfs(int root, emp::vector<bool>& visited, emp::vector<emp::vector<double>>& interactions, emp::vector<int>& component){
    visited[root] = true;
    component.push_back(root);

    for(size_t i = 0; i < interactions[root].size(); i++){
      //Check for incoming and outgoing edges
      if((interactions[root][i] != 0 && (!visited[i])) || (interactions[i][root] != 0 && (!visited[i]))){
        dfs(i, visited, interactions, component);
      }
    }
  }

  //Species cannot be a part of a community they have no interaction with
  //Find the connected components of the interaction matrix, and determine sub-communites
  emp::vector<emp::vector<int>> findSubCommunities(){
    emp::vector<emp::vector<double> > interactions = GetInteractions();
    emp::vector<bool> visited(interactions.size(), false);
    emp::vector<emp::vector<int>> components;

    for (size_t i = 0; i < interactions.size(); i++) {
      if (!visited[i]) {
        emp::vector<int> component;
        dfs(i, visited, interactions, component);
        components.push_back(component);
      }
    }
    return components;
  }
  

  world_t StochasticModel(int num_updates, bool repro, double prob_clear, double seeding_prob) {
    world_t model_world;
    world_t next_model_world;
    if(repro)
      worldType = "Repro";
    else
      worldType = "Soup";
    model_world.resize(config->WORLD_X() * config->WORLD_Y());
    for (emp::vector<double> & v : model_world) {
      v.resize(N_TYPES, 0);
    }
    next_model_world.resize(config->WORLD_X() * config->WORLD_Y());
    for (emp::vector<double> & v : next_model_world) {
      v.resize(N_TYPES, 0);
    }

    for(int i = 0; i < num_updates; i++) {
      //handle in cell growth
      for (size_t pos = 0; pos < model_world.size(); pos++) {
        DoGrowth(pos, model_world, next_model_world);
      }
      //handle abiotic parameters and group repro
      //adj will be empty for each pos, since there is no spatial structure
      //diffusion will be zero for the same reason
      emp::vector<int> adj = {};
      double diff = 0.0;
      for (size_t pos = 0; pos < model_world.size(); pos++) {
        DoRepro(pos, adj, model_world, next_model_world, config->SEEDING_PROB(), config->PROB_CLEAR(), diff, repro);
      }

      if(i%10 == 0){
        curr_update2 = i;
        stochasticWorldState = next_model_world;
        stochastic_data_file.Update(curr_update2);
      }

      std::swap(model_world, next_model_world);

      stochasticWorldState.clear();
    }

    return model_world;
  }

  //This function should be called to create a stable copy of the world for graph calculations
  world_t stableUpdate(int num_updates, world_t custom_world, int world_x, int world_y){
    world_t stable_world;
    world_t next_stable_world;

    stable_world.resize(world_x * world_y);
    for (emp::vector<double> & v : next_stable_world) {
      v.resize(N_TYPES, 0);
    }

    next_stable_world.resize(world_x * world_y);
    for (emp::vector<double> & v : next_stable_world) {
      v.resize(N_TYPES, 0);
    }

    //copy current world into stable world
    for(size_t i = 0; i < custom_world.size(); i++){
      stable_world[i].assign(custom_world[i].begin(), custom_world[i].end());
    }

    for(int i = 0; i < num_updates; i++){
      // Handle population growth for each cell
      for (size_t pos = 0; pos < custom_world.size(); pos++) {
        DoGrowth(pos, stable_world, next_stable_world);
      }

      std::swap(stable_world, next_stable_world);
    }

    return stable_world;
  }

  // Handle the process of running the program through
  // all time steps
  void Run() {
    //If there are any sub communities, we should know about them
    if(subCommunities.size() != 1){
      std::cout << "Sub-communities detected!: " << "\n";
      for (const auto& community : subCommunities) {
        std::cout << "sub_community:" << " ";
        for (const auto& species : community) {
          std::cout << species << " ";
        }
        std::cout << std::endl;
      }
    }

    // Call update the specified number of times
    for (int i = 0; i < config->UPDATES(); i++) {
      Update(i);
    }
    
    world_t stable_world = stableUpdate(10000, world, config->WORLD_X(), config->WORLD_Y());
    
    std::map<std::string, double> finalCommunities = getFinalCommunities(stable_world);
    std::map<std::string, double> assemblyFinalCommunities;
    std::map<std::string, double> adaptiveFinalCommunities;
    // Average n stochastic worlds
    int n = 10;
    for(int i = 0; i < n; i++){
      world_t assemblyModel = StochasticModel(1000, false, config->PROB_CLEAR(), config->SEEDING_PROB());
      world_t stableAssemblyModel = stableUpdate(10000, assemblyModel, config->WORLD_X(), config->WORLD_Y());
      world_t adaptiveModel = StochasticModel(1000, true, config->PROB_CLEAR(), config->SEEDING_PROB());
      world_t stableAdaptiveModel = stableUpdate(10000, adaptiveModel, config->WORLD_X(), config->WORLD_Y());
      std::map<std::string, double> assemblyCommunities = getFinalCommunities(stableAssemblyModel);
      std::map<std::string, double> adaptiveCommunities = getFinalCommunities(stableAdaptiveModel);
      for(auto& [node, proportion] : assemblyCommunities){
        if(assemblyFinalCommunities.find(node) == assemblyFinalCommunities.end()){
          assemblyFinalCommunities.insert({node, proportion});
        }
        else{
          assemblyFinalCommunities[node] += proportion;
        }
      }
      for(auto& [node, proportion] : adaptiveCommunities){
        if(adaptiveFinalCommunities.find(node) == adaptiveFinalCommunities.end()){
          adaptiveFinalCommunities.insert({node, proportion});
        }
        else{
          adaptiveFinalCommunities[node] += proportion;
        }
      }
    }
    // Actually average the summed proportions
    for (auto & comm : assemblyFinalCommunities) comm.second = comm.second/double(n);
    for (auto & comm : adaptiveFinalCommunities) comm.second = comm.second/double(n);

    std::ofstream FinalCommunitiesFile("stochastic_scores.txt");
    FinalCommunitiesFile << "Community Proportion AssemblyProportion AdaptiveProportion" << std::endl;
    for(auto& [node, proportion] : finalCommunities)
    {
      FinalCommunitiesFile << node << " " << proportion << " " << assemblyFinalCommunities[node] << " " << adaptiveFinalCommunities[node] << std::endl;
    }
    FinalCommunitiesFile.close();

    std::cout << "World" << std::endl;
    for(auto& [node, proportion] : finalCommunities)
    {
      std::cout << node << " " << proportion << std::endl;
    }

    std::cout << "Assembly" << std::endl;
    for(auto& [node, proportion] : assemblyFinalCommunities)
    {
      std::cout << node << " " << proportion << std::endl;
    }
    std::cout << "Adaptive" << std::endl;
    for(auto& [node, proportion] : adaptiveFinalCommunities)
    {
      std::cout << node << " " << proportion << std::endl;
    }
    
    //Print out final state if in verbose mode
    if(config->V()){
      std::cout << "World Vectors:" << std::endl;
      for (auto & v : world) {
        std::cout << emp::to_string(v) << std::endl;
      }
      std::cout << "Stable World Vectors:" << std::endl;
      for (auto & v : stable_world) {
        std::cout << emp::to_string(v) << std::endl;
      }
    }
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

  // Handles population growth of each type within a cell
  void DoGrowth(size_t pos, world_t & curr_world, world_t & next_world) {
    //For each species i
    for (int i = 0; i < N_TYPES; i++) {
      double modifier = 0;
      for (int j = 0; j < N_TYPES; j++) {
        // Sum up growth rate modifier for type i 
        modifier += interactions[i][j] * curr_world[pos][j];
      }

      // Logistic Growth
      double new_pop = modifier*curr_world[pos][i] * (1 - (curr_world[pos][i]/MAX_POP)); // * ((double)(MAX_POP - pos[i])/MAX_POP));
      // Population size cannot be negative
      next_world[pos][i] = std::max(curr_world[pos][i] + new_pop, 0.0);
      // Population size capped at MAX_POP
      next_world[pos][i] = std::min(next_world[pos][i], MAX_POP);
    }
  }

  // The probability of group reproduction is proportional to 
  // the biomass of the community 
  void doGroupRepro(size_t pos, world_t & w, world_t & next_w) {
    //Get these values once so they can be reused 
    int max_pop = config->MAX_POP();
    int types = config->N_TYPES();
    double dilution = config->REPRO_DILUTION();

    //Need to do GR in a random order, so the last sub-community does not have more repro power
    emp::Shuffle(rnd, subCommunities);
    for (const auto& community : subCommunities) {
      double pop = 0;
      for (const auto& species : community) {
        pop += w[pos][species];
      }
      double ratio = pop/(max_pop*types);
      // If group repro 
      if (rnd.P(ratio)){
        // Get a random cell that is not this cell to group level reproduce into 
        size_t new_pos;
        do {
          new_pos = rnd.GetInt((int)w.size() - 1);
        } while (new_pos == pos);
        for (int i = 0; i < N_TYPES; i++) {
          // Add a portion (configured by REPRO_DILUTION) of the quantity of the type
          // in the focal cell to the cell we're replicating into
          // NOTE: An important decision here is whether to clear the cell first.
          // We have chosen to, but can revisit that choice
          next_w[new_pos][i] = w[pos][i] * dilution;
        }
      }
    }
  }

  // Handle movement of biomass between cells
  void DoRepro(size_t pos, emp::vector<int> & adj, world_t & curr_world, world_t & next_world, double seed_prob, double prob_clear, double diffusion, bool GROUP_REPRO) {
    // Handle group repro for adaptive world
    if(GROUP_REPRO){
      doGroupRepro(pos, curr_world, next_world);
    }

    // Each cell has a chance of being cleared on every time step
    if (rnd.P(prob_clear)) {
      for (int i = 0; i < N_TYPES; i++) {
        next_world[pos][i] = 0;
      }
    }

    // Handle diffusion
    // Only diffuse in the real world
    // The adj vector should be empty for stochastic worlds, which do not have spatial structure
    if(adj.size() > 0){
      for (int direction : adj) {
        for (int i = 0; i < N_TYPES; i++) {
          // Calculate amount diffusing
          double avail = curr_world[pos][i] * diffusion;
          
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
        next_world[pos][i] -= curr_world[pos][i] * diffusion;
        // We can't have negative population sizes
        next_world[pos][i] = std::max(next_world[pos][i], 0.0);
      }
    }
    //Seed in
    if (rnd.P(seed_prob)) {
      int species = rnd.GetInt(N_TYPES);
      next_world[pos][species]++;
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

  std::map<std::string, double> getFinalCommunities(world_t stable_world) {
    double size = stable_world.size();
    std::map<std::string, double> finalCommunities;
    for(emp::vector<double> cell: stable_world){
      std::string temp = "";
      for(double species: cell){
      //for(int i = cell.size(); i --> 0;){
        //double species = cell[i];
        /*if(species > 0.001){
          temp.append("1");
        }
        else{
          temp.append("0");
        }*/
        temp.append(std::to_string(species)+" ");
      }
      //std::reverse(temp.begin(), temp.end());
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