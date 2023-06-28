#pragma once

#include <iostream>
#include <string>
#include <math.h>
#include "emp/Evolve/World.hpp"
#include "emp/math/distances.hpp"
#include "emp/math/random_utils.hpp"
#include "emp/datastructs/Graph.hpp"
#include "emp/bits/BitArray.hpp"
#include "emp/data/DataNode.hpp"
#include "emp/datastructs/map_utils.hpp"
#include "chemical-ecology/config_setup.hpp"
#include "chemical-ecology/SpatialStructure.hpp"

// TERMINOLOGY NOTES:
//
// - Cell = location in the world grid
// - Type = a type of thing that can live in a cell
//          (e.g. a chemical or a species). Usually
//          called a species in practice
// - Community = a set of types that interact and can
//                coexist in the world

// A class to handle running the simple ecology
class AEcoWorld {
  private:
    // The world is a vector of cells (where each cell is
    // represented as a vector ints representing the count of
    // each type in each cell).

    // Although the world is stored as a flat vector of cells
    // it represents a grid of cells
    using world_t = emp::vector< emp::vector<double> >;

    // The matrix of interactions between types
    // The diagonal of this matrix represents the
    // intrinsic growth rate (r) of each type
    emp::vector< emp::vector<double> > interactions;

    chemical_ecology::SpatialStructure spatial_structure;

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
    emp::DataNode<double, emp::data::Stats> biomass_node;
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

    data_file.AddVar(curr_update, "Time", "Time");
    data_file.AddFun((std::function<std::string()>)[this](){return emp::to_string(worldState);}, "worldState", "world state");
    data_file.SetTimingRepeat(10);
    data_file.PrintHeaderKeys();

    stochastic_data_file.AddVar(curr_update2, "Time", "Time");
    //Takes all communities in the current world, and prints them to a file along with their stochastic world proportions
    stochastic_data_file.AddFun((std::function<std::string()>)[this](){return emp::to_string(stochasticWorldState);}, "stochasticWorldState", "stochastic world state");
    stochastic_data_file.AddVar(worldType, "worldType", "world type");
    stochastic_data_file.SetTimingRepeat(10);
    stochastic_data_file.PrintHeaderKeys();

    // Make sure you get sub communities after setting up the matrix
    // otherwise the matrix will be empty when this is called
    subCommunities = findSubCommunities();
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

    world_t stable_world = stableUpdate(world, config->WORLD_X(), config->WORLD_Y());

    std::map<std::string, double> finalCommunities = getFinalCommunities(stable_world);
    std::map<std::string, double> assemblyFinalCommunities;
    std::map<std::string, double> adaptiveFinalCommunities;
    // Run n stochastic worlds
    int n = 10;
    for(int i = 0; i < n; i++){
      world_t assemblyModel = StochasticModel(config->UPDATES(), false, config->PROB_CLEAR(), config->SEEDING_PROB(), i);
      world_t stableAssemblyModel = stableUpdate(assemblyModel, config->WORLD_X(), config->WORLD_Y());
      world_t adaptiveModel = StochasticModel(config->UPDATES(), true, config->PROB_CLEAR(), config->SEEDING_PROB(), i);
      world_t stableAdaptiveModel = stableUpdate(adaptiveModel, config->WORLD_X(), config->WORLD_Y());
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
    // Average the summed proportions
    for (auto & comm : assemblyFinalCommunities) comm.second = comm.second/double(n);
    for (auto & comm : adaptiveFinalCommunities) comm.second = comm.second/double(n);

    // std::ofstream FinalCommunitiesFile("stochastic_scores.txt");
    // FinalCommunitiesFile << "Community Proportion AssemblyProportion AdaptiveProportion" << std::endl;
    // for(auto& [node, proportion] : finalCommunities)
    // {
    //   FinalCommunitiesFile << node << " " << proportion << " " << assemblyFinalCommunities[node] << " " << adaptiveFinalCommunities[node] << std::endl;
    // }
    // FinalCommunitiesFile.close();

    std::cout << "World" << std::endl;
    for(auto& [node, proportion] : finalCommunities)
    {
      std::cout << "Community: " << node << " Proportion: " << proportion << std::endl;
    }

    std::cout << "Assembly" << std::endl;
    for(auto& [node, proportion] : assemblyFinalCommunities)
    {
      std::cout << "Community: " << node << " Proportion: " << proportion << std::endl;
    }
    std::cout << "Adaptive" << std::endl;
    for(auto& [node, proportion] : adaptiveFinalCommunities)
    {
      std::cout << "Community: " << node << " Proportion: " << proportion << std::endl;
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
    worldState.clear();
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

  //This function should be called to create a stable copy of the world
  world_t stableUpdate(world_t custom_world, int world_x, int world_y, int max_updates=10000){
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

    for(int i = 0; i < max_updates; i++){
      // Handle population growth for each cell
      for (size_t pos = 0; pos < custom_world.size(); pos++) {
        DoGrowth(pos, stable_world, next_stable_world);
      }

      // We can stop iterating early if the world has already stabilized
      double delta = 0;
      double epsilon = .0001;
      for (size_t j = 0; j < stable_world.size(); j++){
        delta += emp::EuclideanDistance(stable_world[j], next_stable_world[j]);
      }
      // If the change from one world to the next is very small, return early
      if(delta < epsilon){
        for (size_t pos = 0; pos < stable_world.size(); pos++) {
          for (size_t s = 0; s < stable_world[pos].size(); s++) {
            if (stable_world[pos][s] < 1) {
              stable_world[pos][s] = 0.0;
            }
          }
        }
        return stable_world;
      }
      std::swap(stable_world, next_stable_world);
    }
    std::cout << "\n Max number of stable updates reached" << std::endl;
    for (size_t pos = 0; pos < stable_world.size(); pos++) {
      for (size_t s = 0; s < stable_world[pos].size(); s++) {
        if (stable_world[pos][s] < 1) {
          stable_world[pos][s] = 0.0;
        }
      }
    }
    return stable_world;
  }

  world_t StochasticModel(int num_updates, bool repro, double prob_clear, double seeding_prob, int iter) {
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

      //Only plot the first runs of each
      if(i%10 == 0 && iter == 0){
        curr_update2 = i;
        stochasticWorldState = next_model_world;
        stochastic_data_file.Update(curr_update2);
      }

      std::swap(model_world, next_model_world);

      if(iter == 0){
        stochasticWorldState.clear();
      }
    }

    return model_world;
  }

  std::map<std::string, double> getFinalCommunities(world_t stable_world) {
    double size = stable_world.size();
    std::map<std::string, double> finalCommunities;
    for(emp::vector<double> cell: stable_world){
      std::string comm = "";
      for(size_t i = 0; i < cell.size(); i++){
        std::string count = std::to_string(cell[i]);
        // Formatting to make the strings more legible
        count.erase(count.find_last_not_of('0') + 1, std::string::npos);
        count.erase(count.find_last_not_of('.') + 1, std::string::npos);
        // Check if the species is present, and if it is make sure it has interactions
        if(count.compare("0") != 0){
          int sc = -1;
          for(size_t j = 0; j < subCommunities.size(); j++){
            for(int x : subCommunities[j]){
              if(x == i){
                sc = j;
              }
            }
          }
          bool interacts = false;
          for(int m : subCommunities[sc]){
            if(cell[m] != 0 && m != i){
              comm.append(count + " ");
              interacts = true;
              break;
            }
          }
          if(!interacts){
            comm.append("PWI ");
          }
        }
        // If there are zero present
        else{
          comm.append("0 ");
        }
      }
      if(finalCommunities.find(comm) == finalCommunities.end()){
          finalCommunities.insert({comm, 1});
      }
      else{
        finalCommunities[comm] += 1;
      }
    }
    // Get proportion
    for(auto& [key, val] : finalCommunities){
      val = val/size;
    }
    return finalCommunities;
  }

  // Helper functions not related to the running of the worlds go down here

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

  //Depth first search allows us to recursively find sub communities (connected components in the interaction matrix)
  void dfs(int root, emp::vector<bool>& visited, emp::vector<emp::vector<double>>& interactions, emp::vector<int>& component){
    visited[root] = true;
    component.push_back(root);

    for(size_t i = 0; i < interactions[root].size(); i++){
      //Both incoming and outgoing edges count as interactions
      if((interactions[root][i] != 0 && (!visited[i])) || (interactions[i][root] != 0 && (!visited[i]))){
        dfs(i, visited, interactions, component);
      }
    }
  }

  // Calculate the growth rate (one measurement of fitness)
  // for a given cell
  double CalcGrowthRate(size_t pos, world_t & curr_world) {
    return doCalcGrowthRate(curr_world[pos]);
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

  // Getter for current update/time step
  int GetTime() {
    return curr_update;
  }

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
};