#pragma once

#include <iostream>
#include <string>
#include <cmath>
#include <sys/stat.h>
#include <algorithm>
#include <functional>
#include <limits>

#include "emp/Evolve/World.hpp"
#include "emp/math/distances.hpp"
#include "emp/math/random_utils.hpp"
#include "emp/datastructs/Graph.hpp"
#include "emp/bits/BitVector.hpp"
#include "emp/data/DataNode.hpp"
#include "emp/base/Ptr.hpp"
#include "emp/datastructs/map_utils.hpp"
#include "emp/datastructs/map_utils.hpp"

#include "chemical-ecology/SpatialStructure.hpp"
#include "chemical-ecology/CommunityStructure.hpp"
#include "chemical-ecology/RecordedCommunitySummarizer.hpp"
#include "chemical-ecology/RecordedCommunitySet.hpp"
#include "chemical-ecology/Config.hpp"
#include "chemical-ecology/utils/graph_utils.hpp"
#include "chemical-ecology/InteractionMatrix.hpp"

namespace chemical_ecology {

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
public:
  // The world is a vector of cells (where each cell is
  // represented as a vector ints representing the count of
  // each type in each cell).
  // Although the world is stored as a flat vector of cells
  // it represents a grid of cells
  using world_t = emp::vector< emp::vector<double> >;
  using config_t = Config;

  // Used to help output recorded community summary sets
  template<typename SUMMARY_SET_KEY_T>
  struct RecordedCommunitySetInfo {
    const RecordedCommunitySet<SUMMARY_SET_KEY_T>& summary_set;
    std::string source;
    bool stabilized;      // Was this community "stabilized" (i.e., GenStablizedWorld)
    size_t updates;       // How many updates was this community run for?

    RecordedCommunitySetInfo(
      const RecordedCommunitySet<SUMMARY_SET_KEY_T>& comm_summary_set,
      const std::string& set_source,
      bool set_stabilized,
      size_t set_updates
    ) :
      summary_set(comm_summary_set),
      source(set_source),
      stabilized(set_stabilized),
      updates(set_updates)
    { }
  };

private:

  // The matrix of interactions between types
  // The diagonal of this matrix represents the
  // intrinsic growth rate (r) of each type
  InteractionMatrix interactions;

  SpatialStructure diffusion_spatial_structure;
  SpatialStructure group_repro_spatial_structure;
  size_t world_size = 0;

  // A random number generator for all our random number
  // generating needs
  emp::Random rnd;

  // All configuration information is stored in config
  emp::Ptr<chemical_ecology::Config> config = nullptr;

  // These values are set in the config but we have local
  // copies for efficiency in accessing their values
  size_t N_TYPES;
  double MAX_POP;
  int curr_update;
  int curr_update2;  // NOTE (@AML): if not using negative values, switch to size_t, rename to be more informative

  // Initialize vector that keeps track of grid
  world_t world;

  // Manages community structure (determined by interaction matrix)
  CommunityStructure community_structure;
  emp::vector<size_t> subcommunity_group_repro_schedule;

  // Used to track activation order of positions in the world
  emp::vector<size_t> position_activation_order;

  emp::Ptr<RecordedCommunitySummarizer> community_summarizer_raw;         // Uses raw species counts.
  emp::Ptr<RecordedCommunitySummarizer> community_summarizer_pwip;        // Keeps all species present with valid interaction paths to other present species

  // Set up data tracking
  std::string output_dir;
  emp::Ptr<emp::DataFile> data_file = nullptr;
  emp::Ptr<emp::DataFile> stochastic_data_file = nullptr;

  world_t worldState;
  world_t stochasticWorldState;
  std::string worldType;

  // Configures spatial structure based on world configuration.
  void SetupSpatialStructure();
  void SetupSpatialStructure_ToroidalGrid(SpatialStructure& spatial_structure);
  void SetupSpatialStructure_WellMixed(SpatialStructure& spatial_structure);
  void SetupSpatialStructure_Load(
    SpatialStructure& spatial_structure,
    const std::string& file_path,
    const std::string& load_mode
  );

  void SetupCommunitySummarizers();

  void AnalyzeWorldCommunities(
    const std::string& output_file_path,
    const RecordedCommunitySet<emp::vector<double>>& world_communities,
    const RecordedCommunitySet<emp::vector<double>>& assembly_communities,
    const RecordedCommunitySet<emp::vector<double>>& adaptive_communities
  );

  // Output a snapshot of how the world is configured
  void SnapshotConfig();

  // Output a snapshot of the community structure identified from interaction matrix
  void SnapshotSubCommunities();

  // Output a snapshot of the community structure
  template<typename SUMMARY_SET_KEY_T>
  void SnapshotRecordedCommunitySets(
    const std::string& file_path,
    const emp::vector<RecordedCommunitySetInfo<SUMMARY_SET_KEY_T>>& community_sets
  );

public:

  AEcoWorld() = default;

  ~AEcoWorld() {
    if (data_file != nullptr) data_file.Delete();
    if (stochastic_data_file != nullptr) stochastic_data_file.Delete();
    if (community_summarizer_raw != nullptr) community_summarizer_raw.Delete();
    if (community_summarizer_pwip != nullptr) community_summarizer_pwip.Delete();
  }

  // Setup the world according to the specified configuration
  void Setup(config_t& cfg) {

    // Store cfg for future reference
    config = &cfg;

    // Set local config variables based on given configuration
    N_TYPES = config->N_TYPES();
    MAX_POP = double(config->MAX_POP());

    // Set seed to configured value for reproducibility
    // NOTE (@AML): Make sure to be using updated version of Empirical with patch for ResetSeed function!
    rnd.ResetSeed(config->SEED());

    // Setup spatial structure (configures world size)
    world_size = 0; // World size not valid until after setting up spatial structure
    SetupSpatialStructure();

    // world vector needs a spot for each cell in the grid
    world.resize(
      world_size,
      emp::vector<double>(N_TYPES, 0.0)
    );

    // Initialize world vector
    for (emp::vector<double>& v : world) {
      for (double& count : v) {
        // The quantity of each type in each cell is either 0 or 1
        // The probability of it being 1 is controlled by SEEDING_PROB
        count = (double)rnd.P(config->SEEDING_PROB());
      }
    }

    // Initialize activation order
    position_activation_order.resize(world.size(), 0);
    std::iota(
      position_activation_order.begin(),
      position_activation_order.end(),
      0
    );

    // Setup interaction matrix based on the method
    // specified in the configuration file
    interactions.SetInteractsFun(
      [](const auto& mat, size_t from, size_t to) -> bool {
        return (mat[from][to] != 0) || (mat[to][from] != 0);
      }
    );
    if (config->INTERACTION_SOURCE() == "") {
      SetupRandomInteractions();
    } else {
      interactions.LoadInteractions(
        config->INTERACTION_SOURCE(),
        N_TYPES
      );
    }

    // Make sure you get sub communities after setting up the matrix
    community_structure.SetStructure(
      interactions
    );
    // Initialize schedule for sub-community group repro
    subcommunity_group_repro_schedule.resize(
      community_structure.GetNumSubCommunities(),
      0
    );
    std::iota(
      subcommunity_group_repro_schedule.begin(),
      subcommunity_group_repro_schedule.end(),
      0
    );

    // Configure community summarizers (used to report summaries of recorded communities)
    // NOTE: should be called after community structure has been configured
    SetupCommunitySummarizers();

    // Configure output directory path, create directory
    output_dir = config->OUTPUT_DIR();
    mkdir(output_dir.c_str(), ACCESSPERMS);
    if(output_dir.back() != '/') {
        output_dir += '/';
    }

    // NOTE (@AML): Could merge data_file and stochastic_data_file => same information being printed for each
    data_file = emp::NewPtr<emp::DataFile>(output_dir + "a-eco_data.csv");
    data_file->AddVar(curr_update, "Time", "Time");
    data_file->AddFun((std::function<std::string()>)[this](){return emp::to_string(worldState);}, "worldState", "world state");
    data_file->SetTimingRepeat(config->OUTPUT_RESOLUTION());
    data_file->PrintHeaderKeys();

    stochastic_data_file = emp::NewPtr<emp::DataFile>(output_dir + "a-eco_model_data.csv");
    stochastic_data_file->AddVar(curr_update2, "Time", "Time");
    // Takes all communities in the current world, and prints them to a file along with their stochastic world proportions
    stochastic_data_file->AddFun((std::function<std::string()>)[this](){return emp::to_string(stochasticWorldState);}, "stochasticWorldState", "stochastic world state");
    stochastic_data_file->AddVar(worldType, "worldType", "world type");
    stochastic_data_file->SetTimingRepeat(config->OUTPUT_RESOLUTION());
    stochastic_data_file->PrintHeaderKeys();

    // Output a snapshot of identified subcommunities
    SnapshotSubCommunities();

    // Output a snapshot of the run configuration
    SnapshotConfig();
  }

  // Handle the process of running the program through
  // all time steps
  void Run() {

    // Call update the specified number of times
    for (size_t i = 0; i < config->UPDATES(); i++) {
      Update(i);
    }

    // Run world forward without diffusion
    world_t stable_world(GenStabilizedWorld(world, config->CELL_STABILIZATION_UPDATES()));

    // NOTE (@AML): Not in total love with this; could possibly use another iteration after chatting
    //   about future functionality that would be useful to have
    RecordedCommunitySet<emp::vector<double>>::summary_key_fun_t recorded_comm_key_fun = [](
      const RecordedCommunitySummary& summary
    ) -> const auto& {
      return summary.counts;
    };

    // Set of recorded communities with "raw" counts
    RecordedCommunitySet<emp::vector<double>> recorded_communities_world_raw(recorded_comm_key_fun);
    RecordedCommunitySet<emp::vector<double>> recorded_communities_assembly_raw(recorded_comm_key_fun);
    RecordedCommunitySet<emp::vector<double>> recorded_communities_adaptive_raw(recorded_comm_key_fun);

    // Set of recorded communities with counts that have had PNI-species zeroed-out
    RecordedCommunitySet<emp::vector<double>> recorded_communities_world_pwip(recorded_comm_key_fun);
    RecordedCommunitySet<emp::vector<double>> recorded_communities_assembly_pwip(recorded_comm_key_fun);
    RecordedCommunitySet<emp::vector<double>> recorded_communities_adaptive_pwip(recorded_comm_key_fun);

    // Add world summaries
    recorded_communities_world_raw.Add(
      community_summarizer_raw->SummarizeAll(stable_world)
    );
    recorded_communities_world_pwip.Add(
      community_summarizer_pwip->SummarizeAll(stable_world)
    );

    // Run n stochastic worlds
    for (size_t i = 0; i < config->STOCHASTIC_ANALYSIS_REPS(); ++i) {
      const bool record_analysis_state = (i == 0); // NOTE: any reason to not output all?
      // Run stochastic assembly model
      // - NOTE (@AML): I'd be tempted to splitting the adaptive and assembly models into their own function.
      //   con: repeated code; pro: simpler parameters, can then have separate implementations down the line if necessary
      world_t assemblyModel = StochasticModel(config->UPDATES(), false, config->PROB_CLEAR(), config->SEEDING_PROB(), record_analysis_state);
      world_t stableAssemblyModel = GenStabilizedWorld(assemblyModel, config->CELL_STABILIZATION_UPDATES());

      // Run stochastic adaptive model
      world_t adaptiveModel = StochasticModel(config->UPDATES(), true, config->PROB_CLEAR(), config->SEEDING_PROB(), record_analysis_state);
      world_t stableAdaptiveModel = GenStabilizedWorld(adaptiveModel, config->CELL_STABILIZATION_UPDATES());

      // Add summarized recorded communities to sets
      recorded_communities_assembly_raw.Add(
        community_summarizer_raw->SummarizeAll(stableAssemblyModel)
      );
      recorded_communities_assembly_pwip.Add(
        community_summarizer_pwip->SummarizeAll(stableAssemblyModel)
      );
      recorded_communities_adaptive_raw.Add(
        community_summarizer_raw->SummarizeAll(stableAdaptiveModel)
      );
      recorded_communities_adaptive_pwip.Add(
        community_summarizer_pwip->SummarizeAll(stableAdaptiveModel)
      );

    }

    // NOTE (@AML): Slightly clunky way to tie together things
    // Snapshot raw recorded community summaries
    SnapshotRecordedCommunitySets</*SUMMARY_SET_KEY_T=*/emp::vector<double>>(
      output_dir + "recorded_communities_raw.csv",
      {
        {recorded_communities_world_raw, "world", true, config->UPDATES()},
        {recorded_communities_assembly_raw, "assembly", true, config->UPDATES()},
        {recorded_communities_adaptive_raw, "adaptive", true, config->UPDATES()}
      }
    );

    // Snapshot summaries where "present-no-interactions" species have been removed
    SnapshotRecordedCommunitySets</*SUMMARY_SET_KEY_T=*/emp::vector<double>>(
      output_dir + "recorded_communities_pwip.csv",
      {
        {recorded_communities_world_pwip, "world", true, config->UPDATES()},
        {recorded_communities_assembly_pwip, "assembly", true, config->UPDATES()},
        {recorded_communities_adaptive_pwip, "adaptive", true, config->UPDATES()}
      }
    );

    // Summarize world
    AnalyzeWorldCommunities(
      output_dir + "world_summary_pwip.csv",
      recorded_communities_world_pwip,
      recorded_communities_assembly_pwip,
      recorded_communities_adaptive_pwip
    );

    // ---- Commandline summary output ---
    // NOTE (@AML): Can add updated commandline output back if folks want!
    //               - Didn't add it because it can be unweildy for many model
    //                 paramterizations and looking at the output files is easy enough

    //Print out final state if in verbose mode
    if (config->V()) {
      std::cout << "World Vectors:" << std::endl;
      for (const auto& v : world) {
        std::cout << emp::to_string(v) << std::endl;
      }
      std::cout << "Stable World Vectors:" << std::endl;
      for (const auto& v : stable_world) {
        std::cout << emp::to_string(v) << std::endl;
      }
    }
  }

  // Handle an individual time step
  // ud = which time step we're on
  void Update(size_t ud) {

    // Create a new world object to store the values
    // for the next time step
    world_t next_world(
      world_size,
      emp::vector<double>(N_TYPES, 0)
    );

    // Handle population growth for each cell
    for (size_t pos = 0; pos < world.size(); pos++) {
      DoGrowth(pos, world, next_world); // @AML: BOOKMARK
    }

    // We need to handle cell updates in a random order
    // so that cells at the end of the world do not have an
    // advantage when group repro triggers
    emp::Shuffle(rnd, position_activation_order);

    // Handle everything that allows biomass to move from
    // one cell to another. Do so for each cell
    for (size_t i = 0; i < world.size(); ++i) {

      const size_t pos = position_activation_order[i];
      // Actually call function that handles between-cell
      // movement
      // (1) Do group reproduction?
      if (config->GROUP_REPRO()) {
        DoGroupRepro(pos, world, next_world);
      }
      // (2) Do cell clearing
      DoClearing(pos, world, next_world, config->PROB_CLEAR());
      // (3) Do diffusion
      DoDiffusion(pos, world, next_world, config->DIFFUSION());
      // (4) Do seeding
      DoSeeding(pos, world, next_world, config->SEEDING_PROB());

    }

    // Update our time tracker
    curr_update = ud;

    // Give data_file the opportunity to write to the file
    worldState = next_world;
    data_file->Update(curr_update);

    // We're done calculating the type counts for the next
    // time step. We can now swap our counts for the next
    // time step into the main world variable
    std::swap(world, next_world);

    // Clean-up data trackers
    worldState.clear();
  }

  // Handles population growth of each type within a cell
  void DoGrowth(size_t pos, const world_t& curr_world, world_t& next_world) {
    //For each species i
    for (size_t i = 0; i < N_TYPES; i++) {
      double modifier = 0;
      for (size_t j = 0; j < N_TYPES; j++) {
        // Sum up growth rate modifier for type i
        // NOTE (@AML): does directionality [i][j] [j][i] matter here? (i.e., are interaction graphs directed or undirected?)
        // Updated species I, species I changing based on interaction[i][j]
        // Effect species j has on species i is inter[i][j]
        modifier += interactions[i][j] * curr_world[pos][j];
      }
      const double cur_count = curr_world[pos][i];
      // Logistic Growth
      const double new_pop = modifier * cur_count * (1 - (cur_count/MAX_POP)); // * ((double)(MAX_POP - pos[i])/MAX_POP));
      // Population size cannot be negative
      next_world[pos][i] = std::max(cur_count + new_pop, 0.0);
      // Population size capped at MAX_POP
      next_world[pos][i] = std::min(next_world[pos][i], MAX_POP);
      emp_assert(next_world[pos][i] <= MAX_POP);
    }
  }

  // The probability of group reproduction is proportional to
  // the biomass of the community
  void DoGroupRepro(size_t pos, const world_t& w, world_t& next_w) {
    // Get these values once so they can be reused
    const int max_pop = config->MAX_POP();
    const size_t types = config->N_TYPES();
    const double dilution = config->REPRO_DILUTION();
    emp_assert(max_pop * types > 0);

    // Get neighbors
    const auto& neighbors = group_repro_spatial_structure.GetNeighbors(pos);
    const size_t num_neighbors = neighbors.size();
    // If no neighbors, no valid destination to reproduce into
    if (num_neighbors == 0) return;

    // Need to do GR in a random order, so the last sub-community does not have more repro power
    // emp::Shuffle(rnd, subCommunities);
    // for (const auto& community : subCommunities) {
    emp::Shuffle(rnd, subcommunity_group_repro_schedule);
    for (size_t schedule_i = 0; schedule_i < subcommunity_group_repro_schedule.size(); ++schedule_i) {
      const size_t community_id = subcommunity_group_repro_schedule[schedule_i];
      const auto& community = community_structure.GetSubCommunity(community_id);
      // Compute population size of this community
      double pop = 0;
      for (const auto& species : community) {
        pop += w[pos][species];
      }
      const double ratio = pop / (max_pop * types);

      // If group repro
      if (rnd.P(ratio)) {

        // Get a random neighboring cell to reproduce into
        const auto rnd_neighbor = group_repro_spatial_structure.GetRandomNeighbor(
          rnd,
          pos
        ); // Returned if not neighbors, this should always be valid.
        emp_assert(rnd_neighbor);
        const size_t new_pos = rnd_neighbor.value();

        for (size_t i = 0; i < N_TYPES; i++) {
          // Add a portion (configured by REPRO_DILUTION) of the quantity of the type
          // in the focal cell to the cell we're replicating into
          // NOTE: An important decision here is whether to clear the cell first.
          // We have chosen to, but can revisit that choice
          next_w[new_pos][i] = w[pos][i] * dilution;
        }
      }

    }
  }

  void DoClearing(size_t pos, const world_t& curr_world, world_t& next_world, double prob_clear) {
    // Each cell has a chance of being cleared on every time step
    if (rnd.P(prob_clear)) {
      for (size_t i = 0; i < N_TYPES; i++) {
        next_world[pos][i] = 0;
      }
    }
  }

  void DoDiffusion(size_t pos, const world_t& curr_world, world_t& next_world, double diffusion) {
    // Handle diffusion
    // Only diffuse in the real world
    // The adj vector should be empty for stochastic worlds, which do not have spatial structure
    const auto& neighbors = diffusion_spatial_structure.GetNeighbors(pos);
    const size_t num_neighbors = neighbors.size();
    if (num_neighbors > 0) {
      // Diffuse to neighbors
      for (size_t neighbor : neighbors) {
        for (size_t i = 0; i < N_TYPES; ++i) {

          // Calculate amount diffusing
          double avail = curr_world[pos][i] * diffusion;

          // Evenly distribute the diffusion
          next_world[neighbor][i] +=  avail / num_neighbors;

          // make sure final value is legal number
          next_world[neighbor][i] = std::min(next_world[neighbor][i], MAX_POP);
          next_world[neighbor][i] = std::max(next_world[neighbor][i], 0.0);
        }
      }

      // Subtract diffusion
      for (size_t i = 0; i < N_TYPES; ++i) {
        next_world[pos][i] -= curr_world[pos][i] * diffusion;
        // We can't have negative population sizes
        next_world[pos][i] = std::max(next_world[pos][i], 0.0);
      }
    }
  }

  void DoSeeding(size_t pos, const world_t& curr_world, world_t& next_world, double seed_prob) {
    // Seed in
    if (rnd.P(seed_prob)) {
      size_t species = rnd.GetUInt(N_TYPES);
      next_world[pos][species]++;
      // NOTE (@AML): want to respect max pop here?
      next_world[pos][species] = std::min(next_world[pos][species], MAX_POP);
    }
  }

  // This function should be called to create a stable copy of the world
  world_t GenStabilizedWorld(const world_t& custom_world, size_t max_updates=10000){

    // Track current and next state of world
    world_t stable_world(
      world_size,
      emp::vector<double>(N_TYPES, 0.0)
    );
    world_t next_stable_world(stable_world);

    // Copy current world into stable world
    std::copy(
      custom_world.begin(),
      custom_world.end(),
      stable_world.begin()
    );
    emp_assert(custom_world == stable_world);

    for (size_t i = 0; i < max_updates; i++) {

      // Handle population growth for each cell
      for (size_t pos = 0; pos < stable_world.size(); pos++) {
        DoGrowth(pos, stable_world, next_stable_world);
      }

      // We can stop iterating early if the world has already stabilized
      double delta = 0;
      const double epsilon = config->CELL_STABILIZATION_EPSILON();
      for (size_t j = 0; j < stable_world.size(); j++) {
        delta += emp::EuclideanDistance(stable_world[j], next_stable_world[j]);
      }
      // If the change from one world to the next is very small, return early
      if (delta < epsilon) {
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

  // NOTE: The base model is also stochastic model ==> this is a confusing name.
  world_t StochasticModel(
    int num_updates,
    bool repro,
    double prob_clear,
    double seeding_prob,
    bool record_world_state=false
  ) {

    worldType = (repro) ? "Repro" : "Soup";

    // Track current and next stochastic model worlds.
    world_t model_world(
      world_size,
      emp::vector<double>(N_TYPES, 0.0)
    );
    world_t next_model_world(model_world);

    for (int i = 0; i < num_updates; i++) {
      // handle in cell growth
      for (size_t pos = 0; pos < model_world.size(); pos++) {
        DoGrowth(pos, model_world, next_model_world);
      }
      // Handle abiotic parameters and group repro
      // There is no spatial structure / no diffusion.
      for (size_t pos = 0; pos < model_world.size(); pos++) {
        // ORIGINAL DoRepro call:
        //   DoRepro(pos, adj, model_world, next_model_world, config->SEEDING_PROB(), config->PROB_CLEAR(), diff, repro);
        // (1) Group repro
        if (repro) {
          DoGroupRepro(pos, model_world, next_model_world);
        }
        // (2) clearing
        DoClearing(pos, model_world, next_model_world, config->PROB_CLEAR());
        // (3) seeding
        DoSeeding(pos, model_world, next_model_world, config->SEEDING_PROB());
      }

      // Record world state
      const bool final_update = ((i + 1) == num_updates);
      const bool res_update =  !(bool)(i % config->OUTPUT_RESOLUTION());
      if (record_world_state && (final_update || res_update)) {
        curr_update2 = i;
        stochasticWorldState = next_model_world;
        stochastic_data_file->Update(); // Don't need to call with update #, already performing timing check.
      }

      std::swap(model_world, next_model_world);

      if (record_world_state) {
        stochasticWorldState.clear();
      }

    }

    return model_world;
  }

  const world_t& GetWorld() const { return world; }

  size_t GetWorldSize() const { return world_size; }

  // Calculate the growth rate (one measurement of fitness)
  // for a given cell
  double CalcGrowthRate(size_t pos, const world_t& curr_world) {
    return doCalcGrowthRate(curr_world[pos]);
  }

  double doCalcGrowthRate(const emp::vector<double>& community){
    double growth_rate = 0;
    for (size_t i = 0; i < N_TYPES; i++) {
      double modifier = 0;
      for (size_t j = 0; j < N_TYPES; j++) {
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
  const InteractionMatrix& GetInteractions() const {
    return interactions;
  }

  // TODO - cleanup after switching to InteractionMatrix class
  // Set the interaction matrix to contain random values
  // (with probabilities determined by configs)
  void SetupRandomInteractions() {
    // interaction matrix will ultimately have size
    interactions.RandomizeInteractions(
      rnd,
      N_TYPES,
      config->PROB_INTERACTION(),
      config->INTERACTION_MAGNITUDE()
      // NOTE (@AML): Self-interactions? Looks like they were enabled in old implementation.
    );
  }

}; // End AEcoWorld class definition

void AEcoWorld::SetupSpatialStructure() {
  // Configure diffusion spatial structure
  if (config->DIFFUSION_SPATIAL_STRUCTURE() == "toroidal-grid") {
    SetupSpatialStructure_ToroidalGrid(diffusion_spatial_structure);
  } else if (config->DIFFUSION_SPATIAL_STRUCTURE() == "well-mixed") {
    SetupSpatialStructure_WellMixed(diffusion_spatial_structure);
  } else if (config->DIFFUSION_SPATIAL_STRUCTURE() == "load") {
    SetupSpatialStructure_Load(
      diffusion_spatial_structure,
      config->DIFFUSION_SPATIAL_STRUCTURE_FILE(),
      config->DIFFUSION_SPATIAL_STRUCTURE_LOAD_MODE()
    );
  } else {
    std::cout << "Unknown diffusion spatial structure: " << config->DIFFUSION_SPATIAL_STRUCTURE() << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(-1);
  }

  // Configure group repro spatial structure
  if (config->GROUP_REPRO_SPATIAL_STRUCTURE() == "toroidal-grid") {
    SetupSpatialStructure_ToroidalGrid(group_repro_spatial_structure);
  } else if (config->GROUP_REPRO_SPATIAL_STRUCTURE() == "well-mixed") {
    SetupSpatialStructure_WellMixed(group_repro_spatial_structure);
  } else if (config->GROUP_REPRO_SPATIAL_STRUCTURE() == "load") {
    SetupSpatialStructure_Load(
      group_repro_spatial_structure,
      config->GROUP_REPRO_SPATIAL_STRUCTURE_FILE(),
      config->GROUP_REPRO_SPATIAL_STRUCTURE_LOAD_MODE()
    );
  } else {
    std::cout << "Unknown group repro spatial structure: " << config->GROUP_REPRO_SPATIAL_STRUCTURE() << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(-1);
  }

  // Group repro and diffusion worlds must have the same number of positions
  emp_assert(
    group_repro_spatial_structure.GetNumPositions() == diffusion_spatial_structure.GetNumPositions(),
    "Group repro and diffusion spatial structures must have an equivalent number of positions"
  );

  if (group_repro_spatial_structure.GetNumPositions() != diffusion_spatial_structure.GetNumPositions()) {
    std::cout << "Group repro and diffusion spatial structures do not have the same number of positions." << std::endl;
    std::cout << "  Exiting." << std::endl;
    exit(-1);
  }

  world_size = diffusion_spatial_structure.GetNumPositions();

}

// Configures spatial structure as 2d toroidal grid
void AEcoWorld::SetupSpatialStructure_ToroidalGrid(SpatialStructure& spatial_structure) {
  ConfigureToroidalGrid(
    spatial_structure,
    config->WORLD_WIDTH(),
    config->WORLD_HEIGHT()
  );
}

// Configures spatial structure to be fully connected
void AEcoWorld::SetupSpatialStructure_WellMixed(SpatialStructure& spatial_structure) {
  ConfigureFullyConnected(
    spatial_structure,
    config->WORLD_WIDTH() * config->WORLD_HEIGHT()
  );
}

// Loads spatial structure from file
void AEcoWorld::SetupSpatialStructure_Load(
  SpatialStructure& spatial_structure,
  const std::string& file_path,
  const std::string& load_mode
) {
  if (load_mode == "edges") {
    spatial_structure.LoadStructureFromEdgeCSV(file_path);
  } else if (load_mode == "matrix") {
    spatial_structure.LoadStructureFromMatrix(file_path);
  } else {
    std::cout << "Unknown spatial structure load mode: " << load_mode << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(-1);
  }
}

// Configures community summerizers
void AEcoWorld::SetupCommunitySummarizers() {
  emp_assert(community_summarizer_raw == nullptr);
  emp_assert(community_summarizer_pwip == nullptr);

  // Reports raw counts (with decimal components truncated)
  community_summarizer_raw = emp::NewPtr<RecordedCommunitySummarizer>(
    community_structure,
    [](double count) -> bool { return count >= 1.0; }
  );

  // Removes species present that have no interactions with other present species
  community_summarizer_pwip = emp::NewPtr<RecordedCommunitySummarizer>(
    community_structure,
    [](double count) -> bool { return count >= 1.0; },
    emp::vector<RecordedCommunitySummarizer::summary_update_fun_t>{
      KeepPresentWithInteractionPath
    }
  );
}

void AEcoWorld::AnalyzeWorldCommunities(
  const std::string& output_file_path,
  const RecordedCommunitySet<emp::vector<double>>& world_communities,
  const RecordedCommunitySet<emp::vector<double>>& assembly_communities,
  const RecordedCommunitySet<emp::vector<double>>& adaptive_communities
) {
  // For each community in the world:
  // - summary information (some repeated in snapshot)
  // - community
  // - found in adaptive
  // - found in assembly
  // - prop(adaptive) / prop(assembly)
  size_t community_id = 0;
  emp::DataFile community_summary_file(output_file_path);

  // proportion
  community_summary_file.AddFun<double>(
    [&community_id, &world_communities]() -> double {
      emp_assert(emp::Sum(world_communities.GetCommunityCounts()) > 0);
      return (double)world_communities.GetCommunityCount(community_id) / (double)emp::Sum(world_communities.GetCommunityCounts());
    },
    "proportion",
    "Proportion of cells where this particular community was found"
  );

  // num_present_species
  community_summary_file.AddFun<size_t>(
    [&community_id, &world_communities]() -> size_t {
      const auto& summary = world_communities.GetCommunitySummary(community_id);
      return summary.GetNumSpeciesPresent();
    },
    "num_present_species"
  );

  // num_possible_species
  community_summary_file.AddFun<size_t>(
    [&community_id, &world_communities]() -> size_t {
      const auto& summary = world_communities.GetCommunitySummary(community_id);
      return summary.counts.size();
    },
    "num_possible_species"
  );


  // species_counts
  community_summary_file.AddFun<std::string>(
    [&community_id, &world_communities]() -> std::string {
      const auto& summary = world_communities.GetCommunitySummary(community_id);
      return emp::to_string(summary.counts);
    },
    "species_counts"
  );

  // present_species_ids
  community_summary_file.AddFun<std::string>(
    [&community_id, &world_communities]() -> std::string {
      const auto& summary = world_communities.GetCommunitySummary(community_id);
      return emp::to_string(summary.present_species_ids);
    },
    "present_species_ids"
  );

  // present_species
  community_summary_file.AddFun<std::string>(
    [&community_id, &world_communities]() -> std::string {
      const auto& summary = world_communities.GetCommunitySummary(community_id);
      return emp::to_string(summary.present);
    },
    "present_species"
  );

  // found in assembly
  community_summary_file.AddFun<bool>(
    [&community_id, &world_communities, &assembly_communities]() -> bool {
      const auto& world_summary = world_communities.GetCommunitySummary(community_id);
      return assembly_communities.Has(world_summary);
    },
    "found_in_assembly"
  );

  community_summary_file.AddFun<bool>(
    [&community_id, &world_communities, &adaptive_communities]() -> bool {
      const auto& world_summary = world_communities.GetCommunitySummary(community_id);
      return adaptive_communities.Has(world_summary);
    },
    "found_in_adaptive"
  );

  community_summary_file.AddFun<double>(
    [&community_id, &world_communities, &adaptive_communities]() -> double {
      const auto& world_summary = world_communities.GetCommunitySummary(community_id);
      emp_assert(emp::Sum(adaptive_communities.GetCommunityCounts()) > 0);
      const auto adaptive_id = adaptive_communities.GetCommunityID(world_summary);
      return (adaptive_id) ? adaptive_communities.GetCommunityProportion(adaptive_id.value()) : 0.0;

    },
    "adaptive_proportion"
  );

  community_summary_file.AddFun<double>(
    [&community_id, &world_communities, &assembly_communities]() -> double {
      const auto& world_summary = world_communities.GetCommunitySummary(community_id);
      emp_assert(emp::Sum(assembly_communities.GetCommunityCounts()) > 0);
      const auto assembly_id = assembly_communities.GetCommunityID(world_summary);
      return (assembly_id) ? assembly_communities.GetCommunityProportion(assembly_id.value()) : 0.0;
    },
    "assembly_proportion"
  );

  community_summary_file.AddFun<std::string>(
    [&community_id, &world_communities, &assembly_communities, &adaptive_communities]() -> std::string {
      const auto& world_summary = world_communities.GetCommunitySummary(community_id);
      emp_assert(emp::Sum(assembly_communities.GetCommunityCounts()) > 0);
      emp_assert(emp::Sum(adaptive_communities.GetCommunityCounts()) > 0);
      const auto assembly_id = assembly_communities.GetCommunityID(world_summary);
      const auto adaptive_id = adaptive_communities.GetCommunityID(world_summary);
      const double assembly_prop = (assembly_id) ? assembly_communities.GetCommunityProportion(assembly_id.value()) : 0.0;
      const double adaptive_prop = (adaptive_id) ? adaptive_communities.GetCommunityProportion(adaptive_id.value()) : 0.0;
      return (assembly_prop != 0) ?
        emp::to_string(adaptive_prop / assembly_prop) :
        "ERR";
    },
    "adaptive_assembly_ratio"
  );

  community_summary_file.PrintHeaderKeys();
  for (community_id = 0; community_id < world_communities.GetSize(); ++community_id) {
    community_summary_file.Update();
  }


}

void AEcoWorld::SnapshotConfig() {
  emp::DataFile snapshot_file(output_dir + "run_config.csv");
  std::function<std::string(void)> get_param;
  std::function<std::string(void)> get_value;

  snapshot_file.AddFun<std::string>(
    [&get_param]() { return get_param(); },
    "parameter"
  );
  snapshot_file.AddFun<std::string>(
    [&get_value]() { return get_value(); },
    "value"
  );
  snapshot_file.PrintHeaderKeys();

  // Snapshot everything from config file
  for (const auto& entry : *config) {
    get_param = [&entry]() { return entry.first; };
    get_value = [&entry]() { return emp::to_string(entry.second->GetValue()); };
    snapshot_file.Update();
  }

  // Snapshot misc. other details
  emp::vector<std::pair<std::string, std::string>> misc_params = {
    std::make_pair("world_size", emp::to_string(world_size))
    // Can add more param-value pairs here that are not included in config object
  };

  for (const auto& entry : misc_params) {
    get_param = [&entry]() { return entry.first; };
    get_value = [&entry]() { return entry.second; };
    snapshot_file.Update();
  }
}

// TODO - does this need to live in AEcoWorld? Or can it be moved?
void AEcoWorld::SnapshotSubCommunities() {
  emp::DataFile subcommunities_file(output_dir + "subcommunities.csv");
  size_t cur_subcomm_id = 0;

  // subcommunity_id
  subcommunities_file.AddFun<size_t>(
    [&cur_subcomm_id, this]() -> size_t {
      return cur_subcomm_id;
    },
    "subcommunity_id"
  );

  // total_subcommunities
  subcommunities_file.AddFun<size_t>(
    [&cur_subcomm_id, this]() -> size_t {
      return community_structure.GetNumSubCommunities();
    },
    "total_subcommunities"
  );

  // total species
  subcommunities_file.AddFun<size_t>(
    [&cur_subcomm_id, this]() -> size_t {
      return community_structure.GetNumSpecies();
    },
    "total_species"
  );

  // num_subcommunity_species
  subcommunities_file.AddFun<size_t>(
    [&cur_subcomm_id, this]() -> size_t {
      return community_structure.GetNumMembers(cur_subcomm_id);
    },
    "num_subcommunity_species"
  );

  // prop_species_in_subcommunity
  subcommunities_file.AddFun<double>(
    [&cur_subcomm_id, this]() -> double {
      return community_structure.GetNumMembers(cur_subcomm_id) / community_structure.GetNumSpecies();
    },
    "prop_species_in_subcommunity"
  );

  // species composition
  subcommunities_file.AddFun<std::string>(
    [&cur_subcomm_id, this]() -> std::string {
      // Make copy (to sort for at a glance readability)
      auto subcomm = community_structure.GetSubCommunity(cur_subcomm_id);
      std::sort(
        subcomm.begin(),
        subcomm.end()
      );
      return emp::to_string(subcomm);
    },
    "species_composition",
    "Species IDs that makeup this subcommunity"
  );

  // species presence / absence
  subcommunities_file.AddFun<std::string>(
    [&cur_subcomm_id, this]() -> std::string {
      return emp::to_string(community_structure.GetFingerprints()[cur_subcomm_id]);
    },
    "species_presnt",
    "Binary vector indicating presence/absence"
  );

  // interaction_matrix
  subcommunities_file.AddFun<std::string>(
    [&cur_subcomm_id, this]() -> std::string {
      return config->INTERACTION_SOURCE();
    },
    "source_interaction_matrix"
  );

  subcommunities_file.PrintHeaderKeys();
  for (cur_subcomm_id = 0; cur_subcomm_id < community_structure.GetNumSubCommunities(); ++cur_subcomm_id) {
    subcommunities_file.Update();
  }
}


template<typename SUMMARY_SET_KEY_T>
void AEcoWorld::SnapshotRecordedCommunitySets(
  const std::string& file_path,
  const emp::vector<RecordedCommunitySetInfo<SUMMARY_SET_KEY_T>>& community_sets
) {
  emp::DataFile recorded_community_file(file_path);
  size_t cur_set_id = 0;
  size_t cur_summary_id = 0;

  // source
  recorded_community_file.AddFun<std::string>(
    [&cur_set_id, &cur_summary_id, &community_sets]() -> std::string {
      return community_sets[cur_set_id].source;
    },
    "source"
  );

  // proportion
  recorded_community_file.AddFun<double>(
    [&cur_set_id, &cur_summary_id, &community_sets]() -> double {
      const auto& community_set = community_sets[cur_set_id].summary_set;
      emp_assert(emp::Sum(community_set.GetCommunityCounts()) > 0);
      return (double)community_set.GetCommunityCount(cur_summary_id) / (double)emp::Sum(community_set.GetCommunityCounts());
    },
    "proportion",
    "Proportion of cells where this particular community was found"
  );

  // stabilized
  recorded_community_file.AddFun<bool>(
    [&cur_set_id, &cur_summary_id, &community_sets]() -> bool {
      return community_sets[cur_set_id].stabilized;
    },
    "stabilized"
  );

  // updates
  recorded_community_file.AddFun<size_t>(
    [&cur_set_id, &cur_summary_id, &community_sets]() -> size_t {
      return community_sets[cur_set_id].updates;
    },
    "updates"
  );

  // --- Summary information ---
  // num_present_species
  recorded_community_file.AddFun<size_t>(
    [&cur_set_id, &cur_summary_id, &community_sets, this]() -> size_t {
      const auto& summary = community_sets[cur_set_id].summary_set.GetCommunitySummary(cur_summary_id);
      return summary.GetNumSpeciesPresent();
    },
    "num_present_species"
  );

  recorded_community_file.AddFun<size_t>(
    [&cur_set_id, &cur_summary_id, &community_sets, this]() -> size_t {
      const auto& summary = community_sets[cur_set_id].summary_set.GetCommunitySummary(cur_summary_id);
      return summary.present_with_interaction_path.CountOnes();
    },
    "num_present_with_interaction_path"
  );

  recorded_community_file.AddFun<size_t>(
    [&cur_set_id, &cur_summary_id, &community_sets, this]() -> size_t {
      const auto& summary = community_sets[cur_set_id].summary_set.GetCommunitySummary(cur_summary_id);
      return summary.present_with_other_subcommunity_members.CountOnes();
    },
    "num_present_with_other_subcommunity_members"
  );

  // num_possible_species
  recorded_community_file.AddFun<size_t>(
    [&cur_set_id, &cur_summary_id, &community_sets, this]() -> size_t {
      const auto& summary = community_sets[cur_set_id].summary_set.GetCommunitySummary(cur_summary_id);
      return summary.counts.size();
    },
    "num_possible_species"
  );

  // counts
  recorded_community_file.AddFun<std::string>(
    [&cur_set_id, &cur_summary_id, &community_sets, this]() -> std::string {
      const auto& summary = community_sets[cur_set_id].summary_set.GetCommunitySummary(cur_summary_id);
      return emp::to_string(summary.counts);
    },
    "species_counts"
  );

  // present_species_ids
  recorded_community_file.AddFun<std::string>(
    [&cur_set_id, &cur_summary_id, &community_sets, this]() -> std::string {
      const auto& summary = community_sets[cur_set_id].summary_set.GetCommunitySummary(cur_summary_id);
      return emp::to_string(summary.present_species_ids);
    },
    "present_species_ids"
  );

  // present
  recorded_community_file.AddFun<std::string>(
    [&cur_set_id, &cur_summary_id, &community_sets, this]() -> std::string {
      const auto& summary = community_sets[cur_set_id].summary_set.GetCommunitySummary(cur_summary_id);
      return emp::to_string(summary.present);
    },
    "present_species"
  );

  // present_no_interactions
  recorded_community_file.AddFun<std::string>(
    [&cur_set_id, &cur_summary_id, &community_sets, this]() -> std::string {
      const auto& summary = community_sets[cur_set_id].summary_set.GetCommunitySummary(cur_summary_id);
      return emp::to_string(summary.present_with_interaction_path);
    },
    "species_present_with_interaction_path"
  );

  recorded_community_file.AddFun<std::string>(
    [&cur_set_id, &cur_summary_id, &community_sets, this]() -> std::string {
      const auto& summary = community_sets[cur_set_id].summary_set.GetCommunitySummary(cur_summary_id);
      return emp::to_string(summary.present_with_other_subcommunity_members);
    },
    "species_present_with_other_subcommunity_members"
  );

  recorded_community_file.AddFun<size_t>(
    [&cur_set_id, &cur_summary_id, &community_sets, this]() -> size_t {
      const auto& summary = community_sets[cur_set_id].summary_set.GetCommunitySummary(cur_summary_id);
      return summary.complete_subcommunities_present.size();
    },
    "num_complete_subcommunities_present"
  );

  recorded_community_file.AddFun<size_t>(
    [&cur_set_id, &cur_summary_id, &community_sets, this]() -> size_t {
      const auto& summary = community_sets[cur_set_id].summary_set.GetCommunitySummary(cur_summary_id);
      return summary.partial_subcommunities_present.size();
    },
    "num_partial_subcommunities_present"
  );

  // complete_subcommunities_present
  recorded_community_file.AddFun<std::string>(
    [&cur_set_id, &cur_summary_id, &community_sets, this]() -> std::string {
      const auto& summary = community_sets[cur_set_id].summary_set.GetCommunitySummary(cur_summary_id);
      return emp::to_string(summary.complete_subcommunities_present);
    },
    "complete_subcommunities_present"
  );

  // partial_subcommunities_present
  recorded_community_file.AddFun<std::string>(
    [&cur_set_id, &cur_summary_id, &community_sets, this]() -> std::string {
      const auto& summary = community_sets[cur_set_id].summary_set.GetCommunitySummary(cur_summary_id);
      return emp::to_string(summary.partial_subcommunities_present);
    },
    "partial_subcommunities_present"
  );

  // proportion_subcommunity_present
  recorded_community_file.AddFun<std::string>(
    [&cur_set_id, &cur_summary_id, &community_sets, this]() -> std::string {
      const auto& summary = community_sets[cur_set_id].summary_set.GetCommunitySummary(cur_summary_id);
      return emp::to_string(summary.proportion_subcommunity_present);
    },
    "proportion_subcommunity_present"
  );

  recorded_community_file.PrintHeaderKeys();
  // Write one line per summary, per community set
  for (cur_set_id = 0; cur_set_id < community_sets.size(); ++cur_set_id) {
    const auto& summary_set = community_sets[cur_set_id].summary_set;
    for (cur_summary_id = 0; cur_summary_id < summary_set.GetSize(); ++cur_summary_id) {
      recorded_community_file.Update();
    }
  }

}

} // End chemical_ecology namespace