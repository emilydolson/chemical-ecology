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

/// Class that helps to manage world community summary files
class WorldCommunitySummaryFile {
public:
  using community_set_t = RecordedCommunitySet<emp::vector<double>>;
protected:
  emp::DataFile summary_file;
  size_t community_id = 0;
  size_t cur_update = 0;
  // Unowned pointers
  emp::Ptr<const community_set_t> cur_world_communities;
  emp::Ptr<const community_set_t> cur_assembly_communities;
  emp::Ptr<const community_set_t> cur_adaptive_communities;

  void Setup() {
    // current update
    summary_file.AddVar(cur_update, "update", "Model update");

    // proportion
    summary_file.AddFun<double>(
      [this]() -> double {
        emp_assert(emp::Sum(cur_world_communities->GetCommunityCounts()) > 0);
        return (double)cur_world_communities->GetCommunityCount(community_id) / (double)emp::Sum(cur_world_communities->GetCommunityCounts());
      },
      "proportion",
      "Proportion of cells where this particular community was found"
    );

    // num_present_species
    summary_file.AddFun<size_t>(
      [this]() -> size_t {
        const auto& summary = cur_world_communities->GetCommunitySummary(community_id);
        return summary.GetNumSpeciesPresent();
      },
      "num_present_species"
    );

    // num_possible_species
    summary_file.AddFun<size_t>(
      [this]() -> size_t {
        const auto& summary = cur_world_communities->GetCommunitySummary(community_id);
        return summary.counts.size();
      },
      "num_possible_species"
    );

    // species_counts
    summary_file.AddFun<std::string>(
      [this]() -> std::string {
        const auto& summary = cur_world_communities->GetCommunitySummary(community_id);
        return emp::to_string(summary.counts);
      },
      "species_counts"
    );

    // present_species_ids
    summary_file.AddFun<std::string>(
      [this]() -> std::string {
        const auto& summary = cur_world_communities->GetCommunitySummary(community_id);
        return emp::to_string(summary.present_species_ids);
      },
      "present_species_ids"
    );

    // present_species
    summary_file.AddFun<std::string>(
      [this]() -> std::string {
        const auto& summary = cur_world_communities->GetCommunitySummary(community_id);
        return emp::to_string(summary.present);
      },
      "present_species"
    );

    // found in assembly
    summary_file.AddFun<bool>(
      [this]() -> bool {
        const auto& world_summary = cur_world_communities->GetCommunitySummary(community_id);
        return cur_assembly_communities->Has(world_summary);
      },
      "found_in_assembly"
    );

    summary_file.AddFun<bool>(
      [this]() -> bool {
        const auto& world_summary = cur_world_communities->GetCommunitySummary(community_id);
        return cur_adaptive_communities->Has(world_summary);
      },
      "found_in_adaptive"
    );

    summary_file.AddFun<double>(
      [this]() -> double {
        const auto& world_summary = cur_world_communities->GetCommunitySummary(community_id);
        emp_assert(emp::Sum(cur_adaptive_communities->GetCommunityCounts()) > 0);
        const auto adaptive_id = cur_adaptive_communities->GetCommunityID(world_summary);
        return (adaptive_id) ? cur_adaptive_communities->GetCommunityProportion(adaptive_id.value()) : 0.0;
      },
      "adaptive_proportion"
    );

    summary_file.AddFun<double>(
      [this]() -> double {
        const auto& world_summary = cur_world_communities->GetCommunitySummary(community_id);
        emp_assert(emp::Sum(cur_assembly_communities->GetCommunityCounts()) > 0);
        const auto assembly_id = cur_assembly_communities->GetCommunityID(world_summary);
        return (assembly_id) ?
          cur_assembly_communities->GetCommunityProportion(assembly_id.value()) :
          0.0;
      },
      "assembly_proportion"
    );

    summary_file.AddFun<std::string>(
      [this]() -> std::string {
        const auto& world_summary = cur_world_communities->GetCommunitySummary(community_id);
        emp_assert(emp::Sum(cur_assembly_communities->GetCommunityCounts()) > 0);
        emp_assert(emp::Sum(cur_adaptive_communities->GetCommunityCounts()) > 0);
        const auto assembly_id = cur_assembly_communities->GetCommunityID(world_summary);
        const auto adaptive_id = cur_adaptive_communities->GetCommunityID(world_summary);
        const double assembly_prop = (assembly_id) ?
          cur_assembly_communities->GetCommunityProportion(assembly_id.value()) :
          0.0;
        const double adaptive_prop = (adaptive_id) ?
          cur_adaptive_communities->GetCommunityProportion(adaptive_id.value()) :
          0.0;
        return (assembly_prop != 0) ?
          emp::to_string(adaptive_prop / assembly_prop) :
          "ERR";
      },
      "adaptive_assembly_ratio"
    );

    summary_file.PrintHeaderKeys();
  }

public:
  WorldCommunitySummaryFile(const std::string& file_path) :
    summary_file(file_path)
  {
    Setup();
  }

  void Update(
    size_t update,
    const community_set_t& world_communities,
    const community_set_t& assembly_communities,
    const community_set_t& adaptive_communities
  ) {
    // Capture given community sets in member variables (giving data file functions access)
    cur_world_communities = &world_communities;
    cur_assembly_communities = &assembly_communities;
    cur_adaptive_communities = &adaptive_communities;
    cur_update = update;

    // Update summary file for every recorded world community
    for (community_id = 0; community_id < world_communities.GetSize(); ++community_id) {
      summary_file.Update();
    }

    // Null out unowned pointers
    cur_world_communities = nullptr;
    cur_assembly_communities = nullptr;
    cur_adaptive_communities = nullptr;
  }
};

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

  size_t world_update;
  size_t analysis_update; // Update inside of "analysis"
  size_t stochastic_rep;

  // Initialize vector that keeps track of grid
  world_t world;

  // Manages community structure (determined by interaction matrix)
  CommunityStructure community_structure;
  emp::vector<size_t> subcommunity_group_repro_schedule;

  // Used to track activation order of positions in the world
  emp::vector<size_t> position_activation_order;

  emp::Ptr<RecordedCommunitySummarizer> community_summarizer_raw;         // Uses raw species counts.
  emp::Ptr<RecordedCommunitySummarizer> community_summarizer_pwip;        // Keeps all species present with valid interaction paths to other present species

  RecordedCommunitySet<emp::vector<double>>::summary_key_fun_t recorded_comm_key_fun;
  emp::Ptr<RecordedCommunitySet<emp::vector<double>>> recorded_communities_assembly_raw;
  emp::Ptr<RecordedCommunitySet<emp::vector<double>>> recorded_communities_adaptive_raw;
  emp::Ptr<RecordedCommunitySet<emp::vector<double>>> recorded_communities_assembly_pwip;
  emp::Ptr<RecordedCommunitySet<emp::vector<double>>> recorded_communities_adaptive_pwip;

  // Set up data tracking
  std::string output_dir;
  emp::Ptr<emp::DataFile> data_file = nullptr;
  emp::Ptr<emp::DataFile> assembly_data_file = nullptr;
  emp::Ptr<emp::DataFile> adaptive_data_file = nullptr;

  emp::Ptr<WorldCommunitySummaryFile> world_community_summary_pwip_file = nullptr; // Summarizes results from world community analysis

  world_t worldState;
  world_t assemblyWorldState;
  world_t adaptiveWorldState;

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
    bool output_snapshots = false
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
    if (assembly_data_file != nullptr) assembly_data_file.Delete();
    if (adaptive_data_file != nullptr) adaptive_data_file.Delete();
    if (community_summarizer_raw != nullptr) community_summarizer_raw.Delete();
    if (community_summarizer_pwip != nullptr) community_summarizer_pwip.Delete();
    if (world_community_summary_pwip_file != nullptr) world_community_summary_pwip_file.Delete();

    if (recorded_communities_assembly_raw != nullptr) recorded_communities_assembly_raw.Delete();
    if (recorded_communities_adaptive_raw != nullptr) recorded_communities_adaptive_raw.Delete();
    if (recorded_communities_assembly_pwip != nullptr) recorded_communities_assembly_pwip.Delete();
    if (recorded_communities_adaptive_pwip != nullptr) recorded_communities_adaptive_pwip.Delete();
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
    data_file->AddVar(world_update, "Time", "Time");
    data_file->AddFun<std::string>(
      [this]() -> std::string { return emp::to_string(worldState); },
      "worldState",
      "world state"
    );
    data_file->SetTimingRepeat(config->OUTPUT_RESOLUTION());
    data_file->PrintHeaderKeys();

    // Assembly file
    assembly_data_file = emp::NewPtr<emp::DataFile>(output_dir + "a-eco_assembly_model_data.csv");
    assembly_data_file->AddVar(stochastic_rep, "replicate", "Replicate of model");
    assembly_data_file->AddVar(analysis_update, "Time", "Time");
    assembly_data_file->AddFun<std::string>(
      [this]() -> std::string { return emp::to_string(assemblyWorldState); },
      "assemblyWorldState",
      "assembly world state"
    );
    assembly_data_file->SetTimingRepeat(config->OUTPUT_RESOLUTION());
    assembly_data_file->PrintHeaderKeys();

    // Adaptive file
    adaptive_data_file = emp::NewPtr<emp::DataFile>(output_dir + "a-eco_adaptive_model_data.csv");
    adaptive_data_file->AddVar(stochastic_rep, "replicate", "Replicate of model");
    adaptive_data_file->AddVar(analysis_update, "Time", "Time");
    adaptive_data_file->AddFun<std::string>(
      [this]() -> std::string { return emp::to_string(adaptiveWorldState); },
      "adaptiveWorldState",
      "adaptive world state"
    );
    adaptive_data_file->SetTimingRepeat(config->OUTPUT_RESOLUTION());
    adaptive_data_file->PrintHeaderKeys();

    // Setup world summary file
    // NOTE (@AML): This should happen wherever we decide to configure the set of
    //              community summary methods. We might want different summary files
    //              for each summary method.
    world_community_summary_pwip_file = emp::NewPtr<WorldCommunitySummaryFile>(
      output_dir + "world_summary_pwip.csv"
    );

    // Output a snapshot of identified subcommunities
    SnapshotSubCommunities();

    // Output a snapshot of the run configuration
    SnapshotConfig();
  }

  // Handle the process of running the program through
  // all time steps
  void Run() {

    // Run N replicates of the adaptive model and assembly model.
    // Save recorded summaries to use when summarizing the world.
    for (size_t i = 0; i < config->STOCHASTIC_ANALYSIS_REPS(); ++i) {
      // Adaptive and assembly data tracking
      stochastic_rep = i;

      // Run assembly model
      world_t assemblyModel = AssemblyModel(config->UPDATES(), config->PROB_CLEAR(), config->SEEDING_PROB());
      world_t stableAssemblyModel = GenStabilizedWorld(assemblyModel, config->CELL_STABILIZATION_UPDATES());

      // Run adaptive model
      world_t adaptiveModel = AdaptiveModel(config->UPDATES(), config->PROB_CLEAR(), config->SEEDING_PROB());
      world_t stableAdaptiveModel = GenStabilizedWorld(adaptiveModel, config->CELL_STABILIZATION_UPDATES());

      // Add summarized recorded communities to sets
      recorded_communities_assembly_raw->Add(
        community_summarizer_raw->SummarizeAll(stableAssemblyModel)
      );
      recorded_communities_assembly_pwip->Add(
        community_summarizer_pwip->SummarizeAll(stableAssemblyModel)
      );
      recorded_communities_adaptive_raw->Add(
        community_summarizer_raw->SummarizeAll(stableAdaptiveModel)
      );
      recorded_communities_adaptive_pwip->Add(
        community_summarizer_pwip->SummarizeAll(stableAdaptiveModel)
      );
    }

    // Call update the specified number of times
    for (world_update = 0; world_update <= config->UPDATES(); ++world_update) {
      Update();
    }

    //Print out final state if in verbose mode
    if (config->V()) {
      std::cout << "World Vectors:" << std::endl;
      for (const auto& v : world) {
        std::cout << emp::to_string(v) << std::endl;
      }
    }
  }

  // Handle an individual time step
  // ud = which time step we're on
  void Update() {
    // Create a new world object to store the values
    // for the next time step
    world_t next_world(
      world_size,
      emp::vector<double>(N_TYPES, 0)
    );

    // Handle population growth for each cell
    for (size_t pos = 0; pos < world.size(); pos++) {
      DoGrowth(pos, world, next_world);
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

    // Give data_file the opportunity to write to the file
    // TODO - can we eliminate worldState variable?
    // E.g., if we move after swap, world should be up-to-date
    worldState = next_world;
    data_file->Update(world_update);

    // We're done calculating the type counts for the next
    // time step. We can now swap our counts for the next
    // time step into the main world variable
    std::swap(world, next_world);

    // NOTE - world member variable should be accurate for end of update
    // On final update if next update would equal or exceed updates param
    const bool is_final_update = world_update >= config->UPDATES();
    const bool world_comm_analysis_interval = ((world_update % config->OUTPUT_RESOLUTION()) == 0) || is_final_update;
    if (world_comm_analysis_interval) {
      AnalyzeWorldCommunities(
        /*output_snapshots = */ is_final_update
      );
    }

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
  world_t GenStabilizedWorld(const world_t& custom_world, size_t max_updates=10000) {

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

  world_t AssemblyModel(
    int num_updates,
    double prob_clear,
    double seeding_prob
  ) {
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
        // (1) clearing
        DoClearing(pos, model_world, next_model_world, prob_clear);
        // (2) seeding
        DoSeeding(pos, model_world, next_model_world, seeding_prob);
      }

      // Record world state
      const bool final_update = (i == num_updates);
      const bool res_update =  !(bool)(i % config->OUTPUT_RESOLUTION());
      if (config->RECORD_ASSEMBLY_MODEL() && (final_update || res_update)) {
        analysis_update = i;
        assemblyWorldState = next_model_world;
        assembly_data_file->Update();
      }

      std::swap(model_world, next_model_world);

      if (config->RECORD_ASSEMBLY_MODEL()) {
        assemblyWorldState.clear();
      }

    }

    return model_world;
  }

  world_t AdaptiveModel(
    int num_updates,
    double prob_clear,
    double seeding_prob
  ) {
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
        DoGroupRepro(pos, model_world, next_model_world);
        // (2) clearing
        DoClearing(pos, model_world, next_model_world, prob_clear);
        // (3) seeding
        DoSeeding(pos, model_world, next_model_world, seeding_prob);
      }

      // Record world state
      const bool final_update = (i == num_updates);
      const bool res_update =  !(bool)(i % config->OUTPUT_RESOLUTION());
      if (config->RECORD_ADAPTIVE_MODEL() && (final_update || res_update)) {
        analysis_update = i;
        adaptiveWorldState = next_model_world;
        adaptive_data_file->Update();
      }

      std::swap(model_world, next_model_world);

      if (config->RECORD_ADAPTIVE_MODEL()) {
        adaptiveWorldState.clear();
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
  int GetWorldTime() {
    return world_update;
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
  emp_assert(recorded_communities_assembly_raw == nullptr);
  emp_assert(recorded_communities_adaptive_raw == nullptr);
  emp_assert(recorded_communities_assembly_pwip == nullptr);
  emp_assert(recorded_communities_adaptive_pwip == nullptr);

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

  // Configure recorded community sets for adaptive / assembly models
  recorded_comm_key_fun = [](
    const RecordedCommunitySummary& summary
  ) -> const auto& {
    return summary.counts;
  };

  recorded_communities_assembly_raw = emp::NewPtr<RecordedCommunitySet<emp::vector<double>>>(recorded_comm_key_fun);
  recorded_communities_adaptive_raw = emp::NewPtr<RecordedCommunitySet<emp::vector<double>>>(recorded_comm_key_fun);
  recorded_communities_assembly_pwip = emp::NewPtr<RecordedCommunitySet<emp::vector<double>>>(recorded_comm_key_fun);
  recorded_communities_adaptive_pwip = emp::NewPtr<RecordedCommunitySet<emp::vector<double>>>(recorded_comm_key_fun);
}

void AEcoWorld::AnalyzeWorldCommunities(
  bool output_snapshots
) {
    // Run world forward without diffusion
    world_t stable_world(
      GenStabilizedWorld(world, config->CELL_STABILIZATION_UPDATES())
    );

    // NOTE (@AML): Not in total love with this; could possibly use another iteration after chatting
    //   about future functionality that would be useful to have
    // - Should all be managed together (to make it easier to add different summary strategies)
    // RecordedCommunitySet<emp::vector<double>>::summary_key_fun_t recorded_comm_key_fun = [](
    //   const RecordedCommunitySummary& summary
    // ) -> const auto& {
    //   return summary.counts;
    // };

    // Set of recorded communities with "raw" counts
    RecordedCommunitySet<emp::vector<double>> recorded_communities_world_raw(recorded_comm_key_fun);
    // RecordedCommunitySet<emp::vector<double>> recorded_communities_assembly_raw(recorded_comm_key_fun);
    // RecordedCommunitySet<emp::vector<double>> recorded_communities_adaptive_raw(recorded_comm_key_fun);

    // Set of recorded communities with counts that have had PNI-species zeroed-out
    RecordedCommunitySet<emp::vector<double>> recorded_communities_world_pwip(recorded_comm_key_fun);
    // RecordedCommunitySet<emp::vector<double>> recorded_communities_assembly_pwip(recorded_comm_key_fun);
    // RecordedCommunitySet<emp::vector<double>> recorded_communities_adaptive_pwip(recorded_comm_key_fun);

    // Add world summaries
    recorded_communities_world_raw.Add(
      community_summarizer_raw->SummarizeAll(stable_world)
    );
    recorded_communities_world_pwip.Add(
      community_summarizer_pwip->SummarizeAll(stable_world)
    );

    // Update world community file
    world_community_summary_pwip_file->Update(
      world_update,
      recorded_communities_world_pwip,
      *recorded_communities_assembly_pwip,
      *recorded_communities_adaptive_pwip
    );

    if (output_snapshots) {
      // NOTE (@AML): Slightly clunky way to tie together things
      // Snapshot raw recorded community summaries
      SnapshotRecordedCommunitySets</*SUMMARY_SET_KEY_T=*/emp::vector<double>>(
        output_dir + "recorded_communities_raw_" + emp::to_string(world_update) + ".csv",
        {
          {recorded_communities_world_raw, "world", true, config->UPDATES()},
          {*recorded_communities_assembly_raw, "assembly", true, config->UPDATES()},
          {*recorded_communities_adaptive_raw, "adaptive", true, config->UPDATES()}
        }
      );

      // Snapshot summaries where "present-no-interactions" species have been removed
      SnapshotRecordedCommunitySets</*SUMMARY_SET_KEY_T=*/emp::vector<double>>(
        output_dir + "recorded_communities_pwip_" + emp::to_string(world_update) + ".csv",
        {
          {recorded_communities_world_pwip, "world", true, config->UPDATES()},
          {*recorded_communities_assembly_pwip, "assembly", true, config->UPDATES()},
          {*recorded_communities_adaptive_pwip, "adaptive", true, config->UPDATES()}
        }
      );
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