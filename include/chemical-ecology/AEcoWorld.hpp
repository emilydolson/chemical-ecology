#pragma once

#include <iostream>
#include <string>
#include <cmath>
#include <sys/stat.h>
#include <algorithm>

#include "emp/Evolve/World.hpp"
#include "emp/math/distances.hpp"
#include "emp/math/random_utils.hpp"
#include "emp/datastructs/Graph.hpp"
#include "emp/bits/BitArray.hpp"
#include "emp/data/DataNode.hpp"
#include "emp/base/Ptr.hpp"
#include "emp/datastructs/map_utils.hpp"

#include "emp/datastructs/map_utils.hpp"
#include "chemical-ecology/config_setup.hpp"
#include "chemical-ecology/SpatialStructure.hpp"
#include "chemical-ecology/utils/graph_utils.hpp"

namespace chemical_ecology {

// TERMINOLOGY NOTES:
//
// - Cell = location in the world grid
// - Type = a type of thing that can live in a cell
//          (e.g. a chemical or a species). Usually
//          called a species in practice
// - Community = a set of types that interact and can
//                coexist in the world

// Extracts and manages community structure from given interaction matrix
// TODO - write tests
class CommunityStructure {
public:
  using interaction_mat_t = emp::vector<emp::vector<double>>;
protected:

  emp::vector<emp::vector<size_t>> subcommunities;
  emp::vector<emp::BitVector> subcommunity_fingerprints; // Each bit vector represents presence / absence of species
  emp::vector<size_t> species_to_subcommunity_id;
  size_t num_species;

  emp::vector<emp::BitVector> species_shares_subcommunity_with; // For each species, what other species does it share a subcommunity with
  emp::vector<emp::BitVector> species_interacts_with;    // For each species, what other species does it interact with?

public:
  CommunityStructure() = default;
  CommunityStructure(
    const interaction_mat_t& interaction_matrix,
    std::function<bool(
      const interaction_mat_t&,
      size_t,
      size_t
    )> interacts
  ) {
    SetStructure(interaction_matrix, interacts);
  }

  // Configure subcommunity structure given an interaction matrix
  // Identifies the set of isolated subcommunities present in the given interaction matrix.
  void SetStructure(
    const interaction_mat_t& interaction_matrix,
    std::function<bool(
      const interaction_mat_t&,
      size_t,
      size_t
    )> interacts
  ) {
    Clear();
    num_species = interaction_matrix.size();

    // Identify all subcommunities (connected components) in interaction matrix
    subcommunities = utils::FindConnectedComponents<double>(
      interaction_matrix,
      interacts
    );

    // Extract subcommunity presence/absence fingerprints and
    //   map each species to its subcommunity id
    subcommunity_fingerprints.resize(subcommunities.size(), {num_species, false});
    species_to_subcommunity_id.resize(num_species, (size_t)-1);
    for (size_t comm_i = 0; comm_i < subcommunities.size(); ++comm_i) {
      const auto& sub_comm = subcommunities[comm_i];
      for (size_t species : sub_comm) {
        subcommunity_fingerprints[comm_i][species] = true;
        // NOTE: Each species should be part of *one* community.
        emp_assert(species_to_subcommunity_id[species] == (size_t)-1);
        species_to_subcommunity_id[species] = comm_i;
      }
    }

    // For each species, cache what set of species it has interactions with
    //   Extract directly from interaction matrix to more easily capture self-interactions
    //   which would be indistinguishable if just looking at the subcommunity membership
    // QUESTION: if something interacts with itself but nothing else, is it present without interactions?
    // Additionally, cache what set of species a species shares a subcommunity with
    species_interacts_with.resize(num_species, {num_species, false});
    species_shares_subcommunity_with.resize(num_species, {num_species, false});
    for (size_t spec_id = 0; spec_id < num_species; ++spec_id) {
      // Record interactions
      auto& interacts_with = species_interacts_with[spec_id];
      for (size_t other_id = 0; other_id < num_species; ++other_id) {
        interacts_with[other_id] = interacts(interaction_matrix, spec_id, other_id);
      }
      // Record subcommunity co-members
      auto& shares_subcommunity = species_shares_subcommunity_with[spec_id];
      const auto& sub_comm = subcommunities[GetSubCommunityID(spec_id)];
      for (size_t other_id : sub_comm) {
        // QUESTION: should a species be counted as sharing a subcommunity with itself?
        shares_subcommunity[other_id] = true;
      }
    }
  }

  // Get number of species represented in community structure
  size_t GetNumSpecies() const { return num_species; }

  // Get number of isolated subcommunities present in community structure
  size_t GetNumSubCommunities() const { return subcommunities.size(); }

  bool SpeciesInteractsWith(size_t focal_species_id, size_t interacts_with_id) const {
    return species_interacts_with[focal_species_id][interacts_with_id];
  }

  bool SpeciesSharesSubCommunityWith(size_t focal_species_id, size_t shares_with_id) const {
    return species_shares_subcommunity_with[focal_species_id][shares_with_id];
  }

  // Get list of all isolated subcommunities
  const emp::vector<emp::vector<size_t>>& GetSubCommunities() const { return subcommunities; }

  // Get a particular isolated subcommunity
  const emp::vector<size_t>& GetSubCommunity(size_t community_id) const {
    emp_assert(community_id < GetNumSubCommunities());
    return subcommunities[community_id];
  }

  // Get presence/absence fingerpring for all isolated subcommunities
  const emp::vector<emp::BitVector>& GetFingerprints() const { return subcommunity_fingerprints; }

  // Get the subcommunity ID to which the given species belongs
  size_t GetSubCommunityID(size_t species_id) const {
    emp_assert(species_id < species_to_subcommunity_id.size());
    return species_to_subcommunity_id[species_id];
  }

  // Clear all stored structure
  void Clear() {
    subcommunities.clear();
    subcommunity_fingerprints.clear();
    species_to_subcommunity_id.clear();
    num_species = 0;
    species_interacts_with.clear();
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
  using interaction_mat_t = emp::vector< emp::vector<double> >;
  using config_t = Config;

  // Community summary information (extracted from cell of a world)
  // NOTE: depending on how sophisticated this ends up, should get moved to own file & promoted to a proper class
  struct CommunityInfo {
    emp::vector<size_t> present_species_ids;
    emp::BitVector present; // species present/absent fingerprint
    emp::BitVector present_no_interactions; // species present without interactions with *OTHER* species
    emp::vector<size_t> counts; // species counts (uniquely identifies this community)

    void SetCommunity(
      const emp::vector<double>& member_counts,
      const CommunityStructure& community_structure
    ) {
      const size_t num_members = member_counts.size();
      // Reset current community info
      Reset(num_members);
      emp_assert(counts.size() == member_counts.size());

      // Process counts
      for (size_t mem_i = 0; mem_i < num_members; ++mem_i) {
        // TODO: change this to configurable functor?
        // Drops fractional component of number
        counts[mem_i] = std::trunc(member_counts[mem_i]);
        // Fingerprint present/absence
        // TODO: change this to a configurable functor?
        present[mem_i] = counts[mem_i] >= 1;
        if (present[mem_i]) {
          present_species_ids.emplace_back(mem_i);
        }
      }

      // Identify members that are present with no interactions
      for (size_t mem_i : present_species_ids) {
        emp_assert(present[mem_i]);
        const size_t member_comm_id = community_structure.GetSubCommunityID(mem_i);
        const auto& subcommunity = community_structure.GetSubCommunity(member_comm_id);
        emp_assert(emp::Has(subcommunity, mem_i));

        // Does this member species have other species present that share a community?
        // - For each other species in this species' subcommunity, are any present?
        bool interacts = false;
        for (size_t other_id : subcommunity) {
          if (other_id != mem_i && present[other_id]) {
            interacts = true;
            break;
          }
        }
        present_no_interactions[mem_i] = !interacts;
      }
    }

    void Reset(size_t num_members=0) {
      (present.Resize(num_members)).Clear();
      (present_no_interactions.Resize(num_members)).Clear();
      present_species_ids.clear();
      counts.clear();
      counts.resize(num_members, 0);
    }

    // Operator< necessary for std::map
    // Uses lexicographic ordering (so not necessarily meaningful in context of numeric comparisons)
    // NOTE: O(N) operation, where 9 is size of counts
    bool operator<(const CommunityInfo& in) const {
      emp_assert(in.counts.size() == counts.size(), "Expected communities to have same number of potential species");
      return counts < in.counts;
    }

  };

private:


  // The matrix of interactions between types
  // The diagonal of this matrix represents the
  // intrinsic growth rate (r) of each type
  interaction_mat_t interactions;

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
  int curr_update2;

  // Initialize vector that keeps track of grid
  world_t world;

  // Manages community structure (determined by interaction matrix)
  CommunityStructure community_structure;
  emp::vector<size_t> subcommunity_group_repro_schedule;

  // Used to track activation order of positions in the world
  emp::vector<size_t> position_activation_order;

  // Set up data tracking
  std::string output_dir;
  emp::Ptr<emp::DataFile> data_file = nullptr;
  emp::Ptr<emp::DataFile> stochastic_data_file = nullptr;

  emp::DataNode<double, emp::data::Stats> biomass_node;
  emp::vector<double> fittest;
  emp::vector<double> dominant;
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

  // Output a snapshot of how the world is configured
  void SnapshotConfig();

public:

  AEcoWorld() = default;

  ~AEcoWorld() {
    if (data_file != nullptr) data_file.Delete();
    if (stochastic_data_file != nullptr) stochastic_data_file.Delete();
  }

  // Setup the world according to the specified configuration
  void Setup(config_t& cfg) {

    // Store cfg for future reference
    config = &cfg;

    // Set local config variables based on given configuration
    N_TYPES = config->N_TYPES();
    MAX_POP = double(config->MAX_POP());

    // Set seed to configured value for reproducibility
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
    if (config->INTERACTION_SOURCE() == "") {
      SetupRandomInteractions();
    } else {
      LoadInteractionMatrix(config->INTERACTION_SOURCE());
    }

    // Configure output directory path, create directory
    output_dir = config->OUTPUT_DIR();
    mkdir(output_dir.c_str(), ACCESSPERMS);
    if(output_dir.back() != '/') {
        output_dir += '/';
    }

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

    // Make sure you get sub communities after setting up the matrix
    community_structure.SetStructure(
      interactions,
      [](const interaction_mat_t& matrix, size_t from, size_t to) -> bool {
        return (matrix[from][to] != 0) || (matrix[to][from] != 0);
      }
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

    // Output a snapshot of the run configuration
    SnapshotConfig();
  }

  // Handle the process of running the program through
  // all time steps
  void Run() {

    // If there are any sub communities, we should know about them
    if (community_structure.GetNumSubCommunities() != 1) {
      std::cout << "Multiple sub-communities detected!: " << "\n";
      for (const auto& community : community_structure.GetSubCommunities()) {
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

    // Run world forward without diffusion
    world_t stable_world(stableUpdate(world));
    // Identify final communities in world
    std::map<std::string, double> finalCommunities = getFinalCommunities(stable_world);
    std::map<emp::BitVector, double> world_community_fingerprints = IdentifyWorldCommunities(stable_world);
    std::map<CommunityInfo, double> world_community_props;
    world_community_props[CommunityInfo()] = 0.5;

    std::map<std::string, double> assemblyFinalCommunities;
    std::map<std::string, double> adaptiveFinalCommunities;

    emp::vector<std::map<emp::BitVector, double>> assembly_fingerprints(config->STOCHASTIC_ANALYSIS_REPS());
    emp::vector<std::map<emp::BitVector, double>> adaptive_fingerprints(config->STOCHASTIC_ANALYSIS_REPS());

    // Run n stochastic worlds
    for (size_t i = 0; i < config->STOCHASTIC_ANALYSIS_REPS(); i++) {
      const bool record_analysis_state = (i == 0);
      world_t assemblyModel = StochasticModel(config->UPDATES(), false, config->PROB_CLEAR(), config->SEEDING_PROB(), record_analysis_state);
      world_t stableAssemblyModel = stableUpdate(assemblyModel);
      world_t adaptiveModel = StochasticModel(config->UPDATES(), true, config->PROB_CLEAR(), config->SEEDING_PROB(), record_analysis_state);
      world_t stableAdaptiveModel = stableUpdate(adaptiveModel);

      std::map<std::string, double> assemblyCommunities = getFinalCommunities(stableAssemblyModel);
      assembly_fingerprints[i] = IdentifyWorldCommunities(stableAssemblyModel);

      std::map<std::string, double> adaptiveCommunities = getFinalCommunities(stableAdaptiveModel);
      adaptive_fingerprints[i] = IdentifyWorldCommunities(stableAdaptiveModel);

      for (auto& [node, proportion] : assemblyCommunities) {
        if (assemblyFinalCommunities.find(node) == assemblyFinalCommunities.end()) {
          assemblyFinalCommunities.insert({node, proportion});
        } else {
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
    for (auto & comm : assemblyFinalCommunities) comm.second = comm.second/double(config->STOCHASTIC_ANALYSIS_REPS());
    for (auto & comm : adaptiveFinalCommunities) comm.second = comm.second/double(config->STOCHASTIC_ANALYSIS_REPS());

    // Average over analysis replicates (for assembly community fingerprint proportions)
    std::map<emp::BitVector, double> assembly_community_fingerprints_final;
    for (const auto& communities : assembly_fingerprints) {
      for (const auto& comm_prop : communities) {
        const auto& community = comm_prop.first;
        const double proportion = comm_prop.second;
        if (!emp::Has(assembly_community_fingerprints_final, community)) {
          assembly_community_fingerprints_final[community] = 0.0;
        }
        assembly_community_fingerprints_final[community] += proportion / (double)config->STOCHASTIC_ANALYSIS_REPS();
      }
    }

    // Average over analysis replicates (for adaptive community fingerprint proportions)
    std::map<emp::BitVector, double> adaptive_community_fingerprints_final;
    for (const auto& communities : adaptive_fingerprints) {
      for (const auto& comm_prop : communities) {
        const auto& community = comm_prop.first;
        const double proportion = comm_prop.second;
        if (!emp::Has(adaptive_community_fingerprints_final, community)) {
          adaptive_community_fingerprints_final[community] = 0.0;
        }
        adaptive_community_fingerprints_final[community] += proportion / (double)config->STOCHASTIC_ANALYSIS_REPS();
      }
    }


    // TODO - output this information in a datafile!
    std::cout << "World" << std::endl;
    for (auto& [node, proportion] : finalCommunities) {
      std::cout << "  Community: " << node << " Proportion: " << proportion << std::endl;
    }

    std::cout << "World fingerprints" << std::endl;
    for (auto& [node, proportion] : world_community_fingerprints) {
      std::cout << "  Community: " << node << " Proportion: " << proportion << std::endl;
    }

    std::cout << "Assembly" << std::endl;
    for (auto& [node, proportion] : assemblyFinalCommunities) {
      std::cout << "  Community: " << node << " Proportion: " << proportion << std::endl;
    }
    std::cout << "Assembly fingerprints" << std::endl;
    for (auto& [node, proportion] : assembly_community_fingerprints_final) {
      std::cout << "  Community: " << node << " Proportion: " << proportion << std::endl;
    }

    std::cout << "Adaptive" << std::endl;
    for (auto& [node, proportion] : adaptiveFinalCommunities) {
      std::cout << "Community: " << node << " Proportion: " << proportion << std::endl;
    }
    std::cout << "Adaptive fingerprints" << std::endl;
    for (auto& [node, proportion] : adaptive_community_fingerprints_final) {
      std::cout << "  Community: " << node << " Proportion: " << proportion << std::endl;
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

      int pos = (int)position_activation_order[i];

      // Actually call function that handles between-cell
      // movement
      // ORIGINAL CALL:
      //  - DoRepro(pos, world, next_world, config->SEEDING_PROB(), config->PROB_CLEAR(), config->DIFFUSION(), false);
      // (1) Do group reproduction?
      //  - Note: group repro was never possible here
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
    }
  }

  // This function should be called to create a stable copy of the world
  world_t stableUpdate(const world_t& custom_world, int max_updates=10000){

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

    for (int i = 0; i < max_updates; i++) {

      // Handle population growth for each cell
      for (size_t pos = 0; pos < custom_world.size(); pos++) {
        DoGrowth(pos, stable_world, next_stable_world);
      }

      // We can stop iterating early if the world has already stabilized
      double delta = 0;
      double epsilon = .0001;
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



  std::map<std::string, double> getFinalCommunities(const world_t& stable_world) {
    double size = stable_world.size();
    std::map<std::string, double> finalCommunities;
    const auto& subCommunities = community_structure.GetSubCommunities();
    for (const emp::vector<double>& cell: stable_world) {
      std::string comm = "";
      for (size_t i = 0; i < cell.size(); i++) {
        std::string count = std::to_string(cell[i]);
        // Formatting to make the strings more legible
        count.erase(count.find_last_not_of('0') + 1, std::string::npos);
        count.erase(count.find_last_not_of('.') + 1, std::string::npos);
        // Check if the species is present, and if it is make sure it has interactions
        if (count.compare("0") != 0) {
          int sc = -1;
          for (size_t j = 0; j < subCommunities.size(); j++) {
            for (size_t x : subCommunities[j]) {
              if (x == i) {
                sc = j;
              }
            }
          }
          bool interacts = false;
          for (size_t m : subCommunities[sc]) {
            if (cell[m] != 0 && m != i) {
              comm.append(count + " ");
              interacts = true;
              break;
            }
          }
          if (!interacts) {
            comm.append("PWI ");
          }
        } else {
          // If there are zero present
          comm.append("0 ");
        }
      }
      if (finalCommunities.find(comm) == finalCommunities.end()) {
          finalCommunities.insert({comm, 1});
      } else {
        finalCommunities[comm] += 1;
      }
    }
    // Get proportion
    for(auto& [key, val] : finalCommunities){
      val = val/size;
    }
    return finalCommunities;
  }

  // Returns mapping from community fingerprint to proportion found in given world.
  std::map<emp::BitVector, double> IdentifyWorldCommunities(
    const world_t& in_world,
    std::function<bool(double)> is_present = [](double count) -> bool { return count >= 1.0; }
  ) {
    std::map<emp::BitVector, double> communities;
    for (const emp::vector<double>& cell : in_world) {
      emp::BitVector community(N_TYPES, false);
      emp_assert(N_TYPES == cell.size());
      // Fingerprint this cell
      for (size_t spec_i = 0; spec_i < cell.size(); ++spec_i) {
        community[spec_i] = is_present(cell[spec_i]);
      }
      // If this is the first time we've seen this community, make note
      if (!emp::Has(communities, community)) {
        communities[community] = 0.0;
      }
      communities[community] += 1;
    }
    // Convert community counts to proportions
    const double world_size = in_world.size();
    emp_assert(world_size > 0);
    for (auto& [key, val] : communities) {
      val = val / world_size;
    }
    return communities;
  }

  // Helper functions not related to the running of the worlds go down here

  // Species cannot be a part of a community they have no interaction with
  // Find the connected components of the interaction matrix, and determine sub-communites
  emp::vector<emp::vector<size_t>> FindSubCommunities() {
    // Subcommunities are connected components of the interaction matrix.
    // "Connectivity" should include connection in either direction
    return utils::FindConnectedComponents<double>(
      interactions,
      [](const interaction_mat_t& matrix, size_t from, size_t to) -> bool {
        return (matrix[from][to] != 0) || (matrix[to][from] != 0);
      }
    );
  }

  // Calculate the growth rate (one measurement of fitness)
  // for a given cell
  double CalcGrowthRate(size_t pos, world_t & curr_world) {
    return doCalcGrowthRate(curr_world[pos]);
  }

  double doCalcGrowthRate(emp::vector<double> community){
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

    for (size_t i = 0; i < N_TYPES; i++) {
      interactions[i].resize(N_TYPES);
      for (size_t j = 0; j < N_TYPES; j++) {
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
    for (size_t i = 0; i < N_TYPES; i++) {
      interactions[i].resize(N_TYPES);
      for (size_t j = 0; j < N_TYPES; j++) {
        interactions[i][j] = interaction_data[i][j];
      }
    }
  }

  // Store the current interaction matrix in a file
  void WriteInteractionMatrix(std::string filename) {
    emp::File outfile;

    for (size_t i = 0; i < N_TYPES; i++) {
      outfile += emp::join(interactions[i], ",");
    }

    outfile.Write(filename);
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

} // End chemical_ecology namespace