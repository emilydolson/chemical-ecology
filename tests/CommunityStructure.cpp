#define CATCH_CONFIG_MAIN

#include "Catch/single_include/catch2/catch.hpp"

#include <iostream>
#include <functional>
#include <string>
#include <set>

#include "chemical-ecology/CommunityStructure.hpp"
#include "chemical-ecology/InteractionMatrix.hpp"

#include "emp/base/vector.hpp"
#include "emp/math/Random.hpp"
#include "emp/datastructs/set_utils.hpp"
#include "emp/datastructs/vector_utils.hpp"

std::function<bool(const emp::vector<emp::vector<double>>&, size_t, size_t)> interacts_fun = [](
  const emp::vector<emp::vector<double>>& matrix,
  size_t from,
  size_t to
) -> bool {
  return (matrix[from][to] != 0) || (matrix[to][from] != 0);
};


void VerifyStructure(const chemical_ecology::CommunityStructure& comm) {
  // subcommunities
  // subcommunity_fingerprints
  // species_to_subcommunity_id
  // num_species
  // species_shares_subcommunity_with
  // species_interacts_with
  // Localize member variables for convenience
  const size_t num_species = comm.GetNumSpecies();
  const size_t num_sub_comms = comm.GetNumSubCommunities();
  const auto& subcommunities = comm.GetSubCommunities();
  const auto& subcommunity_fingerprints = comm.GetFingerprints();
  const auto& species_to_subcommunity_id = comm.GetSpeciesToSubCommunityMapping();
  const auto& species_shares_subcommunity_with = comm.GetSharedSubCommunityInfo();
  const auto& species_interacts_with = comm.GetSpeciesInteractsWith();

  // (1) Check num_species accuracy
  // valid = valid && species_to_subcommunity_id.size() == num_species;
  // valid = valid && species_shares_subcommunity_with.size() == num_species;
  // valid = valid && species_interacts_with.size() == num_species;

  REQUIRE(species_to_subcommunity_id.size() == num_species);
  REQUIRE(species_shares_subcommunity_with.size() == num_species);
  REQUIRE(species_interacts_with.size() == num_species);

  // Number of members across subcommunities should sum to num_species
  size_t spec_cnt = 0;
  for (size_t i = 0; i < subcommunities.size(); ++i) {
    spec_cnt += subcommunities[i].size();
    REQUIRE(subcommunity_fingerprints[i].GetSize() == num_species);
  }
  REQUIRE(spec_cnt == num_species);
  for (const auto& share : species_shares_subcommunity_with) {
    REQUIRE(share.GetSize() == num_species);
  }
  for (const auto& interact : species_interacts_with) {
    REQUIRE(interact.GetSize() == num_species);
  }

  // (2) Verify number of subcommunities
  REQUIRE(num_sub_comms == subcommunities.size());
  REQUIRE(num_sub_comms == subcommunity_fingerprints.size());

  // (3) Verify species' subcommunity ids
  for (size_t i = 0; i < num_species; ++i) {
    const size_t subcomm_id = species_to_subcommunity_id[i];
    REQUIRE(subcommunity_fingerprints[subcomm_id][i]);
    REQUIRE(emp::Has(subcommunities[subcomm_id], i));
  }

  // (4) Verify that each species id is represented once (and only once) in subcommunities
  std::set<size_t> represented;
  for (const auto& subcomm : subcommunities) {
    for (size_t spec_id : subcomm) {
      REQUIRE(!emp::Has(represented, spec_id));
      represented.insert(spec_id);
    }
  }
  REQUIRE(represented.size() == num_species);

  // (5) Verify that interactions are represented in subcommunity
  for (size_t spec_id = 0; spec_id < num_species; ++spec_id) {
    const auto& interacts_with = species_interacts_with[spec_id];
    const auto& shares_with = species_shares_subcommunity_with[spec_id];
    const size_t subcomm_id = species_to_subcommunity_id[spec_id];
    for (size_t other_id = 0; other_id < num_species; ++other_id) {
      if (interacts_with[other_id]) {
        REQUIRE(emp::Has(subcommunities[subcomm_id], other_id));
        REQUIRE(shares_with[other_id]);
      }
      if (shares_with[other_id]) {
        REQUIRE(emp::Has(subcommunities[subcomm_id], other_id));
      }
    }
  }

}

void VerifyStructure(const std::string& matrix_path) {

  chemical_ecology::InteractionMatrix interaction_matrix(
    interacts_fun
  );
  interaction_matrix.LoadInteractions(matrix_path);


  chemical_ecology::CommunityStructure comm_struct(
    interaction_matrix
  );

  VerifyStructure(comm_struct);
}

TEST_CASE("Should be able to default-initialize CommunityStructure object") {
  chemical_ecology::CommunityStructure comm_struct;
}

TEST_CASE("CommunityStructure should correctly reflect no-connection interaction matrix") {
  const std::string mat_path = "data/graph-no-connections.dat";

  chemical_ecology::InteractionMatrix interaction_matrix(
    interacts_fun
  );

  interaction_matrix.LoadInteractions(
    mat_path,
    3
  );

  chemical_ecology::CommunityStructure comm_struct(
    interaction_matrix
  );

  // Check number of species / subcommunities as expected
  REQUIRE(comm_struct.GetNumSpecies() == 3);
  REQUIRE(comm_struct.GetNumSubCommunities() == 3);

  // Check that species interactions are as expected
  for (size_t spec_id = 0; spec_id < 3; ++spec_id) {
    for (size_t other_id = 0; other_id < 3; ++other_id) {
      const bool valid_interactions = (!comm_struct.SpeciesInteractsWith(spec_id, other_id)) || spec_id == other_id;
      REQUIRE( valid_interactions );
      const bool valid_shared = (!comm_struct.SpeciesSharesSubCommunityWith(spec_id, other_id)) || spec_id == other_id;
      REQUIRE( valid_shared );
    }
  }

  // Check that each species is represented once and only once
  std::set<size_t> represented;
  for (const auto& subcomm : comm_struct.GetSubCommunities()) {
    REQUIRE(subcomm.size() == 1);
    for (size_t id : subcomm) {
      REQUIRE(!emp::Has(represented, id));
      represented.insert(id);
    }
  }
  REQUIRE(represented.size() == 3);

  // Check fingerprints
  for (const auto& present : comm_struct.GetFingerprints()) {
    REQUIRE(present.size() == 3);
    REQUIRE(present.CountOnes() == 1);
  }

  VerifyStructure(comm_struct);

}

TEST_CASE("Verify structure for misc loaded matrices") {
  SECTION("data/403_matrix.dat") {
    VerifyStructure("data/403_matrix.dat");
  }
  SECTION("data/class1.dat") {
    VerifyStructure("data/class1.dat");
  }
  SECTION("data/class2.dat") {
    VerifyStructure("data/class2.dat");
  }
  SECTION("data/class3.dat") {
    VerifyStructure("data/class3.dat");
  }
  SECTION("data/class4.dat") {
    VerifyStructure("data/class4.dat");
  }
  SECTION("data/interaction_matrix.dat") {
    VerifyStructure("data/interaction_matrix.dat");
  }
  SECTION("data/newPOC.dat") {
    VerifyStructure("data/newPOC.dat");
  }
  SECTION("data/newPOC.dat") {
    VerifyStructure("data/newPOC.dat");
  }
}