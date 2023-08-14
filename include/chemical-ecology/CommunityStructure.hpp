#pragma once

#include <functional>
#include <algorithm>
#include <iostream>

#include "emp/base/vector.hpp"
#include "emp/bits/BitVector.hpp"
#include "emp/datastructs/vector_utils.hpp"

#include "chemical-ecology/utils/graph_utils.hpp"
#include "chemical-ecology/InteractionMatrix.hpp"

namespace chemical_ecology {

// Extracts and manages community structure from given interaction matrix
class CommunityStructure {
public:
  using interaction_mat_t = emp::vector<emp::vector<double>>;

protected:
  emp::vector<emp::vector<size_t>> subcommunities;
  // NOTE (@AML): Move away from or embrace fingerprint terminology?
  emp::vector<emp::BitVector> subcommunity_fingerprints; // Each bit vector represents presence / absence of species
  emp::vector<size_t> species_to_subcommunity_id;
  size_t num_species;

  emp::vector<emp::BitVector> species_shares_subcommunity_with; // For each species, what other species does it share a subcommunity with
  emp::vector<emp::BitVector> species_interacts_with;    // For each species, what other species does it interact with?

public:
  CommunityStructure() = default;

  CommunityStructure(
    const InteractionMatrix& interaction_matrix
  ) {
    SetStructure(interaction_matrix);
  }

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

  void SetStructure(const InteractionMatrix& interaction_matrix) {
    SetStructure(
      interaction_matrix.GetInteractions(),
      [&interaction_matrix](
        const interaction_mat_t& mat,
        size_t from,
        size_t to
      ) -> bool {
        return interaction_matrix.Interacts(from, to);
      }
    );
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
      interacts,
      true
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

  // Get number of species part of specified subcommunity
  size_t GetNumMembers(size_t subcommunity_id) const { return subcommunities[subcommunity_id].size(); }

  bool SpeciesInteractsWith(size_t focal_species_id, size_t interacts_with_id) const {
    return species_interacts_with[focal_species_id][interacts_with_id];
  }

  bool SpeciesSharesSubCommunityWith(size_t focal_species_id, size_t shares_with_id) const {
    return species_shares_subcommunity_with[focal_species_id][shares_with_id];
  }

  const emp::vector<emp::BitVector>& GetSharedSubCommunityInfo() const {
    return species_shares_subcommunity_with;
  }

  const emp::vector<emp::BitVector>& GetSpeciesInteractsWith() const {
    return species_interacts_with;
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

  const emp::BitVector& GetSubCommunityPresent(size_t community_id) const {
    return subcommunity_fingerprints[community_id];
  }

  const emp::vector<size_t>& GetSpeciesToSubCommunityMapping() const {
    return species_to_subcommunity_id;
  }

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
    // species_interacts_with.clear();
  }

};

bool PathExists(
  const CommunityStructure& community_structure,
  size_t from,
  size_t to,
  const emp::BitVector& limit_path_to
) {

  if (limit_path_to.None()) return false;

  std::deque<size_t> next;
  std::unordered_set<size_t> discovered;
  next.emplace_back((size_t)limit_path_to.FindOne());
  discovered.emplace(next.back());
  emp_assert(next.back() < limit_path_to.size());

  while (next.size() > 0) {
    // Current node at front of queue
    const size_t cur_node = next.front();
    emp_assert(emp::Has(discovered, cur_node));
    emp_assert(limit_path_to[cur_node]);
    // Queue all connected nodes that have not yet been discovered and are allowed by 'limit_path_to'
    for (size_t node_id = 0; node_id < community_structure.GetNumSpecies(); ++node_id) {
      if (!emp::Has(discovered, node_id) && limit_path_to[node_id] && community_structure.SpeciesInteractsWith(cur_node, node_id)) {
        if (node_id == to) return true;
        discovered.emplace(node_id);
        next.emplace_back(node_id);
      }
    }

    // Pop front of queue
    next.pop_front();
  }

  // Completed search without discovering target
  return false;
}

} // End chemical_ecology namespace