#pragma once

#include <functional>
#include <algorithm>
#include <iostream>

#include "emp/base/vector.hpp"
#include "emp/bits/BitVector.hpp"

namespace chemical_ecology {

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

// Community summary information (extracted from cell of a world)
// NOTE: depending on how sophisticated this ends up, should get moved to own file & promoted to a proper class
// TODO - write Print functions
struct CommunitySummary {
  emp::vector<size_t> counts; // species counts (uniquely identifies this community)
  emp::vector<size_t> present_species_ids;
  emp::BitVector present; // species present/absent fingerprint
  // NOTE: as per slack conversation, might be worth renaming 'present_no_interactions'
  emp::BitVector present_no_interactions; // species present without interactions with *OTHER* species

  CommunitySummary() = default;

  CommunitySummary(
    const emp::vector<double>& member_counts,
    std::function<bool(double)> is_present,
    const CommunityStructure& community_structure
  ) {
    SummarizeCommunity(member_counts, is_present, community_structure);
  }

  // Summarize given community (extracted from a world)
  // Process / summarize member counts
  // Identifies which species are present
  // Identifies which species are present but have no co-subcommunity members also present
  void SummarizeCommunity(
    const emp::vector<double>& member_counts,
    std::function<bool(double)> is_present,
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

    // Identify members that are present with no interactions (with *OTHER* members)
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
  bool operator<(const CommunitySummary& in) const {
    emp_assert(in.counts.size() == counts.size(), "Expected communities to have same number of potential species");
    return counts < in.counts;
  }

};

} // End chemical_ecology namespace