#pragma once

#include <functional>
#include <map>
#include <iostream>
#include <string>
#include <algorithm>
#include <optional>

#include "emp/bits/BitVector.hpp"
#include "emp/base/vector.hpp"
#include "emp/datastructs/vector_utils.hpp"
#include "emp/tools/string_utils.hpp"

#include "chemical-ecology/CommunityStructure.hpp"

// This file defines:
// - RecordedCommunitySummary: stores summary information about a recorded community
// - RecordedCommunitySet: manages a set of summary information instances

namespace chemical_ecology {

// Manages a set of RecordedCommunitySummary instances
// TODO - test!
// NOTE (@AML): Set of existing accessors probably not complete. Added as needed by world.
// NOTE (@AML): If you remove something from the set, size_t ids are no longer guaranteed to
//              be the same before / after the remove
template<typename SUMMARY_KEY_T>
class RecordedCommunitySet {
public:
  using summary_key_fun_t = std::function<const SUMMARY_KEY_T& (const RecordedCommunitySummary&)>;
protected:

  // REMINDER - summary IDs are not necessarily stable between adds! (Remove operations can change IDs)
  //  - I.e., don't give out internal summary IDs!
  emp::vector<RecordedCommunitySummary> summary_set;  // Stores recorded community summaries, *one* summary for each type added
  emp::vector<size_t> community_counts;               // How many of each recorded community "type" has been recorded?
  std::map<SUMMARY_KEY_T, size_t> summary_id_map;     // Map summary keys to position in summaries vector
  summary_key_fun_t get_summary_key_fun;              // Given a summary, extracts component that uniquely identifies the summary
                                                      // Determines what recorded communities should be considered identical

  std::optional<size_t> GetCommunityID(const SUMMARY_KEY_T& key) const {
    return (emp::Has(summary_id_map, key)) ? summary_id_map[key] : std::nullopt;
  }

  // Removes summary with given ID, adjusts other IDs as necessary
  void Remove(size_t summary_id) {
    emp_assert(summary_id < summary_set.size());
    size_t back_id = summary_set.size() - 1;
    if (summary_id != back_id) {
      // Swap summary to be removed with summary with last id (to avoid changing more than one other ID)
      const auto back_key = get_summary_key_fun(summary_set[back_id]);
      // Change the ID assocaited with the last item (it is about to be swapped forward)
      summary_id_map[back_key] = summary_id;
      // Swap the last summary in the set forward with the summary to be removed
      std::swap(summary_set[summary_id], summary_set[back_id]);
      std::swap(community_counts[summary_id], community_counts[back_id]);
    }
    // Summary to be removed should be last item in the vector
    summary_set.pop_back();
    community_counts.pop_back();
  }

public:

  RecordedCommunitySet(
    summary_key_fun_t get_summary_key
  ) :
    get_summary_key_fun(get_summary_key)
  { }

  void Clear() {
    summary_set.clear();
    community_counts.clear();
    summary_id_map.clear();
  }

  size_t GetSize() const {
    emp_assert(summary_id_map.size() == summary_set.size());
    emp_assert(community_counts.size() == summary_set.size());
    return summary_set.size();
  }

  // Get numeric ID assigned to given summary
  std::optional<size_t> GetCommunityID(const RecordedCommunitySummary& summary) const {
    const auto key = get_summary_key_fun(summary);
    return (emp::Has(summary_id_map, key)) ? std::optional<size_t>{summary_id_map.at(key)} : std::nullopt;
  }

  bool Has(const RecordedCommunitySummary& summary) const {
    return emp::Has(summary_id_map, get_summary_key_fun(summary));
  }

  const emp::vector<size_t>& GetCommunityCounts() const { return community_counts; }

  size_t GetCommunityCount(size_t id) const {
    return community_counts[id];
  }

  double GetCommunityProportion(size_t id) const {
    emp_assert(emp::Sum(community_counts) > 0);
    return (double)community_counts[id] / (double)emp::Sum(community_counts);
  }

  const RecordedCommunitySummary& GetCommunitySummary(size_t id) const {
    return summary_set[id];
  }

  void Add(const RecordedCommunitySummary& summary, size_t community_count=1) {
    const auto summary_key = get_summary_key_fun(summary);
    size_t summary_id = 0;

    // If we haven't encounted this community type before, add to set.
    // Otherwise, get the correct id
    if (!emp::Has(summary_id_map, summary_key)) {
      summary_id = summary_set.size();
      summary_id_map[summary_key] = summary_id;
      summary_set.emplace_back(summary);
      community_counts.emplace_back(0);
    } else {
      summary_id = summary_id_map[summary_key];
    }

    community_counts[summary_id] += community_count;
  }

  // Add multiple summaries at once (each entry will increment associated community count by 1)
  void Add(const emp::vector<RecordedCommunitySummary>& summaries) {
    for (const auto& summary : summaries) {
      Add(summary, 1);
    }
  }

  // Remove 'remove_count' number of recorded communities of specified type
  void Remove(const RecordedCommunitySummary& summary, size_t remove_count) {
    const auto summary_key = get_summary_key_fun(summary);
    emp_assert(emp::Has(summary_id_map, summary_key), "Summary not in set.");

    const size_t summary_id = summary_id_map[summary_key];
    const size_t prev_count = community_counts[summary_id];
    // If asked to remove more than we have, delete the community all-together
    if (remove_count >= prev_count) {
      Remove(summary_id);
    } else {
      community_counts[summary_id] -= remove_count;
    }
  }

  // Completely removes summary from set (regardless of current summary count)
  // - Note that this is a somewhat expensive operation
  void Remove(const RecordedCommunitySummary& summary) {
    const auto summary_key = get_summary_key_fun(summary);
    emp_assert(emp::Has(summary_id_map, summary_key), "Summary not in set.");
    const size_t summary_id = summary_id_map[summary_key];
    Remove(summary_id);
  }

}; // End of RecordedCommunitySet definition

} // End of chemical_ecology namespace