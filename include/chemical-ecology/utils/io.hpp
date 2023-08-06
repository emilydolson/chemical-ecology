#pragma once

#include <optional>

#include "emp/base/vector.hpp"
#include "emp/io/File.hpp"

// Collection of stand-alone IO functions

namespace chemical_ecology::utils {

// Load an interaction matrix from the specified file
emp::vector<emp::vector<double>> LoadInteractionMatrix(
  const std::string& filename,
  std::optional<size_t> expected_size=std::nullopt
) {
  emp::File infile(filename);
  emp::vector<emp::vector<double>> interaction_data = infile.ToData<double>();

  const size_t N_TYPES = (expected_size) ? expected_size.value() : interaction_data.size();

  emp::vector<emp::vector<double>> interactions(N_TYPES, emp::vector<double>(N_TYPES, 0.0));

  // NOTE (@AML): this doesn't look like it does anything extra to "interaction_data" ==> could just return interaction data?
  for (size_t i = 0; i < N_TYPES; i++) {
    if (interaction_data.size() == 0) continue; // Skip over empty lines
    std::copy(
      interaction_data[i].begin(),
      interaction_data[i].end(),
      interactions[i].begin()
    );
    emp_assert(interactions[i].size() == N_TYPES);
  }

  return interactions;
}

// Store the current interaction matrix in a file
void WriteInteractionMatrix(std::string filename, const emp::vector<emp::vector<double>>& interactions) {
  emp::File outfile;
  const size_t N_TYPES = interactions.size();

  for (size_t i = 0; i < N_TYPES; i++) {
    outfile += emp::join(interactions[i], ",");
  }

  outfile.Write(filename);
}

} // End chemical_ecology::utils namespace