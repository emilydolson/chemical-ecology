#pragma once

#include <functional>
#include <string>
#include <iostream>
#include <algorithm>

#include "emp/base/vector.hpp"
#include "emp/io/File.hpp"

namespace chemical_ecology {

// NOTE (@AML): This InteractionMatrix class is optimized for "ReadOnly" operations.
// I.e., configure it (by loading or setting the matrix) and details about the configuration
//       get cached for fast access. As such, the current implementation does not allow for on-the-fly
//       edits to the interaction matrix. If that's an important functionality, we could redesign this class a bit,
//       or I might setup a "mutable" version that can generate this read-only version.
// NOTE (@AML): The rest of the code hasn't necessarily been fully refactored to make use of
//              this class!
class InteractionMatrix {
public:
  using matrix_t = emp::vector< emp::vector<double> >;
  using interacts_fun_t =  std::function<bool(const matrix_t&, size_t, size_t)>;

protected:
  matrix_t interactions;         // Interaction matrix with weighted connections
  interacts_fun_t interacts_fun = [](const matrix_t& mat, size_t from, size_t to) -> bool {
    return (mat[from][to] != 0) || (mat[to][from] != 0);
  }; // Determines if two types should be considered "interacting"



public:
  InteractionMatrix() = default;

  InteractionMatrix(interacts_fun_t interacts) {
    SetInteractsFun(interacts);
  }

  const emp::vector<double>& operator[](size_t idx) const {
    return interactions[idx];
  }

  size_t GetNumTypes() const {
    return interactions.size();
  }

  const matrix_t& GetInteractions() const {
    return interactions;
  }

  bool Interacts(size_t from, size_t to) const {
    return interacts_fun(interactions, from, to);
  }

  void SetInteractsFun(interacts_fun_t interacts) {
    interacts_fun = interacts;
  }

  void SetInteractions(
    const matrix_t& in_interactions
  ) {
    interactions = in_interactions;
  }

  void LoadInteractions(
    const std::string& filename,
    std::optional<size_t> expected_size=std::nullopt
  ) {
    emp::File infile(filename);
    emp::vector<emp::vector<double>> interaction_data = infile.ToData<double>();

    const size_t num_types = (expected_size) ? expected_size.value() : interaction_data.size();

    interactions.resize(num_types, emp::vector<double>(num_types, 0.0));

    // NOTE (@AML): this doesn't look like it does anything extra to "interaction_data" ==> could just return interaction data?
    for (size_t i = 0; i < num_types; i++) {
      if (interaction_data.size() == 0) continue; // Skip over empty lines
      std::copy(
        interaction_data[i].begin(),
        interaction_data[i].end(),
        interactions[i].begin()
      );
      emp_assert(interactions[i].size() == num_types);
    }
  }

  void RandomizeInteractions(
    emp::Random& rnd,
    size_t num_types,
    double prob_interaction,
    double interaction_magnitude,
    bool self_interactions=true
  ) {
    interactions.resize(num_types, emp::vector<double>(num_types, 0.0));

    for (size_t i = 0; i < num_types; i++) {
      for (size_t j = 0; j < num_types; j++) {
        // If self interactions aren't allowed, leave as 0
        if (!self_interactions && (i == j)) continue;
        // PROB_INTERACTION determines probability of there being
        // an interaction between a given pair of types.
        // Controls sparsity of interaction matrix
        interactions[i][j] = (rnd.P(prob_interaction)) ?
          rnd.GetDouble(interaction_magnitude * -1, interaction_magnitude) :
          interactions[i][j];
      }
    }
  }

  // Store the current interaction matrix in a file
  void WriteInteractionMatrix(
    const std::string& filename
  ) const {
    emp::File outfile;
    const size_t num_types = interactions.size();

    for (size_t i = 0; i < num_types; i++) {
      emp_assert(interactions[i].size() == num_types);
      outfile += emp::join(interactions[i], ",");
    }

    outfile.Write(filename);
  }


};

} // End chemical_ecology namespace