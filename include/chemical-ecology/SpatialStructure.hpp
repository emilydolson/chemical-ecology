#pragma once

// This file contains the SpatialStructure class that defines and manages
// spatial structure for the AEcoWorld.
//
// The class is designed with the following trade-offs:
// - Efficient random neighbor selection
// - Efficient neighbor checking

//

#include <iostream>
#include <algorithm>
#include <unordered_set>
#include "emp/base/vector.hpp"
#include "emp/math/Random.hpp"
#include "emp/math/random_utils.hpp"

namespace chemical_ecology {

// Functionality:

// Internals:
// - Track spatial structure as
//   - vector of vectors (for easy random access)
//   - vector of maps

class SpatialStructure {
public:

protected:

  // Mapping of source positions ==> destination positions.
  // Destination IDs are kept in sorted order.
  emp::vector< emp::vector<size_t> > ordered_connections;

  // Matrix specifying whether two positions are connected ([from][to]).
  emp::vector< emp::vector<bool> > connection_matrix;

  // Internal verification that ordered_connections and connection_matrix are
  // consistent. Used for internal asserts in debug mode.
  bool VerifyConnectionConsistency() const {
    const size_t num_positions = connection_matrix.size();
    // Check that ordered_connections has same number of positions as matrix
    if (num_positions != ordered_connections.size()) {
      return false;
    }
    // Check consistency in connections across representations
    for (size_t from = 0; from < num_positions; ++from) {
      emp_assert(connection_matrix[from].size() == num_positions);
      for (size_t to = 0; to < num_positions; ++to) {
        const bool connected = connection_matrix[from][to];
        const auto& neighbors = ordered_connections[from];
        emp_assert(std::is_sorted(neighbors.begin(), neighbors.end()));
        const bool is_neighbor = std::binary_search(
          neighbors.begin(),
          neighbors.end(),
          to
        );
        if (connected != is_neighbor) {
          return false;
        }
      }
    }
    return true;
  }

public:

  // Configure spatial structure from mapping of "from" positions to "to" positions
  void SetStructure(const emp::vector< emp::vector<size_t> >& in_struct) {
    // Configure ordered connections (copy over and sort connections)
    const size_t num_positions = in_struct.size();
    ordered_connections = in_struct;
    for (emp::vector<size_t>& neighbors : ordered_connections) {
      std::sort(
        neighbors.begin(),
        neighbors.end()
      );
    }
    // Configure connection matrix
    connection_matrix.clear();
    connection_matrix.resize(
      num_positions,
      emp::vector<bool>(num_positions, false)
    );
    for (size_t from = 0; from < num_positions; ++from) {
      const auto& neighbors = ordered_connections[from];
      for (size_t to : neighbors) {
        connection_matrix[from][to] = true;
      }
    }
    emp_assert(VerifyConnectionConsistency());
  }

  // Configure spatial structure from a connection matrix (maps [from][to])
  void SetStructure(const emp::vector< emp::vector<bool> >& in_struct) {
    // Configure connection matrix (copy from parameter)
    const size_t num_positions = in_struct.size();
    connection_matrix = in_struct;
    // Configure ordered connections
    ordered_connections.clear();
    ordered_connections.resize(
      num_positions,
      emp::vector<size_t>(0)
    );
    for (size_t from = 0; from < num_positions; ++from) {
      emp_assert(connection_matrix.size() == num_positions, "Connection matrix must be square");
      for (size_t to = 0; to < num_positions; ++to) {
        const bool connected = connection_matrix[from][to];
        if (connected) {
          ordered_connections[from].emplace_back(to);
        }
      }
    }
    emp_assert(VerifyConnectionConsistency());
  }

  void Connect(size_t from, size_t to) {
    // TODO
  }

  void Disconnect(size_t from, size_t to) {
    // TODO
  }

  bool IsConnected(size_t from, size_t to) {
    // TODO
    return false;
  }

  // Get the total number of positions in the spatial structure
  size_t GetNumPositions() {
    return connection_matrix.size();
  }

  // Get an ordered list of neighbors for given position
  const emp::vector<size_t>& GetNeighbors(size_t pos) {
    emp_assert(pos < ordered_connections.size());
    return ordered_connections[pos];
  }

  size_t GetRandomNeighbor(emp::Random& rnd, size_t pos) {
    // TODO
    return 0;
  }

  void LoadStructure(/* TODO */) {
    // TODO
  }

  void Print(std::ostream& os, bool as_mapping = true) const {
    if (as_mapping) {
      PrintConnectionMapping(os);
    } else {
      PrintConnectionMatrix(os);
    }
  }

  void PrintConnectionMatrix(std::ostream& os) const {
    // TODO
  }

  // Print mapping of "from" positions to "to" positions
  void PrintConnectionMapping(std::ostream& os) const {
    emp_assert(VerifyConnectionConsistency());
    const size_t num_positions = ordered_connections.size();
    for (size_t from = 0; from < num_positions; ++from) {
      const auto& neighbors = ordered_connections[from];
      os << from << ":";
      for (size_t i = 0; i < neighbors.size(); ++i) {
        if (i) os << ",";
        os << " " << neighbors[i];
      }
      os << std::endl;
    }
  }

}; // End SpatialStructure class

} // End chemical_ecology namespace