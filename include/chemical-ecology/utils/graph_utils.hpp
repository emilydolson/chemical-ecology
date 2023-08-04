#pragma once

#include <functional>
#include <unordered_set>
#include <algorithm>

#include "emp/base/vector.hpp"

namespace chemical_ecology::utils {


// Stateful depth-first traversal.
// Allows repeated traversals from different starting points while maintaining
// a shared visited list.
// Operates on graphs in *matrix* form!
template<typename T>
class DepthFirstTraversal {
public:
  using matrix_t = emp::vector< emp::vector<T> >;
  using connected_fun_t = std::function<bool(const matrix_t&, size_t from, size_t to)>;

protected:
  const matrix_t& connection_matrix;  // Square matrix of connections
  connected_fun_t is_connected;       // Function that defines whether things are connected
  emp::vector<bool> visited;
  std::unordered_set<size_t> unvisited;

  // Recursively traverse graph starting from root, adding traversed nodes to component vector
  void Traverse(size_t root, emp::vector<size_t>& component) {
    // Mark root as visited, add to current component
    visited[root] = true;
    component.emplace_back(root);
    unvisited.erase(root);

    for (size_t i = 0; i < connection_matrix[root].size(); ++i) {
      emp_assert(
        connection_matrix.size() == connection_matrix[root].size(),
        "Connection matrix must be square."
      );
      const bool connected = is_connected(connection_matrix, root, i);
      if (connected && !visited[i]) {
        Traverse(i, component);
      }
    }

  }

public:
  // Initialize traversal with connections and connection function
  DepthFirstTraversal(
    const matrix_t& connections,
    connected_fun_t is_connected_fun
  ) :
    connection_matrix(connections),
    is_connected(is_connected_fun),
    visited(connection_matrix.size(), false)
  {
    // Populate unvisited
    for (size_t i = 0; i < visited.size(); ++i) {
      unvisited.insert(i);
    }
  }

  // Performs traversal from a given starting node
  emp::vector<size_t> operator()(size_t root) {
    emp::vector<size_t> component;
    Traverse(root, component);
    return component;
  }

  // Get visited nodes
  const emp::vector<bool>& GetVisited() const {
    return visited;
  }

  // Get nodes not yet visited
  const std::unordered_set<size_t>& GetUnvisited() const {
    return unvisited;
  }

  // Has node ID been visited?
  bool Visited(size_t id) const {
    emp_assert(id < visited.size());
    return visited[id];
  }

  // Reset visited state
  void Reset() {
    visited.clear();
    visited.resize(connection_matrix.size(), false);
    unvisited.clear();
    for (size_t i = 0; i < visited.size(); ++i) {
      unvisited.insert(i);
    }
  }

};

// Find connected components in a graph (represented as a matrix)
template<typename T>
emp::vector< emp::vector<size_t> > FindConnectedComponents(
  const emp::vector< emp::vector<T> >& connection_matrix,
  std::function<bool(
    const emp::vector< emp::vector<T> >&,
    size_t,
    size_t)
  > is_connected,
  bool sort_node_ids=false
) {

  emp::vector< emp::vector<size_t> > components;
  DepthFirstTraversal<T> dfs(connection_matrix, is_connected);
  while (dfs.GetUnvisited().size() > 0) {
    const size_t next_root = *(dfs.GetUnvisited().begin()); // Grab arbitrary value from unvisited
    components.emplace_back(
      dfs(next_root)
    );
    // dfs makes no guarantees about id ordering
    if (sort_node_ids) {
      std::sort(
        components.back().begin(),
        components.back().end()
      );
    }
    emp_assert(components.back().size() > 0, "Most recent component should always have at least one node");
  }


  return components;
}

} // End chemical_ecology::utils namespace