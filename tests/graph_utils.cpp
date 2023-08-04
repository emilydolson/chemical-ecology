#define CATCH_CONFIG_MAIN

#include "Catch/single_include/catch2/catch.hpp"

#include <iostream>
#include <functional>

#include "emp/base/vector.hpp"
#include "emp/math/Random.hpp"

#include "chemical-ecology/utils/graph_utils.hpp"
#include "chemical-ecology/SpatialStructure.hpp"

TEST_CASE("DepthFirstTraversal should perform a depth first traversal starting at a given root node") {
  SECTION("DepthFirstTraversal should return one graph-sized component when given a toroidal grid") {
    using mat_t = emp::vector< emp::vector<bool> >;
    using conn_fun_t = std::function<bool(const mat_t&, size_t, size_t)>;

    chemical_ecology::SpatialStructure structure;
    chemical_ecology::ConfigureToroidalGrid(structure, 3, 3);

    conn_fun_t conn_fun = [](const mat_t& mat, size_t from, size_t to) -> bool {
      return mat[from][to];
    };

    auto components = chemical_ecology::utils::FindConnectedComponents(
      structure.GetConnectionMatrix(),
      conn_fun
    );

    REQUIRE(components.size() == 1);
    REQUIRE(components[0].size() == structure.GetNumPositions());
  }

  SECTION("DepthFirstTraversal should return many size-one components for a fully disconnected graph") {
    using mat_t = emp::vector< emp::vector<bool> >;
    using conn_fun_t = std::function<bool(const mat_t&, size_t, size_t)>;

    const std::string mat_path = "data/graph-no-connections.dat";
    chemical_ecology::SpatialStructure structure;
    structure.LoadStructureFromMatrix(mat_path);

    conn_fun_t conn_fun = [](const mat_t& mat, size_t from, size_t to) -> bool {
      return mat[from][to];
    };

    auto components = chemical_ecology::utils::FindConnectedComponents(
      structure.GetConnectionMatrix(),
      conn_fun
    );

    REQUIRE(components.size() == structure.GetNumPositions());
    for (size_t comp_i = 0; comp_i < components.size(); ++comp_i) {
      REQUIRE(components[comp_i].size() == 1);
    }
  }

}