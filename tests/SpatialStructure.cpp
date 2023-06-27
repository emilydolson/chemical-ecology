#define CATCH_CONFIG_MAIN

#include "Catch/single_include/catch2/catch.hpp"

#include <iostream>

#include "chemical-ecology/SpatialStructure.hpp"

#include "emp/base/vector.hpp"

TEST_CASE("Test SetStructure(emp::vector<emp::vector<size_t>)&)") {



  chemical_ecology::SpatialStructure ring_structure;
  emp::vector< emp::vector<size_t> > ring_map = {
    /* 0 -> */ {1},
    /* 1 -> */ {2},
    /* 2 -> */ {0}
  };
  ring_structure.SetStructure(ring_map);
  ring_structure.Print(std::cout);

  emp::vector< emp::vector<bool> > ring_matrix = {
    { 0, 1, 0 },
    { 0, 0, 1 },
    { 1, 0, 0 }
  };
  ring_structure.SetStructure(ring_matrix);
  ring_structure.Print(std::cout);

  // TODO - test connections

}