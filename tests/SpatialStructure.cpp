#define CATCH_CONFIG_MAIN

#include "Catch/single_include/catch2/catch.hpp"

#include <iostream>

#include "chemical-ecology/SpatialStructure.hpp"

#include "emp/base/vector.hpp"

TEST_CASE("Can define spatial structure from a connection mapping") {
  chemical_ecology::SpatialStructure ring_structure;

  emp::vector< emp::vector<size_t> > ring_map = {
    /* 0 -> */ {1},
    /* 1 -> */ {2},
    /* 2 -> */ {0}
  };
  ring_structure.SetStructure(ring_map);
  // ring_structure.Print(std::cout);

  REQUIRE(!ring_structure.IsConnected(0, 0));
  REQUIRE(ring_structure.IsConnected(0, 1));
  REQUIRE(!ring_structure.IsConnected(0, 2));

  REQUIRE(!ring_structure.IsConnected(1, 0));
  REQUIRE(!ring_structure.IsConnected(1, 1));
  REQUIRE(ring_structure.IsConnected(1, 2));

  REQUIRE(ring_structure.IsConnected(2, 0));
  REQUIRE(!ring_structure.IsConnected(2, 1));
  REQUIRE(!ring_structure.IsConnected(2, 2));

}

TEST_CASE("Can define spatial structure from a connection matrix") {
  chemical_ecology::SpatialStructure ring_structure;

  emp::vector< emp::vector<bool> > ring_matrix = {
    { 0, 1, 0 },
    { 0, 0, 1 },
    { 1, 0, 0 }
  };
  ring_structure.SetStructure(ring_matrix);
  // ring_structure.Print(std::cout);

  REQUIRE(!ring_structure.IsConnected(0, 0));
  REQUIRE(ring_structure.IsConnected(0, 1));
  REQUIRE(!ring_structure.IsConnected(0, 2));

  REQUIRE(!ring_structure.IsConnected(1, 0));
  REQUIRE(!ring_structure.IsConnected(1, 1));
  REQUIRE(ring_structure.IsConnected(1, 2));

  REQUIRE(ring_structure.IsConnected(2, 0));
  REQUIRE(!ring_structure.IsConnected(2, 1));
  REQUIRE(!ring_structure.IsConnected(2, 2));
}


TEST_CASE("Can make new connections and remove existing connections") {
  chemical_ecology::SpatialStructure ring_structure;
  emp::vector< emp::vector<bool> > ring_matrix = {
    { 0, 1, 0 },
    { 0, 0, 1 },
    { 1, 0, 0 }
  };
  ring_structure.SetStructure(ring_matrix);
  // Make a connection
  ring_structure.Connect(1, 0);
  // ring_structure.Print(std::cout);

  REQUIRE(!ring_structure.IsConnected(0, 0));
  REQUIRE(ring_structure.IsConnected(0, 1));
  REQUIRE(!ring_structure.IsConnected(0, 2));

  REQUIRE(ring_structure.IsConnected(1, 0));
  REQUIRE(!ring_structure.IsConnected(1, 1));
  REQUIRE(ring_structure.IsConnected(1, 2));

  REQUIRE(ring_structure.IsConnected(2, 0));
  REQUIRE(!ring_structure.IsConnected(2, 1));
  REQUIRE(!ring_structure.IsConnected(2, 2));

  // Remove a connection
  ring_structure.Disconnect(1, 0);
  // ring_structure.Print(std::cout);

  REQUIRE(!ring_structure.IsConnected(0, 0));
  REQUIRE(ring_structure.IsConnected(0, 1));
  REQUIRE(!ring_structure.IsConnected(0, 2));

  REQUIRE(!ring_structure.IsConnected(1, 0));
  REQUIRE(!ring_structure.IsConnected(1, 1));
  REQUIRE(ring_structure.IsConnected(1, 2));

  REQUIRE(ring_structure.IsConnected(2, 0));
  REQUIRE(!ring_structure.IsConnected(2, 1));
  REQUIRE(!ring_structure.IsConnected(2, 2));
}