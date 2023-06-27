#define CATCH_CONFIG_MAIN

#include "Catch/single_include/catch2/catch.hpp"

#include <iostream>

#include "chemical-ecology/SpatialStructure.hpp"

#include "emp/base/vector.hpp"
#include "emp/math/Random.hpp"

TEST_CASE("Can define spatial structure from a connection mapping") {
  chemical_ecology::SpatialStructure ring_structure;

  emp::vector< emp::vector<size_t> > ring_map = {
    /* 0 -> */ {1},
    /* 1 -> */ {2},
    /* 2 -> */ {0}
  };
  ring_structure.SetStructure(ring_map);
  // ring_structure.Print(std::cout);
  // ring_structure.Print(std::cout, false);

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

TEST_CASE("Can get a random neighbor") {
  chemical_ecology::SpatialStructure structure;
  emp::Random rnd(5);

  structure.SetStructure(
    emp::vector< emp::vector<size_t> >{
      /* 0 -> */ {1, 2, 3},
      /* 1 -> */ {4},
      /* 2 -> */ {},
      /* 3 -> */ {},
      /* 4 -> */ {}
    }
  );

  for (size_t i = 0; i < 100; ++i) {
    const size_t neighbor = structure.GetRandomNeighbor(rnd, 0).value();
    REQUIRE(emp::Has(structure.GetNeighbors(0), neighbor));
  }

  auto result = structure.GetRandomNeighbor(rnd, 1);
  REQUIRE(result.value() == 4);

  result = structure.GetRandomNeighbor(rnd, 2);
  REQUIRE(!result);

}

TEST_CASE("Can read structure from csv that lists edges") {
  const std::string csv_path = "data/spatial-structure-edges.csv";

  chemical_ecology::SpatialStructure structure;

  structure.LoadStructureFromEdgeCSV(csv_path);
  // structure.Print(std::cout);
  REQUIRE(structure.GetNumPositions() == 5);
  REQUIRE(structure.GetNeighbors(0) == emp::vector<size_t>{1, 2, 3});
  REQUIRE(structure.GetNeighbors(1) == emp::vector<size_t>{0, 2});
  REQUIRE(structure.GetNeighbors(2) == emp::vector<size_t>{0});
  REQUIRE(structure.GetNeighbors(3) == emp::vector<size_t>{});
  REQUIRE(structure.GetNeighbors(4) == emp::vector<size_t>{});

}

TEST_CASE("Can read structure from matrix file") {
  const std::string mat_path = "data/spatial-structure-matrix.dat";

  chemical_ecology::SpatialStructure structure;

  structure.LoadStructureFromMatrix(mat_path);
  // structure.Print(std::cout);
  // structure.Print(std::cout, false);
  REQUIRE(structure.GetNumPositions() == 5);
  REQUIRE(structure.GetNeighbors(0) == emp::vector<size_t>{1, 2, 3});
  REQUIRE(structure.GetNeighbors(1) == emp::vector<size_t>{0, 2});
  REQUIRE(structure.GetNeighbors(2) == emp::vector<size_t>{0});
  REQUIRE(structure.GetNeighbors(3) == emp::vector<size_t>{});
  REQUIRE(structure.GetNeighbors(4) == emp::vector<size_t>{});

}

