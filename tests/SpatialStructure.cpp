#define CATCH_CONFIG_MAIN

#include "Catch/single_include/catch2/catch.hpp"

#include <iostream>

#include "chemical-ecology/SpatialStructure.hpp"

TEST_CASE("Test SetStructure") {

  chemical_ecology::SpatialStructure spatial_structure;

  spatial_structure.Print(std::cout);
}