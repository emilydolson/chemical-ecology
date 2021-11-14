#pragma once

#include "emp/config/config.hpp"

namespace chemical_ecology {

  EMP_BUILD_CONFIG(Config,
    GROUP(GLOBAL_SETTINGS, "Global settings"),
    VALUE(SEED, int, -1, "Seed for a simulation"),
    VALUE(N_TYPES, int, 100, "Number of types"),
    VALUE(WORLD_X, int, 10, "Width world"),
    VALUE(WORLD_Y, int, 10, "Height world"),
    VALUE(UPDATES, int, 1000, "Number of time steps to run for"),
    VALUE(MAX_POP, int, 10000, "Maximum population size for one type in one cell"),
    VALUE(INTERACTION_MAGNITUDE, double, 1, "Range of interaction intensities (from negative of this value to positive of this value"),
    VALUE(DIFFUSION, double, .01, "Proportion of each population that diffuses to adjacent cells"),
    VALUE(SEEDING_PROB, double, .01, "Probability that a member of a given species will randomly be introduced to a cell"),
    VALUE(REPRO_THRESHOLD, double, 5, "Total population required for reproduction; expressed as a multiplier on MAX_POP"),

    // GROUP(OTHER_SETTINGS, "Miscellaneous settings"),
    // VALUE(LUNCH_ORDER, std::string, "ham on five", "What's for lunch today"),
    // VALUE(HOLD_MAYO, bool, true, "Whether or not to hold the mayo"),
    // VALUE(SUPER_SECRET, bool, true, "It's a hidden switch"),
  );
}
