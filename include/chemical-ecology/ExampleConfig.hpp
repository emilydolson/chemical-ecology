#pragma once

#include "emp/config/config.hpp"

namespace chemical_ecology {

  EMP_BUILD_CONFIG(Config,
    GROUP(GLOBAL_SETTINGS, "Global settings"),
    VALUE(SEED, int, -1, "Seed for a simulation"),
    VALUE(N_TYPES, int, 9, "Number of types"),
    VALUE(WORLD_X, int, 10, "Width world"),
    VALUE(WORLD_Y, int, 10, "Height world"),
    VALUE(UPDATES, int, 10000, "Number of time steps to run for"),
    VALUE(MAX_POP, int, 10000, "Maximum population size for one type in one cell"),
    VALUE(DIFFUSION, double, 1, "Proportion of each population that diffuses to adjacent cells"),
    VALUE(SEEDING_PROB, double, .001, "Probability that a member of a given species will randomly be introduced to a cell"),
    VALUE(REPRO_THRESHOLD, double, 1000000, "Total population required for reproduction; expressed as a multiplier on MAX_POP"),
    VALUE(REPRO_DILUTION, double, .1, "Proportion of contents to propogate on reproduction"),
    VALUE(PROB_CLEAR, double, .25, "Probability of cell being cleared out"),    

    GROUP(INTERACTION_SETTINGS, "Interaction Sources"),
    VALUE(INTERACTION_SOURCE, std::string, "config/proof_of_concept_interaction_matrix.dat", "Where to load interaction matrix from; empty string will generate randomly"),    
    VALUE(INTERACTION_MAGNITUDE, double, 1, "Range of interaction intensities (from negative of this value to positive of this value"),
    VALUE(PROB_INTERACTION, double, .1, "Probability of there being an interaction between any given pair of species"),

  );
}
