#pragma once

#include <string>

#include "emp/config/config.hpp"

namespace chemical_ecology {

  EMP_BUILD_CONFIG(Config,
    GROUP(GLOBAL_SETTINGS, "Global settings"),
    VALUE(SEED, int, -1, "Seed for a simulation"),
    VALUE(N_TYPES, int, 9, "Number of types"),
    VALUE(WORLD_X, int, 10, "Width world"),
    VALUE(WORLD_Y, int, 10, "Height world"),
    VALUE(UPDATES, int, 1000, "Number of time steps to run for"),
    VALUE(MAX_POP, int, 10000, "Maximum population size for one type in one cell"),
    VALUE(DIFFUSION, double, 1, "Proportion of each population that diffuses to adjacent cells"),
    VALUE(SEEDING_PROB, double, .001, "Probability that a member of a given species will randomly be introduced to a cell"),
    VALUE(PROB_CLEAR, double, .25, "Probability of cell being cleared out"),
    VALUE(REPRO_DILUTION, double, .1, "Proportion of contents to propogate on reproduction"),
    //Group repro enabled for all adaptive stochastic worlds under current architecture
    //VALUE(GROUP_REPRO, bool, false, "True if proportional group level reproduction is enabled this run"),
    VALUE(V, bool, false, "True if running in verbose mode (Will print out all world vectors)"),

    GROUP(INTERACTION_SETTINGS, "Interaction Sources"),
    VALUE(INTERACTION_SOURCE, std::string, "config/proof_of_concept_interaction_matrix.dat", "Where to load interaction matrix from; empty string will generate randomly"),
    VALUE(INTERACTION_MAGNITUDE, double, 1, "Range of interaction intensities (from negative of this value to positive of this value"),
    VALUE(PROB_INTERACTION, double, .1, "Probability of there being an interaction between any given pair of species"),

    GROUP(SPATIAL_STRUCTURE_SETTINGS, "Settings related to connectivity of communities in the world"),
    VALUE(SPATIAL_STRUCTURE, std::string, "default", "Specifies spatial structure to use. Options:\n  'default' (toroidal grid)\n  'well-mixed' (all connected to all)\n  'load' (loads spatial structure from specified file)"),
    VALUE(SPATIAL_STRUCTURE_LOAD_MODE, std::string, "edges", "Species loading mode for spatial structure. Options:\n  'edges'\n  'matrix'"),
    VALUE(SPATIAL_STRUCTURE_FILE, std::string, "spatial-structure-edges.csv", "File to load spatial structure from. File format must be consistent with specified SPATIAL_STRUCTURE_LOAD_MODE"),

    GROUP(OUTPUT_SETTINGS, "Settings related to data output"),
    VALUE(OUTPUT_DIR, std::string, "./output/", "What directory are we dumping data?"),
    VALUE(OUTPUT_RESOLUTION, size_t, 10, "How often should we output data?")

  );
}
