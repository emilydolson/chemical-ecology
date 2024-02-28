#pragma once

#include <string>

#include "emp/config/config.hpp"

namespace chemical_ecology {

  EMP_BUILD_CONFIG(Config,
    GROUP(GLOBAL_SETTINGS, "Global settings"),
    VALUE(SEED, int, -1, "Seed for a simulation"),
    VALUE(N_TYPES, size_t, 9, "Number of types"),
    VALUE(UPDATES, size_t, 1000, "Number of time steps to run for"),
    VALUE(MAX_POP, int, 10000, "Maximum population size for one type in one cell"),
    VALUE(DIFFUSION, double, .25, "Proportion of each population that diffuses to adjacent cells"),
    VALUE(SEEDING_PROB, double, .25, "Probability that a member of a given species will randomly be introduced to a cell"),
    VALUE(PROB_CLEAR, double, .1, "Probability of cell being cleared out"),
    VALUE(REPRO_DILUTION, double, .1, "Proportion of contents to propogate on reproduction"),
    //Group repro enabled for all adaptive stochastic worlds under current architecture
    VALUE(GROUP_REPRO, bool, false, "True if proportional group level reproduction is enabled this run"),
    VALUE(V, bool, false, "True if running in verbose mode (Will print out all world vectors)"),
    VALUE(THRESHOLD_VALUE, double, 10.0, "Anything below threshold value will be rounded down for ranked threshold analysis"),

    GROUP(INTERACTION_SETTINGS, "Interaction Sources"),
    VALUE(INTERACTION_SOURCE, std::string, "", "Where to load interaction matrix from; empty string will generate randomly"),
    VALUE(INTERACTION_MAGNITUDE, double, 1, "Range of interaction intensities (from negative of this value to positive of this value"),
    VALUE(PROB_INTERACTION, double, .75, "Probability of there being an interaction between any given pair of species"),

    GROUP(SPATIAL_STRUCTURE_SETTINGS, "Settings related to connectivity of communities in the world"),
    VALUE(DIFFUSION_SPATIAL_STRUCTURE, std::string, "toroidal-grid", "Specifies spatial structure to use. Options:\n  'toroidal-grid' (2D toroidal grid)\n  'well-mixed' (all connected to all)\n  'load' (loads spatial structure from specified file)"),
    VALUE(DIFFUSION_SPATIAL_STRUCTURE_LOAD_MODE, std::string, "edges", "Species loading mode for spatial structure. Options:\n  'edges'\n  'matrix'"),
    VALUE(DIFFUSION_SPATIAL_STRUCTURE_FILE, std::string, "spatial-structure-edges.csv", "File to load spatial structure from. File format must be consistent with specified SPATIAL_STRUCTURE_LOAD_MODE"),
    VALUE(GROUP_REPRO_SPATIAL_STRUCTURE, std::string, "well-mixed", "Specifies spatial structure to use. Options:\n  'toroidal-grid' (2D toroidal grid)\n  'well-mixed' (all connected to all)\n  'load' (loads spatial structure from specified file)"),
    VALUE(GROUP_REPRO_SPATIAL_STRUCTURE_LOAD_MODE, std::string, "edges", "Species loading mode for spatial structure. Options:\n  'edges'\n  'matrix'"),
    VALUE(GROUP_REPRO_SPATIAL_STRUCTURE_FILE, std::string, "spatial-structure-edges.csv", "File to load spatial structure from. File format must be consistent with specified SPATIAL_STRUCTURE_LOAD_MODE"),
    VALUE(WORLD_WIDTH, size_t, 10, "Width of world. Used only for toroidal-grid and well-mixed spatial structure options."),
    VALUE(WORLD_HEIGHT, size_t, 10, "Height of world. Used only for toroidal-grid and well-mixed spatial structure options."),

    GROUP(ANALYSIS_SETTINGS, "Settings related to post-hoc model analyses"),
    VALUE(STOCHASTIC_ANALYSIS_REPS, size_t, 10, "Number of times to run post-hoc stochastic analyses"),
    VALUE(CELL_STABILIZATION_UPDATES, size_t, 10000, "Number of updates to run growth for cell stabilization"),
    VALUE(CELL_STABILIZATION_EPSILON, double, 0.0001, "If cell doesn't change more than this, can break stabilization early"),

    GROUP(OUTPUT_SETTINGS, "Settings related to data output"),
    VALUE(OUTPUT_DIR, std::string, "./output/", "What directory are we dumping data?"),
    VALUE(OUTPUT_RESOLUTION, size_t, 10, "How often should we output data?"),
    VALUE(RECORD_ASSEMBLY_MODEL, bool, false, "Should we output the assembly model updating over time?"),
    VALUE(RECORD_ADAPTIVE_MODEL, bool, false, "Should we output the adaptive model updating over time?"),
    VALUE(RECORD_A_ECO_DATA, bool, false, "Should we output a-eco_data?")
  );
}
