#pragma once

#include <iostream>
#include "emp/Evolve/World.hpp"
#include "chemical-ecology/config_setup.hpp"

struct Particle {
  int type;
  double r = 1;
  int position;
};


class AEcoWorld { //: public emp::World<Particle> {
  private:
    emp::vector<emp::vector<double> > interactions;
    double r = 1;
    emp::Random rnd;
    chemical_ecology::Config * config = nullptr;
    int N_TYPES;
    int MAX_POP;

  public:
  AEcoWorld(){;}
  emp::vector<emp::vector<int> > world;

  emp::vector<emp::vector<double> > GetInteractions() {
    return interactions;
  }

  void Setup(chemical_ecology::Config & cfg) {
    config = &cfg;

    N_TYPES = cfg.N_TYPES();
    MAX_POP = cfg.MAX_POP();

    rnd.ResetSeed(cfg.SEED());

    world.resize(cfg.WORLD_SIZE());
    for (emp::vector<int> & v : world) {
      v.resize(N_TYPES);
      for (int & count : v) {
        count = rnd.GetInt(100);
      }
    }

    interactions.resize(N_TYPES);
    for (int i = 0; i < N_TYPES; i++) {
      interactions[i].resize(N_TYPES);
      for (int j = 0; j < N_TYPES; j++) {
        if (i != j ){
          interactions[i][j] = rnd.GetDouble(cfg.INTERACTION_MAGNITUDE() * -1, cfg.INTERACTION_MAGNITUDE());
        }
      }
      std::cout << emp::to_string(interactions[i]) << std::endl;
    }
    std::cout << std::endl;
  }

  void Update() {
    // for (emp::vector<Particle> pos : world) {
    //   emp::vector<int> counts(interactions.size(), 0);
    //   for (Particle p : pos) {
    //     counts[p.type]++;
    //   }      
    // }

    for (emp::vector<int> & pos : world) {
      emp::vector<int> new_counts(N_TYPES);
      for (int i = 0; i < N_TYPES; i++) {
        double modifier = 0;
        for (int j = 0; j < N_TYPES; j++) {
          modifier += interactions[i][j] * pos[j];
        }
  
        // Logistic growth
        int new_pop = ceil((r+modifier)*pos[i]); // * ((double)(MAX_POP - pos[i])/MAX_POP));
        // std::cout << pos[i] << " " << new_pop << " " << modifier << " " << ((double)(MAX_POP - pos[i])/MAX_POP) <<  std::endl;
        new_counts[i] = std::max(pos[i] + new_pop, 0);
        new_counts[i] = std::min(new_counts[i], MAX_POP);
      }
      // std::cout << emp::to_string(pos) << std::endl;
      std::swap(pos, new_counts);
      std::cout << emp::to_string(pos) << std::endl;
    }

  }

  void Run() {
    for (int i = 0; i < config->UPDATES(); i++) {
      Update();
      std::cout << std::endl;
    }

    for (auto & v : world) {
      std::cout << emp::to_string(v) << std::endl;
    }
  }

};