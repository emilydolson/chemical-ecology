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

    N_TYPES = config->N_TYPES();
    MAX_POP = config->MAX_POP();

    rnd.ResetSeed(config->SEED());

    world.resize(config->WORLD_X() * config->WORLD_Y());
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
          interactions[i][j] = rnd.GetDouble(config->INTERACTION_MAGNITUDE() * -1, config->INTERACTION_MAGNITUDE());
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

    emp::vector<emp::vector<int> > next_world;
    next_world.resize(config->WORLD_X() * config->WORLD_Y());
    for (emp::vector<int> & v : next_world) {
      v.resize(N_TYPES, 0);
    }

    for (int pos = 0; pos < world.size(); pos++) {
      for (int i = 0; i < N_TYPES; i++) {
        double modifier = 0;
        for (int j = 0; j < N_TYPES; j++) {
          modifier += interactions[i][j] * world[pos][j];
        }
  
        // Logistic growth
        int new_pop = ceil((r+modifier)*world[pos][i]); // * ((double)(MAX_POP - pos[i])/MAX_POP));
        // std::cout << pos[i] << " " << new_pop << " " << modifier << " " << ((double)(MAX_POP - pos[i])/MAX_POP) <<  std::endl;
        next_world[pos][i] = std::max(world[pos][i] + new_pop, 0);
        next_world[pos][i] = std::min(next_world[pos][i], MAX_POP);
      }
      // std::cout << emp::to_string(pos) << std::endl;
    }

    for (int pos = 0; pos < world.size(); pos++) {
      int x = pos % config->WORLD_X();
      int y = pos / config->WORLD_Y();
      int left = pos - 1;     
      int right = pos + 1;

      if (x == 0) {
        left += config->WORLD_X(); 
      } else if (x == config->WORLD_X() - 1) {
        right -= config->WORLD_X();
      }
 
      int up = pos - config->WORLD_X();
      int down = pos + config->WORLD_X();

      if (y == 0) {
        up += config->WORLD_X() * (config->WORLD_Y());
      } else if (y == config->WORLD_Y() - 1) {
        down -= config->WORLD_X() * (config->WORLD_Y());
      }

      for (int i = 0; i < N_TYPES; i++) {
        double avail = world[pos][i] * config->DIFFUSION();
        next_world[left][i] +=  avail / 4;
        next_world[right][i] +=  avail / 4;
        next_world[up][i] +=  avail / 4;
        next_world[down][i] +=  avail / 4;
        
        if (rnd.P(config->SEEDING_PROB())) {
          next_world[pos][i]++;
        }
      }
    }

    std::swap(world, next_world);

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