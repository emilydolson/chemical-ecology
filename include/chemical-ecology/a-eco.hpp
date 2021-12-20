#pragma once

#include <iostream>
#include "emp/Evolve/World.hpp"
#include "emp/math/distances.hpp"
#include "emp/data/DataNode.hpp"
#include "chemical-ecology/config_setup.hpp"

struct Particle {
  int type;
  double r = 1;
  int position;
};

struct CellData {
  double fitness;
  double heredity;
  double growth_rate;
  double equilib_growth_rate;
  int count = 1;
};

class AEcoWorld { //: public emp::World<Particle> {
  private:
    using world_t = emp::vector<emp::vector<int> >;
    emp::vector<emp::vector<double> > interactions;
    // emp::vector<double> rs;
    emp::Random rnd;
    chemical_ecology::Config * config = nullptr;
    int N_TYPES;
    int MAX_POP;
    int curr_update;

    emp::DataFile data_file;
    emp::DataNode<double, emp::data::Stats> fitness_node;
    emp::DataNode<double, emp::data::Stats> heredity_node;
    std::map<double, CellData> dom_map;

  public:
  AEcoWorld() : data_file("a-eco_data.csv") {;}
  world_t world;

  emp::vector<emp::vector<double> > GetInteractions() {
    return interactions;
  }

  void SetInteraction(int x, int y, double w) {
    interactions[x][y] = w;
  }

  void SetupRandomInteractions() {
    interactions.resize(N_TYPES);
    // rs.resize(N_TYPES);
    for (int i = 0; i < N_TYPES; i++) {
      // rs[i] = rnd.GetDouble(0,2);
      interactions[i].resize(N_TYPES);
      for (int j = 0; j < N_TYPES; j++) {
        if (rnd.P(config->PROB_INTERACTION())){
          interactions[i][j] = rnd.GetDouble(config->INTERACTION_MAGNITUDE() * -1, config->INTERACTION_MAGNITUDE());
        }
      }
      std::cout << emp::to_string(interactions[i]) << std::endl;
    }
    std::cout << std::endl;
  }

  void LoadInteractionMatrix(std::string filename) {
    emp::File infile(filename);
    emp::vector<emp::vector<double>> interaction_data = infile.ToData<double>();

    interactions.resize(N_TYPES);
    for (int i = 0; i < N_TYPES; i++) {
      interactions[i].resize(N_TYPES);
      for (int j = 0; j < N_TYPES; j++) {
        interactions[i][j] = interaction_data[i][j];
      }
      std::cout << emp::to_string(interactions[i]) << std::endl;
    }
    std::cout << std::endl;
  }

  void WriteInteractionMatrix(std::string filename) {
    emp::File outfile;
    
    for (int i = 0; i < N_TYPES; i++) {
      outfile += emp::join(interactions[i], ",");
    }
    
    outfile.Write(filename);
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

    if (config->INTERACTION_SOURCE() == "") {
      SetupRandomInteractions();
    } else {
      LoadInteractionMatrix(config->INTERACTION_SOURCE());
    }

    data_file.AddVar(curr_update, "Time", "Time");
    data_file.AddStats(fitness_node, "Fitness", "Fitness", true);
    data_file.AddStats(heredity_node, "Heredity", "Heredity", true);
    data_file.AddFun((std::function<int()>)[this](){return dom_map[fitness_node.GetMax()].count;}, "dominant_count", "dominant_count");
    data_file.AddFun((std::function<int()>)[this](){return dom_map[fitness_node.GetMax()].heredity;}, "dominant_heredity", "dominant_heredity");
    data_file.AddPreFun([this](){CalcAllFitness();});
    data_file.SetTimingRepeat(10);
    data_file.PrintHeaderKeys();

  }

  void Update(int ud) {

    world_t next_world;
    next_world.resize(config->WORLD_X() * config->WORLD_Y());
    for (emp::vector<int> & v : next_world) {
      v.resize(N_TYPES, 0);
    }

    for (size_t pos = 0; pos < world.size(); pos++) {
      DoGrowth(pos, world, next_world);
    }

    for (int pos = 0; pos < (int)world.size(); pos++) {
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

      emp::vector<int> adj = {up, down, left, right};
      DoRepro(pos, adj, world, next_world, config->SEEDING_PROB(), config->PROB_CLEAR());
      curr_update = ud;
    }
    data_file.Update(curr_update);
    std::swap(world, next_world);
    dom_map.clear();

  }

  void Run() {
    for (int i = 0; i < config->UPDATES(); i++) {
      Update(i);
      std::cout << i << std::endl;
    }

    for (auto & v : world) {
      std::cout << emp::to_string(v) << std::endl;
    }

    WriteInteractionMatrix("interaction_matrix.dat");
  }

  double CalcGrowthRate(size_t pos, world_t & curr_world) {
    double growth_rate = 0;
    for (int i = 0; i < N_TYPES; i++) {
      double modifier = 0;
      for (int j = 0; j < N_TYPES; j++) {
        modifier += interactions[i][j] * curr_world[pos][j];
      }

      growth_rate += ceil((interactions[i][i]+modifier)*((double)curr_world[pos][i]/(double)MAX_POP)); // * ((double)(MAX_POP - pos[i])/MAX_POP));
    }
    return growth_rate;
  }

  void DoGrowth(size_t pos, world_t & curr_world, world_t & next_world) {
    for (int i = 0; i < N_TYPES; i++) {
      double modifier = 0;
      for (int j = 0; j < N_TYPES; j++) {
        modifier += interactions[i][j] * curr_world[pos][j];
      }

      // Logistic growth
      int new_pop = ceil((interactions[i][i]+modifier)*curr_world[pos][i]); // * ((double)(MAX_POP - pos[i])/MAX_POP));
      // std::cout << pos[i] << " " << new_pop << " " << modifier << " " << ((double)(MAX_POP - pos[i])/MAX_POP) <<  std::endl;
      next_world[pos][i] = std::max(curr_world[pos][i] + new_pop, 0);
      next_world[pos][i] = std::min(next_world[pos][i], MAX_POP);
    }
    // std::cout << emp::to_string(pos) << std::endl;
  }

  bool IsReproReady(size_t pos, world_t & w) {
    return emp::Sum(w[pos]) > config->MAX_POP()*(config->REPRO_THRESHOLD());
  }

  void DoRepro(size_t pos, emp::vector<int> & adj, world_t & curr_world, world_t & next_world, double seed_prob, double prob_clear) {

    if (IsReproReady(pos, curr_world)) {
      int direction = adj[rnd.GetInt(adj.size())];
      for (int i = 0; i < N_TYPES; i++) {
        next_world[direction][i] += curr_world[pos][i] * config->REPRO_DILUTION();
      }
    }

    if (rnd.P(prob_clear)) {
      for (int i = 0; i < N_TYPES; i++) {
        next_world[pos][i] = 0;
      }
    }

    for (int direction : adj) {
      for (int i = 0; i < N_TYPES; i++) {
        double avail = curr_world[pos][i] * config->DIFFUSION();
        next_world[direction][i] +=  avail / 4;

        next_world[direction][i] = std::min(next_world[direction][i], MAX_POP);
        next_world[direction][i] = std::max(next_world[direction][i], 0);


        if (rnd.P(seed_prob)) {
          next_world[pos][i]++;
        }
      }
    }
  }

  void CalcAllFitness() {
    for (size_t pos = 0; pos < world.size(); pos++) {
      CellData res = GetFitness(pos);
      if (emp::Has(dom_map, res.equilib_growth_rate)) {
        dom_map[res.equilib_growth_rate].count++;
      } else {
        dom_map[res.equilib_growth_rate] = res;
      }
    }
  }

  CellData GetFitness(size_t orig_pos) {
    CellData data;
    int num_cells = 5;
    int time_limit = 250;
    world_t test_world;
    test_world.resize(num_cells);
    for (emp::vector<int> & v : test_world) {
      v.resize(config->N_TYPES(), 0);
    }

    test_world[0] = world[orig_pos];

    world_t next_world;
    next_world.resize(num_cells);
    for (emp::vector<int> & v : next_world) {
      v.resize(config->N_TYPES(), 0);
    }

    for (int i = 0; i < 20; i++) {
      DoGrowth(0, test_world, next_world);
      std::swap(test_world, next_world);
    }

    emp::vector<int> starting_point = test_world[0];
    data.growth_rate = CalcGrowthRate(orig_pos, world);
    data.equilib_growth_rate = CalcGrowthRate(0, test_world);

    fitness_node.Add(data.equilib_growth_rate);

    int time = 0;
    while (!IsReproReady(num_cells - 1, test_world) && time < time_limit) {

      for (int pos = 0; pos < num_cells; pos++) {
        DoGrowth(pos, test_world, next_world);
      }

      for (int pos = 0; pos < num_cells; pos++) {
        emp::vector<int> adj;
        if (pos > 0) {
          adj.push_back(pos - 1);
        }
        if (pos < num_cells - 1) {
          adj.push_back(pos + 1);
        }

        DoRepro(pos, adj, test_world, next_world, 0, 0);
      }

      std::swap(test_world, next_world);
      time++;
    }

    data.fitness = time;
    data.heredity = 1 - emp::EuclideanDistance(test_world[num_cells - 1], starting_point) / ((double)MAX_POP*N_TYPES);

    heredity_node.Add(data.heredity);

    return data;
  }

};