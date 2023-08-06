#include <iostream>
#include <map>
#include <fstream>
#include <string>

#include "emp/base/vector.hpp"
#include "emp/datastructs/Graph.hpp"

#include "chemical-ecology/AEcoWorld.hpp"
#include "chemical-ecology/utils/config_setup.hpp"

void printPageRank(AEcoWorld world, emp::Graph g, ofstream& File)
{
  std::map<std::string, float> map = world.calculatePageRank(g);
  std::map<std::string, float>::iterator it = map.begin();
  while (it != map.end())
  {
      File << it->first << " " << it->second << std::endl;
      ++it;
  }
}

int main(int argc, char* argv[])
{
  // Read custom assembly graph
  std::ifstream graph_file("custom_graph.txt");
  std::string file_line;
  std::set<std::string> nodes;
  std::vector<std::vector<std::string>> node_map;
  if (graph_file.is_open())
  {
    while (graph_file)
    {
      std::getline(graph_file, file_line);
      std::string node0 = file_line.substr(0, file_line.find(" "));
      std::string node1 = file_line.substr(file_line.find(" ")+1, file_line.size());

      std::vector v = {node0, node1};
      node_map.push_back(v);
      nodes.insert(node0);
      nodes.insert(node1);
    }
  }

  // Turn graph into Empirical graph
  emp::Graph g(nodes.size());
  for(auto x : nodes)
  {
    std::size_t xint = std::stoi(x);
    g.SetLabel(xint, x);
  }

  for(int i = 0; i < node_map.size(); i++)
  {
    std::string node0 = node_map[i][0];
    std::string node1 = node_map[i][1];
    std::size_t node0int = std::stoi(node0);
    std::size_t node1int = std::stoi(node1);
    std::cout << node0int << " " << node1int << std::endl;
    g.AddEdge(node0int, node1int);
  }

  AEcoWorld world;
  ofstream PageRankFile("custom_pagerank.txt");
  printPageRank(world, g, PageRankFile);
  PageRankFile.close();

  return 0;
}