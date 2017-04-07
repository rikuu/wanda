// Copyright 2017 Riku Walve

#include <vector>
#include <string>
#include <iostream>

#include "interval.h"
#include "graph.h"

std::vector<std::vector<interval_t>> compute_unitigs(const graph_t &graph) {
  std::vector<std::vector<interval_t>> paths;

  const std::vector<interval_t> kmers = graph.distinct_kmers(0);
  for (size_t i = 0; i < kmers.size(); i++) {
    const interval_t node = kmers[i];
    std::cout << node.left << ", " << node.right << std::endl;

    const std::vector<interval_t> outgoing = graph.outgoing(node);
    if (outgoing.size() == 1 && graph.indegree(node) == 1) {
      std::vector<interval_t> path;
      path.push_back(node);

      interval_t n = outgoing[0];
      std::vector<interval_t> out = graph.outgoing(n);
      while (out.size() == 1) {
        path.push_back(n);
        n = out[0];
        out = graph.outgoing(n);
      }

      paths.push_back(path);
    }
  }

  return paths;
}

int main(int argc, char* argv[]) {
  if (argc != 3) {
    std::cerr << "Usage: " << argv[0] << " <graph prefix> <unitigs>" << std::endl;
    return 1;
  }

  const std::string prefix = argv[1];
  const std::string out = argv[2];

  // Load graph
  const graph_t graph = graph_t::load(prefix, 2);
  const size_t k = 2;

  const std::vector<interval_t> kmers = graph.distinct_kmers(0);
  for (size_t i = 0; i < kmers.size(); i++) {
    std::cout << graph.label(kmers[i]) << std::endl;
  }

  const auto unitigs = compute_unitigs(graph);
  for (size_t i = 0; i < unitigs.size(); i++) {
    std::string unitig = graph.label(unitigs[i][0]);
    for (size_t j = 1; j < unitigs[i].size(); j++) {
      unitig += graph.label(unitigs[i][j]).substr(k-1);
    }
    std::cout << unitig << std::endl;
  }

  return 0;
}
