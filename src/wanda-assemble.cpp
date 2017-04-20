// Copyright 2017 Riku Walve

#include <vector>
#include <string>
#include <iostream>

#include "interval.h"
#include "graph.h"

#define frequency(n) ((n.right - n.left) + 1)

void print_path(const graph_t &graph, const std::vector<interval_t> &path) {
  if (path.size() == 0) return;

  std::string unitig = graph.label(path[0]);
  const size_t k = unitig.length();
  for (size_t i = 1; i < path.size(); i++) {
    unitig += graph.label(path[i]).substr(k-1);
  }

  std::cout << unitig << std::endl;
}

void compute_unitigs(const graph_t &graph, const size_t solid, const size_t min_length) {
  const std::vector<interval_t> kmers = graph.distinct_kmers(solid);
  std::cerr << "[V::" << __func__ << "]: " << kmers.size() << " nodes" << std::endl;

  sdsl::bit_vector visited = sdsl::bit_vector(graph.rank(kmers.back()) + 1, false);

  size_t unitig_count = 0;
  for (size_t i = 0; i < kmers.size(); i++) {
    const interval_t node = kmers[i];
    const size_t rank = graph.rank(node);
    if (visited[rank]) continue;
    visited[rank] = true;

    const std::vector<interval_t> outgoing = graph.outgoing(node);

    #ifdef DEBUG
      std::cerr << "[D::" << __func__ << "]: " <<
        "(" << node.left << ", " << node.right << "): " <<
        outgoing.size() << ", " << graph.indegree(node) << ", " << frequency(node) << std::endl;

      for (size_t j = 0; j < outgoing.size(); j++) {
        std::cerr << "\t (" << outgoing[j].left << ", " << outgoing[j].right << ")" << std::endl;
      }
      std::cerr << std::endl;
    #endif

    if (outgoing.size() == 1 && graph.indegree(node) == 1) {
      std::vector<interval_t> path;
      path.push_back(node);

      interval_t n = outgoing[0];
      std::vector<interval_t> out = graph.outgoing(n);
      while (out.size() == 1 && frequency(n) >= solid) {
        const size_t r = graph.rank(n);
        if (visited[r]) break;
        visited[r] = true;

        path.push_back(n);
        n = out[0];
        out = graph.outgoing(n);
      }

      if ((graph.k() + path.size() - 1) >= min_length) {
        std::cout << ">contig" << unitig_count << std::endl;
        print_path(graph, path);
        unitig_count++;
      }
    }
  }

  std::cerr << "[V::" << __func__ << "]: " << unitig_count << " unitigs" << std::endl;
}

int main(int argc, char* argv[]) {
  if (argc != 5) {
    std::cerr << "Usage: " << argv[0] << " <graph prefix> <k> <s> <min length>" << std::endl;
    return 1;
  }

  const std::string prefix = argv[1];
  const size_t k = std::stoi(argv[2]);
  const size_t solid = std::stoi(argv[3]);
  const size_t min_length = std::stoi(argv[4]);

  // Load graph
  const graph_t graph = graph_t::load(prefix, k);

  // Compute unitigs
  compute_unitigs(graph, solid, min_length);

  return 0;
}
