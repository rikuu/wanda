// Copyright 2017 Riku Walve

#include <vector>
#include <string>
#include <iostream>

#include "interval.h"
#include "graph.h"

void print_path(const graph_t &graph, const std::vector<interval_t> &path) {
  // TODO: Print with a lock

  if (path.size() == 0) return;

  std::string unitig = graph.label(path[0]);
  for (size_t i = 1; i < path.size(); i++) {
    unitig = graph.label(path[i])[0] + unitig;
  }

  std::cout << unitig << std::endl;
}

void compute_unitigs(const graph_t &graph, const size_t solid, const size_t min_length) {
  // TODO: Compute unitigs starting from nodes with rank in [i, j] in parallel

  const std::vector<interval_t> kmers = graph.distinct_kmers(solid);
  std::cerr << "[V::" << __func__ << "]: " << kmers.size() << " nodes" << std::endl;

  sdsl::bit_vector visited = sdsl::bit_vector(graph.rank(kmers.back()) + 1, false);

  size_t unitig_count = 0;
  for (size_t i = 0; i < kmers.size(); i++) {
    const interval_t node = kmers[i];
    const size_t rank = graph.rank(node);
    if (visited[rank]) continue;
    visited[rank] = true;

    const std::vector<interval_t> incoming = graph.incoming(node, solid);

    #ifdef DEBUG
      std::cerr << "[D::" << __func__ << "]: " <<
        "(" << node.left << ", " << node.right << ") = (" << graph.label(node) << "): " <<
        "in: " << incoming.size() << ", out: " << graph.outdegree(node, solid) << ", f: " << frequency(node) << std::endl;

      for (size_t j = 0; j < incoming.size(); j++) {
        std::cerr << "\t (" << incoming[j].left << ", " << incoming[j].right << ")" << std::endl;
      }
      std::cerr << std::endl;
    #endif

    // Maximal unitigs start from nodes with outdegree > 1
    if (incoming.size() == 1 && graph.outdegree(node, solid) > 1) {
      std::vector<interval_t> path;
      path.push_back(node);

      // Traverse graph backwards until a non-unary node or the starting node
      // is reached
      interval_t n = incoming[0];
      std::vector<interval_t> in = graph.incoming(n, solid);
      while (in.size() == 1 && n != node) {
        // Mark all visited nodes
        const size_t r = graph.rank(n);
        visited[r] = true;

        // Add node to path and get next node
        path.push_back(n);
        n = in[0];
        in = graph.incoming(n, solid);
      }

      // k + |v| - 1 = |path|
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
  if (argc != 4) {
    std::cerr << "Usage: " << argv[0] << " <graph prefix> <s> <min length>" << std::endl;
    return 1;
  }

  const std::string prefix = argv[1];
  const size_t solid = std::stoi(argv[2]);
  const size_t min_length = std::stoi(argv[3]);

  // Load graph
  const graph_t graph = graph_t::load(prefix);

  // Compute unitigs
  compute_unitigs(graph, solid, min_length);

  return 0;
}
