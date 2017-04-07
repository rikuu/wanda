// Copyright 2017 Riku Walve

#include <vector>
#include <string>
#include <iostream>

#include "graph.h"

int main(int argc, char* argv[]) {
  if (argc != 4) {
    std::cerr << "Usage: " << argv[0] << " <stream> <k> <graph prefix>" << std::endl;
    return 1;
  }

  const std::string in = argv[1];
  const size_t k = std::stoi(argv[2]);
  const std::string prefix = argv[3];

  // Construct graph
  const graph_t graph(in, k);

  // Save graph to file
  graph.store_to_file(prefix);

  return 0;
}
