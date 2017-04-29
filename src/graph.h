// Copyright 2017 Riku Walve

#ifndef WANDA_GRAPH_H_
#define WANDA_GRAPH_H_

#include <vector>
#include <string>

#include <sdsl/bit_vectors.hpp>

#include "index.h"
#include "interval.h"

#define MARKER '$'

#define frequency(n) ((n.right - n.left) + 1)

class graph_t {
public:
  graph_t(const std::string &kernel_filename, const size_t k) :
      m_k(k), m_index(index_t(kernel_filename)) {
    m_buffer = new char[k];

    m_first = build_first(kernel_filename, k);
    m_first_ss = sdsl::select_support_rrr<1, 127>(&m_first);
    m_first_rs = sdsl::rank_support_rrr<1, 127>(&m_first);
  }

  graph_t(const size_t k, const index_t index, const sdsl::rrr_vector<127> first) :
      m_k(k), m_index(index), m_first(first) {
    m_buffer = new char[k];

    m_first_ss = sdsl::select_support_rrr<1, 127>(&m_first);
    m_first_rs = sdsl::rank_support_rrr<1, 127>(&m_first);
  }

  // Copy constructor
  graph_t(const graph_t& graph) :
      m_k(graph.m_k), m_index(graph.m_index), m_first(graph.m_first) {
    m_buffer = new char[m_k];
  }

  // Move constructor
  graph_t(graph_t&& graph) noexcept :
      m_k(graph.m_k), m_index(graph.m_index), m_first(graph.m_first) {
    m_buffer = new char[m_k];

    // Invalidate other graph here
    // sdsl::util::clear(graph.m_first);
  }

  // Copy assignment operator
  graph_t& operator=(const graph_t& graph) {
    graph_t tmp(graph);
    *this = std::move(tmp);
    return *this;
  }

  // Move assignment operator
  graph_t& operator=(graph_t&& graph) noexcept {
    m_k = graph.m_k;
    m_index = graph.m_index;
    m_first = graph.m_first;

    // Invalidate other graph here

    return *this;
  }

  ~graph_t() noexcept {
    delete[] m_buffer;
  }

  // Loads a graph from a file
  static graph_t load(const std::string &base) {
    index_t index = index_t::load(base);

    size_t k;
    sdsl::rrr_vector<127> first;

    std::ifstream in(base + ".first");
    if (in.good()) {
      in.read(reinterpret_cast<char*>(&k), sizeof(k));
      first.load(in);
    } else {
      std::cerr << "[E::" << __func__ << "]: Unable to read \"" << base + ".first\"!" << std::endl;
      exit(1);
    }

    #ifdef DEBUG
      std::cerr << "[D::" << __func__ << "]: k = " << k << std::endl;
      std::cerr << "[D::" << __func__ << "]: first = ";
      for (size_t i = 0; i < first.size(); i++)
        std::cerr << first[i];
      std::cerr << std::endl;
    #endif

    return graph_t(k, index, first);
  }

  // Stores the graph to a file
  void store_to_file(const std::string &base) const {
    m_index.store_to_file(base);

    std::ofstream out(base + ".first");
    if (out.good()) {
      out.write(reinterpret_cast<const char*>(&m_k), sizeof(m_k));
      m_first.serialize(out);
    } else {
      std::cerr << "[E::" << __func__ << "]: Unable to write to \"" << base + ".first\"!" << std::endl;
    }
  }

  void change_k(const size_t k) {
    if (k == m_k) {
      return;
    }

    std::cerr << "[E::" << __func__ << "]: Not yet implemented!" << std::endl;
    exit(1);

    m_k = k;

    // TODO: Implement build_first() without text
    m_first = build_first("", k);
    m_first_ss = sdsl::select_support_rrr<1, 127>(&m_first);
    m_first_rs = sdsl::rank_support_rrr<1, 127>(&m_first);
  }

  inline size_t size() const {
    return m_index.size();
  }

  // All distinct nodes with a minimum frequency
  std::vector<interval_t> distinct_kmers(const size_t solid = 0) const;

  // Returns the label of a node (i.e. the "content" of the corresponding kmer)
  std::string label(const interval_t &node) const;

  // Returns all the nodes which have an outgoing edge to a node
  std::vector<interval_t> incoming(const interval_t &node, const size_t solid = 0) const;

  // Returns all the nodes which have an incoming edge from a node
  std::vector<interval_t> outgoing(const interval_t &node, const size_t solid = 0) const;

  // The in-degree of a node
  inline size_t outdegree(const interval_t &node, const size_t solid = 0) const {
    return outgoing(node, solid).size();
  }

  // The out-degree of a node
  inline size_t indegree(const interval_t &node, const size_t solid = 0) const {
    return incoming(node, solid).size();
    // const std::vector<uint8_t> symbols = m_index.interval_symbols(node.left, node.right);
    //
    // size_t count = 0;
    // for (size_t i = 0; i < symbols.size(); i++) {
    //   if (symbols[i] != 0 && symbols[i] != MARKER) {
    //     count++;
    //   }
    // }
    //
    // return count;
  }

  // Returns all occurrences of a kmer in the text
  std::vector<size_t> occurrences(const interval_t &node) const {
    std::vector<size_t> occurrences;
    for (size_t j = node.left; j <= node.right; j++) {
      occurrences.push_back(m_index.sa(j));
    }
    return occurrences;
  }

  inline size_t rank(const interval_t &node) const {
    return m_first_rs.rank(node.left);
  }

  inline size_t k() const {
    return m_k;
  }

private:
  // Follows an edge in the graph from a node to a node
  interval_t follow_edge(const interval_t &node, uint8_t c) const;

  static sdsl::rrr_vector<127> build_first(const std::string &filename, const size_t k);

private:
  size_t m_k;

  // FM-m_index
  index_t m_index;

  // Bitvector marking starting positions for k-mers
  sdsl::rrr_vector<127> m_first;
  sdsl::select_support_rrr<1, 127> m_first_ss;
  sdsl::rank_support_rrr<1, 127> m_first_rs;

  // Buffer for kmer labels
  char *m_buffer;
};

#endif
