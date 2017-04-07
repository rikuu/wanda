
#include <vector>
#include <string>

#include <sdsl/bit_vectors.hpp>
#include <sdsl/lcp.hpp>

#include "index.h"
#include "interval.h"
#include "graph.h"

// TODO: Use a more space efficient construction based on FM-index
// TODO: Possible to remove "false-positive" nodes/edges, i.e. containing '$'?
void graph_t::build_first(const std::string &kernel_filename) {
  sdsl::lcp_wt<> lcp;
  sdsl::construct(lcp, kernel_filename.c_str(), 1);

  sdsl::bit_vector first = sdsl::bit_vector(m_index.size() - 1, false);
  for (size_t i = 0; i < m_index.size()-1; i++) {
    if (lcp[i+2] < m_k) {
      first[i] = true;
    }
  }

  sdsl::util::clear(lcp);

  m_first = sdsl::rrr_vector<127>(first);
}

std::string graph_t::label(const interval_t &node) const {
  interval_t interval = node;
  uint8_t c = '\0';
  for (size_t i = 0; i < m_k; i++) {
    std::cout << i << ", " << m_k << ", " << static_cast<char>(c) << ", " << interval.left << ", " << interval.right << std::endl;

    interval = m_index.inverse_lf(interval, &c);
    m_buffer[m_k - i - 1] = static_cast<char>(c);

    if (c == '$') {
      return "";
    }
  }

  return std::string(m_buffer, m_k);
}

std::vector<interval_t> graph_t::distinct_kmers(const size_t frequency) const {
  std::vector<interval_t> kmers;
  for (size_t i = 1; i <= m_first_rs.rank(m_first.size()); i++) {
    const size_t start = m_first_ss.select(i);
    const size_t end = m_first_ss.select(i + 1) - 1;

    if (end - start >= frequency) {
      kmers.push_back(interval_t(start, end));
    }
  }
  return kmers;
}

interval_t graph_t::follow_edge(const interval_t &node, const uint8_t c) const {
  // First find the inteval e corresponding to c1 .. ck+1
  const interval_t e = m_index.extend(node, c);

  // std::cout << "n: " << node.left << ", " << node.right << std::endl;
  // std::cout << "c: " << c << ", " << c1 << std::endl;
  // std::cout << "c1: " << start_e << ", " << end_e << std::endl;

  if (m_first[e.left] && (e.right == m_index.size()+1 || m_first[e.right+1]))
    return interval_t(e.left, e.right);

  // Take interval inside it corresponding to c2 .. ck+1
  const size_t start = m_first[e.left] ? e.left : m_first_ss.select(m_first_rs.rank(e.left));
  const size_t end = (e.right == m_index.size()) ? e.right : m_first_ss.select(m_first_rs.rank(e.right) + 1) - 1;

  // std::cout << "c2: " << start << ", " << end << std::endl;

  return interval_t(start, end);
}

std::vector<interval_t> graph_t::outgoing(const interval_t &node) const {
  const std::vector<uint8_t> edges = m_index.interval_symbols(node.left, node.right);

  std::vector<interval_t> nodes;
  for (size_t i = 0; i < edges.size(); i++) {
    nodes.push_back(follow_edge(node, edges[i]));
  }

  return nodes;
}

std::vector<interval_t> graph_t::incoming(const interval_t &node) const {
  // std::cout << std::endl;
  std::vector<interval_t> edges;
  for (size_t i = node.left; i <= node.right; i++) {
    const size_t ilf = m_index.inverse_lf(i);
    if (ilf == 0) continue;

    const size_t rank = m_first_rs.rank(ilf+1);

    const size_t left = m_first_ss.select(rank);
    const size_t right = m_first_ss.select(rank + 1) - 1;

    // std::cout << i << ", " << index[i] << ", " << c_array[index[i]]
    //   << ", "<< ilf << ", " << rank <<
    //   ", " << left << ", " << right << std::endl;

    const interval_t n = interval_t(left, right);
    bool in = false;
    for (size_t j = 0; j < edges.size(); j++) {
      if (edges[j] == n) in = true;
    }
    if (!in) edges.push_back(n);

    // i += right - ilf + 1;
  }

  return edges;
}
