// Copyright 2017 Riku Walve

#ifndef WANDA_INDEX_H_
#define WANDA_INDEX_H_

#include <vector>

#include <sdsl/bit_vectors.hpp>
#include <sdsl/wavelet_trees.hpp>

#include "interval.h"

#define SA_SAMPLE_DENSITY 32

class index_t {
public:
  index_t(const std::string &kernel_filename);

  index_t(const sdsl::wt_huff<sdsl::rrr_vector<127> > &tree, const sdsl::int_vector<> &sa_samples) :
      m_tree(tree), m_sa_samples(sa_samples) {
    build_c_array();
  }

  static index_t load(const std::string &base) {
    sdsl::wt_huff<sdsl::rrr_vector<127> > tree;
    sdsl::int_vector<> sa_samples;

    sdsl::load_from_file(tree, base + ".bwt");
    sdsl::load_from_file(sa_samples, base + ".sa");

    #ifdef DEBUG
      std::cerr << "[D::" << __func__ << "]: ";
      for (size_t i = 0; i < tree.size(); i++)
        std::cerr << tree[i];
      std::cerr << std::endl;

      std::cerr << "[D::" << __func__ << "]: ";
      for (size_t i = 0; i < tree.size(); i++)
        std::cerr << i;
      std::cerr << std::endl;
    #endif

    return index_t(tree, sa_samples);
  }

  void store_to_file(const std::string &base) const {
    sdsl::store_to_file(m_tree, base + ".bwt");
    sdsl::store_to_file(m_sa_samples, base + ".sa");
  }

  inline size_t size() const {
    return m_tree.size();
  }

  size_t sa(const size_t position) const {
    size_t i = position, offset = 0;
    while ((i % SA_SAMPLE_DENSITY) != 0) {
      i = lf(i);
      offset++;
    }

    const size_t result = m_sa_samples[i / SA_SAMPLE_DENSITY] + offset;
    return (result < m_tree.size()) ? result : result - m_tree.size();
  }

  std::vector<uint8_t> interval_symbols(const size_t left, const size_t right) const {
    if (left == right) {
      std::vector<uint8_t> alphabet = { m_tree[left] };
      return alphabet;
    }

    sdsl::int_vector_size_type extensions;
    std::vector<uint8_t> alphabet(m_tree.sigma);
    std::vector<uint64_t> ranks_i(m_tree.sigma);
    std::vector<uint64_t> ranks_j(m_tree.sigma);

    m_tree.interval_symbols(left, right + 1, extensions, alphabet, ranks_i, ranks_j);

    while (alphabet.size() > extensions)
      alphabet.pop_back();

    return alphabet;
  }

  interval_t extend(const interval_t &interval, const uint8_t c) const {
    const size_t c1 = m_c_array[c];
    const size_t left = interval.left > 0 ? c1 + m_tree.rank(interval.left - 1, c) + 1 : c1 + 1;
    const size_t right = c1 + m_tree.rank(interval.right, c);
    return interval_t(left, right);
  }

  inline size_t lf(const size_t i) const {
    const uint8_t c = m_tree[i];

    if (m_c_array[c] == 0) {
      return 0;
    }

    return m_c_array[c] + m_tree.rank(i, c);
  }

  size_t inverse_lf(const size_t i, uint8_t *_c = nullptr) const {
    uint8_t c = '\0';
    for (size_t j = 1; j < m_c_array.size(); j++) {
      if (m_c_array[j] > 0 && m_c_array[j] <= i)
        c = static_cast<uint8_t>(j);
    }
    if (_c != nullptr) *_c = c;

    if (m_c_array[c] == 0)
      return 0;

    return m_tree.select(i - m_c_array[c] + 1, c);
  }

  interval_t inverse_lf(const interval_t &interval, uint8_t *_c = nullptr) const {
    uint8_t c;
    const size_t start = inverse_lf(interval.left, &c);

    if (m_c_array[c] == 0)
      return interval_t(0, 0);

    #ifdef DEBUG
      std::cout << "[D::" << __func__ << "]: " <<
        "(" << interval.left << ", " << interval.right << "), " <<
        static_cast<char>(c) << ", " << m_c_array[c] << std::endl;
    #endif

    const size_t end = m_tree.select(interval.right - m_c_array[c] + 1, c);

    if (_c != nullptr) *_c = c;
    return interval_t(start, end);
  }

private:
  void build_c_array() {
    std::vector<uint8_t> alphabet = this->interval_symbols(0, m_tree.size()-1);
    std::sort(alphabet.begin(), alphabet.end());
    assert(alphabet.size() != 0);

    std::vector<size_t> counts(256);
    counts[alphabet[0]] = 0;
    for (size_t i = 1; i < alphabet.size(); i++) {
      const size_t count_prev = m_tree.rank(m_tree.size(), alphabet[i-1]) - m_tree.rank(0, alphabet[i-1]);
      counts[alphabet[i]] = counts[alphabet[i-1]] + count_prev;
    }

    m_c_array = counts;
  }

private:
  sdsl::wt_huff<sdsl::rrr_vector<127> > m_tree;
  sdsl::int_vector<> m_sa_samples;
  std::vector<size_t> m_c_array;
};

#endif
