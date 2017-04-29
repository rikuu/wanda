// Copyright 2017 Riku Walve

#include <vector>

#include <sdsl/bit_vectors.hpp>
#include <sdsl/wavelet_trees.hpp>

#include "index.h"
#include "interval.h"

static inline size_t filelength(FILE * fp) {
  fseek(fp, 0, SEEK_END);
  size_t file_len = static_cast<size_t>(ftell(fp));
  fseek(fp, 0, SEEK_SET);
  return file_len;
}

// Construct BWT from input file and sa5 file
static void stream_bwt(const std::string &input, const std::string &suffix,
    const std::string &bwt, sdsl::int_vector<> *samples) {
  // Load input file to memory
  FILE *in = fopen(input.c_str(), "r");
  const size_t n = filelength(in);
  char *in_buffer = new char[n];
  fread(in_buffer, sizeof(char), n, in);
  fclose(in);

  // Open sa5 file for reading
  FILE *suf = fopen(suffix.c_str(), "rb");
  uint8_t *buffer = new uint8_t[5];

  std::ofstream bwt_stream(bwt.c_str());
  for (size_t i = 0; i < n; i++) {
    // TODO: Read 5*N characters at a time
    fread(buffer, sizeof(char), 5, suf);

    // Shift the char into position
    int64_t sa = static_cast<int64_t>(buffer[0]);
    sa |= (static_cast<int64_t>(buffer[1])) << 8;
    sa |= (static_cast<int64_t>(buffer[2])) << 16;
    sa |= (static_cast<int64_t>(buffer[3])) << 24;
    sa |= (static_cast<int64_t>(buffer[4])) << 32;

    // TODO: Not load to memory, but still random accesses?
    bwt_stream << (in_buffer[sa == 0 ? n - 1 : sa - 1]);

    if ((i % SA_SAMPLE_DENSITY) == 0) {
      (*samples)[i / SA_SAMPLE_DENSITY] = sa;
    }
  }

  bwt_stream.close();

  fclose(suf);

  delete[] buffer;
  delete[] in_buffer;
}

index_t::index_t(const std::string &kernel_filename) {
  const std::string suffix_filename = kernel_filename + ".sa5";
  const std::string bwt_filename = kernel_filename + ".bwt";

  // Construct suffix array with pSAscan
  // TODO: pSAscan public API
  const std::string sa_construct_command = std::string(PROJECT_ROOT) + std::string("/ext/pSAscan/src/psascan ") + kernel_filename;
  std::cout << "We wil call:\n" << sa_construct_command << std::endl;
  if (system(sa_construct_command.c_str())) {
    std::cout << "Command failed. " << std::endl;
    exit(-1);
  }

  // !!!
  FILE *in = fopen(kernel_filename.c_str(), "r");
  const size_t n = filelength(in);
  fclose(in);

  const size_t num_of_samples = n / SA_SAMPLE_DENSITY;
  m_sa_samples = sdsl::int_vector<>(num_of_samples + 1, 0, 64);

  stream_bwt(kernel_filename, suffix_filename, bwt_filename, &m_sa_samples);
  sdsl::construct(m_tree, bwt_filename, 1);

  build_c_array();

  // Delete SA and BWT files
  const std::string delete_command = "rm -f " + suffix_filename + " " + bwt_filename;
  std::cout << "We wil call:\n" <<  delete_command << std::endl;
  if (system(delete_command.c_str())) {
    std::cout << "Command failed. " << std::endl;
    exit(-1);
  }
}
