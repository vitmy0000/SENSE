#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <iomanip>

#include <cstdint>
#include <cstdlib>
#include <cassert>

using std::string;
using std::vector;

class Kmer {
 public:
  Kmer(const string& seq_read, const std::size_t k) :
      seq_read_{seq_read}, seq_len_{seq_read.size()} {
    uint64_t index = 0;
    for (std::size_t i = 0; i < k - 1; ++i) {
      index <<= 2u;
      index |= char_map[static_cast<std::size_t>(seq_read[i])];
    }
    uint64_t mask = (1 << 2 * k) - 1;
    std::map<uint16_t, std::size_t> index_2_cnt;
    for (std::size_t i = 0; i < seq_len_ - k + 1; ++i) {
      index <<= 2u;
      index |= char_map[static_cast<std::size_t>(seq_read[k - 1 + i])];
      index &= mask;
      ++index_2_cnt[index];
    }
    // fetch non-zero entries only
    index_array_.resize(index_2_cnt.size());
    cnt_array_.resize(index_2_cnt.size());
    std::size_t increment = 0;
    for (const auto& kv : index_2_cnt) {
      index_array_[increment] = kv.first;
      cnt_array_[increment] = kv.second;
      ++increment;
    }
  }
  vector<uint64_t> get_index_array() const {
    return index_array_;
  }
  vector<std::size_t> get_cnt_array() const {
    return cnt_array_;
  }
  std::size_t get_seq_len() const {
    return seq_len_;
  }
 private:
  const string seq_read_;
  const std::size_t seq_len_;
  vector<uint64_t> index_array_;
  vector<std::size_t> cnt_array_;
  static uint8_t char_map[];
};

uint8_t Kmer::char_map[] = {
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
};

float kmer_dist(const Kmer& kmer0,
                const Kmer& kmer1,
                const std::size_t k) {
  std::size_t i0 = 0;
  std::size_t i1 = 0;
  const vector<uint64_t>& index_array0 = kmer0.get_index_array();
  const vector<uint64_t>& index_array1 = kmer1.get_index_array();
  const vector<std::size_t>& cnt_array0 = kmer0.get_cnt_array();
  const vector<std::size_t>& cnt_array1 = kmer1.get_cnt_array();
  std::size_t intersect_cnt = 0;
  while (i0 < index_array0.size() && i1 < index_array1.size()) {
    if (index_array0[i0] < index_array1[i1]) {
      ++i0;
    } else if (index_array0[i0] > index_array1[i1]) {
      ++i1;
    } else {
      intersect_cnt += std::min(cnt_array0[i0], cnt_array1[i1]);
      ++i0;
      ++i1;
    }
  }
  return 1 - static_cast<float>(intersect_cnt) /
      (std::min(kmer0.get_seq_len(), kmer1.get_seq_len()) - k + 1);
}

int main(int argc, char *argv[]) {
  std::ifstream fa_file(argv[1]);
  std::ofstream out_file(argv[2]);
  const std::size_t K = strtoul(argv[3], nullptr, 10);

  // load data
  vector<string> headers;
  vector<string> seq_reads;
  if (fa_file.is_open()) {
    std::size_t cnt = 0;
    string line;
    string header;
    string seq_read;
    while (getline(fa_file, line)) {
      if (cnt % 2 == 0) {
        headers.push_back(line.substr(1));
      } else {
        seq_reads.push_back(line);
      }
      ++cnt;
    }
    fa_file.close();
  } else {
    std::cout << "unable to open file" << std::endl;
  }
  assert(headers.size() == seq_reads.size());

  // convert all seqs
  vector<Kmer> kmers;
  for (string& seq_read : seq_reads) {
    kmers.push_back(Kmer(seq_read, K));
  }
  assert(seq_reads.size() == kmers.size());
  std::cout << "Converstion done, make sure k is smaller than 32" << std::endl;

  // pair_wise distance
  if (out_file.is_open()) {
    std::cout << std::setprecision(4);
    for (decltype(headers.size()) i = 0, sz = headers.size(); i != sz; ++i) {
      for (decltype(headers.size()) j = 0, sz = headers.size(); j != sz; ++j) {
        if (i < j) {
          out_file << headers[i] << "-" << headers[j]
                   << "\t" << kmer_dist(kmers[i], kmers[j], K) << "\n";
        }
      }
    }
  } else {
    std::cout << "unable to open file" << std::endl;
  }

  return 0;
}
