#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <numeric>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <cassert>
#include <cstdlib>

#include "cxxopts.hh"

enum class Direction {
  noDir, fromM, fromX, fromY
};

std::tuple<std::string, std::string, int>
nw_align(const std::string& seq0, const std::string& seq1) {
  const int matchScore = 5;
  const int missScore = -4;
  const int gapOpenScore = -10;
  const int gapExtendScore = -1;

  auto nrow = seq0.size() + 1;
  auto ncol = seq1.size() + 1;
  // score matrix of best alignment of x[1..i] and y[1..j] ending with a char-char match or mismatch
  std::vector<int> M(nrow * ncol);
  // score matrix of best alignment of x[1..i] and y[1..j] ending with a space in x
  std::vector<int> X(nrow * ncol);
  // score matrix of best alignment of x[1..i] and y[1..j] ending with a space in y
  std::vector<int> Y(nrow * ncol);
  // direction matrix for backtracking
  std::vector<Direction> dirM(nrow * ncol);
  std::vector<Direction> dirX(nrow * ncol);
  std::vector<Direction> dirY(nrow * ncol);

  for (decltype(nrow) i = 1; i < nrow; ++i) {
    M[i * ncol + 0] = -999999999;
    X[i * ncol + 0] = -999999999;
    Y[i * ncol + 0] = gapOpenScore + i * gapExtendScore;
  }
  for (decltype(ncol) j = 1; j < ncol; ++j) {
    M[0 * ncol + j] = -999999999;
    X[0 * ncol + j] = gapOpenScore + j * gapExtendScore;
    Y[0 * ncol + j] = -999999999;
  }
  for (decltype(nrow) i = 1; i < nrow; ++i) {
    for (decltype(ncol) j = 1; j < ncol; ++j) {
      int addScore = 0;
      if (seq0[i - 1] == seq1[j - 1]) {
        addScore = matchScore;
      } else {
        addScore = missScore;
      }
      // M
      int scoreMM = M[(i - 1) * ncol + (j - 1)] + addScore;
      int scoreMX = X[(i - 1) * ncol + (j - 1)] + addScore;
      int scoreMY = Y[(i - 1) * ncol + (j - 1)] + addScore;
      if (scoreMM >= scoreMX && scoreMM >= scoreMY){
        M[i * ncol + j] = scoreMM;
        dirM[i * ncol + j] = Direction::fromM;
      } else if (scoreMX >= scoreMM && scoreMX >= scoreMY) {
        M[i * ncol + j] = scoreMX;
        dirM[i * ncol + j] = Direction::fromX;
      } else if (scoreMY >= scoreMM && scoreMY >= scoreMX) {
        M[i * ncol + j] = scoreMY;
        dirM[i * ncol + j] = Direction::fromY;
      }
      // X
      int scoreXM = M[i * ncol + (j - 1)] + gapOpenScore + gapExtendScore;
      int scoreXX = X[i * ncol + (j - 1)] + gapExtendScore;
      if (scoreXM >= scoreXX) {
        X[i * ncol + j] = scoreXM;
        dirX[i * ncol + j] = Direction::fromM;
      } else {
        X[i * ncol + j] = scoreXX;
        dirX[i * ncol + j] = Direction::fromX;
      }
      // Y
      int scoreYM = M[(i - 1) * ncol + j] + gapOpenScore + gapExtendScore;
      int scoreYY = Y[(i - 1) * ncol + j] + gapExtendScore;
      if (scoreYM >= scoreYY) {
        Y[i * ncol + j] = scoreYM;
        dirY[i * ncol + j] = Direction::fromM;
      } else {
        Y[i * ncol + j] = scoreYY;
        dirY[i * ncol + j] = Direction::fromY;
      }
    }
  }
  // backtracking
  int final_socre = 0;
  int i = nrow - 1;
  int j = ncol - 1;
  std::string aligned_seq0;
  std::string aligned_seq1;
  std::vector<Direction>* dir_matrix_ptr = nullptr;
  if (X[i * ncol + j] >= M[i * ncol + j] &&
      X[i * ncol + j] >= Y[i * ncol + j]) {
    dir_matrix_ptr = &dirX;
    final_socre = X[i * ncol + j];
  } else if (Y[i * ncol + j] >= M[i * ncol + j] &&
             Y[i * ncol + j] >= X[i * ncol + j]) {
    dir_matrix_ptr = &dirY;
    final_socre = Y[i * ncol + j];
  } else {
    dir_matrix_ptr = &dirM;
    final_socre = M[i * ncol + j];
  }
  while (i > 0 && j > 0) {
    Direction dir = (*dir_matrix_ptr)[i * ncol + j];
    if (dir_matrix_ptr == &dirM) {
      aligned_seq0.push_back(seq0[i - 1]);
      aligned_seq1.push_back(seq1[j - 1]);
      i--;
      j--;
    } else if (dir_matrix_ptr == &dirX) {
      aligned_seq0.push_back('-');
      aligned_seq1.push_back(seq1[j - 1]);
      j--;
    } else if (dir_matrix_ptr == &dirY) {
      aligned_seq0.push_back(seq0[i - 1]);
      aligned_seq1.push_back('-');
      i--;
    }
    switch (dir) {
      case Direction::fromM:
        dir_matrix_ptr = &dirM;
        break;
      case Direction::fromX:
        dir_matrix_ptr = &dirX;
        break;
      case Direction::fromY:
        dir_matrix_ptr = &dirY;
        break;
      default:
        break;
    }
  }
  std::reverse(std::begin(aligned_seq0), std::end(aligned_seq0));
  std::reverse(std::begin(aligned_seq1), std::end(aligned_seq1));
  return std::make_tuple(aligned_seq0, aligned_seq1, final_socre);
}

float nw_dist(const std::string& aligned_seq0,
              const std::string& aligned_seq1) {
  assert(aligned_seq0.size() == aligned_seq1.size());
  std::size_t match_cnt = 0;
  std::size_t miss_cnt = 0;
  for (decltype(aligned_seq0.size()) i = 0; i != aligned_seq0.size(); ++i) {
    if (aligned_seq0[i] == aligned_seq1[i]) {
      ++match_cnt;
    } else {
      ++miss_cnt;
    }
  }
  return static_cast<float>(miss_cnt) / (miss_cnt + match_cnt);
}

int main(int argc, char** argv) {
  cxxopts::Options options("select_train",
                           "A program for selecting training pairs.");
  options.add_options()
    ("f,fasta_file", "input fasta file", cxxopts::value<std::string>())
    ("s,seq_ids_file", "output seq ids file", cxxopts::value<std::string>())
    ("p,pairs_file", "output pairs file", cxxopts::value<std::string>())
    ("d,dists_file", "output dists file", cxxopts::value<std::string>())
    ("a,abundance_threshold", "abundance threshold", cxxopts::value<std::size_t>())
    ("t,target_num_landmarks", "target number of landmarks", cxxopts::value<std::size_t>())
    ("n,num_random_sample", "number of random sample", cxxopts::value<std::size_t>())
    ;
  options.parse(argc, argv);

  std::ifstream fasta_file(options["fasta_file"].as<std::string>());
  std::ofstream seq_ids_file(options["seq_ids_file"].as<std::string>());

  // read input and derep
  std::unordered_map<std::string, std::vector<std::string>> derep_map;
  if (fasta_file.is_open()) {
    std::size_t line_num = 0;
    std::string line, seq_id, seq_read;
    while (std::getline(fasta_file, line)) {
      if (line_num % 2 == 0) {
        seq_id = line.substr(1);
      } else {
        seq_read = line;
        derep_map[seq_read].push_back(seq_id);
      }
      ++line_num;
    }
    fasta_file.close();
  } else {
    std::cout << "unable to open file" << std::endl;
  }
  std::cout << "Reading finished!" << std::endl;
  std::cout << "Number of unique sequences: " << derep_map.size() << std::endl;

  // select abundant seqs
  std::size_t abundance_threshold =
      options["abundance_threshold"].as<std::size_t>();
  std::vector<std::string> abundant_seq_ids;
  std::vector<std::string> abundant_seq_reads;
  if (seq_ids_file.is_open()) {
    for (const auto& kv : derep_map) {
      const std::string& seq_read = kv.first;
      const std::vector<std::string>& seq_ids = kv.second;
      std::copy(seq_ids.begin(), seq_ids.end(),
                std::ostream_iterator<std::string>(seq_ids_file, "\t"));
      seq_ids_file << "\n";
      if (seq_ids.size() >= abundance_threshold) {
        abundant_seq_ids.push_back(seq_ids[0]);
        abundant_seq_reads.push_back(seq_read);
      }
    }
    seq_ids_file.close();
  } else {
    std::cout << "unable to open file" << std::endl;
  }
  assert(abundant_seq_ids.size() == abundant_seq_reads.size());
  derep_map.clear();
  std::cout << "Number of abundant sequences: "
            << abundant_seq_ids.size() << std::endl;

  // adaptive landmark selection; save pairs and dists
  std::srand(0);
  std::size_t landmark_index = std::rand() % abundant_seq_reads.size();
  std::size_t num_landmarks = 0;
  std::ofstream pairs_file(options["pairs_file"].as<std::string>());
  std::ofstream dists_file(options["dists_file"].as<std::string>());
  const std::size_t target_num_landmarks =
      options["target_num_landmarks"].as<std::size_t>();
  if (pairs_file.is_open() && dists_file.is_open()) {
    while (num_landmarks < target_num_landmarks) {
      std::cout << "landmark # " << num_landmarks << std::endl;
      std::string landmark_id = abundant_seq_ids[landmark_index];
      std::string landmark_read = abundant_seq_reads[landmark_index];

      std::size_t num_random_sample =
          options["num_random_sample"].as<std::size_t>();
      assert(num_random_sample <= abundant_seq_reads.size());
      std::unordered_set<std::size_t> selected_idx;
      std::vector<std::size_t> idxs;
      // random sample and choose the furthest as next landmark
      std::vector<float> min_dists(num_random_sample, 1.0f);
      while (selected_idx.size() != num_random_sample) {
        std::size_t idx = std::rand() % abundant_seq_reads.size();
        auto find_res = selected_idx.find(idx);
        if (find_res == selected_idx.end()) {
          selected_idx.insert(idx);
          idxs.push_back(idx);
        } else {
          continue;
        }
        std::string seq_id = abundant_seq_ids[idx];
        std::string seq_read = abundant_seq_reads[idx];
        const auto& align_res = nw_align(seq_read, landmark_read);
        float dist = nw_dist(std::get<0>(align_res), std::get<1>(align_res));
        min_dists[idxs.size() - 1] = std::min(dist, min_dists[idxs.size() - 1]);
        pairs_file << ">" << landmark_id << "-" << seq_id << "\n"
                   << landmark_read << "\n"
                   << ">" << landmark_id << "-" << seq_id << "\n"
                   << seq_read << "\n";
        dists_file << landmark_id << "-" << seq_id << "\t"
                   << std::setprecision(4) << dist << "\n";
      }
      std::iota(idxs.begin(), idxs.end(), 0);
      std::sort(idxs.begin(), idxs.end(),
          [&min_dists](size_t i1, size_t i2) {return min_dists[i1] < min_dists[i2];});
      landmark_index = idxs[std::rand() % (idxs.size()/2) + (idxs.size()/2)];
      ++num_landmarks;
    }
    pairs_file.close();
    dists_file.close();
  } else {
    std::cout << "unable to open file" << std::endl;
  }

  return EXIT_SUCCESS;
}
