#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>

using std::string;
using std::vector;

enum class Direction {
  noDir, fromM, fromX, fromY
};

std::tuple<string, string, int> nw_align(const string& seq0, const string& seq1) {
  const int matchScore = 5;
  const int missScore = -4;
  const int gapOpenScore = -10;
  const int gapExtendScore = -1;

  auto nrow = seq0.size() + 1;
  auto ncol = seq1.size() + 1;
  // score matrix of best alignment of x[1..i] and y[1..j] ending with a char-char match or mismatch
  vector<int> M(nrow * ncol);
  // score matrix of best alignment of x[1..i] and y[1..j] ending with a space in x
  vector<int> X(nrow * ncol);
  // score matrix of best alignment of x[1..i] and y[1..j] ending with a space in y
  vector<int> Y(nrow * ncol);
  // direction matrix for backtracking
  vector<Direction> dirM(nrow * ncol);
  vector<Direction> dirX(nrow * ncol);
  vector<Direction> dirY(nrow * ncol);

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
  string aligned_seq0;
  string aligned_seq1;
  vector<Direction>* dir_matrix_ptr = nullptr;
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

int main(int argc, char *argv[]) {
  std::ifstream fa_file(argv[1]);
  std::ofstream out_file(argv[2]);
  if (fa_file.is_open() && out_file.is_open()) {
    std::size_t cnt = 0;
    string line;
    string header;
    string seq0, seq1;
    while (getline(fa_file, line)) {
      if (cnt % 4 == 0) {
        header = line.substr(1);
      } else if (cnt % 4 == 1) {
        seq0 = line;
      } else if (cnt % 4 == 3) {
        seq1 = line;
        auto res = nw_align(seq0, seq1);
        // write output
        out_file << header << std::endl;
        out_file << std::get<0>(res) << std::endl;
        out_file << std::get<1>(res) << std::endl;
        out_file << "score: " << std::get<2>(res) << std::endl;
      }
      ++cnt;
    }
    fa_file.close();
    out_file.close();
  } else {
    std::cout << "unable to open file" << std::endl;
  }
  return 0;
}
