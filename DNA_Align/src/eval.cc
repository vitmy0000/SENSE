#include <iostream>
#include <fstream>
#include <vector>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

int main(int argc, char** argv) {
  std::ifstream embedding_file(argv[1]);
  std::ofstream dist_file(argv[2]);
  std::size_t N = boost::lexical_cast<std::size_t>(argv[3]);
  std::size_t dim = boost::lexical_cast<std::size_t>(argv[4]);
  if (embedding_file.is_open() && dist_file.is_open()) {
    std::string line;
    std::getline(embedding_file, line);
    std::vector<std::string> str_vec;
    boost::split(str_vec, line, boost::is_any_of(","));
    std::vector<double> vec;
    for (const auto& s : str_vec) {
      vec.push_back(boost::lexical_cast<double>(s));
    }
    std::cout << "read done!" << std::endl;
    std::vector<double> dists;
    for (std::size_t i = 0; i != N; ++i) {
      for (std::size_t j = 0; j != N; ++j) {
        if (i < j) {
          double max_sum = 0.0;
          double min_sum = 0.0;
          for (std::size_t k = 0; k != dim; ++k) {
            double v_i = vec[i * dim + k];
            double v_j = vec[j * dim + k];
            if (v_i > v_j) {
              max_sum += v_i;
              min_sum += v_j;
            } else {
              min_sum += v_i;
              max_sum += v_j;
            }
          }
          dists.push_back(1 - min_sum / max_sum);
        }
      }
    }
    for (const auto& d : dists) {
      dist_file << d << "\n";
    }
    std::cout << "dist done!" << std::endl;
  } else {
    std::cout << "unable to open file" << std::endl;
  }
}

