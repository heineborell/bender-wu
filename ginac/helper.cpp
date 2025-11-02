#include "helper.h"

void printArray(std::vector<ex> &arr) {
  for (ex &item : arr) {
    std::cout << item << '\n';
  }
}

void printDict(std::vector<std::vector<ex>> &dict) {
  std::size_t rows{dict.size()};
  std::size_t columns{dict[0].size()};
  std::cout << "the dictionary is " << rows - 1 << 'x' << columns - 1 << '\n';
  for (std::size_t i{1}; i < rows; ++i) {
    for (std::size_t j{1}; j < columns; ++j) {
      std::cout << '(' << i << ',' << j << "): " << dict[i][j] << '\n';
    }
  }
}
