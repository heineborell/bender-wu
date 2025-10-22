#include "helper.h"

void printDict(std::vector<std::vector<RCP<const Basic>>> &dict) {
  std::size_t rows{dict.size()};
  std::size_t columns{dict[0].size()};
  std::cout << "the dictionary is " << rows - 1 << 'x' << columns - 1 << '\n';
  for (std::size_t i{1}; i < rows; ++i) {
    for (std::size_t j{1}; j < columns; ++j) {
      if (dict[i][j].get()) // check if its nullptr (since these are smart gotta
                            // use .get)
        std::cout << '(' << i << ',' << j << "): " << *dict[i][j] << '\n';
      else
        std::cout << '(' << i << ',' << j << "): " << "nullptr" << '\n';
    }
  }
}

void printArray(std::vector<RCP<const Basic>> &arr) {
  for (const RCP<const Basic> &item : arr) {
    std::cout << *item << '\n';
  }
}
