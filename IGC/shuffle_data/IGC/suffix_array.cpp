#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <limits>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <algorithm>

using namespace std;

// g++ -std=c++0x suffix_array.cpp -o suffix_array_sort

struct SortBySuffixCMP {
  explicit SortBySuffixCMP(const vector<pair<string, string> >& _proteins)
      : proteins(_proteins) {
  }
  bool operator()(const pair<uint32_t, uint32_t>& p1,
                  const pair<uint32_t, uint32_t>& p2) {
    uint32_t i = p1.second, j = p2.second, k = 0;
    while (i < proteins[p1.first].second.size()
        && j < proteins[p2.first].second.size() && k < 500) {
      if (proteins[p1.first].second[i] < proteins[p2.first].second[j]) {
        return true;
      } else if (proteins[p1.first].second[i] > proteins[p2.first].second[j]) {
        return false;
      }
      i++;
      j++;
      k++;
    }
    if (k == 500)
      return false;
    if (i == proteins[p1.first].second.size()
        && j != proteins[p2.first].second.size()) {
      return true;
    } else if (i == proteins[p1.first].second.size()
        && j == proteins[p2.first].second.size()) {
      return false;
    } else if (i != proteins[p1.first].second.size()) {
      return false;
    }
  }

  const vector<pair<string, string> >& proteins;
};

int main(int argc, const char **argv) {
  ifstream fin("/home/rcf-40/haifengc/panfs/github/IGC/data/IGC.pep");
  // ifstream fin("test_peptides.txt");
  string name, peptide;
  vector < pair<string, string> > proteins;
  uint64_t total_length = 0;
  uint32_t max_length = 0, min_length = std::numeric_limits < uint32_t
      > ::max();
  while (getline(fin, name)) {
    getline(fin, peptide);
    if (peptide[0] == '>') {
      cout << "error" << endl;
    }
    proteins.push_back(make_pair(name, peptide));
    total_length += peptide.size();
    if (peptide.size() > max_length) {
      max_length = peptide.size();
    }
    if (peptide.size() < min_length) {
      min_length = peptide.size();
    }
  }

  uint64_t k = 0;
  vector < pair<uint32_t, uint32_t> > positions(total_length);
  for (uint32_t i = 0; i < proteins.size(); ++i) {
    for (uint32_t j = 0; j < proteins[i].second.size(); ++j) {
      positions[k++] = make_pair(i, j);
    }
  }
  cout << "The number of positions is " << total_length << "." << endl;
  cout << "The number of peptides is " << proteins.size() << "." << endl;
  cout << "Max length is " << max_length << "." << endl;
  cout << "Min length is " << min_length << "." << endl;
  clock_t start_t = clock();
  sort(positions.begin(), positions.end(), SortBySuffixCMP(proteins));
  cout << "Sorting takes " << (double(clock() - start_t) / CLOCKS_PER_SEC)
      << " seconds." << endl;

  ofstream fout("peptides_suffix_array_sort.txt");
  for (uint64_t i = 0; i < positions.size(); ++i) {
    fout << positions[i].first << "\t" << positions[i].second << "\t";
    fout << endl;
  }
  fout.close();

  start_t = clock();
  fout.open("peptides_sort.txt");
  for (uint64_t i = 0; i < positions.size(); ++i) {
    fout << positions[i].first << "\t" << positions[i].second << "\t";
    for (uint32_t j = positions[i].second, k = 0;
        j < proteins[positions[i].first].second.size() && k < 10; ++j, ++k) {
      fout << proteins[positions[i].first].second[j];
    }
    fout << endl;
  }
  fout.close();
  cout << "Writing takes " << (double(clock() - start_t) / CLOCKS_PER_SEC)
      << " seconds." << endl;

  return 0;
}
