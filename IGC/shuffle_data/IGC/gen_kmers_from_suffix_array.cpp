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

// g++ -std=c++0x gen_kmers_from_suffix_array.cpp -o gen_kmers

int main(int argc, const char **argv) {
  ifstream fin("/home/rcf-40/haifengc/panfs/github/IGC/data/IGC.pep");
  // ifstream fin("test_peptides.txt");
  string name, peptide;
  vector < pair<string, string> > proteins;
  uint64_t total_length = 0;
  uint32_t max_length = 0;
  uint32_t min_length = std::numeric_limits < uint32_t > ::max();
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
  fin.close();

  uint32_t kmer_length = 0;
  sscanf(argv[1], "%u", &kmer_length);
  uint32_t id = 0;
  string cur;
  pair<uint32_t, uint32_t> cur_pos;
  uint32_t cnt = 0;
  uint32_t proID, proPos;
  ifstream fsuffix("peptides_suffix_array_sort.txt");
  string output_file = "kmers.gen";
  output_file += argv[1];
  output_file += ".txt";
  ofstream fout(output_file.c_str());
  while (fsuffix >> proID >> proPos) {
    if (proPos + kmer_length > proteins[proID].second.size()) {
      continue;
    }
    string tmp = proteins[proID].second.substr(proPos, kmer_length);
    if (tmp == cur) {
      cnt++;
    } else {
      if (cnt != 0)
        fout << proID << "\t" << proPos << "\t" << cnt << endl;
      cur = tmp;
      cnt = 1;
      cur_pos = make_pair(proID, proPos);
    }
  }
  if (cnt != 0) {
    fout << cur << "\t" << cnt << endl;
  }
  fout.close();

  return 0;
}
