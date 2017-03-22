#include "option.hpp"

#include "sdk.hpp"

#include <stdio.h>
#include <ctime>

#include <vector>
#include <string>
#include <fstream>
#include <iostream>

#include "distance_matrix.hpp"

using namespace std;

string RandomGenerateKmer(const int& k_mer) {
  string rand_seed;
  for (int i = 0; i < k_mer; ++i) {
    int r = rand() % 20;
    rand_seed += AA20[r];
  }
  return rand_seed;
}

int SimilarityScore(const string& s1, const string& s2) {
  int s = 0;
  for (size_t i = 0; i < s1.size(); ++i) {
    s += BLOSUM62[base[s1[i] - 'A']][base[s2[i] - 'B']];
  }
  return s;
}

int DistanceScore(const string& s1, const string& s2,
                  const vector<vector<int> >& distance) {
  int d = 0;
  for (size_t i = 0; i < s1.size(); ++i) {
    d += distance[base[s1[i] - 'A']][base[s2[i] - 'B']];
  }
  return d;
}

int main(int argc, const char* argv[]) {
  InitProgram(argc, argv);
  srand(time(NULL));
  for (int k_mer = 2; k_mer <= 200; k_mer++) {
    vector<string> seeds_1;
    vector<string> seeds_2;
    for (int i = 0; i < 1000000; ++i) {
      seeds_1.push_back(RandomGenerateKmer(k_mer));
      seeds_2.push_back(RandomGenerateKmer(k_mer));
    }

    vector<vector<int> > distance;
    SimilarityMatrix2DistanceMatrix(distance);
    char fs_chr[100], fd_chr[100];
    sprintf(fs_chr, "smimilarty_%d.txt", k_mer);
    sprintf(fd_chr, "distance_%d.txt", k_mer);
    ofstream fs(fs_chr);
    ofstream fd(fd_chr);
    for (size_t i = 0; i < seeds_1.size(); ++i) {
      fs << SimilarityScore(seeds_1[i], seeds_2[i]) << endl;
      fd << DistanceScore(seeds_1[i], seeds_2[i], distance) << endl;
    }

    fs.close();
    fd.close();
  }

  return 0;
}
