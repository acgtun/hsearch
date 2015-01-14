#include "sdk.hpp"

#include <ctime>

#include <vector>
#include <string>
#include <fstream>
#include <iostream>

#include "distance_matrix.hpp"

#include "k_nearest_neighbor_similarity.hpp"
#include "k_nearest_neighbor_distance.hpp"

using namespace_k_nearest_neighbor_similarity::KNeighborsSimilarity;
using namespace_k_nearest_neighbor_distance::KNeighborsDistance;

using namespace std;

#define MAX_NUM_CANDIDATE 10

string RandomGenerateKmer() {
  string rand_seed;
  for (int i = 0; i < HASHAALEN; ++i) {
    int r = rand() % 20;
    rand_seed += AA20[r];
  }
  return rand_seed;
}

void FindSimilarityCandidate(KNeighborsSimilarity& k_similarity,
                             const string& rand_seed,
                             vector<string>& candidate) {
  k_similarity.UpdateWeight(rand_seed.c_str());
  bool bPath = true;
  int pathID = 0;
  int rawScore = 0;
  char strCandPep[HASHAALEN + 1];
  while (bPath && pathID < MAX_NUM_CANDIDATE) {
    bPath = k_similarity.FindCandidate(strCandPep, pathID, rawScore);
    candidate.push_back(strCandPep);
  }
}

void FindDistanceCandidate(KNeighborsDistance& k_distance,
                           const string& rand_seed, vector<string>& candidate) {
  k_distance.UpdateWeight(rand_seed.c_str());
  bool bPath = true;
  int pathID = 0;
  double rawScore = 0;
  char strCandPep[HASHAALEN + 1];
  while (bPath && pathID < MAX_NUM_CANDIDATE) {
    bPath = k_distance.FindCandidate(strCandPep, pathID, rawScore);
    candidate.push_back(strCandPep);
  }
}

int main(int argc, const char* argv[]) {
  vector<vector<double> > distance;
  GetMetricDistance(distance);

  KNeighborsSimilarity k_similarity;
  KNeighborsDistance k_distance(&distance);

  srand(time(NULL));

  ofstream fout("test.out");
  for (int i = 0; i < 1000; ++i) {
    string rand_seed = RandomGenerateKmer();
    vector<string> candidate_similarity;
    vector<string> candidate_distance;
    FindSimilarityCandidate(k_similarity, rand_seed, candidate_similarity);
    FindDistanceCandidate(k_distance, rand_seed, candidate_distance);
    fout << "test " << i + 1 << endl;
    fout << "--" << rand_seed << endl;
    size_t p = 0;
    while (p < candidate_similarity.size() && p < candidate_distance.size()) {
      fout << candidate_similarity[p] << " " << candidate_distance[p] << endl;
    }
  }

  return 0;
}
