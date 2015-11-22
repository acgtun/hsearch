#include "sdk.hpp"

#include <stdio.h>
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

//#define MAX_NUM_CANDIDATE 5

string RandomGenerateKmer() {
  string rand_seed;
  for (int i = 0; i < HASHAALEN; ++i) {
    int r = rand() % 20;
    rand_seed += AA20[r];
  }
  return rand_seed;
}

void FindSimilarityCandidate(KNeighborsSimilarity& k_similarity,
                             const string& rand_seed, vector<string>& candidate,
                             const int& MAX_NUM_CANDIDATE) {
  k_similarity.UpdateWeight(rand_seed.c_str());
  bool bPath = true;
  int pathID = 0;
  int rawScore = 0;
  char strCandPep[HASHAALEN + 1];
  while (bPath && pathID < MAX_NUM_CANDIDATE) {
    bPath = k_similarity.FindCandidate(strCandPep, pathID, rawScore);
    candidate.push_back(strCandPep);
    pathID++;
  }
}

void FindDistanceCandidate(KNeighborsDistance& k_distance,
                           const string& rand_seed, vector<string>& candidate,
                           const int& MAX_NUM_CANDIDATE) {
  k_distance.UpdateWeight(rand_seed.c_str());
  bool bPath = true;
  int pathID = 0;
  double rawScore = 0;
  char strCandPep[HASHAALEN + 1];
  while (bPath && pathID < MAX_NUM_CANDIDATE) {
    bPath = k_distance.FindCandidate(strCandPep, pathID, rawScore);
    candidate.push_back(strCandPep);
    pathID++;
  }
}

void evaluate(const vector<vector<double> >& distance,
              const vector<string>& seeds, const int& MAX_NUM_CANDIDATE) {
  KNeighborsSimilarity k_similarity;
  KNeighborsDistance k_distance(&distance);
  ofstream fout("test.out");
  int total = 0, cnt = 0;
  for (size_t i = 0; i < seeds.size(); ++i) {
    vector<string> candidate_similarity;
    vector<string> candidate_distance;
    FindSimilarityCandidate(k_similarity, seeds[i], candidate_similarity,
                            MAX_NUM_CANDIDATE);
    FindDistanceCandidate(k_distance, seeds[i], candidate_distance,
                          MAX_NUM_CANDIDATE);
    fout << "test " << i + 1 << endl;
    fout << "--" << seeds[i] << endl;
    fout << "size " << candidate_similarity.size() << " "
         << candidate_distance.size() << endl;
    size_t p = 0;
    while (p < candidate_similarity.size() && p < candidate_distance.size()) {
      fout << candidate_similarity[p] << " " << candidate_distance[p] << endl;
      p++;
    }

    for (size_t p = 0; p < candidate_similarity.size(); ++p) {
      total++;
      for (size_t q = 0; q < candidate_distance.size(); ++q) {
        if (candidate_similarity[p] == candidate_distance[q]) {
          cnt++;
          break;
        }
      }
    }
  }

  cout << "total = " << total << " " << cnt << " " << cnt / (double) total
       << endl;
}

int main(int argc, const char* argv[]) {

  srand(time(NULL));
  vector<string> seeds;
  for (int i = 0; i < 100000; ++i) {
    seeds.push_back(RandomGenerateKmer());
  }

  vector<vector<double> > distance;
  GetMetricDistanceother(distance);
  int candidates[] = { 1, 5, 10, 20, 30, 50, 80, 100, 150, 200, 300, 500, 800,
      1000, 2000, 5000 };
  vector<int> MAX_NUM_CANDIDATE(candidates,
                                candidates + sizeof(candidates) / sizeof(int));
  for (size_t i = 0; i < MAX_NUM_CANDIDATE.size(); ++i) {
    cout << "MAX_NUM_CANDIDATE = " << MAX_NUM_CANDIDATE[i] << endl;
    evaluate(distance, seeds, MAX_NUM_CANDIDATE[i]);
  }
  //////////////////////
  /* cout << "new method end..............." << endl;

   //  for (double e = 1.1; e < 100; e += 0.1) {
   double e = 2.76;
   for (double lamda = -0.5; lamda < 0; lamda += 0.01) {
   cout << "e lamda = " << e << " " << lamda << endl;
   vector<vector<double> > distance;
   if (!GetMetricDistance(distance, lamda, e))
   continue;

   for (int i = 0; i < 20; ++i) {
   for (int j = 0; j < 20; ++j) {
   //    printf("%.3lf ", distance[i][j]);
   }
   cout << endl;
   }

   evaluate(distance, seeds);
   }
   // }
   */

  return 0;
}
