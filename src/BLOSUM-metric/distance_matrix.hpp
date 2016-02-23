#include "sdk.hpp"
#include "bio_util.hpp"

#include <cmath>
#include <vector>
#include <iostream>

using std::vector;
using std::cerr;
using std::cout;
using std::endl;

void SimilarityMatrix2DistanceMatrix(vector<vector<int> >& distance) {
  distance.resize(20);
  for (int i = 0; i < 20; ++i) {
    distance[i].resize(20);
    for (int j = 0; j < 20; ++j) {
      distance[i][j] = BLOSUM62[i][i] + BLOSUM62[j][j] - 2 * BLOSUM62[i][j];
    }
  }

  for (int i = 0; i < 20; ++i) {
    for (int j = 0; j < 20; ++j) {
      cout << BLOSUM62[i][j] << " ";
    }
    cout << endl;
  }
  cout << "..........................." << endl;
  for (int i = 0; i < 20; ++i) {
    for (int j = 0; j < 20; ++j) {
      cout << distance[i][j] << " ";
    }
    cout << endl;
  }

  int cnt = 0;
  for (int i = 0; i < 20; i++) {
    for (int j = 0; j < 20; j++) {
      for (int k = 0; k < 20; k++) {
        if (distance[i][j] + distance[j][k] < distance[i][k]) {
          cnt++;
        }
      }
    }
  }
  if (cnt == 0) {
    cout << "Triangle Inequality Holds" << endl;
  } else {
    cout << "Triangle Inequality Doesn't Hold" << endl;
  }
}

