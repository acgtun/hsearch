
#include "sdk.hpp"
#include "bio_util.hpp"

#include <cmath>
#include <vector>
#include <iostream>

using std::vector;
using std::cerr;
using std::cout;
using std::endl;

void GetMetricDistance(vector<vector<double> >& distance) {
  double lamda = -0.17;
  int cnt = 0;
  for (int i = 0; i < 20; i++) {
    for (int j = 0; j < 20; j++) {
      for (int k = 0; k < 20; k++) {
        if (exp(lamda * BLOSUM62[i][j]) + exp(lamda * BLOSUM62[j][k])
            < exp(lamda * BLOSUM62[i][k])) {
          cnt++;
        } else {
        }
      }
    }
  }
  if (cnt > 0) {
    cerr << "error~!~" << endl;
    exit (EXIT_FAILURE);
  }

  distance.resize(20);
  for (int i = 0; i < 20; i++) {
    distance[i].resize(20);
    for (int j = 0; j < 20; j++) {
      if (i == j) {
        distance[i][j] = 0.000;
      } else {
        distance[i][j] = exp(lamda * BLOSUM62[i][j]);
      }
    }
  }
}
