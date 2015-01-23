#include "sdk.hpp"
#include "bio_util.hpp"

#include <cmath>
#include <vector>
#include <iostream>

using std::vector;
using std::cerr;
using std::cout;
using std::endl;

bool GetMetricDistance2(vector<vector<double> >& distance, const double& lamda,
                       const double& e) {

  distance.resize(20);
  for (int i = 0; i < 20; i++) {
    distance[i].resize(20);
    for (int j = 0; j < 20; j++) {
      if (i == j) {
        distance[i][j] = 0.000;
      } else {
        distance[i][j] = 1 - BLOSUM62[i][j];
      }
    }
  }

  int cnt = 0;
  for (int i = 0; i < 20; i++) {
    for (int j = 0; j < 20; j++) {
      for (int k = 0; k < 20; k++) {
        if (distance[i][j] + distance[j][k] < distance[i][k]) {
          cnt++;
        } else {
        }
      }
    }
  }
  if (cnt > 0) {
    cerr << "error~!~" << endl;
    //exit (EXIT_FAILURE);
    return false;
  }

  return true;
}

bool GetMetricDistance(vector<vector<double> >& distance, const double& lamda,
                        const double& e) {
  int cnt = 0;
  for (int i = 0; i < 20; i++) {
    for (int j = 0; j < 20; j++) {
      for (int k = 0; k < 20; k++) {
        if (pow(e, (lamda * BLOSUM62[i][j])) + pow(e, (lamda * BLOSUM62[j][k]))
            < pow(e, (lamda * BLOSUM62[i][k]))) {
          cnt++;
        } else {
        }
      }
    }
  }
  if (cnt > 0) {
    cerr << "error~!~" << endl;
    //exit (EXIT_FAILURE);
    return false;
  }

  distance.resize(20);
  for (int i = 0; i < 20; i++) {
    distance[i].resize(20);
    for (int j = 0; j < 20; j++) {
      if (i == j) {
        distance[i][j] = 0.000;
      } else {
        distance[i][j] = pow(e, (lamda * BLOSUM62[i][j]));
      }
    }
  }

  return true;
}
