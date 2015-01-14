#include <vector>
#include <iostream>


#include "distance_matrix.hpp"

#include "k_nearest_neighbor_similarity.hpp"
#include "k_nearest_neighbor_distance.hpp"

using namespace_k_nearest_neighbor_similarity::KNeighborsSimilarity;
using namespace_k_nearest_neighbor_distance::KNeighborsDistance;


using namespace std;

int main(int argc, const char* argv[]) {
  vector<vector<double> > distance;
  GetMetricDistance(distance);

  KNeighborsSimilarity k_similarity;
  KNeighborsDistance k_distance(&distance);

  return 0;
}
