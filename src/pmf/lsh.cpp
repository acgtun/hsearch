#include "lsh.hpp"

#include <iostream>
using namespace std;

double LSH::DotProduct(const vector<double>& px, const vector<double>& py) {
  double dot_product = 0;
  for (uint32_t i = 0; i < m_dimension; ++i) {
    dot_product += px[i] * py[i];
  }

  return dot_product;
}

int LSH::HashBucketIndex(const vector<double>& point,
                         const uint32_t& hash_K_id) {
  return static_cast<int>((DotProduct(a[hash_K_id], point) + b[hash_K_id])
      / m_bucket_width);
}
