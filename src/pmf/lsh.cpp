#include "lsh.hpp"

double LSH::DotProduct(const vector<double>& point,
                       const uint32_t& hash_K_id) const {
  double dot_product = 0;
  for (uint32_t i = 0; i < m_dimension; ++i) {
    dot_product += point[i] * a[hash_K_id][i];
  }

  return dot_product;
}

int LSH::HashBucketIndex(const vector<double>& point,
                         const uint32_t& hash_K_id) const {
  double val = DotProduct(point, hash_K_id) + b[hash_K_id];
  return floor(val / hash_W);
}

string LSH::HashKey(const vector<double>& point) const {
  string hash_value;
  for (uint32_t k = 0; k < hash_K; ++k) {
    int bucket = HashBucketIndex(point, k);
    hash_value += to_string(bucket);
  }
  return hash_value;
}
