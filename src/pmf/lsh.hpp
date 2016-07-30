#ifndef LSH_H
#define LSH_H

#pragma once

#include "util.hpp"

class LSH {
 public:
  LSH(const uint32_t& dimension, const uint32_t& _hash_K = 4,
      const double& _hash_W = 1.0)
      : generator(rd()),
        m_dimension(dimension),
        hash_K(_hash_K),
        hash_W(_hash_W),
        m_normal(normal_distribution<double>(0.0, 1.0)),
        m_uniform_width(uniform_real_distribution<double>(0, _hash_W)),
        a(hash_K, vector<double>(dimension, 0)),
        b(hash_K, 0.0) {
    for (uint32_t k = 0; k < hash_K; ++k) {
      for (uint32_t i = 0; i < dimension; ++i) {
        a[k][i] = m_normal(generator);
      }
      b[k] = m_uniform_width(generator);
    }
  }

  double DotProduct(const vector<double>& point,
                    const uint32_t& hash_K_id) const;
  int HashBucketIndex(const vector<double>& point,
                      const uint32_t& hash_K_id) const;
  string HashKey(const vector<double>& point) const;

 private:
  random_device rd;
  default_random_engine generator;
  uint32_t m_dimension;
  uint32_t hash_K;
  double hash_W;
  normal_distribution<double> m_normal;
  uniform_real_distribution<double> m_uniform_width;
  vector<vector<double> > a;
  vector<double> b;
};

#endif // LSH_H
