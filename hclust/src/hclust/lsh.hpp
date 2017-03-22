#ifndef _LSH_H
#define _LSH_H

#pragma once

#include "util.hpp"

class LSH {
 public:
  LSH(const uint32_t& dimension, const uint32_t& _hash_K = 4,
      const double& _hash_W = 1.0)
      : m_dimension(dimension),
        hash_K(_hash_K),
        hash_W(_hash_W),
        m_normal(normal_distribution<double>(0.0, 1.0)),
        m_uniform_width(uniform_real_distribution<double>(0, _hash_W)),
        a(hash_K, vector<double>(dimension, 0)),
        b(hash_K, 0.0) {
    random_device rd;
    default_random_engine generator(rd());
    for (uint32_t k = 0; k < hash_K; ++k) {
      for (uint32_t i = 0; i < dimension; ++i) {
        a[k][i] = m_normal(generator);
        //cout << a[k][i] << " ";
      }
      //cout << endl;
      b[k] = m_uniform_width(generator);
      //cout << b[k] << endl;
      //cout << "------------------" << endl;
    }
  }

  double DotProduct(const vector<double>& point,
                    const uint32_t& hash_K_id) const {
    double dot_product = 0;
    for (uint32_t i = 0; i < m_dimension; ++i) {
      dot_product += point[i] * a[hash_K_id][i];
      //cout << i << " "  << m_dimension << " x " << point[i] << " " << a[hash_K_id][i] << " "  <<   point[i] * a[hash_K_id][i] << endl;
    }
    //cout << "dot_product = " << dot_product << endl;
    return dot_product;
  }

  int HashBucketIndex(const vector<double>& point,
                      const uint32_t& hash_K_id) const {
    double val = DotProduct(point, hash_K_id) + b[hash_K_id];
    //cout << "val: " << val << endl;
    return floor(val / hash_W);
  }

  string HashKey(const vector<double>& point) const {
    string hash_value;
    for (uint32_t k = 0; k < hash_K; ++k) {
      int bucket = HashBucketIndex(point, k);
      hash_value += to_string(bucket);
    }
    //cout << "h " << hash_value << endl;
    return hash_value;
  }

 private:
  uint32_t m_dimension;
  uint32_t hash_K;
  double hash_W;
  normal_distribution<double> m_normal;
  uniform_real_distribution<double> m_uniform_width;
  vector<vector<double> > a;
  vector<double> b;
};

#endif // _LSH_H
