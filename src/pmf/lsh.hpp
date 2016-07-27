// export CC=/usr/usc/gnu/gcc/4.8.1/bin/gcc
// export CXX=/usr/usc/gnu/gcc/4.8.1/bin/g++
// source /usr/usc/gnu/gcc/4.8.1/setup.sh

#ifndef LSH_H
#define LSH_H

#pragma once

#include <cmath>
#include <cstdint>

#include <vector>
#include <random>
#include <iostream>

using namespace std;

class LSH {
 public:
  LSH(const uint32_t& dimension, const double& bucket_width = 1.0,
      const uint32_t& hash_K = 4)
      : generator(rd()),
        m_dimension(dimension),
        m_bucket_width(bucket_width),
        m_normal(normal_distribution<double>(0.0, 1.0)),
        m_uniform_width(uniform_real_distribution<double>(0, bucket_width)),
        a(hash_K, vector<double>(dimension, 0)),
        b(hash_K, 0.0) {
    for (uint32_t k = 0; k < hash_K; ++k) {
      for (uint32_t i = 0; i < dimension; ++i) {
        a[k][i] = m_normal(generator);
        cout << a[k][i] << " ";
      }
      cout << endl;
      b[k] = m_uniform_width(generator);
      cout << "b=" << b[k] << endl;
    }
  }

  int HashBucketIndex(const vector<double>& point, const uint32_t& hash_K_id);
  double DotProduct(const vector<double>& px, const vector<double>& py);

 private:
  random_device rd;
  default_random_engine generator;
  uint32_t m_dimension;
  double m_bucket_width;
  normal_distribution<double> m_normal;
  uniform_real_distribution<double> m_uniform_width;
  vector<vector<double> > a;
  vector<double> b;
};

#endif // LSH_H
