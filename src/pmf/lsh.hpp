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

using namespace std;

class LSH {
 public:
  LSH(const uint32_t& dimension, const double& bucket_width = 1.0)
      : generator(rd()),
        m_dimension(dimension),
        m_bucket_width(bucket_width),
        m_normal(normal_distribution<double>(0.0, 1.0)),
        m_uniform_width(uniform_real_distribution<double>(0, bucket_width)),
        a(dimension, 0),
        b(m_uniform_width(generator)) {
    for (uint32_t i = 0; i < dimension; ++i) {
      a[i] = m_normal(generator);
    }
  }

  int HashBucketIndex(const vector<double>& point);
  double DotProduct(const vector<double>& px, const vector<double>& py);
  double DotProductPlusRandom(const vector<double>& point);

 private:
  random_device rd;
  default_random_engine generator;
  uint32_t m_dimension;
  double m_bucket_width;
  normal_distribution<double> m_normal;
  uniform_real_distribution<double> m_uniform_width;
  vector<double> a;
  double b;
};

#endif // LSH_H
