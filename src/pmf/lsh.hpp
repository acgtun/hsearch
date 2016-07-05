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
  LSH(const uint32_t& dimension)
      : generator(rd()),
        m_dimension(dimension),
        bucket_width(1.0),
        m_uniform(uniform_real_distribution<double>(-10, 10)),
        m_uniform_width(uniform_real_distribution<double>(0, bucket_width)),
        unit_vector(dimension, 0) {
    for (uint32_t i = 0; i < dimension; ++i) {
      unit_vector[i] = m_uniform(generator);
    }
    double length = DotProduct(unit_vector, unit_vector);
    length = sqrt(length);

    for (uint32_t i = 0; i < dimension; ++i) {
      unit_vector[i] = unit_vector[i] / length;
    }
  }

  int HashBucketIndex(const vector<double>& point);
  double DotProduct(const vector<double>& px, const vector<double>& py);

 private:
  random_device rd;
  default_random_engine generator;
  uint32_t m_dimension;
  double bucket_width;
  uniform_real_distribution<double> m_uniform;
  uniform_real_distribution<double> m_uniform_width;
  vector<double> unit_vector;
};

#endif // LSH_H
