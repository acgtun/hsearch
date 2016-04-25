// export CC=/usr/usc/gnu/gcc/4.8.1/bin/gcc
// export CXX=/usr/usc/gnu/gcc/4.8.1/bin/g++
// source /usr/usc/gnu/gcc/4.8.1/setup.sh

#ifndef KERNEL_LSH_H
#define KERNEL_LSH_H

#pragma once
#include <vector>
#include <random>

#include <cmath>
#include <cstdint>

using std::vector;

//struct feature_node {
//  uint32_t index;
//  double value;
//  feature_node() {
//  }
//
//  feature_node(const uint32_t& _index, const double& _value)
//      : index(_index),
//        value(_value) {
//  }
//};

class KLSH {
 public:
  KLSH(const uint32_t& max_feat_num, const uint32_t& hash_bit_num,
       const double& sigma);
  uint64_t GetHashValue(const vector<double>& pfeat);

 private:
  uint32_t m_max_feat_num;
  uint32_t m_hash_bit_num;
  double m_sigma;

  std::normal_distribution<double> m_normal;
  std::uniform_real_distribution<double> m_uniform_1;
  std::uniform_real_distribution<double> m_uniform_pi;

  // std::vector is not used because operator "=" is not allowed in vector
  vector<vector<double> > m_project_w;
  vector<double> m_project_t;
  vector<double> m_project_b;

  std::default_random_engine generator;
};

#endif // KERNEL_LSH_H
