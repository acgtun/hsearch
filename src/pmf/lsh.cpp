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

int LSH::HashBucketIndex(const vector<double>& point) {
  //double b = m_uniform_width(generator);
  double len = DotProduct(a, point);
  //cerr << len + b << endl;
  //cout << "bucket_width " << m_bucket_width << endl;
  return static_cast<int>((len + b) / m_bucket_width);
}

double LSH::DotProductPlusRandom(const vector<double>& point) {
  //double b = m_uniform_width(generator);
  double len = DotProduct(a, point);
  return len + b;
}
