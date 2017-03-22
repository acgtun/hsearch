#include "lsh.hpp"

template<class T>
inline uint32_t sign(T val) {
  return val >= 0 ? 1 : 0;
}

double Dot(const vector<double>& px, const vector<double>& py) {
  double sum = 0;
  for (uint32_t i = 0; i < px.size(); ++i) {
    sum += px[i] * py[i];
  }

  return sum;
}

KLSH::KLSH(const uint32_t& max_feat_num, const uint32_t& hash_bit_num,
           const double& sigma)
    : m_max_feat_num(max_feat_num),
      m_hash_bit_num(hash_bit_num),
      m_sigma(sigma),
      m_normal(std::normal_distribution<double>(0.0, sigma * sigma)),
      m_uniform_1(std::uniform_real_distribution<double>(-1.0, 1.0)),
      m_uniform_pi(std::uniform_real_distribution<double>(0.0, 2.0 * M_PI)),
      m_project_w(m_hash_bit_num, vector<double>(m_max_feat_num)),
      m_project_t(m_hash_bit_num),
      m_project_b(m_hash_bit_num) {
  for (uint32_t i = 0; i < m_hash_bit_num; ++i) {
    // init project_t, project_b
    m_project_t[i] = m_uniform_1(generator);
    m_project_b[i] = m_uniform_pi(generator);
    // init project_w (m_max_feat_num + 1, 1 more space for index -1)
    for (uint32_t j = 0; j < m_max_feat_num; ++j) {
      // svm index begin from 1
      m_project_w[i][j] = m_normal(generator);
    }
  }
}

uint64_t KLSH::GetHashValue(const vector<double>& pfeat) {
  double sum = 0;
  uint64_t hash_value = 0;
  for (uint32_t i = 0; i < m_hash_bit_num; ++i) {
    sum = static_cast<double>(Dot(pfeat, m_project_w[i]) + m_project_b[i]);
    hash_value |= static_cast<uint64_t>(sign(std::cos(sum) + m_project_t[i]))
        << i;
  }
  return hash_value;
}
