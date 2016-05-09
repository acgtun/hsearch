#ifndef PCLUSTER_UNION_FIND_HPP_
#define PCLUSTER_UNION_FIND_HPP_

#include "util.hpp"

class UnionFind {
 public:
  UnionFind(const vector<uint32_t>& _protienIDS, const ProteinDB& _proteinDB);
  ~UnionFind();

  int FindRoot(uint32_t& x);
  void JoinUnion(uint32_t& x, uint32_t& y);
  void ProteinClustering();

  unordered_map<uint32_t, uint32_t> root;

  const vector<uint32_t>& protienIDS;
  const ProteinDB& proteinDB;
};

#endif /* PCLUSTER_UNION_FIND_HPP_ */
