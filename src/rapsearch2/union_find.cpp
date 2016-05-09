#include "union_find.hpp"

UnionFind::UnionFind(const vector<uint32_t>& _protienIDS,
                     const ProteinDB& _proteinDB)
    : protienIDS(_protienIDS),
      proteinDB(_proteinDB) {
  for (uint32_t i = 0; i < protienIDS.size(); ++i) {
    root.insert(make_pair(protienIDS[i], protienIDS[i]));
  }
}

UnionFind::~UnionFind() {

}

int UnionFind::FindRoot(uint32_t& x) {
  int t = x;
  while (t != root[t]) {
    t = root[t];
  }

  while (x != root[x]) {
    int tmp = root[x];
    root[x] = t;
    x = tmp;
  }

  return t;
}

void UnionFind::JoinUnion(uint32_t& x, uint32_t& y) {
  root[x] = root[y];
}

void UnionFind::ProteinClustering() {
  //Build Prortein Index

  // Search
  for (uint32_t i = 0; i < protienIDS.size(); ++i) {

  }

}
