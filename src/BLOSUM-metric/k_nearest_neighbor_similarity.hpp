#ifndef K_NEAREST_NEIGHBOR_SIMILARITY_H_
#define K_NEAREST_NEIGHBOR_SIMILARITY_H_

#include "sdk.hpp"
#include "bio_util.hpp"

#include <vector>
#include <queue>
#include <algorithm>

using std::vector;

namespace namespace_k_nearest_neighbor_similarity {

#define NODE_NUM 5002
#define MAX_K  5001

class AdjLinkDP {
 public:
  int v, c, nxt;
};

class QueNodeDP {
 public:
  int c, tc, pre;
  bool operator<(const QueNodeDP & b) const {
    return c < b.c;
  }
};

class KNeighborsSimilarity {
 private:
  std::priority_queue<QueNodeDP> cand[NODE_NUM];
  int *dp[NODE_NUM], *pre[NODE_NUM];
  struct AdjLinkDP adj[100010], invAdj[100010];
  int adjN;
  int invAdjN;

  int nNumNode;
  vector<int> vNodeLabel;
  vector<vector<uint32_t> > vLinkList;
  vector<uint32_t> vLinkNodeSize;
  vector<vector<int> > vEdgeWeight;

 public:
  KNeighborsSimilarity() {
    for (uint32_t i = 0; i < NODE_NUM; ++i) {
      dp[i] = new int[MAX_K];
      pre[i] = new int[MAX_K];
    }
    adjN = 0;
    invAdjN = 0;
    nNumNode = 0;

    vNodeLabel.resize(NODE_NUM);
    vLinkList.resize(NODE_NUM);
    vEdgeWeight.resize(NODE_NUM);
    vLinkNodeSize.resize(NODE_NUM);
    for (uint32_t i = 0; i < vLinkList.size(); i++) {
      vLinkList[i].resize(70);
      vEdgeWeight[i].resize(70);
    }
    BuildDAG();
  }
  ~KNeighborsSimilarity() {
    //INFO("Release K longest path memory...");
    for (uint32_t i = 0; i < NODE_NUM; ++i) {
      delete[] dp[i];
      delete[] pre[i];
    }
  }

 private:
  void insert(AdjLinkDP* adj, int &adjN, int a, int b, int c);
  int Query(int u, int k);
  void ShowPath(int n, int k, vector<uint32_t>& pathTmp);
  int GetWeight(const char& queryAA, const char& AA);
  void BuildDAG();

 public:
  void UpdateWeight(const char* querySeq);
  bool FindCandidate(char* strCandPep, const int& i, int& rawScore);
};

}  // namespace_k_nearest_neighbor_similarity
#endif /* K_NEAREST_NEIGHBOR_SIMILARITY_H_ */
