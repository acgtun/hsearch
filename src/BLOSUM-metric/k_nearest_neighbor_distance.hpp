#ifndef K_NEAREST_NEIGHBOR_DISTANCE_H_
#define K_NEAREST_NEIGHBOR_DISTANCE_H_

#include "sdk.hpp"
#include "bio_util.hpp"

#include <vector>
#include <queue>
#include <algorithm>

using namespace std;

namespace namespace_k_nearest_neighbor_distance {

#define NODE_NUM 5002
#define MAX_K  5001

class AdjLinkDP {
 public:
  int v, nxt;
  double c;
};

class QueNodeDP {
 public:
  double c, tc;
  int pre;
  bool operator<(const QueNodeDP & b) const {
    return c > b.c;
  }
};

class KNeighborsDistance {
 private:
  std::priority_queue<QueNodeDP> cand[NODE_NUM];
  double *dp[NODE_NUM];
  int *pre[NODE_NUM];
  int dp0[NODE_NUM];
  struct AdjLinkDP adj[100010], invAdj[100010];
  int adjN;
  int invAdjN;

  int nNumNode;
  vector<int> vNodeLabel;
  vector<vector<uint32_t> > vLinkList;
  vector<uint32_t> vLinkNodeSize;
  vector<vector<double> > vEdgeWeight;
  const vector<vector<double> >* distance_matrix;

 public:
  KNeighborsDistance(const vector<vector<double> >* _distance_matrix) {
    distance_matrix = _distance_matrix;
    for (uint32_t i = 0; i < NODE_NUM; ++i) {
      dp[i] = new double[MAX_K];
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
  ~KNeighborsDistance() {
    //INFO("Release K longest path memory...");
    for (uint32_t i = 0; i < NODE_NUM; ++i) {
      delete[] dp[i];
      delete[] pre[i];
    }
  }

 private:
  void insert(AdjLinkDP* adj, int &adjN, int a, int b, double c);
  double Query(int u, int k);
  void ShowPath(int n, int k, vector<uint32_t>& pathTmp);
  double GetWeight(const char& queryAA, const char &AA);
  void BuildDAG();

 public:
  void UpdateWeight(const char* querySeq);
  bool FindCandidate(char* strCandPep, const int& i, double& rawScore);
};

}  // namespace namespace_k_nearest_neighbor_distance
#endif /* K_NEAREST_NEIGHBOR_DISTANCE_H_ */
