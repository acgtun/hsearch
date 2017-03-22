#include "k_nearest_neighbor_similarity.hpp"

namespace namespace_k_nearest_neighbor_similarity {

void KNeighborsSimilarity::BuildDAG() {
  /* construct nodes, there are totally 20*len + 2 nodes, the first node is source node,
   * and the last node is the destination node */
  int len = HASHAALEN, nodeID = 0;
  /* set nodes label */
  vNodeLabel[nodeID++] = 0;
  for (int i = 0; i < len; i++) {
    for (int j = 0; j < 20; j++) {
      vNodeLabel[nodeID++] = j + 1;
    }
  }
  vNodeLabel[nodeID++] = 0;

  /* set edges */
  nNumNode = nodeID;
  int size = nodeID - 20 - 1, base = 0;
  for (int i = 0; i < size; i++) {
    if (i % 20 == 1)
      base += 20;
    for (int j = 1; j <= 20; j++) {
      vLinkList[i][j - 1] = base + j;
    }
    vLinkNodeSize[i] = 20;
  }
  for (int i = size; i < nNumNode - 1; i++) {
    vLinkList[i][0] = nNumNode - 1;
    vLinkNodeSize[i] = 1;
  }
  vLinkNodeSize[nNumNode - 1] = 0;
}

void KNeighborsSimilarity::insert(AdjLinkDP * adj, int &adjN, int a, int b,
                                  int c) {
  adj[adjN].nxt = adj[a].nxt;
  adj[adjN].c = c;
  adj[adjN].v = b;
  adj[a].nxt = adjN;
  adjN++;
}

int KNeighborsSimilarity::Query(int u, int k) {
  if (dp[u][0] >= k)
    return dp[u][k];
  if (dp[u][0] == 0) {
    if (invAdj[u].nxt == -1) {
      if (u != 0)
        return INT_MIN;
      QueNodeDP tmp;
      tmp.c = 0;
      tmp.pre = -1;
      tmp.tc = 0;
      cand[u].push(tmp);
      dp[u][0]++;
      dp[u][1] = cand[u].top().c;
      pre[u][0]++;
      pre[u][1] = cand[u].top().pre;
      return dp[u][1];
    } else {
      for (int i = invAdj[u].nxt; i != -1; i = invAdj[i].nxt) {
        int v = invAdj[i].v;
        QueNodeDP tmp;
        int res = Query(v, 1);
        if (res == INT_MIN)
          continue;
        tmp.c = res + invAdj[i].c;
        tmp.pre = v + 1 * 10000;
        tmp.tc = invAdj[i].c;
        cand[u].push(tmp);
      }
      if (cand[u].size() == 0)
        return -1;
      dp[u][0]++;
      dp[u][1] = cand[u].top().c;
      pre[u][0]++;
      pre[u][1] = cand[u].top().pre;
      return dp[u][1];
    }
  } else {
    if (cand[u].empty() || invAdj[u].nxt == -1)
      return INT_MIN;
    QueNodeDP cur = cand[u].top(), tmp;
    cand[u].pop();
    int res = Query(cur.pre % 10000, cur.pre / 10000 + 1);
    if (res != INT_MIN) {
      tmp.c = res + cur.tc;
      tmp.pre = cur.pre % 10000 + (cur.pre / 10000 + 1) * 10000;
      tmp.tc = cur.tc;
      cand[u].push(tmp);
    }
    if (!cand[u].empty()) {
      dp[u][0]++;
      dp[u][dp[u][0]] = cand[u].top().c;
      pre[u][0]++;
      pre[u][pre[u][0]] = cand[u].top().pre;
      return dp[u][dp[u][0]];
    } else
      return INT_MIN;
  }
}

void KNeighborsSimilarity::ShowPath(int n, int k, vector<uint32_t> & pathTmp) {
  if (pre[n][k] == -1) {
    pathTmp.push_back(n);
    return;
  }
  ShowPath(pre[n][k] % 10000, pre[n][k] / 10000, pathTmp);
  pathTmp.push_back(n);
}

bool KNeighborsSimilarity::FindCandidate(char * strCandPep, const int & i,
                                         int & rawScore) {
  rawScore = Query(nNumNode - 1, i + 1);
  if (rawScore != INT_MIN) {
    vector < uint32_t > pathTmp;
    ShowPath(nNumNode - 1, i + 1, pathTmp);
    for (uint32_t i = 1; i < pathTmp.size() - 1; i++) {
      strCandPep[i - 1] = NODELABEL[vNodeLabel[pathTmp[i]]];
    }
    strCandPep[HASHAALEN] = 0;
    return true;
  }
  return false;

}

int KNeighborsSimilarity::GetWeight(const char & queryAA, const char & AA) {
  if (queryAA == ' ')
    return 0;
  return BLOSUM62[base[AA - 'A']][base[queryAA - 'A']];
}

void KNeighborsSimilarity::UpdateWeight(const char * querySeq) {
  /* set edges' weight */
  for (int i = 0; i < nNumNode; i++) {
    for (uint32_t j = 0; j < vLinkNodeSize[i]; j++) {
      vEdgeWeight[i][j] = GetWeight(NODELABEL[vNodeLabel[vLinkList[i][j]]],
                                    querySeq[(i + 19) / 20]);
    }
  }

  for (int i = 0; i < nNumNode; i++) {
    adj[i].nxt = -1;
    invAdj[i].nxt = -1;
  }
  adjN = invAdjN = nNumNode;

  for (int i = 0; i < nNumNode; i++) {
    for (uint32_t j = 0; j < vLinkNodeSize[i]; j++) {
      int a, b, c;
      a = i;
      b = vLinkList[i][j];
      c = vEdgeWeight[a][j];
      insert(adj, adjN, a, b, c);
      insert(invAdj, invAdjN, b, a, c);
    }
  }
  for (int i = 0; i < nNumNode; i++) {
    dp[i][0] = 0;
    pre[i][0] = 0;
    while (!cand[i].empty())
      cand[i].pop();
  }
}

}  // namespace namespace_k_nearest_neighbor_similarity
