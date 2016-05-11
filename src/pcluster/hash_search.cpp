#include <cmath>
#include <cstring>
#include <cstdio>
#include <ctime>
#include <climits>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <algorithm>
#include <functional>
#include <iterator>
#include <numeric>

#include "merge_unit.hpp"
#include "hash_search.hpp"
#include "weight.hpp"
#include "aa.hpp"
#include "n2a.hpp"

using namespace std;

const ushort ONEBYTE = 15;
const ushort TWOBYTE = 255;
const ushort THRBYTE = 4095;
const ushort FOUBYTE = 65535;

CHashSearch::CHashSearch() {
  // for any letter which is not in the 20 aa
  m_uMask = 10;
  m_uSeg = 8;
  m_unMer = 6;
  m_bFast = true;
  m_unTotalIdx = static_cast<uint>(pow(10.0, int(m_unMer)));

  fill_n(m_aChar2Code, 256, (m_uMask << 4));
  fill_n(m_aCode2Char, 256, m_uMask);
  fill_n(m_aCode2Ten, 256, m_uMask);
  // read group info in aa.h and build mapping array
  // default use murphy10s
  for (int i = 0; i < 500; ++i) {
    if ('\0' == murphy10s[i][0]) {
      break;
    }

    char* p = murphy10s[i];
    for (uint j = 0; j < strlen(p); ++j) {
      // use new format to fix alignment break in SEGed region
      // the first four bits: group id
      // then one bit: SEGed? 1 : 0
      // then three bits: offset in a group
      uint unIdx = (i << 4) + j + 1;
      m_aChar2Code[p[j]] = unIdx;
      m_aChar2Code[p[j] + 32] = unIdx;	// add lower-case character
      m_aCode2Char[unIdx] = p[j];
      m_aCode2Ten[unIdx] = i;

      unIdx |= m_uSeg;	// set SEGed bit
      m_aChar2Code[p[j] + 128] = unIdx;  // SEGed letter
    }
  }

  // build substitute matrix for compressed char
  fill_n((int*) m_aSubMatrix, 256 * 256, -5);
  for (uint i = 0; i < strlen(aAAAlph); ++i) {
    for (uint j = 0; j < strlen(aAAAlph); ++j) {
      if (i < 20 && j < 20) {
        m_aSubMatrix[m_aChar2Code[aAAAlph[i]]][m_aChar2Code[aAAAlph[j]]] =
            blosum62[i][j];
      }
    }
  }
  m_aChar2Code['.'] = (10 << 4);
  m_aCode2Char[(10 << 4)] = '.';
  m_aCode2Ten[(10 << 4)] = m_uMask;
  m_aCode2Char['-'] = '-';

  m_pBlastSig = NULL;
  m_pComptor = NULL;

  int LONGQUERY = 4096;
  m_vTrace.assign(LONGQUERY, vector<char>(LONGQUERY));
  m_vETrace.assign(LONGQUERY, vector<char>(LONGQUERY));
  m_vDTrace.assign(LONGQUERY, vector<char>(LONGQUERY));

  m_unTotalSeeds = 0;
  //m_unTotalQuery = 0;
  m_unTotalSubj = 0;

  m_llOutCum = 0;
  //m_llM8Cum = 0;
  //m_nSeqBase = 0;

  lnDTotalSeqs = 0;
  lnDTotalAa = 0;

  // for test on gap extension
  m_unGapExt = 0;

  // hssp
  m_vCriteria.assign(100, 0);
  for (int i = 1; i < 100; ++i) {
    float f = 290.15 * pow(i, -0.562);
    f = f * i / 100;
    m_vCriteria[i] = (int) ceil(f);
  }
  of.open("chen.out");
}

struct CompDbObj {
  CompDbObj(VUCHAR& vSeqs, VUINT& vLens, uint& unMer)
      : m_vSeqs(vSeqs),
        m_vLens(vLens),
        m_unMer(unMer) {
  }

  bool operator()(const uint pos1, const uint pos2) const {
    // init paras
    int nIdx1 = pos1 >> 11;
    int nLen1 = m_vLens[nIdx1 + 1] - m_vLens[nIdx1];
    int nOff1 = (pos1 & 0x7ff) + m_unMer;
    uchar* p1 = &m_vSeqs[m_vLens[nIdx1]] + nOff1;

    int nIdx2 = pos2 >> 11;
    int nLen2 = m_vLens[nIdx2 + 1] - m_vLens[nIdx2];
    int nOff2 = (pos2 & 0x7ff) + m_unMer;
    uchar* p2 = &m_vSeqs[m_vLens[nIdx2]] + nOff2;

    // comp
    int nDiff = 0;
    if (nLen1 - nOff1 >= 4 && nLen2 - nOff2 >= 4) {
      nDiff = 4;
    } else {
      nDiff =
          (nLen1 - nOff1) < (nLen2 - nOff2) ? (nLen1 - nOff1) : (nLen2 - nOff2);
    }
    for (int i = 0; i < nDiff; ++i) {
      if ((*(p1 + i) >> 4) != (*(p2 + i) >> 4)) {
        return (*(p1 + i) >> 4) < (*(p2 + i) >> 4);
      }
    }
    return ((nLen1 - nOff1) < (nLen2 - nOff2));
  }

  VUCHAR& m_vSeqs;
  VUINT& m_vLens;
  uint& m_unMer;
};

void CHashSearch::BuildProteinsIndex(const vector<uint32_t>& protienIDS,
                                     const ProteinDB& proteinDB) {
  vDHash.assign(m_unTotalIdx, VUINT());  // all k-mer of database
  vDFreq.assign(strlen(murphy10r), 0);
  vDWordCnts.assign(m_unTotalIdx, 0);
  unDMedian = 0;
  lnDTotalAa = 0;
  lnDTotalSeqs = protienIDS.size();

  vDLens.push_back(0);
  for (size_t i = 0; i < protienIDS.size(); ++i) {
    vDNames.push_back(proteinDB.pro_names[i]);
    for (size_t j = 0; j < proteinDB.pro_seqs[i].size(); ++j) {
      vDSeqs.push_back(proteinDB.pro_seqs[i][j]);
    }
    vDLens.push_back(vDSeqs.size());
  }

  // char to code
  lnDTotalAa += Encode(vDSeqs, vDFreq);
  fprintf(stderr, "[THE TOTAL NUMBER OF AA IN THIS GROUP IS %lu]\n",
          lnDTotalAa);

  for (uint i = 0; i < vDLens.size() - 1; ++i) {
    // -1 or no, I need to think about it
    for (uint j = vDLens[i]; j < vDLens[i + 1] - m_unMer; ++j) {
      int nIdx = Tran2Ten(vDSeqs, j);
      if (-1 != nIdx) {
        // the left 21 bits denotes the index,
        // the right 11 bits denotes the starting position of the seed
        vDHash[nIdx].push_back((i << 11) | (j - vDLens[i]));
      }
    }
  }

  vDComp.assign(m_unTotalIdx, VUSHORT());
  for (uint i = 0; i < vDHash.size(); ++i) {
    sort(vDHash[i].begin(), vDHash[i].end(),
         CompDbObj(vDSeqs, vDLens, m_unMer));
    for (uint j = 0; j < vDHash[i].size(); ++j) {
      uint pos1 = vDHash[i][j];
      int nIdx1 = pos1 >> 11;
      int nLen1 = vDLens[nIdx1 + 1] - vDLens[nIdx1];
      int nOff1 = (pos1 & 0x7ff) + m_unMer;
      uchar* p1 = &vDSeqs[vDLens[nIdx1]] + nOff1;

      int m = nLen1 - nOff1;
      int n = 0;
      ushort nSuff = 0;
      //string sSuff;
      if (m >= 4) {
        nSuff |= (m_aCode2Ten[p1[n]]) << 12;
        nSuff |= (m_aCode2Ten[p1[++n]]) << 8;
        nSuff |= (m_aCode2Ten[p1[++n]]) << 4;
        nSuff |= (m_aCode2Ten[p1[++n]]);
      } else if (m == 3) {
        nSuff |= (m_aCode2Ten[p1[n]]) << 12;
        nSuff |= (m_aCode2Ten[p1[++n]]) << 8;
        nSuff |= (m_aCode2Ten[p1[++n]]) << 4;
        nSuff |= ONEBYTE;
      } else if (m == 2) {
        nSuff |= (m_aCode2Ten[p1[n]]) << 12;
        nSuff |= (m_aCode2Ten[p1[++n]]) << 8;
        nSuff |= TWOBYTE;
      } else if (m == 1) {
        nSuff |= (m_aCode2Ten[p1[n]]) << 12;
        nSuff |= THRBYTE;
      } else if (m == 0) {
        nSuff |= FOUBYTE;
      }
      vDComp[i].push_back(nSuff);
    }
  }

  for (uint i = 0; i < vDHash.size(); ++i) {
    vDWordCnts[i] += vDHash[i].size();
  }

  //serialize
  for (uint i = 0; i < vDFreq.size(); ++i) {
    vDFreq[i] /= lnDTotalAa;
  }

  sort(vDWordCnts.begin(), vDWordCnts.end());
  unDMedian = vDWordCnts[m_unTotalIdx / 2];
}

void CHashSearch::ProteinSearching(const vector<uint32_t>& proteinIDS,
                                   const ProteinDB& proteinDB, bool bEvalue,
                                   bool bLogE, double dThr, int nMaxOut,
                                   int nMaxM8, int nQueryTypeq,
                                   bool bPrintEmpty, bool bGapExt, bool bAcc,
                                   bool bHssp, int nMinLen, bool bXml,
                                   uint unDSize, uint unQSize, uint unMer) {
  m_bEvalue = bEvalue;
  m_bLogE = bLogE;
  if (m_bEvalue == true) {
    m_pComptor = new CompEval();
  } else {
    m_pComptor = new CompBits();
  }
  m_dThr = dThr;
  if (nMaxOut == -1) {
    m_nMaxOut = LLONG_MAX;
  } else {
    m_nMaxOut = abs(nMaxOut);
  }
  if (nMaxM8 == -1) {
    m_nMaxM8 = LLONG_MAX;
  } else {
    m_nMaxM8 = abs(nMaxM8);
  }

  m_bPrintEmpty = false;
  m_bGapExt = bGapExt;
  m_bAcc = bAcc;
  m_bHssp = bHssp;
  m_nMinLen = nMinLen;

  if (true == m_bFast) {
    m_unMutSeedLen = 10;
    m_vMutation.push_back(static_cast<uint>(pow(10.0, int(m_unMer - 4 - 1))));
    m_vMutation.push_back(static_cast<uint>(pow(10.0, int(m_unMer - 5 - 1))));
    m_vMutation.push_back(static_cast<uint>(pow(10.0, int(m_unMer - 3 - 1))));
    if (m_unMer > 6) {
      m_vMutation.push_back(static_cast<uint>(pow(10.0, int(m_unMer - 6 - 1))));
    }
  } else {
    m_unMutSeedLen = 9;
    m_vMutation.push_back(static_cast<uint>(pow(10.0, int(m_unMer - 3 - 1))));
    m_vMutation.push_back(static_cast<uint>(pow(10.0, int(m_unMer - 5 - 1))));
  }

  // set BlastStat
  InitAlignPara();
  // construct query package
  CDbPckg Db(vDHash, vDSeqs, vDLens, vDNames, vDComp, vDFreq, vDWordCnts,
             unDMedian);
  for (size_t i = 0; i < proteinIDS.size(); ++i) {
    vector<uchar> vQSeqs;
    vector<uint> vQLens;
    VNAMES vQNames;

    vQLens.push_back(0);
    vQNames.push_back(proteinDB.pro_names[i]);
    of << "Query: ";
    for (size_t j = 0; j < proteinDB.pro_seqs[i].size(); ++j) {
      vQSeqs.push_back(proteinDB.pro_seqs[i][j]);
      of << proteinDB.pro_seqs[i][j];
    }
    of << endl;
    vQLens.push_back(vQSeqs.size());
    Encode(vQSeqs);

    // construct query package
    CQrPckg Query(vQSeqs, vQLens, vQNames);
    int nSeqNum = 1;
    m_vOutIdx.assign(nSeqNum, CIndex());

    for (int nn = 0; nn < nSeqNum; ++nn) {
      m_vOutIdx[nn].m_llBeg = 0;
      m_vOutIdx[nn].m_nSize = 0;
    }
    Searching(Query, Db);
  }
}

void CHashSearch::Searching(CQrPckg& Query, CDbPckg& Db) {
  uint32_t L = Query.m_vLens[1] - Query.m_vLens[0];
  int nFoundHit = 0;
  MRESULT mRes;

  int nQrIdx = 0;

  // original length of query
  uint unQLen = Query.m_vLens[nQrIdx + 1] - Query.m_vLens[nQrIdx];
  if (unQLen < m_unMer) {
    return;
  }

  int nQOriLen = unQLen;

  // set up BlastStat
  m_vpBlastSig->blastComputeLengthAdjustmentComp(unQLen);

  uchar* pQ = &Query.m_vSeqs[0] + Query.m_vLens[nQrIdx];
  CAlnPckg QrAln(pQ, unQLen, 0);

  // build invalid index position
  vector<char> vValid(unQLen, 0);
  for (uint i = 0; i < unQLen; ++i) {
    vValid[i] = m_aCode2Ten[pQ[i]];
    if (0 != (m_uSeg & pQ[i])) {
      pQ[i] &= ~m_uSeg;
      vValid[i] = m_uMask;
    }
  }

  // for consistence with swift
  uint unPrvSdLen = 6;
  for (uint i = 0; i < unQLen - m_unMer; ++i) {
    uint unCnt = 0;
    // pick seed length
    uint unQSeedBeg = QrAln.m_unSeedBeg = i;
    int nSeed = Tran2Ten(QrAln, vValid);
    if (-1 == nSeed) {
      continue;
    }
    uint unLocalSeed = 0;
    uint unIdx = 0;
    // get local seed length
    if (m_bAcc == false) {
      uint unIncr = 0;
      double dFold = 0.0;
      int nLeft = unQLen - unQSeedBeg - m_unMer;
      uint unRng = nLeft + 1 >= 3 ? 3 : nLeft + 1;
      uint unFreq = Db.m_vWordCnts[nSeed];
      if (unFreq <= Db.m_unMedian) {
        unLocalSeed = m_unMer;
      } else {
        double dExpFreq = unFreq;
        for (unIncr = 1; unIncr < unRng; unIncr++) {
          if ((unIdx = vValid[unQSeedBeg + m_unMer + unIncr - 1]) != m_uMask) {
            dFold = Db.m_vFreq[unIdx];
          } else {
            //dFold = 1.0 / strlen(murphy10r);
            break;
          }
          dExpFreq *= dFold;
          if (dExpFreq <= Db.m_unMedian) {
            break;
          }
        }
        unLocalSeed = m_unMer + unIncr;
      }

      // if there is a unacceptable char, give up this seed
      if (m_uMask == unIdx) {
        continue;
      }

      if (unLocalSeed < unPrvSdLen - 1) {
        unLocalSeed = unPrvSdLen - 1;
      }

      // no enough letters
      if (unQSeedBeg + unLocalSeed > unQLen) {
        continue;
      }
    } else {
      unLocalSeed = 10;
      // no enough letters
      if (unQSeedBeg + unLocalSeed > unQLen) {
        continue;
      }
      for (uint i = m_unMer; i < unLocalSeed; ++i) {
        if ((unIdx = vValid[unQSeedBeg + i]) == m_uMask) {
          break;
        }
      }
      if (m_uMask == unIdx) {
        continue;
      }
    }

    //#################################################
    //  after 6 AA
    vector<uchar> vExtra(pQ + unQSeedBeg + m_unMer,
                         pQ + unQSeedBeg + unLocalSeed);
    for (uint idx = 0; idx < vExtra.size(); ++idx) {
      vExtra[idx] = m_aCode2Ten[vExtra[idx]];
    }

    if (!Db.m_vHash[nSeed].empty()) {
      int nCnt = ExtendSeq2Set(nSeed, unLocalSeed, vExtra, nQrIdx, QrAln,
                               nQOriLen, vValid, Db.m_vHash[nSeed], Db,
                               Query.m_vNames, Db.m_vNames, mRes);

      if (nCnt > 0) {
        unPrvSdLen = unLocalSeed;
      } else {
        unPrvSdLen = m_unMer;
      }
    }

    if (m_bAcc == true) {
      continue;
    }

    // mutation, pos: 4, 5, 3, (6)
    // check whether or not the length is enough
    if (unQLen < unQSeedBeg + m_unMutSeedLen) {
      continue;
    }
    // check non-aa char
    for (uint j = unQSeedBeg + unLocalSeed; j < unQSeedBeg + m_unMutSeedLen;
        ++j) {
      if ((unIdx = vValid[j]) == m_uMask) {
        break;
      }
    }
    if (m_uMask == unIdx) {
      continue;
    }

    vExtra.assign(pQ + unQSeedBeg + m_unMer, pQ + unQSeedBeg + m_unMutSeedLen);
    for (uint idx = 0; idx < vExtra.size(); ++idx) {
      vExtra[idx] = m_aCode2Ten[vExtra[idx]];
    }

    for (uint m = 0; m < m_vMutation.size(); ++m) {
      int nVal = (nSeed / m_vMutation[m]) % 10;
      int nMutIdx = nSeed - nVal * m_vMutation[m];
      for (int n = 0; n < 10; ++n) {
        if (nMutIdx == nSeed) {
          nMutIdx += m_vMutation[m];
          continue;
        }

        if (Db.m_vHash[nMutIdx].empty()) {
          nMutIdx += m_vMutation[m];
          continue;
        }

        int nCnt = ExtendSeq2Set(nMutIdx, m_unMutSeedLen, vExtra, nQrIdx, QrAln,
                                 nQOriLen, vValid, Db.m_vHash[nMutIdx], Db,
                                 Query.m_vNames, Db.m_vNames, mRes);
        nFoundHit += nCnt;
        unCnt += nCnt;

        nMutIdx += m_vMutation[m];
      }
    }

    /*********************************************************/
    if (6 == m_unMer && m_unMutSeedLen > m_unMer) {
      // mutation pos 6
      // think it as a mutation at pos 5, then do set intersection with current seed set
      //int nBase = (nQrIdx % 100000) * 10;
      int nNextNum = m_aCode2Ten[QrAln.m_pSeq[QrAln.m_unSeedBeg + m_unMer]];
      for (int i = 0; i < 10; ++i) {
        if (i == nNextNum) {
          // if the mutation is equal to the original next char
          continue;
        }

        // change the 6th position
        vExtra[0] = i;
        int nCnt = ExtendSeq2Set(nSeed, m_unMutSeedLen, vExtra, nQrIdx, QrAln,
                                 nQOriLen, vValid, Db.m_vHash[nSeed], Db,
                                 Query.m_vNames, Db.m_vNames, mRes);
        nFoundHit += nCnt;
        unCnt += nCnt;
      }
    }
  }

  PrintRes(mRes, Query, Db);
}

struct CompShortLow {
  bool operator()(const ushort& s1, const ushort& s2) {
    int nLen1 = 4;
    if ((s1 & ONEBYTE) == ONEBYTE) {
      --nLen1;
    }
    if ((s1 & TWOBYTE) == TWOBYTE) {
      --nLen1;
    }
    if ((s1 & THRBYTE) == THRBYTE) {
      --nLen1;
    }
    if ((s1 & FOUBYTE) == FOUBYTE) {
      --nLen1;
    }

    int nLen2 = 4;
    if ((s2 & ONEBYTE) == ONEBYTE) {
      --nLen2;
    }
    if ((s2 & TWOBYTE) == TWOBYTE) {
      --nLen2;
    }
    if ((s2 & THRBYTE) == THRBYTE) {
      --nLen2;
    }
    if ((s2 & FOUBYTE) == FOUBYTE) {
      --nLen2;
    }

    int nLen = nLen1 < nLen2 ? nLen1 : nLen2;
    if (0 == nLen) {
      return nLen1 < nLen2;
    }
    bool b = ((s1 >> ((4 - nLen) << 2)) == (s2 >> ((4 - nLen) << 2)));
    if (true == b) {
      return nLen1 < nLen2;
    } else {
      return ((s1 >> ((4 - nLen) << 2)) < (s2 >> ((4 - nLen) << 2)));
    }
  }
};

struct CompShortUp {
  bool operator()(const ushort& s1, const ushort& s2) {
    int nLen1 = 4;
    if ((s1 & ONEBYTE) == ONEBYTE) {
      --nLen1;
    }
    if ((s1 & TWOBYTE) == TWOBYTE) {
      --nLen1;
    }
    if ((s1 & THRBYTE) == THRBYTE) {
      --nLen1;
    }
    if ((s1 & FOUBYTE) == FOUBYTE) {
      --nLen1;
    }

    int nLen2 = 4;
    if ((s2 & ONEBYTE) == ONEBYTE) {
      --nLen2;
    }
    if ((s2 & TWOBYTE) == TWOBYTE) {
      --nLen2;
    }
    if ((s2 & THRBYTE) == THRBYTE) {
      --nLen2;
    }
    if ((s2 & FOUBYTE) == FOUBYTE) {
      --nLen2;
    }

    int nLen = nLen1 < nLen2 ? nLen1 : nLen2;
    if (0 == nLen) {
      return nLen1 < nLen2;
    }
    bool b = ((s1 >> ((4 - nLen) << 2)) == (s2 >> ((4 - nLen) << 2)));
    if (true == b) {
      return false;
      //return nLen1<nLen2;
    } else {
      return ((s1 >> ((4 - nLen) << 2)) < (s2 >> ((4 - nLen) << 2)));
    }
  }
};

int CHashSearch::ExtendSeq2Set(int nSeed, uint unLocalSeedLen,
                               vector<uchar>& vExtra, int nQSeqIdx,
                               CAlnPckg& QrAln, int nQOriLen,
                               vector<char>& vValid, VUINT& vDSet, CDbPckg& Db,
                               VNAMES& vQNames, VNAMES& vDNames,
                               MRESULT& mRes) {
  // find a proper range for the comparisons
  int nSt = 0;
  int nEd = 0;
  if (unLocalSeedLen > m_unMer) {
    ushort nExtra = 0;
    for (uint i = 0; i < vExtra.size(); ++i) {
      nExtra |= (vExtra[i]) << (12 - (i << 2));
    }
    for (int i = vExtra.size(); i < 4; ++i) {
      nExtra |= (ONEBYTE) << (12 - (i << 2));
    }

    VUSHORT::iterator itShort = lower_bound(Db.m_vComp[nSeed].begin(),
                                            Db.m_vComp[nSeed].end(), nExtra,
                                            CompShortLow());
    nSt = itShort - Db.m_vComp[nSeed].begin();

    if (static_cast<int>(vDSet.size()) == nSt) {
      return 0;
    }

    ushort s1 = Db.m_vComp[nSeed][nSt];
    ushort s2 = nExtra;
    int nLen1 = 4;
    if ((s1 & ONEBYTE) == ONEBYTE) {
      --nLen1;
    }
    if ((s1 & TWOBYTE) == TWOBYTE) {
      --nLen1;
    }
    if ((s1 & THRBYTE) == THRBYTE) {
      --nLen1;
    }
    if ((s1 & FOUBYTE) == FOUBYTE) {
      --nLen1;
    }

    int nLen2 = 4;
    if ((s2 & ONEBYTE) == ONEBYTE) {
      --nLen2;
    }
    if ((s2 & TWOBYTE) == TWOBYTE) {
      --nLen2;
    }
    if ((s2 & THRBYTE) == THRBYTE) {
      --nLen2;
    }
    if ((s2 & FOUBYTE) == FOUBYTE) {
      --nLen2;
    }

    int nLen = nLen1 < nLen2 ? nLen1 : nLen2;
    if (0 == nLen) {
      return 0;
    }
    bool b = ((s1 >> ((4 - nLen) << 2)) == (s2 >> ((4 - nLen) << 2)));
    if (true != b) {
      return 0;
    }

    itShort = upper_bound(Db.m_vComp[nSeed].begin(), Db.m_vComp[nSeed].end(),
                          nExtra, CompShortUp());
    nEd = itShort - Db.m_vComp[nSeed].begin();
  } else {
    nSt = 0;
    nEd = vDSet.size();
  }

  //sequence extension
  STAlnmnt stAlnmnt;
  for (int j = nSt; j < nEd; ++j) {
    if (vDNames[vDSet[j] >> 11] == "7719.ENSCINP00000006706") {
      //int zya = 0;
    }

    uint unDLen, unDSeedBeg;
    uchar* pD = GetSeq(Db.m_vSeqs, Db.m_vLens, Db.m_vNames, vDSet[j], unDLen,
                       unDSeedBeg);

    // no enough letters
    if (unDLen < unDSeedBeg + unLocalSeedLen) {
      continue;
    }
    // have the same previous char, so don't extend them at this time
    if (4 != vExtra.size()	// if this is a mutation case, do not allow to ignore it
    && m_bAcc == false && 0 != QrAln.m_unSeedBeg && 0 != unDSeedBeg
        && m_aCode2Ten[QrAln.m_pSeq[QrAln.m_unSeedBeg - 1]]
            == m_aCode2Ten[pD[unDSeedBeg - 1]]
        && m_uMask != vValid[QrAln.m_unSeedBeg - 1]) {
      continue;
    }

    CAlnPckg DbAln(pD, unDLen, unDSeedBeg);

    ResetResult(stAlnmnt);
    // for debug

    uint unLocalCopy = unLocalSeedLen;
    uchar* pQAlign = QrAln.m_pSeq + QrAln.m_unSeedBeg;
    uchar* pDAlign = DbAln.m_pSeq + DbAln.m_unSeedBeg;

    int ii = 0;
    for (; ii < unLocalCopy; ++ii) {
      stAlnmnt.nScore += m_aSubMatrix[pQAlign[ii]][pDAlign[ii]];
      if (pQAlign[ii] == pDAlign[ii]) {
        ++stAlnmnt.nMatch;
      }
    }
    // forward maximal extension
    uint unTempLen =
        QrAln.m_unLen - QrAln.m_unSeedBeg < DbAln.m_unLen - DbAln.m_unSeedBeg ?
            QrAln.m_unLen - QrAln.m_unSeedBeg :
            DbAln.m_unLen - DbAln.m_unSeedBeg;
    while (ii < unTempLen
        && m_aCode2Ten[pQAlign[ii]] == m_aCode2Ten[pDAlign[ii]]) {
      ++unLocalCopy;
      stAlnmnt.nScore += m_aSubMatrix[pQAlign[ii]][pDAlign[ii]];
      if (pQAlign[ii] == pDAlign[ii]) {
        ++stAlnmnt.nMatch;
      }
      ++ii;
    }
    // backward maximal extension
    uint unQSeed = QrAln.m_unSeedBeg;
    int nRange =
        QrAln.m_unSeedBeg < DbAln.m_unSeedBeg ?
            QrAln.m_unSeedBeg : DbAln.m_unSeedBeg;
    nRange = -nRange;
    ii = -1;
    while (ii >= nRange && m_aCode2Ten[pQAlign[ii]] == m_aCode2Ten[pDAlign[ii]]) {
      ++unLocalCopy;
      stAlnmnt.nScore += m_aSubMatrix[pQAlign[ii]][pDAlign[ii]];
      if (pQAlign[ii] == pDAlign[ii]) {
        ++stAlnmnt.nMatch;
      }
      --QrAln.m_unSeedBeg;
      --DbAln.m_unSeedBeg;
      --ii;
    }

    if (stAlnmnt.nScore >= UngapExtSCut && stAlnmnt.nMatch >= MinMatch4Exp) {
      AlignSeqs(nSeed, QrAln, DbAln, unLocalCopy, stAlnmnt);

      int nDSeqIdx = vDSet[j] >> 11;
      CalRes(nQSeqIdx, QrAln.m_pSeq, nQOriLen, QrAln.m_unSeedBeg, nDSeqIdx,
             DbAln.m_pSeq, DbAln.m_unSeedBeg, Db, unLocalCopy, stAlnmnt, mRes);
    }

    if (unQSeed != QrAln.m_unSeedBeg) {
      QrAln.m_unSeedBeg = unQSeed;
    }
  }

  return (nEd - nSt);
}

bool CHashSearch::AlignSeqs(int nSeed, CAlnPckg& QrAln, CAlnPckg& DbAln,
                            uint& unSeedLen, STAlnmnt& stAlnmnt) {
  uchar* pQAlign = QrAln.m_pSeq + QrAln.m_unSeedBeg;
  uchar* pDAlign = DbAln.m_pSeq + DbAln.m_unSeedBeg;
  int ext_f = 0;
  int ext_b = 0;
  int match_f = 0;
  int match_b = 0;

  int nScore0 = stAlnmnt.nScore;

  int nQLeft = 0;
  int nDLeft = 0;
  nQLeft = QrAln.m_unLen - QrAln.m_unSeedBeg - unSeedLen;
  nDLeft = DbAln.m_unLen - DbAln.m_unSeedBeg - unSeedLen;
  //if (nQLeft > 0 && nDLeft > 0)
  {
    pQAlign = QrAln.m_pSeq + QrAln.m_unSeedBeg + unSeedLen;
    pDAlign = DbAln.m_pSeq + DbAln.m_unSeedBeg + unSeedLen;
    stAlnmnt.nScore += AlignFwd(pQAlign, pDAlign, nQLeft, nDLeft, &ext_f,
                                &match_f, nScore0);
    stAlnmnt.nMatch += match_f;
    stAlnmnt.nQFwd += ext_f;
    stAlnmnt.nDFwd += ext_f;
  }

  nQLeft = QrAln.m_unSeedBeg - 1;
  nDLeft = DbAln.m_unSeedBeg - 1;
  //if (nQLeft > 0 && nDLeft > 0)
  {
    pQAlign = QrAln.m_pSeq;
    pDAlign = DbAln.m_pSeq;
    stAlnmnt.nScore += AlignBwd(pQAlign, pDAlign, nQLeft, nDLeft, &ext_b,
                                &match_b, nScore0);
    stAlnmnt.nMatch += match_b;
    stAlnmnt.nQBwd += ext_b;
    stAlnmnt.nDBwd += ext_b;
  }

  int hsplen = unSeedLen + ext_f + ext_b;

  stAlnmnt.vMode.push_back('s');
  stAlnmnt.vLen.push_back(hsplen);

  if (stAlnmnt.nScore < GapExtSCut) {
    //cout << "less than cut off" << endl;
    return true;
  }

  if (m_bAcc == false && true == m_bGapExt) {
    ++m_unGapExt;
    //cout << "1 hit ..." << endl;
    // vector for alignment path
    vector<char> vMode;
    vector<short> vLen;

    // forward gapped alignment
    int exta = 0;
    int extb = 0;
    int gap = 0;
    int nQOff = QrAln.m_unSeedBeg + unSeedLen + stAlnmnt.nQFwd;
    int nDOff = DbAln.m_unSeedBeg + unSeedLen + stAlnmnt.nDFwd;
    nQLeft = QrAln.m_unLen - nQOff;
    nDLeft = DbAln.m_unLen - nDOff;
    if (nQLeft > 2 && nDLeft > 2) {
      pQAlign = QrAln.m_pSeq + nQOff;
      pDAlign = DbAln.m_pSeq + nDOff;
      int n1 = AlignGapped(pQAlign, pDAlign, nQLeft, nDLeft, &exta, &extb,
                           &match_f, &gap, vMode, vLen);
      if (n1 > 0) {
        stAlnmnt.nScore += n1;
        stAlnmnt.nMatch += match_f;
        stAlnmnt.nQFwd += exta;
        stAlnmnt.nDFwd += extb;
        stAlnmnt.vMode.insert(stAlnmnt.vMode.end(), vMode.rbegin(),
                              vMode.rend());
        stAlnmnt.vLen.insert(stAlnmnt.vLen.end(), vLen.rbegin(), vLen.rend());
      }
    }

    // backward gapped alignment
    nQOff = QrAln.m_unSeedBeg - stAlnmnt.nQBwd;
    nDOff = DbAln.m_unSeedBeg - stAlnmnt.nDBwd;
    if (nQOff > 2 && nDOff > 2) {
      vector<uchar> vQ;
      vQ.reserve(nQOff);
      for (int i = nQOff - 1; i >= 0; --i) {
        vQ.push_back(QrAln.m_pSeq[i]);
      }
      pQAlign = &vQ[0];
      vector<uchar> vD;
      vD.reserve(nDOff);
      for (int i = nDOff - 1; i >= 0; --i) {
        vD.push_back(DbAln.m_pSeq[i]);
      }
      pDAlign = &vD[0];
      int n2 = AlignGapped(pQAlign, pDAlign, nQOff, nDOff, &exta, &extb,
                           &match_b, &gap, vMode, vLen);
      if (n2 > 0) {
        stAlnmnt.nScore += n2;
        stAlnmnt.nMatch += match_b;
        stAlnmnt.nQBwd += exta;
        stAlnmnt.nDBwd += extb;
        stAlnmnt.vMode.insert(stAlnmnt.vMode.begin(), vMode.begin(),
                              vMode.end());
        stAlnmnt.vLen.insert(stAlnmnt.vLen.begin(), vLen.begin(), vLen.end());
      }
    }
  }

  return true;
}

int CHashSearch::AlignFwd(uchar *queryseq, uchar *dataseq, uint len_queryseq,
                          uint len_dataseq, int *extl, int *match, int score0) {
  int i, j, l, s, maxs, ma;

  i = j = l = 0;
  ma = 0;
  maxs = s = score0;
  *extl = 0;
  *match = 0;
  while (i < len_queryseq && j < len_dataseq && s >= MINSCORE
      && s >= maxs - UngapExtDrop) {
    // uncompleted
    s += m_aSubMatrix[queryseq[i]][dataseq[j]];
    if (queryseq[i] == dataseq[j]) {
      ma++;
    }
    l++;
    if (s > maxs) {
      maxs = s;
      *extl = l;
      *match = ma;
    }
    i++;
    j++;
  }
  return maxs - score0;
}

int CHashSearch::AlignBwd(uchar *queryseq, uchar *dataseq, int pos1, int pos2,
                          int *extl, int *match, int score0) {
  int i, j, l, s, maxs, ma;

  i = pos1;
  j = pos2;
  l = 0;
  ma = 0;
  maxs = s = score0;
  *match = *extl = 0;
  while (i >= 0 && j >= 0 && s >= MINSCORE && s >= maxs - UngapExtDrop) {
    //	Skip stop codons
    //	uncompleted
    s += m_aSubMatrix[queryseq[i]][dataseq[j]];
    if (queryseq[i] == dataseq[j]) {
      ma++;
    }
    l++;
    if (s > maxs) {
      maxs = s;
      *extl = l;
      *match = ma;
    }
    i--;
    j--;
  }
  return maxs - score0;
}

int CHashSearch::AlignGapped(uchar *seq1, uchar *seq2, int M, int N, int *ext1,
                             int *ext2, int *match_len, int *gap,
                             vector<char>& vMode, vector<short>& vLen) {
  int i, j;
  int t, s, e, c, d, wa;
  int *CC = new int[N + 1];  //note N + 1
  int *DD = new int[N + 1];
  int g = GapIni;
  int h = GapExt;
  int m = g + h;  //gap-create + gap-extend
  int maxs, E1, E2, match;
  char trace_e, trace_d;
  maxs = E1 = E2 = match = 0;

  //forward-phase
  CC[0] = 0;
  DD[0] = -g;
  t = -g;

  int bb = 1;  //band_begin
  int be = int((GapExtDrop - GapIni) / GapExt);
  int bb_pre, be_pre;
  //these two parameters will be adjusted during the alignment based on the dropoff score

  vector<vector<char> > &trace = m_vTrace;
  vector<vector<char> > &etrace = m_vETrace;
  vector<vector<char> > &dtrace = m_vDTrace;

  // the aligning sequences may be longer than 4096
  bool bModify = false;
  int nMemory = trace.size();
  if (trace.size() - 1 < M) {
    int nSz = M + 1;
    bModify = true;
    trace.clear();
    etrace.clear();
    dtrace.clear();
    trace.assign(nSz, vector<char>(nSz));
    etrace.assign(nSz, vector<char>(nSz));
    dtrace.assign(nSz, vector<char>(nSz));
  }

  trace[0][0] = '0';
  for (j = 1; j <= N && j <= be; j++) {
    CC[j] = t = t - h;  //j - 1 ? or j; when j is used, check score is not the same as alignment score
    DD[j] = CC[j] - g;
    if (j == 1) {
      trace[0][j] = etrace[0][j] = 'E';
    } else {
      trace[0][j] = etrace[0][j] = 'e';
    }
    dtrace[0][j] = 'D';
  }  //global-alignment, with terminal penalty

  MaxGap = 100;
  for (i = 1; i <= M; i++) {
    bb_pre = bb;
    be_pre = be;
    if (be <= bb)
      break;  //band shrinks to zero
    s = CC[bb - 1];
    if (i == 1) {
      trace[i][bb - 1] = dtrace[i][bb - 1] = 'D';
      etrace[i][bb - 1] = 'E';
    } else {
      trace[i][bb - 1] = dtrace[i][bb - 1] = 'd';
      etrace[i][bb - 1] = 'e';
    }
    if (DD[bb - 1] - h > CC[bb - 1] - m) {
      c = DD[bb - 1] - h;
    } else {
      c = CC[bb - 1] - m;
    }
    CC[bb - 1] = DD[bb - 1] = c;  //update it with current row
    e = c - g;
    for (j = bb; j <= be && j <= N; j++) {
      trace_e = 'e';  //insertion extension
      if ((c = c - m) >= (e = e - h)) {
        e = c;
        trace_e = 'E';  //new insertion
      }  //insertion
      trace_d = 'd';  //deletion extension
      if ((c = CC[j] - m) >= (d = DD[j] - h)) {
        d = c;
        trace_d = 'D';  //new deletion
      }  //deletion
         //here   CC[j]==CC[i-1][j]   DD[j]==DD[i-1][j]

      wa = m_aSubMatrix[seq1[i - 1]][seq2[j - 1]];
      //sij[i - 1][j - 1]; //note i - 1, j - 1
      c = s + wa;  //s==CC[i-1][j-1], substitution
      trace[i][j] = 's';  //substitution

      if (e > c) {
        c = e;
        trace[i][j] = trace_e;
      }
      if (d > c) {
        c = d;
        trace[i][j] = trace_d;
      }
      etrace[i][j] = trace_e;
      dtrace[i][j] = trace_d;
      s = CC[j];  //important for next replace
      CC[j] = c;  //CC[i][j]
      DD[j] = d;  //DD[i][j]
      if (c > maxs) {
        E1 = i;
        E2 = j;
        maxs = c;
      }  //local -C
      else if (c < maxs - GapExtDrop && j > E2)  //score drops too much, stop filling this row, note j > E2
          {
        be = j;
        break;
      }
    }
    //after band_e, only allows insertion
    if (be < be_pre)
      continue;
    for (j = be + 1; j <= N; j++) {
      trace_e = 'e';  //insertion extension
      if ((c = c - m) > (e = e - h)) {
        e = c;
        trace_e = 'E';  //new insertion
      }  //insertion
      c = e;
      trace[i][j] = trace_e;
      etrace[i][j] = trace_e;

      s = CC[j];  //important for next replace
      CC[j] = c;  //CC[i][j]
      DD[j] = c - g;
      if (c > maxs) {
        E1 = i;
        E2 = j;
        maxs = c;
      }  //local -C
      else if (c < maxs - GapExtDrop)  //score drops too much, stop filling this row
          {
        be = j;
        break;
      }
    }
    //now infer new bb (starting from E2 going backward)
    for (j = E2; j >= bb; j--) {
      if (CC[j] < maxs - GapExtDrop) {
        bb = j;
        break;
      }
    }
  }

  *ext1 = E1;
  *ext2 = E2;

  delete[] CC;
  delete[] DD;

  //get alignment
  *match_len = 0;
  *gap = 0;

  if (maxs <= 0)
    return maxs;

  if (trace[E1][E2] != 's') {
    printf("E1 %d E2 %d, Not end with substitution %c\n", E1, E2,
           trace[E1][E2]);
    exit(1);
  }

  char mod = trace[E1][E2];
  i = E1;
  j = E2;
  vMode.clear();
  vLen.clear();
  while (mod != '0' && (!(i == 0 && j == 0))) {
    if (vMode.empty() || toupper(mod) != toupper(vMode.back())) {
      vMode.push_back(mod);
      vLen.push_back(0);
    }
    ++vLen.back();

    if (mod == 's') {
      if (seq1[i - 1] == seq2[j - 1])
        *match_len += 1;
      i -= 1;
      j -= 1;
      mod = trace[i][j];
    } else if (mod == 'D' || mod == 'd') {
      i -= 1;
      if (mod == 'D')
        mod = trace[i][j];
      else
        mod = dtrace[i][j];
      *gap += 1;
    } else {
      j -= 1;
      if (mod == 'E')
        mod = trace[i][j];
      else
        mod = etrace[i][j];
      *gap += 1;
    }
    if (i < 0 || j < 0) {
      cout << "This is a bug!" << endl;
      for (int m = 0; m < M; ++m) {
        cout << m_aCode2Char[seq1[m]];
      }
      cout << endl;
      for (int n = 0; n < N; ++n) {
        cout << m_aCode2Char[seq2[n]];
      }
      cout << endl;
      break;
    }
  }

  // reset the size of the buffer
  if (bModify == true) {
    trace.clear();
    etrace.clear();
    dtrace.clear();
    trace.assign(nMemory, vector<char>(nMemory));
    etrace.assign(nMemory, vector<char>(nMemory));
    dtrace.assign(nMemory, vector<char>(nMemory));
  }

  return maxs;
}

void CHashSearch::CalRes(int nQIdx, uchar* pQ, int nQOriLen, uint unQSeedBeg,
                         int nDIdx, uchar* pD, uint unDSeedBeg, CDbPckg& Db,
                         uint unLocalSeedLen, STAlnmnt& stAlnmnt,
                         MRESULT& mRes) {
  double dEValue = 0.0;
  if (m_bLogE == true) {
    dEValue = m_vpBlastSig->rawScore2ExpectLog(stAlnmnt.nScore);
  } else {
    dEValue = m_vpBlastSig->rawScore2Expect(stAlnmnt.nScore);
  }
  double dBits = m_vpBlastSig->rawScore2Bit(stAlnmnt.nScore);

  int nTotGap = 0;
  int nGapOpen = 0;
  int nTotAlnLen = 0;

  for (uint i = 0; i < stAlnmnt.vMode.size(); ++i) {
    nTotAlnLen += stAlnmnt.vLen[i];
    if ('s' != stAlnmnt.vMode[i]) {
      ++nGapOpen;
      nTotGap += stAlnmnt.vLen[i];
    }
  }

  // evalue criteria
  if (m_bHssp == false
      && !(stAlnmnt.nScore > SUMHSP_MINRAWSCORE
          || (m_bEvalue == true && dEValue <= m_dThr)
          || (m_bEvalue == false && dBits >= m_dThr))) {
    return;
  }
  // hssp criteria
  else if (m_bHssp == true
      && (nTotAlnLen < m_nMinLen || stAlnmnt.nMatch < m_vCriteria[nTotAlnLen])) {
    return;
  }

  // compute frame
  //cout << nQIdx << endl;
  int nQSt = 0;
  int nQEd = 0;

  nQSt = unQSeedBeg - stAlnmnt.nQBwd + 1;
  nQEd = unQSeedBeg + unLocalSeedLen + stAlnmnt.nQFwd;

  // print aligned sequences
  uint nAllc = nTotAlnLen > unLocalSeedLen ? nTotAlnLen : unLocalSeedLen;
  VUCHAR vQ;
  vQ.reserve(nAllc);
  VUCHAR vD;
  vD.reserve(nAllc);

  uchar* pQAligned = pQ + unQSeedBeg - stAlnmnt.nQBwd;
  uchar* pDAligned = pD + unDSeedBeg - stAlnmnt.nDBwd;

  if (0 == stAlnmnt.vMode.size()) {
    // only one hits
    vQ.insert(vQ.end(), pQAligned, pQAligned + unLocalSeedLen);
    vD.insert(vD.end(), pDAligned, pDAligned + unLocalSeedLen);
  } else if (1 == stAlnmnt.vMode.size()) {
    // only one hits
    vQ.insert(vQ.end(), pQAligned, pQAligned + stAlnmnt.vLen[0]);
    vD.insert(vD.end(), pDAligned, pDAligned + stAlnmnt.vLen[0]);
  } else {
    for (uint i = 0; i < stAlnmnt.vMode.size(); ++i) {
      char cMode = stAlnmnt.vMode[i];
      if ('s' == cMode) {
        vQ.insert(vQ.end(), pQAligned, pQAligned + stAlnmnt.vLen[i]);
        pQAligned += stAlnmnt.vLen[i];
        vD.insert(vD.end(), pDAligned, pDAligned + stAlnmnt.vLen[i]);
        pDAligned += stAlnmnt.vLen[i];
      } else if ('D' == cMode || 'd' == cMode) {
        vQ.insert(vQ.end(), pQAligned, pQAligned + stAlnmnt.vLen[i]);
        pQAligned += stAlnmnt.vLen[i];
        vD.insert(vD.end(), stAlnmnt.vLen[i], '-');
      } else if ('E' == cMode || 'e' == cMode) {
        vQ.insert(vQ.end(), stAlnmnt.vLen[i], '-');
        vD.insert(vD.end(), pDAligned, pDAligned + stAlnmnt.vLen[i]);
        pDAligned += stAlnmnt.vLen[i];
      }
    }
  }

  string sQ;
  string sD;
  Decode(vQ, sQ);
  Decode(vD, sD);

  string sInfo;
  for (uint i = 0; i < vQ.size(); ++i) {
    if (vQ[i] == vD[i]) {
      sInfo += sQ[i];
    } else if (m_aSubMatrix[vQ[i]][vD[i]] > 0) {
      sInfo += '+';
    } else {
      sInfo += ' ';
    }
  }

  MRESULT::iterator it = mRes.lower_bound(pair<int, int>(nQIdx, nDIdx));
  /****************************************************************/
  // for sum evalue, comment this
  // note: here, the hits are stored according to it's real query index, not 1->6 frame query index
  // store all results
  if (mRes.end() != it && (*it).first.first == nQIdx
      && (*it).first.second == nDIdx && (*it).second.nFrame == 0
      && (*it).second.nQSt == unQSeedBeg - stAlnmnt.nQBwd
      && (*it).second.nDSt == unDSeedBeg - stAlnmnt.nDBwd
      && (*it).second.nQEd == unQSeedBeg + unLocalSeedLen + stAlnmnt.nQFwd - 1
      && (*it).second.nDEd
          == unDSeedBeg + unLocalSeedLen + stAlnmnt.nDFwd - 1) {
    CHitUnit& st = (*it).second;
    if (st.dEValue > dEValue) {
      st.nScore = stAlnmnt.nScore;
      st.dBits = dBits;
      st.dEValue = dEValue;
      st.dIdent = stAlnmnt.nMatch * 100.0 / nTotAlnLen;
      st.nAlnLen = nTotAlnLen;
      st.nMismatch = nTotAlnLen - stAlnmnt.nMatch - nTotGap;
      st.nGapOpen = nGapOpen;
      st.nQBeg = nQSt;
      st.nQEnd = nQEd;
      st.sQ = sQ;
      st.sInfo = sInfo;
      st.sD = sD;
    }
  } else
  /****************************************************************/
  {
    MRESULT::iterator itTmp = mRes.insert(
        it, MRESULT::value_type(pair<int, int>(nQIdx, nDIdx), CHitUnit()));
    CHitUnit& st = (*itTmp).second;
    st.nQrLen = nQOriLen;
    st.nDbIdx = nDIdx;
    st.nDbLen = Db.m_vLens[nDIdx + 1] - Db.m_vLens[nDIdx];
    st.nScore = stAlnmnt.nScore;
    st.dBits = dBits;
    st.dEValue = dEValue;
    st.dIdent = stAlnmnt.nMatch * 100.0 / nTotAlnLen;
    st.nAlnLen = nTotAlnLen;
    st.nMismatch = nTotAlnLen - stAlnmnt.nMatch - nTotGap;
    st.nGapOpen = nGapOpen;
    st.nFrame = 0;
    st.nQSt = unQSeedBeg - stAlnmnt.nQBwd;
    st.nQEd = unQSeedBeg + unLocalSeedLen + stAlnmnt.nQFwd - 1;
    st.nQBeg = nQSt;
    st.nQEnd = nQEd;
    st.nDSt = unDSeedBeg - stAlnmnt.nDBwd;
    st.nDEd = unDSeedBeg + unLocalSeedLen + stAlnmnt.nDFwd - 1;
    st.sQ = sQ;
    st.sInfo = sInfo;
    st.sD = sD;
    //mRes.insert(it, MRESULT::value_type(pair<int, int>(nQIdx/m_nIdxScl, nDIdx), st));
  }
}

struct SetCompObj {
  bool operator()(const uint& p1, const uint& p2) const {
    return ((p1 >> 11) < (p2 >> 11)) && ((p1 & 0x7ff));
  }
} mySetComp;

void CHashSearch::PrintRes(MRESULT& mRes, CQrPckg& Query, CDbPckg& Db) {
  if (mRes.empty()) {
    return;
  }

  MIT it = mRes.begin();
  int nQrIdx = (*it).first.first;
  MRESULT::iterator itFind = mRes.end();
  vector<CHitUnit> vTemp;
  vTemp.reserve(distance(it, itFind));

  // for sum evalue, comment this
  int nDIdx = it->first.second;
  vTemp.push_back(it->second);
  int nSt = 0;
  MRESULT::iterator itTemp = it;
  ++itTemp;
  for (; itTemp != itFind; ++itTemp) {
    if (itTemp->first.second != nDIdx) {
      if (vTemp.size() - nSt > 1) {
        int nLen = Db.m_vLens[nDIdx + 1] - Db.m_vLens[nDIdx];
        SumEvalue(vTemp, nSt, vTemp.size(), nLen);
      }

      nDIdx = itTemp->first.second;
      nSt = vTemp.size();
    }
    vTemp.push_back(itTemp->second);
  }
  // process the last one
  if (vTemp.size() - nSt > 1) {
    int nLen = Db.m_vLens[nDIdx + 1] - Db.m_vLens[nDIdx];
    SumEvalue(vTemp, nSt, vTemp.size(), nLen);
  }

  if (0 == vTemp.size()) {
    return;
  }

  uint nMax = max(m_nMaxOut, m_nMaxM8);
  nMax = min((uint) vTemp.size(), nMax);
  vector<CHitUnit>::iterator itPrint = vTemp.begin() + nMax;
  partial_sort(vTemp.begin(), itPrint, vTemp.end(), ComptorWrapper(m_pComptor));

  int nBegStrAligned = 6;
  vector<CHitUnit>::iterator itSt = vTemp.begin();
  for (; itSt != itPrint; ++itSt) {
    CHitUnit& st = *itSt;
    if ((m_bEvalue == true && st.dEValue > m_dThr)
        || (m_bEvalue == false && st.dBits < m_dThr)) {
      break;
    }
    // note: here, the hits are stored according to it's real query index, not 1->6 frame query index

    st.sInfo.insert(0, 7, ' ');

    ++st.nDSt;
    ++st.nDEd;
    int nFac = 0;
    while (0 <= (st.nDbIdx - nFac - 1)
        && Db.m_vNames[st.nDbIdx] == Db.m_vNames[st.nDbIdx - nFac - 1]) {
      ++nFac;
    }
    st.nDSt = 1848 * nFac + st.nDSt;
    st.nDEd = 1848 * nFac + st.nDEd;

    st.sQName = Query.m_vNames[nQrIdx];
    st.sDName = Db.m_vNames[st.nDbIdx];
  }

  vTemp.resize(itSt - vTemp.begin());
  if (vTemp.size() != 0) {
    // remove possible redundancy in overlapped region
    // there should be only one redundant hit for each overlapped region
    for (uint i = 1; i < vTemp.size(); ++i) {
      if (vTemp[i].nScore == vTemp[i - 1].nScore
          && vTemp[i].sDName == vTemp[i - 1].sDName
          && vTemp[i].sQName == vTemp[i - 1].sQName
          && vTemp[i].nDSt == vTemp[i - 1].nDSt
          && vTemp[i].nDEd == vTemp[i - 1].nDEd
          && vTemp[i].nQBeg == vTemp[i - 1].nQBeg
          && vTemp[i].nQEnd == vTemp[i - 1].nQEnd) {
        vTemp[i].dEValue = 1000000;
      }
    }

    uint nf = 0;
    uint nl = vTemp.size();
    while (nf < nl) {
      if (vTemp[nf].dEValue == 1000000) {
        --nl;
        swap(vTemp[nf], vTemp[nl]);
      }
      ++nf;
    }
    vTemp.resize(nl);
  }

  PrintM8(vTemp);
  mRes.clear();
}

void CHashSearch::SumEvalue(vector<CHitUnit>& v, int nSt, int nEd, int nLen) {
  typedef vector<CHitUnit>::iterator STIT;
  STIT itSt = v.begin() + nSt;
  STIT itEd = v.begin() + nEd;
  // sort by nFrame
  sort(itSt, itEd, CompFrame());
  CHitUnit st;
  st.nFrame = 3;
  STIT itDir = lower_bound(itSt, itEd, st, CompFrame());
  // if there are more than one hit in one direction
  int nDisPos = distance(itSt, itDir);
  int nDisNeg = distance(itDir, itEd);
  if (nDisPos > 1 || nDisNeg > 1) {
    vector<CHitUnit> vRes;
    STIT itStart = itSt;
    STIT itEnd = itDir;
    for (int i = 0; i < 2; ++i) {
      if (distance(itStart, itEnd) == 0) {
        itStart = itEnd;
        itEnd = itEd;
        continue;
      } else if (distance(itStart, itEnd) == 1) {
        if ((m_bEvalue == true && itStart->dEValue <= m_dThr)
            || (m_bEvalue == false && itStart->dBits >= m_dThr)) {
          vRes.push_back(*itStart);
        }
        itStart = itEnd;
        itEnd = itEd;
        continue;
      }
      // sort by score and start position of query
      sort(itStart, itEnd, CompQSt());
      stable_sort(itStart, itEnd, ComptorWrapper(m_pComptor));
      // check overlap and logevalue
      vector<CHitUnit> vNew;
      vNew.push_back(*itStart);
      for (STIT itTemp = itStart + 1; itTemp != itEnd; ++itTemp) {
        int nHalfLen = (itTemp->nQEd - itTemp->nQSt + 1) >> 1;
        int nOverlap = SUMHSP_OVERLAP < nHalfLen ? SUMHSP_OVERLAP : nHalfLen;
        if (itTemp->dEValue >= SUMHSP_MINEVALUE
            && itTemp->nScore <= SUMHSP_MINRAWSCORE) {
          continue;
        }
        bool bNonOvlp = true;
        for (STIT itIso = vNew.begin(); itIso != vNew.end(); ++itIso) {
          if ((itTemp->nQSt <= itIso->nQEd - nOverlap
              && itTemp->nQEd >= itIso->nQSt + nOverlap)
              || (itIso->nQSt <= itTemp->nQEd - nOverlap
                  && itIso->nQEd >= itTemp->nQSt + nOverlap)) {
            bNonOvlp = false;
            break;
          }
        }
        if (true == bNonOvlp) {
          vNew.push_back(*itTemp);
        }
      }
      if (vNew.size() == 1) {
        if ((m_bEvalue == true && vNew[0].dEValue <= m_dThr)
            || (m_bEvalue == false && vNew[0].dBits >= m_dThr)) {
          vRes.push_back(vNew[0]);
        }
        //continue;
      } else {
        // calculate the sum of evalue
        double aRawScore[DEFAULT_SCORE_TOP];
        int nNo = 0;
        for (; nNo < 5 && nNo < vNew.size(); ++nNo) {
          aRawScore[nNo] = vNew[nNo].nScore;
        }

        if (m_bEvalue == true) {
          double dTmp = m_vpBlastSig->sumScore2Expect(nNo, aRawScore, nLen);
          double dSumEvalue = -10000.00;
          if (0 != dTmp) {
            dSumEvalue = log(dTmp) / LOG10;
          }
          if (m_bLogE == false) {
            dSumEvalue = dTmp;
          }
          // modify the logevalue
          if (dSumEvalue < m_dThr) {
            for (uint i = 0; i < vNew.size(); ++i) {
              vNew[i].dEValue = dSumEvalue;
            }
            vRes.insert(vRes.end(), vNew.begin(), vNew.end());
          }
        } else {
          double dTmpScore = m_vpBlastSig->sumScore(nNo, aRawScore, nLen);
          double dTmpBits = m_vpBlastSig->rawScore2Bit(dTmpScore);
          // modify the logevalue
          if (dTmpBits >= m_dThr) {
            for (uint i = 0; i < vNew.size(); ++i) {
              vNew[i].dBits = dTmpBits;
            }
            vRes.insert(vRes.end(), vNew.begin(), vNew.end());
          }
        }
      }
      itStart = itEnd;
      itEnd = itEd;
    }
    // replace
    if (!vRes.empty()) {
      v.erase(itSt, itEd);
      v.insert(v.begin() + nSt, vRes.begin(), vRes.end());
    }
  }
}

void CHashSearch::PrintAln(vector<CHitUnit>& v, ofstream& of) {
  int nPrint = min((long long) v.size(), m_nMaxOut);
  for (int i = 0; i < nPrint; ++i) {
    CHitUnit& c = v[i];

    of << c.sQName << " vs " << c.sDName << " bits=" << c.dBits;
    if (m_bLogE == true) {
      of << " log(E-value)=" << c.dEValue;
    } else {
      of << " E-value=" << c.dEValue;
    }
    of << " identity=" << c.dIdent << "%" << " aln-len=" << c.nAlnLen
       << " mismatch=" << c.nMismatch << " gap-openings=" << c.nGapOpen
       << " nFrame=" << c.nFrame << "\n" << "Query:\t" << c.sQ << "\n"
       << "      \t" << c.sInfo << "\n" << "Sbjct:\t" << c.sD << "\n" << "\n";
  }
}

void CHashSearch::PrintM8(vector<CHitUnit>& v) {
  int nPrint = min((long long) v.size(), m_nMaxM8);
  for (int i = 0; i < nPrint; ++i) {
    CHitUnit& c = v[i];
    of << c.sQName << "\t" << c.sDName << setprecision(1)
       << setiosflags(ios::fixed) << "\t" << c.dIdent << "\t" << c.nAlnLen
       << "\t" << c.nMismatch << "\t" << c.nGapOpen << "\t" << c.nQBeg << "\t"
       << c.nQEnd << "\t" << c.nDSt << "\t" << c.nDEd;
    if (m_bLogE == true) {
      of << setprecision(1) << setiosflags(ios::fixed) << "\t" << c.dEValue;
    } else {
      if (c.dEValue < 0.01) {
        of << setprecision(1) << setiosflags(ios::scientific)
           << setiosflags(ios::fixed) << "\t" << c.dEValue;
        of << resetiosflags(ios::scientific);
      } else if (c.dEValue < 10.0) {
        of << setprecision(2) << setiosflags(ios::fixed) << "\t" << c.dEValue;
      } else {
        of << setprecision(0) << setiosflags(ios::fixed) << "\t" << c.dEValue;
      }
    }
    of << setprecision(1) << setiosflags(ios::fixed) << "\t" << c.dBits << "\n";
  }
}
