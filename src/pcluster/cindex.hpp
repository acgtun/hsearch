#ifndef __CINDEX_H__
#define __CINDEX_H__

class CIndex {
 public:
  CIndex() {
    m_llBeg = 0;
    m_nSize = 0;
  }

 public:
  long long m_llBeg;
  int m_nSize;
};

#endif
