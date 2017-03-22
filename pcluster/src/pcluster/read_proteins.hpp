#ifndef _READONEPROTEIN_H_
#define _READONEPROTEIN_H_

#include "util.hpp"

class ProteinDB {
 public:
  ProteinDB(const string& _file_name)
      : file_name(_file_name) {
    ReadFASTAFile();
  }

  ~ProteinDB() {
#ifdef DEBUG
    for (uint32_t i = 0; i < pro_seqs.size(); ++i) {
      cout << pro_names[i] << endl;
      for (uint32_t j = 0; j < pro_seqs[i].size(); ++j) {
        cout << pro_seqs[i][j];
      }
      cout << endl;
    }
#endif
  }

  void ReadFASTAFile();
  void OutputProtein(const uint32_t& i) const;

  vector<vector<char> > pro_seqs;
  vector<string> pro_names;

  uint32_t num_of_proteins;
  string file_name;
};

#endif /* _READONEPROTEIN_H_ */
