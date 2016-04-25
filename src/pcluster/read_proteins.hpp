#ifndef _READONEPROTEIN_H_
#define _READONEPROTEIN_H_

#include "util.hpp"

#include "smithlab_os.hpp"
#include "OptionParser.hpp"

#define MAX_LINE_LENGTH 5000

class ReadOneProtein {
 public:
  ReadOneProtein(const string& _file_name)
      : file_name(_file_name) {
    fin.open(_file_name.c_str());
    if (!fin) {
      throw SMITHLABException("cannot open input file " + file_name);
    }
    cout << file_name << endl;
    num_of_proteins = 0;
    compelted = false;

    getline(fin, line);
    pre_name = GetProteinName(line);
    protein_ID_Names.insert(make_pair(num_of_proteins, pre_name));
    num_of_proteins++;
  }

  ~ReadOneProtein() {
//    for (unordered_map<uint32_t, string>::iterator it =
//        protein_ID_Names.begin(); it != protein_ID_Names.end(); ++it) {
//      cout << it->first << " " << it->second << endl;
//    }
    fin.close();
  }

  string GetProteinName(const string& input) {
    if (input[0] != '>') {
      throw SMITHLABException(
          "the protein file should be fasta format " + file_name);
    }

    uint32_t space_pos = input.find_first_of(' ');
    if (space_pos == string::npos) {
      return input.substr(1);
    } else {
      return input.substr(1, space_pos - 1);
    }
  }

  vector<char> GetNextProtein() {
    vector<char> seq;
    if (pre_name.size() == 0) {
      compelted = true;
      return seq;
    }

    while (getline(fin, line)) {
      if (line[0] == '>') {
        pre_name = GetProteinName(line);
        protein_ID_Names.insert(make_pair(num_of_proteins, pre_name));
        num_of_proteins++;

        return seq;
      }

      for (uint32_t i = 0; i < line.size(); ++i) {
        if (AA20.find_first_of(line[i]) != string::npos) {
          seq.push_back(toupper(line[i]));
        } else if (isalpha(line[i])) {
          seq.push_back(AA20[rand() % 20]);  //todo
        }
      }
    }

    pre_name.clear();
    compelted = true;
    return seq;
  }

  unordered_map<uint32_t, string> protein_ID_Names;

  string pre_name;
  uint32_t num_of_proteins;
  string file_name;
  ifstream fin;
  string line;
  bool compelted;
};

#endif /* _READONEPROTEIN_H_ */
