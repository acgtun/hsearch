#include "read_proteins.hpp"

#include "smithlab_os.hpp"
#include "OptionParser.hpp"

void ProteinDB::ReadFASTAFile() {
  ifstream fin(file_name.c_str());
  if (!fin) {
    throw SMITHLABException("cannot open input file " + file_name);
  }
  vector<char> sequence;
  string line;
  while (getline(fin, line)) {
    if (line[0] == '>') {
      if (sequence.size() != 0) {
        pro_seqs.push_back(sequence);
        sequence.clear();
      }
      uint32_t space_pos = line.find_first_of(' ');
      if (space_pos == string::npos) {
        pro_names.push_back(line.substr(1));
      } else {
        pro_names.push_back(line.substr(1, space_pos - 1));
      }
      continue;
    }
    for (uint32_t i = 0; i < line.size(); ++i) {
      if (AA20.find_first_of(line[i]) != string::npos) {
        sequence.push_back(toupper(line[i]));
      } else if (isalpha(line[i])) {
        sequence.push_back(AA20[rand() % 20]);  //todo
      }
    }
  }
  if (sequence.size() != 0) {
    pro_seqs.push_back(sequence);
  }

  fin.close();
  num_of_proteins = pro_seqs.size();
}

void ProteinDB::OutputProtein(const uint32_t& i) const {
  for (uint32_t j = 0; j < pro_seqs[i].size(); ++j) {
    cout << pro_seqs[i][j];
  }
  cout << endl;
}

