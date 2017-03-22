/*
max size = 29409
total length = 2470523769
num_of_peptides = 9878647
*/

struct ProteinDB {
  /* protein name */
  vector<string> name;

  /* protein length */
  vector<uint32_t> length;

  /* all proteins are concatenated to one string, and stored
   * in vector<char> sequence. start_index indicates the start position
   * of the chromosome in vector<char> sequence. */
  vector<uint32_t> start_index;

  /* number of chromosomes in the genome */
  uint32_t num_of_proteins;

  /* The total number of proteins */
  uint32_t length_of_proteins;

  /* all proteins are concatenated to one string, and stored in sequence. */
  vector<char> sequence;

  uint32_t ProteinID(const uint32_t& pos) {
    uint32_t l = 0, h = start_index.size() - 1;
    while (l < h) {
      uint32_t m = (l + h + 1) / 2;
      if (pos >= start_index[m])
        l = m;
      else
        h = m - 1;
    }

    return l;
  }

  ProteinDB(const string& protein_file) {
    ifstream fin(protein_file.c_str());
    cout << "Read protein sequences from " << protein_file << endl;
    string peptide;
    srand(time(NULL));
    length_of_proteins = 0;
    num_of_proteins = 0;
    while (getline(fin, peptide)) {
      if (peptide.size() == 0)
        continue;
      if (peptide[0] == '>') {
        name.push_back(peptide.substr(1));
      } else {
        length.push_back(peptide.size());
        start_index.push_back(sequence.size());
        num_of_proteins++;
        length_of_proteins += peptide.size();
        for (uint32_t i = 0; i < peptide.size(); ++i) {
          int AA = base[peptide[i] - 'A'];
          if (AA == -1) {
            AA = rand() % 20;
          }
          sequence.push_back(AA20[AA]);
        }
      }
    }
    start_index.push_back(sequence.size());
    fin.close();
    cout << "number of proteins " << num_of_proteins << endl;
    cout << "total length " << length_of_proteins << endl;
  }
};
