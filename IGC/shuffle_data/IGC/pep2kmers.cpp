#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <unordered_set>

using namespace std;

// g++ -std=c++0x pep2kmers.cpp -o test

int main(int argc, const char **argv) {
  ifstream fin("../data/IGC.pep");
  string peptide;
  size_t num_of_peptides = 0;
  int kmer_length = 10;
  uint64_t total_length = 0;
  uint32_t max_size = 0;
  sscanf(argv[1], "%d", &kmer_length);
  cout << "kmer_length " << kmer_length << endl;
  unordered_set<string> kmers;
  while (getline(fin, peptide)) {
    if (peptide.size() == 0)
      continue;
    if (peptide[0] == '>')
      continue;
    num_of_peptides++;

    if(num_of_peptides % 10000 == 0) cout << "num_of_peptides = " << num_of_peptides << endl;
    total_length += peptide.size();
    max_size = max_size > peptide.size() ? max_size : peptide.size();
    /*for (size_t i = 0; i <= peptide.size() - kmer_length; ++i) {
      kmers.insert(peptide.substr(i, kmer_length));
    }*/
  }
  fin.close();
  cout << "max size = " << max_size << endl;
  cout << "total length = " << total_length << endl;
  cout << "num_of_peptides = " << num_of_peptides << endl;

  //////////////////////////////////
  char output_file[100];
  sprintf(output_file, "kmer.gen_len%d.txt", kmer_length);
  ofstream fout(output_file);
  for (unordered_set<string>::iterator it = kmers.begin(); it != kmers.end();
      ++it) {
    fout << *it << endl;
  }
  fout.close();
  return 0;
}
