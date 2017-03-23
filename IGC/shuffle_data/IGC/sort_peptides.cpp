#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <unordered_set>

using namespace std;

// g++ -std=c++0x pep2kmers.cpp -o sort_peps

bool sortCMP(const pair<string, string>& a, const pair<string, string>& b) {
  return a.second < b.second;
}

int main(int argc, const char **argv) {
  ifstream fin("../data/IGC.pep");
  string name, peptide;
  size_t num_of_peptides = 0;
  vector<pair<string, string> > proteins;

  while (getline(fin, name)) {
    getline(fin, peptide);
    proteins.push_back(make_pair(name, peptide));
    num_of_peptides++;
    if(num_of_peptides % 1000000 == 0) {
      cout << num_of_peptides << endl;
    }
  }

  // sort(proteins.begin(), proteins.end(), sortCMP);
  //ofstream fout("peptides_sort.txt");
  size_t length_of_database = 0;
  size_t num
  for (size_t i = 0; i < proteins.size(); ++i) {
    //fout << proteins[i].first << endl;
   // fout << proteins[i].second << endl;
  }
 // fout.close();

  return 0;
}
