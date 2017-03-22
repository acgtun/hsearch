#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include <set>
#include <random>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>

#include <unordered_set>
#include <unordered_map>

using namespace std;

int main(int argc, const char *argv[]) {
  unordered_map<string, set<string> > motif_proteins;
  //meme
  string line, motif, protein;
  ifstream fin(argv[1]);
  cout << argv[1] << endl;
  fin >> line;
  while(getline(fin, line)) {
    istringstream iss(line);
    iss >> motif >> protein;
    motif_proteins[motif].insert(protein);
  }
  fin.close();
  ////////////////////////////////
  ////////////////////////////////

  // hclust
  fin.open(argv[2]);
  cout << argv[2] << endl;
  double distance;
  unordered_map<string, set<string> > motif_proteins_hclust;
  while(fin >> motif >> protein >> distance) {
    motif_proteins_hclust[motif].insert(protein);
  }
  fin.close();

  set<string> motifs;
  for( unordered_map<string, set<string> >::iterator it = motif_proteins.begin();it != motif_proteins.end();++it) {
    motifs.insert(it->first);
  }
  for( unordered_map<string, set<string> >::iterator it = motif_proteins_hclust.begin();it != motif_proteins_hclust.end();++it) {
    motifs.insert(it->first);
  }
  uint32_t sum1 = 0, sum2 = 0;
  for(set<string>::iterator it = motifs.begin();it != motifs.end();++it) {
    if(motif_proteins.find(*it) != motif_proteins.end()) {
      sum1 += motif_proteins[*it].size();
    }

    if(motif_proteins_hclust.find(*it) != motif_proteins_hclust.end()) {
      sum2 += motif_proteins_hclust[*it].size();
    }
  }
  cout << "ACCURACY: " << sum1 << " " << sum2 << " " << sum2 / (double)sum1 << endl;
  return 0;
}
