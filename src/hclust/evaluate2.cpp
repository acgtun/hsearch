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
#include <algorithm>

#include <unordered_set>
#include <unordered_map>

#include "smithlab_os.hpp"

using namespace std;

struct MOTIFRES {
  MOTIFRES(const string& _motif, const string& _protein, const double& _dis)
      : motif(_motif),
        protein(_protein),
        dis(_dis) {

  }
  string motif;
  string protein;
  double dis;
};

bool sortCMP(const MOTIFRES& a, const MOTIFRES& b) {
  if (a.motif == b.motif) {
    return a.protein < b.protein;
  }
  return a.motif < b.motif;
}

int CompMOTIT(const MOTIFRES& a, const MOTIFRES& b) {

  //is large
  if (a.motif == b.motif && a.protein == b.protein) {
    return 0;
  }

  if (a.motif == b.motif) {
    if (a.protein == b.protein)
      return 0;
    if (a.protein > b.protein)
      return 1;
    return -1;
  }

  if (a.motif > b.motif)
    return 1;
  if (a.motif < b.motif)
    return -1;
}

double weight(const double& dis) {
  if (dis > 49.38) {
    double w = dis / (2 * 49.38);
    if (w > 1)
      return 1;
    return dis / (2 * 49.38);
  }

  return 1 - dis / (2 * 49.38);
}

int main(int argc, const char *argv[]) {
  vector<MOTIFRES> brute_force;
  string line, motif, protein;
  double distance;
  ;
  ifstream fin2(argv[1]);
  cout << argv[1] << endl;
  while (fin2 >> motif >> protein >> distance) {
    brute_force.push_back(MOTIFRES(motif, protein, distance));
  }
  fin2.close();

  ////////////////////////////////
  ////////////////////////////////

 sort(brute_force.begin(), brute_force.end(), sortCMP);
 string out = argv[1];
 out += "sort.txt";
 ofstream fout(out.c_str());
 for(uint32_t i = 0;i < brute_force.size();++i) {
   fout << brute_force[i].motif << "\t" << brute_force[i].protein << "\t" << brute_force[i].dis << endl;
 }
 fout.close();
 return 0;


  string files;
  files = argv[2];
  vector<string> file_names;
  if (isdir(files.c_str())) {
    read_dir(files, file_names);
  } else {
    file_names.push_back(files);
  }
  for (uint32_t k = 0; k < file_names.size(); ++k) {
    cout << file_names[k] << endl;
  }
  for (uint32_t k = 0; k < file_names.size(); ++k) {
    // hclust
    cout << "file " <<  file_names[k] << endl;
    cout << "file " <<  file_names[k] << endl;
    cout << "file " <<  file_names[k] << endl;
    if(file_names[k].back() == '.') {
      continue;
    }
    cout << "file3 " <<  file_names[k] << endl;
    ifstream fin(file_names[k]);
    if(!fin.good()) {
      fin.close();
      continue;
    }
    cout << "file2 " <<  file_names[k] << endl;
    vector<MOTIFRES> hclust;
    while (fin >> motif >> protein >> distance) {
      hclust.push_back(MOTIFRES(motif, protein, distance));
    }
    fin.close();
    sort(hclust.begin(), hclust.end(), sortCMP);

    uint32_t i = 0, j = 0;
    double tp = 0.0, fn = 0.0;
    while (i < brute_force.size() && j < hclust.size()) {
      int cmp = CompMOTIT(brute_force[i], hclust[j]);
      if (cmp == 0) {
        tp += weight(brute_force[i].dis);
        i++;
        j++;
      } else if (cmp == 1) {
        j++;
      } else {
        fn += weight(brute_force[i].dis);
        i++;
      }
    }
    while (i < brute_force.size()) {
      fn += weight(brute_force[i].dis);
      i++;
    }

    cout << "ACCURACY: " << tp << " " << fn << " " << tp / (tp + fn)  << "\t" << file_names[k] << endl;
  }
  return 0;
}
