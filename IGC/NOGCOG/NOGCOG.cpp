#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>

using namespace std;

#define MAX_LINE_LENGTH 1000

int main(int argc, const char **argv) {
  FILE * fin = fopen("./../IGC.annotation_OF.summary", "r");
  string ID;
  string Name;
  uint32_t Length;
  string Completeness;
  string GroupOrigin;
  string PhylumLevel;
  string GenusLevel;
  string KEGGAnnotation;
  string eggNOGAnnotation;
  string Sampleo_ccurence_frequency;
  string Individual_occurence_frequency;
  string KEGG_functional_category;
  string eggNOG_functional_category;

  uint64_t total_length = 0;
  uint64_t un_length = 0;
  int total = 0;
  int un = 0;
  char cline[MAX_LINE_LENGTH];
  ofstream fout("length.txt");
  ofstream fout_t("total_length.txt");
  while (fgets(cline, MAX_LINE_LENGTH, fin)) {
    cline[strlen(cline) - 1] = 0;
    istringstream iss(cline);
    iss >> ID >> Name >> Length >> Completeness >> GroupOrigin >> PhylumLevel
        >> GenusLevel >> KEGGAnnotation >> eggNOGAnnotation;
    total++;
    total_length += Length;
    fout_t << Length << endl;
    if (PhylumLevel == "unknown" && GenusLevel == "unknown"
        && KEGGAnnotation == "unknown" && eggNOGAnnotation == "unknown") {
      un++;
      un_length += Length;
      fout << Length << endl;
    }
  }
  fout_t.close();
  fout.close();
  cout << total << " " << un << endl;
  cout << un / (double) total << endl;
  cout << total_length << " " << un_length << endl;
  cout << un_length / (double) un << endl;

  return 0;
}
