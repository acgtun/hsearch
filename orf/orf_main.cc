#include "./../util/option.h"
#include "orf.h"
#include "fasta_file.h"
#include "./../util/bio_util.h"

int main(int argc, const char* argv[]) {
  InitProgram(argc, argv);
  string query_file;
  Option::GetOption("-q", query_file);
  FastaFile queries(query_file);
  ORF orf;
  ofstream fout((query_file + "_translatedAA.fasta").c_str());
  for (usint32_t i = 0; i < queries.num_of_sequences_; i++) {
    vector < string > translatedAA;
    orf.orf6(queries.sequences_[i], translatedAA);
    for (usint32_t j = 0; j < translatedAA.size(); j++) {
      fout << queries.sequences_names_[i] << "_" << j << endl;
      fout << translatedAA[j] << endl;
    }
  }
  fout.close();

  return 0;
}
