#include "util.hpp"
#include "lsh.hpp"
#include "protein.hpp"

#include "smithlab_os.hpp"
#include "OptionParser.hpp"

#include <stdint.h>
#include <sstream>
#include <unordered_set>
#include <unordered_map>

using namespace std;

uint32_t DIMENSION = 0;
uint32_t KMERLENGTH = 0;

struct Point {
  Point()
      : data(DIMENSION, 0.0) {
  }
  vector<double> data;
  void Output(ofstream& fout) {
    fout << data[0];
    for(uint32_t i = 1;i < DIMENSION;++i) {
      fout << " " << data[i];
    }
    fout << endl;
  }
};


void Protein2Datapoints(const ProteinDB& prodb, const uint32_t& num_of_protein_out, const string& output_file) {
  cout << "protein to data points... " << endl;
  ofstream fout(output_file.c_str());
  Point point;
  uint32_t cnt = 0;
  srand(time(NULL));
  unordered_set<string> proteins;
  for (uint32_t i = 0; i < prodb.num_of_proteins; ++i) {
    cout << i  << " "  << prodb.num_of_proteins << endl;
    if(i >= num_of_protein_out) break;
    for (uint32_t j = 0; j <= prodb.length[i] - KMERLENGTH;  ) {
      uint32_t q = 0, pos = prodb.start_index[i] + j;
      string kmer;
      for (uint32_t k = 0; k < KMERLENGTH; ++k) {
        kmer += prodb.sequence[pos + k];
      }
      if(proteins.find(kmer) != proteins.end()) {
        int r = 30 + rand() % 20;
        j += r;
        continue;
      }
      proteins.insert(kmer);
      for(uint32_t l = 0;l < KMERLENGTH;++l) {
        int AA = base[kmer[l] - 'A'];
        for (size_t p = 0; p < AACoordinateSize; ++p) {
          point.data[q++] = coordinates[AA][p];
        }
      }
      istringstream iss(prodb.name[i]);
      string name;
      iss >> name;
      fout << name << "#" << i << "$" << j <<"@" << kmer << "*" << cnt << endl;
      point.Output(fout);
      pos++;
      cnt++;
      //if(cnt > 1000000) break;
      int r = 30 + rand() % 20;
      j += r;
    }
  }
}

int main(int argc, const char *argv[]) {
  srand (time(NULL));
  try {
    string command = argv[0];
    bool help_info = false;
    for (int i = 1; i < argc; i++) {
      command += " ";
      command += argv[i];
      if (strcmp(argv[i], "-help") == 0 || strcmp(argv[i], "-about") == 0
          || strcmp(argv[i], "-?") == 0) {
        help_info = true;
      }
    }

    if (argc > 1 && help_info == false) {
      /* show the command line one the screen */
      fprintf(stdout, "[WELCOME TO HSEARCH v%s]\n", hclust_version);
      fprintf(stdout, "[%s", argv[0]);
      for (int i = 1; i < argc; i++) {
        fprintf(stdout, " %s", argv[i]);
      }
      fprintf(stdout, "]\n");
    }

    /* protein sequences file */
    string protein_file;

    /* kmer length */
    uint32_t kmer_length = 25;

    /* num_of_protein_out */
    uint32_t num_of_protein_out = 0;


    /* output file */
    string output_file;
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "protein sequences to data points",
        "");
    opt_parse.add_opt("db", 'd', "protein database file", true,
        protein_file);
    opt_parse.add_opt("len", 'l', "kmer length", true,
        kmer_length);
    opt_parse.add_opt("nnn", 'n', "num of proteins out", true,
        num_of_protein_out);
    opt_parse.add_opt("output", 'o', "output file name", true, output_file);

    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      fprintf(stderr, "%s\n", opt_parse.help_message().c_str());
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested()) {
      fprintf(stderr, "%s\n", opt_parse.about_message().c_str());
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      fprintf(stderr, "%s\n", opt_parse.option_missing_message().c_str());
      return EXIT_SUCCESS;
    }
    /****************** END COMMAND LINE OPTIONS *****************/
    clock_t start = clock();
    KMERLENGTH = kmer_length;
    DIMENSION = AACoordinateSize * KMERLENGTH;
    ProteinDB prodb(protein_file);
    Protein2Datapoints(prodb, num_of_protein_out, output_file);
    printf ("It takes %lf seconds\n", (clock() - start)/ (double)CLOCKS_PER_SEC);
  } catch (const SMITHLABException &e) {
    fprintf(stderr, "%s\n", e.what().c_str());
    return EXIT_FAILURE;
  } catch (std::bad_alloc &ba) {
    fprintf(stderr, "ERROR: could not allocate memory\n");
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
