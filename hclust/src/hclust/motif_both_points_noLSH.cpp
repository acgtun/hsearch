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
};

typedef unordered_map<string, vector<uint32_t> > HashTable;

double PairwiseDistance(const Point& a, const Point& b) {
  double dis = 0.0, r = 0.0;
  for (uint32_t i = 0; i < DIMENSION; ++i) {
    r = a.data[i] - b.data[i];
    dis += r * r;
  }
  return sqrt(dis);
}

void Search(const vector<Point>& kmers, const vector<Point>& centers,
            const vector<string>& kmer_names,
            const vector<string>& center_names, const double& hash_R,
            const string& output_file) {
  string nxot = output_file;
  nxot += "notlessthan.txt";
  ofstream fnot(nxot.c_str());
  ofstream fout(output_file.c_str());
  for (uint32_t i = 0; i < centers.size(); ++i) {
    for (uint32_t j = 0; j < kmers.size(); ++j) {
      double dis = PairwiseDistance(kmers[j], centers[i]);
      if (dis > hash_R) {
        fnot << center_names[i] << " " << kmer_names[j] << " " << dis << endl;
        continue;
      }
      fout << center_names[i] << " " << kmer_names[j] << " " << dis << endl;
    }
  }
  fnot.close();
  fout.close();
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

    /* kmer_file */
    string kmer_file; //both input are data points

    /* center file */
    string center_file;

    /* kmer length */
    uint32_t kmer_length = 25;

    /* distance threshold */
    double hash_R = 200;

    /* output file */
    string output_file;
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "cluster kmers to motifs",
        "");
    opt_parse.add_opt("db", 'd', "protein database file", true,
                      kmer_file);
    opt_parse.add_opt("center", 'c', "centers from Pfam database", true,
        center_file);
    opt_parse.add_opt("len", 'l', "kmer length", true,
        kmer_length);
    opt_parse.add_opt("threshold", 'T', "kmer threshold", true,
        hash_R);
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
    KMERLENGTH = kmer_length;
    DIMENSION = AACoordinateSize * KMERLENGTH;
    string line;
    cout << "Read Kmers..." << endl;
    vector<string> kmer_names;
    vector<Point> kmers;
    ifstream fin(kmer_file.c_str());
    while(getline(fin, line)) {
      kmer_names.push_back(line);
      getline(fin, line);
      istringstream iss(line);
      Point point;
      for(uint32_t i = 0;i < DIMENSION;++i) {
        iss >> point.data[i];
      }
      kmers.push_back(point);
    }
    fin.close();

    cout << "Read Centers..." << endl;
    vector<Point> centers;
    vector<string> center_names;
    fin.open(center_file.c_str());
    while(getline(fin, line)) {
      center_names.push_back(line);
      getline(fin, line);
      istringstream iss(line);
      Point point;
      for(uint32_t i = 0;i < DIMENSION;++i) {
        iss >> point.data[i];
      }
      centers.push_back(point);
    }
    fin.close();
    cout << "number of kmers " << kmers.size() << endl;
    cout << "number of centers " << centers.size() << endl;
    clock_t start = clock();
    Search(kmers, centers, kmer_names, center_names, hash_R, output_file);
    printf ("Searching takes %lf seconds\n", (clock() - start)/ (double)CLOCKS_PER_SEC);
  } catch (const SMITHLABException &e) {
    fprintf(stderr, "%s\n", e.what().c_str());
    return EXIT_FAILURE;
  } catch (std::bad_alloc &ba) {
    fprintf(stderr, "ERROR: could not allocate memory\n");
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
