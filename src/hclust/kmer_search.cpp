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
typedef unordered_map<uint32_t, pair<uint32_t, double> > ProteinCenterID;

Point KmerToCoordinates(const string& kmer) {
  Point point;
  size_t k = 0;
  for (size_t i = 0; i < kmer.size(); ++i) {
    int AA = base[kmer[i] - 'A'];
    if (AA == -1) {
      AA = rand() % 20;
    }
    for (size_t j = 0; j < AACoordinateSize; ++j) {
      point.data[k++] = coordinates[AA][j];
    }
  }
  return point;
}

double PairwiseDistance(const Point& a, const Point& b) {
  double dis = 0.0, r = 0.0;
  for (int i = 0; i < DIMENSION; ++i) {
    r = a.data[i] - b.data[i];
    dis += r * r;
  }
  return sqrt(dis);
}

Point PointByPosition(const ProteinDB& prodb, const uint32_t& position) {
  Point point;
  uint32_t q = 0, pos =position;
  for (uint32_t k = 0; k < KMERLENGTH; ++k) {
    int AA = base[prodb.sequence[pos++] - 'A'];
    for (size_t p = 0; p < AACoordinateSize; ++p) {
      point.data[q++] = coordinates[AA][p];
    }
  }
  return point;
}

void BuildLSHTalbe(const ProteinDB& prodb, const LSH& lsh,
                   HashTable& lsh_table) {
  cout << "BuildLSHTalbes... " << endl;
  Point point;
  for (uint32_t i = 0; i < prodb.num_of_proteins; ++i) {
    //cout << i  << " "  << prodb.num_of_proteins << endl;
    for (uint32_t j = 0; j <= prodb.length[i] - KMERLENGTH; ++j) {
      uint32_t q = 0, pos = prodb.start_index[i] + j;
      for (uint32_t k = 0; k < KMERLENGTH; ++k) {
        int AA = base[prodb.sequence[pos] - 'A'];
        for (size_t p = 0; p < AACoordinateSize; ++p) {
          point.data[q++] = coordinates[AA][p];
        }
      }
      string key = lsh.HashKey(point.data);
      lsh_table[key].push_back(pos);
      pos++;
    }
  }
}

void Search(const vector<Point>& centers, const ProteinDB& prodb,
            const uint32_t& hash_K, const uint32_t& hash_L,
            const double& hash_W, const double& hash_R,
            const string& output_file) {
  // matches for proteins, first is center id, second is distance
  ProteinCenterID matches;
  for (uint32_t l = 0; l < hash_L; ++l) {
    clock_t start = clock();
    HashTable lsh_table;
    LSH lsh(DIMENSION, hash_K, hash_W);
    BuildLSHTalbe(prodb, lsh, lsh_table);
    printf("Build LSHTable l=%d takes %lf seconds\n", l,
           (clock() - start) / (double) CLOCKS_PER_SEC);
    start = clock();
    cout << "Searching..." << endl;
    for (uint32_t i = 0; i < centers.size(); ++i) {
      string key = lsh.HashKey(centers[i].data);
      HashTable::iterator it = lsh_table.find(key);
      if (it == lsh_table.end()) {
        continue;
      }
      for (uint32_t j = 0; j < it->second.size(); ++j) {
        double dis = PairwiseDistance(PointByPosition(prodb, it->second[j]),
                                      centers[i]);
        if (dis > hash_R) {
          continue;
        }

        ProteinCenterID::iterator it2 = matches.find(it->second[j]);
        if (it2 == matches.end()) {
          matches.insert(make_pair(it->second[j], make_pair(i, dis)));
        } else {
          if (it2->second.second > dis) {
            it2->second = make_pair(i, dis);
          }
        }
      }
    }
    printf("Clustering l=%d takes %lf seconds\n", l,
           (clock() - start) / (double) CLOCKS_PER_SEC);
    lsh_table.clear();
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

    /* center file */
    string center_file;

    /* kmer length */
    uint32_t kmer_length = 25;

    /* number of random lines for each LSH */
    uint32_t hash_K = 16;

    /* number of hash tables */
    uint32_t hash_L = 32;

    /* bucket width */
    double hash_W = 50;

    /* distance threshold */
    double hash_R = 200;

    /* output file */
    string output_file;
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "cluster kmers to motifs",
        "");
    opt_parse.add_opt("db", 'd', "protein database file", true,
        protein_file);
    opt_parse.add_opt("center", 'c', "centers from Pfam database", true,
        center_file);
    opt_parse.add_opt("len", 'l', "kmer length", true,
        kmer_length);
    opt_parse.add_opt("hash_K", 'K', "number of random lines", true,
        hash_K);
    opt_parse.add_opt("hash_L", 'L', "number of hash tables", true,
        hash_L);
    opt_parse.add_opt("window", 'W', "bucket width", true,
        hash_W);
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
    ProteinDB prodb(protein_file);

    cout << "Read Centers..." << endl;
    string line;
    vector<Point> centers;
    ifstream fin(center_file.c_str());
    while(getline(fin, line)) {
      istringstream iss(line);
      Point point;
      for(uint32_t i = 0;i < DIMENSION;++i) {
        iss >> point.data[i];
      }
      centers.push_back(point);
    }
    fin.close();

    clock_t start = clock();
    Search(centers, prodb, hash_K, hash_L, hash_W, hash_R, output_file);
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
