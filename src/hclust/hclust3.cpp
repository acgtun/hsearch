#include "util.hpp"
#include "lsh.hpp"

#include "smithlab_os.hpp"
#include "OptionParser.hpp"

#include <stdint.h>
#include <unordered_set>
#include <unordered_map>

using namespace std;

uint32_t DIMENSION = 0;

struct Point {
  Point()
      : data(DIMENSION, 0.0) {
  }
  vector<double> data;
};

struct Cluster {
  vector<uint32_t> ids;

  void AddPoint(const uint32_t& id) {
    ids.push_back(id);
  }

  void AddCluster(const Cluster& cluster) {
    for (uint32_t i = 0; i < cluster.ids.size(); ++i) {
      ids.push_back(cluster.ids[i]);
    }
  }
};
Point KmerToCoordinates(const string& kmer);
struct KMER {
  KMER(const string& _name, const string& _seq)
      : name(_name),
        seq(_seq) {
  }
  string name;
  string seq;
  Point point() const {
  return KmerToCoordinates(seq);
}
};

typedef unordered_map<string, vector<uint32_t> > HashTable;

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


void BuildLSHTalbe(const vector<KMER>& kmers, const vector<uint8_t>& merged,
                   const LSH& lsh, HashTable& lsh_table) {
  cout << "BuildLSHTalbes... " << endl;
  for (uint32_t i = 0; i < kmers.size(); ++i) {
    if (merged[i] == 2) {
      continue;
    }
    string key = lsh.HashKey(kmers[i].point().data);
    lsh_table[key].push_back(i);
  }
}

void Clustering(const vector<KMER>& kmers,
                const uint32_t& hash_K,
                const uint32_t& hash_L,
                const double& hash_W,
                const double& hash_R,
                const string& output_file) {
  cout << "Clustering... " << endl;
  vector<uint8_t> merged(kmers.size(), 0);
  // 0 unprocessed
  // 1 center
  // 2 has been added to other cluster
  vector<Cluster> clusters(kmers.size());
  for (uint32_t i = 0; i < kmers.size(); ++i) {
    clusters[i].AddPoint(i);
  }
  clock_t start = clock();
  for (uint32_t l = 0; l < hash_L; ++l) {

    HashTable lsh_table;
    LSH lsh(DIMENSION, hash_K, hash_W);
    BuildLSHTalbe(kmers, merged, lsh, lsh_table);
    for (HashTable::iterator it = lsh_table.begin(); it != lsh_table.end();
        ++it) {
      const vector<uint32_t>& ids = it->second;
      vector<uint32_t> centers;
      for (uint32_t i = 0; i < ids.size(); ++i) {
        if (merged[ids[i]] == 1) {
          centers.push_back(ids[i]);
        }
      }
      for (uint32_t i = 0; i < ids.size(); ++i) {
        if (merged[ids[i]] == 0) {
          for (uint32_t j = 0; j < centers.size(); ++j) {
            if (PairwiseDistance(kmers[ids[i]].point(), kmers[centers[j]].point())
                <= hash_R) {
              clusters[centers[j]].AddPoint(ids[i]);
              merged[centers[j]] = 1;  // to be the real center
              merged[ids[i]] = 2;
              break;
            }
          }
        }
        if (merged[ids[i]] == 0) {
          centers.push_back(ids[i]);
        }
      }
    }
  }
  printf("ClusteringTime takes %lf seconds\n",
         (clock() - start) / (double) CLOCKS_PER_SEC);

  ofstream fout(output_file.c_str());
  uint32_t cluster_id = 0;
  uint32_t num_of_kmers = 0;
  for (uint32_t i = 0; i < kmers.size(); ++i) {
    if (merged[i] == 1 || merged[i] == 0) {
      fout << "#clusterid:" << cluster_id++ << ":size" << clusters[i].ids.size() << endl;
      num_of_kmers += clusters[i].ids.size() ;
      for (uint32_t j = 0; j < clusters[i].ids.size(); ++j) {
        fout << kmers[clusters[i].ids[j]].name << endl;
      }
    }
  }
  cout << "num_of_kmers = " << num_of_kmers << endl;
  fout.close();
}


int main(int argc, const char *argv[]) {
//  for (int i = 0; i < 20; ++i) {
//    printf("{");
//    for (int j = 0; j < 20; ++j) {
//      double dis = 0;
//      for (int k = 0; k < AACoordinateSize; ++k) {
//        dis += (coordinates[i][k] - coordinates[j][k])
//            * (coordinates[i][k] - coordinates[j][k]);
//      }
//      if(j != 19) printf("% 4.6lf, ", dis);
//      else  printf("%3.6lf},\n", dis);
//    }
//  }
//  return 0;
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
      fprintf(stdout, "[WELCOME TO PMF v%s]\n", hclust_version);
      fprintf(stdout, "[%s", argv[0]);
      for (int i = 1; i < argc; i++) {
        fprintf(stdout, " %s", argv[i]);
      }
      fprintf(stdout, "]\n");
    }

    /* kmers file */
    string kmers_file;

    /* kmer length */
    uint32_t len = 25;

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
    opt_parse.add_opt("kmers", 'k', "kmers file", true,
        kmers_file);
    opt_parse.add_opt("len", 'l', "kmer length", true,
        len);
    opt_parse.add_opt("hash_K", 'K', "number of random lines", true,
        hash_K);
    opt_parse.add_opt("hash_L", 'L', "number of hash tables", true,
        hash_L);
    opt_parse.add_opt("window", 'W', "bucket width", true,
        hash_W);
    opt_parse.add_opt("threshold", 'T', "clustering threshold", true,
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
    DIMENSION = AACoordinateSize * len;
    vector<KMER> kmers;
    ifstream fin(kmers_file.c_str());
    string line, name;
    while(fin >> line) {
      if(line[0] == '>') {
        name = line.substr(1);
        fin >> line;
        kmers.push_back(KMER(name, line));
      }
    }
    fin.close();
    printf("The number of kmers is %u\n", kmers.size());
    vector<HashTable> lsh_tables;
    clock_t start = clock();
    Clustering(kmers, hash_K, hash_L, hash_W, hash_R, output_file);
    printf ("Clustering takes %lf seconds\n", (clock() - start)/ (double)CLOCKS_PER_SEC);
  } catch (const SMITHLABException &e) {
    fprintf(stderr, "%s\n", e.what().c_str());
    return EXIT_FAILURE;
  } catch (std::bad_alloc &ba) {
    fprintf(stderr, "ERROR: could not allocate memory\n");
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
