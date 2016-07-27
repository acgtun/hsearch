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


struct KMER {
  KMER(const string& _name, const string& _seq, const Point& _point) :
    name(_name),
    seq(_seq),
    point(_point) {
  }
  string name;
  string seq;
  Point point;
};

struct LSHTable {
  vector<string> hash_keys;
  unordered_map<string, vector<uint32_t> > buckets;
};

void KmerToCoordinates(const string& kmer, Point& point) {
  size_t k = 0;
  for (size_t i = 0; i < kmer.size(); ++i) {
    int AA = base[kmer[i] - 'A'];
    for (size_t j = 0; j < AACoordinateSize; ++j) {
      point.data[k++] = coordinates[AA][j];
    }
  }
}

double PairwiseDistance(const Point& a, const Point& b) {
  double sum = 0;
  for (int i = 0; i < a.data.size(); ++i) {
    sum += (a.data[i] - b.data[i]) * (a.data[i] - b.data[i]);
  }
  return sqrt(sum);
}

void Evaluation(const vector<KMER>& kmers, const vector<Cluster>& clusters,
                const uint32_t& hash_K, const uint32_t& hash_L,
                const double& D_T, const double& gamma) {
  cout << "number of clusters by our methods " << clusters.size() << endl;

  vector<string> ground_truth(kmers.size());
  for (uint32_t i = 0; i < kmers.size(); ++i) {
    size_t pos = kmers[i].name.find("_motif");
    size_t pos_ = kmers[i].name.find_last_of('_');
    string motif_num = kmers[i].name.substr(pos + 6, pos_ - pos - 6);
    ground_truth[i] = motif_num;
  }

  vector<uint32_t> cluster_id(kmers.size(), 0);
  for (uint32_t i = 0; i < clusters.size(); ++i) {
    for (uint32_t j = 0; j < clusters[i].ids.size(); ++j) {
      cluster_id[clusters[i].ids[j]] = i;
    }
  }

  uint64_t tp = 0, tn = 0, fp = 0, fn = 0;
  for (uint32_t i = 0; i < kmers.size(); ++i) {
    if (i % 50000 == 0 && i != 0) {
      cout << i / (double) kmers.size() << " tp: " << tp << " tn: " << tn << " fp: " << fp
      << " fn: " << fn << endl;
      cout << "recall: " << (double) tp / ((double) tp + fn) << endl;
      cout << "precision: " << (double) tp / ((double) tp + fp) << endl;
    }
    // if(i > 50000) break;
    for (uint32_t j = i + 1; j < kmers.size(); ++j) {
      if (ground_truth[i] == ground_truth[j]) {
        if (cluster_id[i] == cluster_id[j])
        tp++;
        else
        fn++;
      } else {
        if (cluster_id[i] == cluster_id[j])
        fp++;
        else
        tn++;
      }
    }
  }
  double recall = (double) tp / ((double) tp + fn);
  double precision = (double) tp / ((double) tp + fp);
  double RI = ((double) tp + tn) / ((double) tp + fp + fn + tn);
  double F1 = 2 * recall * precision / (recall + precision);
  cout << "recall_precision:KLT " << hash_K << "\t"
  << hash_L << "\t" << D_T <<  "\t" << gamma << "\t"
  << recall << "\t"
  << precision << "\t"
  << RI << "\t"
  << F1 << "\t"
  << tp << "\t"
  << tn << "\t"
  << fp << "\t"
  << fn << endl;
}


void BuildLSHTalbes(const vector<KMER>& kmers, const uint32_t& hash_K,
                    const uint32_t& hash_L, const double& D_T, const double gamma,
                    vector<LSHTable>& lsh_tables) {
  cout << hash_K << endl;
  cout << hash_L << endl;
  cout << D_T << endl;
  cout << gamma << endl;
  cout << "BuildLSHTalbes... " << endl;
  Point point;
  vector<vector<uint32_t> > hash_values(kmers.size(),
      vector<uint32_t>(hash_K, 0));
  for(uint32_t l = 0;l < hash_L;++l) {
    for(uint32_t k = 0;k < hash_K;++k) {
      LSH lsh(DIMENSION, gamma);
      for (uint32_t i = 0; i < kmers.size(); ++i) {
        hash_values[i][k] = lsh.HashBucketIndex(kmers[i].point.data);
      }
    }

    LSHTable table;
    char hash_value_chr[100];
    table.hash_keys.resize(kmers.size());
    for (uint32_t i = 0; i < kmers.size(); ++i) {
      string hash_value;
      for (uint32_t j = 0; j < hash_K; ++j) {
        sprintf(hash_value_chr, "%u", hash_values[i][j]);
        hash_value += hash_value_chr;
      }
      table.hash_keys[i] = hash_value;
      table.buckets[hash_value].push_back(i);
    }
    lsh_tables.push_back(table);
  }
}

void Clustering(const vector<KMER>& kmers,
                const uint32_t& hash_K,
                const uint32_t& hash_L,
                const double& D_T,
                const double& gamma,
                const vector<LSHTable>& lsh_tables,
                const string& output_file) {
  cout << hash_K << endl;
  cout << hash_L << endl;
  cout << D_T << endl;
  cout << "Clustering... " << endl;
  vector<Cluster> clusters(kmers.size());
  for (uint32_t i = 0; i < kmers.size(); ++i) {
    clusters[i].AddPoint(i);
  }

  vector<int> merged(kmers.size(), 0);
  for (uint32_t i = 0; i < kmers.size(); ++i) {
    if (merged[i] != 0) {
      continue;
    }

    for (uint32_t l = 0; l < hash_L; ++l) {
      string key = lsh_tables[l].hash_keys[i];
      unordered_map<string, vector<uint32_t> >::const_iterator
        it = lsh_tables[l].buckets.find(key);
      const vector<uint32_t>& ids = it->second;
      for (uint32_t j = 0; j < ids.size(); ++j) {
        if (merged[ids[j]] != 0 || ids[j] == i) {
          continue;
        }
        if (PairwiseDistance(kmers[i].point, kmers[ids[j]].point) <= D_T) {
          clusters[i].AddCluster(clusters[ids[j]]);
          merged[ids[j]] = 1;
        }
      }
    }
  }

  vector<Cluster> new_clusters;
  ofstream fout(output_file.c_str());
  uint32_t cluster_id = 0;
  for (uint32_t i = 0; i < kmers.size(); ++i) {
    if (merged[i] == 0) {
      new_clusters.push_back(clusters[i]);
      fout << "#clusterid:" << cluster_id++ << endl;
      for(uint32_t j = 0;j < clusters[i].ids.size();++j) {
        fout << kmers[clusters[i].ids[j]].name << endl;
      }
    }
  }
  fout.close();

  Evaluation(kmers, new_clusters, hash_K, hash_L, D_T, gamma);
}

int main(int argc, const char *argv[]) {
  srand (time(NULL));try {
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
      fprintf(stdout, "[WELCOME TO PMF v%s]\n", pmf_version);
      fprintf(stdout, "[%s", argv[0]);
      for (int i = 1; i < argc; i++) {
        fprintf(stdout, " %s", argv[i]);
      }
      fprintf(stdout, "]\n");
    }

    /* kmers file */
    string kmers_file;

    /* kmer length */
    uint32_t len = 10;

    /* number of random lines for each LSH */
    uint32_t hash_K = 8;

    /* number of hash functions */
    uint32_t hash_L = 8;

    /* distance threshold */
    double D_T = 50;

    /* bucket width */
    double gamma = 50;

    /* output file */
    string output_file;
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "cluster kmers to motifs",
        "");
    opt_parse.add_opt("kmers", 'k', "kmers file", true,
        kmers_file);
    opt_parse.add_opt("len", 'l', "kmer length", true,
        len);
    opt_parse.add_opt("hash_K", 'K', "number of hash functions", true,
        hash_K);
    opt_parse.add_opt("hash_L", 'L', "number of hash functions", true,
        hash_L);
    opt_parse.add_opt("threshold", 't', "cluster center threshold", true,
        D_T);
    opt_parse.add_opt("gamma", 'g', "bucketwidth", true,
                      gamma);
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
    Point p;
    ifstream fin(kmers_file.c_str());
    string line, name;
    while(fin >> line) {
      if(line[0] == '>') {
        name = line;
        fin >> line;
        //cout << line << endl;
        KmerToCoordinates(line, p);
        kmers.push_back(KMER(name, line, p));
      }
    }
    fin.close();
    printf("The number of kmers is %u\n", kmers.size());
    vector<LSHTable> lsh_tables;
    BuildLSHTalbes(kmers, hash_K, hash_L, D_T, gamma, lsh_tables);
    Clustering(kmers, hash_K, hash_L, D_T, gamma, lsh_tables, output_file);
  } catch (const SMITHLABException &e) {
    fprintf(stderr, "%s\n", e.what().c_str());
    return EXIT_FAILURE;
  } catch (std::bad_alloc &ba) {
    fprintf(stderr, "ERROR: could not allocate memory\n");
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
