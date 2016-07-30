#include <math.h>
#include <stdint.h>

#include "util.hpp"
#include "lsh.hpp"

#include "smithlab_os.hpp"
#include "OptionParser.hpp"

#include <unordered_set>
#include <unordered_map>

using namespace std;

uint32_t DIMENSION = 0;
typedef unordered_map<string, vector<uint32_t> > HashTable;

struct Point {
  Point()
      : data(DIMENSION, 0.0) {
  }
  vector<double> data;
};

struct Cluster {
  void AddPoint(const uint32_t& id) {
    ids.push_back(id);
  }

  void AddCluster(const Cluster& cluster) {
    for (uint32_t i = 0; i < cluster.ids.size(); ++i) {
      ids.push_back(cluster.ids[i]);
    }
  }
  vector<uint32_t> ids;
};

/*struct KMER {
  KMER(const string& _name, const string& _seq, const Point& _point) :
    name(_name),
    seq(_seq),
    point(_point) {
  }
  string name;
  string seq;
  Point point;
};*/

void KmerToCoordinates(const string& kmer, Point& point) {
  size_t k = 0;
  for (size_t i = 0; i < kmer.size(); ++i) {
    int AA = base[kmer[i] - 'A'];
    for (size_t j = 0; j < AACoordinateSize; ++j) {
      point.data[k++] = coordinates[AA][j];
    }
  }
}

Point KmerToCoordinates(const string& kmer) {
  Point point;
  size_t k = 0;
  for (size_t i = 0; i < kmer.size(); ++i) {
    int AA = base[kmer[i] - 'A'];
    for (size_t j = 0; j < AACoordinateSize; ++j) {
      point.data[k++] = coordinates[AA][j];
    }
  }
  return point;
}

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

Point Center(const vector<KMER>& kmers, const vector<uint32_t>& ids) {
  Point center;
  for (uint32_t i = 0; i < ids.size(); ++i) {
    for (uint32_t j = 0; j < DIMENSION; ++j) {
      center.data[j] += kmers[ids[i]].point().data[j];
    }
  }
  for (uint32_t j = 0; j < DIMENSION; ++j) {
    center.data[j] /= DIMENSION;
  }
  return center;
}

double Distance(const Point& a, const Point& b) {
  double dis = 0, r = 0;
  for (int i = 0; i < a.data.size(); ++i) {
    r = (a.data[i] - b.data[i]);
    dis += r * r;
  }
  return sqrt(dis);
}

void PairwiseDistanceSampling(const vector<Point>& points, ofstream& fout) {
  int r1 = 0, r2 = 0;
  for (int i = 0; i < points.size(); i += r1) {
    for (int j = i + 1; j < points.size(); j += r2) {
      r1 = random() % 10;
      r2 = random() % 10;
      if (r1 == 0 || r2 == 0) {
        fout << Distance(points[i], points[j]) << endl;
        r1 += 5;
        r2 += 5;
      }
    }
  }
}

void PairwiseDistanceSampling(const vector<KMER>& kmers,
                              const vector<uint32_t>& ids,
                              const uint32_t& hash_W,
                              pair<uint32_t, uint32_t>& count,
                              pair<uint32_t, uint32_t>& count25,
                              ofstream& fout) {
  for (int i = 0; i < ids.size(); ++i) {
    for (int j = i + 1; j < ids.size(); ++j) {
      double dis = Distance(kmers[ids[i]].point(), kmers[ids[j]].point());
      if (dis <= hash_W)
        count.first++;
      else
        count.second++;

      if (dis <= 25)
        count25.first++;
      else
        count25.second++;
      fout << dis << endl;
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

int ModuloPrime(const int& hash_value) {
  return (hash_value % 251 + 251) % 251;
}


void BuildLSHTalbe(const vector<KMER>& kmers, const LSH& lsh,
                   vector<uint8_t>& merged, HashTable& lsh_table) {
  cout << "BuildLSHTalbes... " << endl;
  for (uint32_t i = 0; i < kmers.size(); ++i) {
    if (merged[i] != 0)
      continue;
    lsh_table[lsh.HashKey(kmers[i].point().data)].push_back(i);
  }
  /*
   ofstream fout("innerdisincluster.txt");
   pair<uint32_t, uint32_t> count;
   pair<uint32_t, uint32_t> count25;
   for(HashTable::iterator it = table.begin();it != table.end();++it) {
   PairwiseDistanceSampling(kmers, it->second, hash_W, count, count25, fout);
   }
   fout.close();
   cout << "seee" << endl;
   cout << count.first / ((double)count.first + count.second) << endl;
   cout << count.second / ((double)count.first + count.second) << endl;
   cout << "see25" << endl;
   cout << count25.first / ((double)count25.first + count25.second) << endl;
   cout << count25.second / ((double)count25.first + count25.second) << endl;
   cout << "complete sampling..." << endl;
   */
}

void ObtainCenters(const vector<KMER>& kmers, const uint32_t& hash_K,
                const uint32_t& hash_L, const uint32_t& hash_W,
                const double& D_T, const double& ratio_D_T,
                const uint32_t& center_S,
                 vector<Point>& queries,
                 const string& output_file) {
  cout << hash_K << endl;
  cout << hash_L << endl;
  cout << hash_W << endl;
  cout << D_T << endl;
  cout << "Clustering... " << endl;
  vector<Cluster> clusters(kmers.size());
  for (uint32_t i = 0; i < kmers.size(); ++i) {
    clusters[i].AddPoint(i);
  }
  double threshold = D_T * ratio_D_T;
  vector<uint8_t> merged(kmers.size(), 0);
  for (uint32_t l = 0; l < hash_L; ++l) {
    clock_t start = clock();
    HashTable lsh_table;
    LSH lsh(DIMENSION, hash_K, hash_W);
    BuildLSHTalbe(kmers, lsh, merged, lsh_table);
    for (uint32_t i = 0; i < kmers.size(); ++i) {
      if (merged[i] != 0) {
        continue;
      }
      string key = lsh.HashKey(kmers[i].point().data);
      HashTable::const_iterator it = lsh_table.find(key);
      const vector<uint32_t>& ids = it->second;
      for (uint32_t j = 0; j < ids.size(); ++j) {
        if (merged[ids[j]] != 0 || ids[j] == i) {
          continue;
        }
        if (PairwiseDistance(kmers[i].point(), kmers[ids[j]].point())
            <= threshold) {
          clusters[i].AddCluster(clusters[ids[j]]);
          merged[ids[j]] = 1;
        }
      }
    }
    printf("Clustering l=%d takes %lf seconds\n", l,
           (clock() - start) / (double) CLOCKS_PER_SEC);
  }

  for (uint32_t i = 0; i < kmers.size(); ++i) {
    if (merged[i] == 0 && clusters[i].ids.size() >= center_S) {
      queries.push_back(Center(kmers, clusters[i].ids));
    }
  }
//////////////////////////////////
  pair<uint32_t, uint32_t> count;
  pair<uint32_t, uint32_t> count25;
  string see = output_file;
  see += "see.dista.txt";
  ofstream fsee(see.c_str());
  int cnt = 0;
  for (uint32_t i = 0; i < clusters.size(); ++i) {
    if(merged[i] != 0 ||  clusters[i].ids.size()  < center_S) continue;
    cnt++;
    PairwiseDistanceSampling(kmers, clusters[i].ids, D_T, count, count25,
                             fsee);
  }
  cout << "cnt = " << cnt << endl;
  cout << "seee" << endl;
  cout << count.first / ((double)count.first + count.second) << endl;
  cout << count.second / ((double)count.first + count.second) << endl;
  cout << "see25" << endl;
  cout << count25.first / ((double)count25.first + count25.second) << endl;
  cout << count25.second / ((double)count25.first + count25.second) << endl;
  fsee.close();
}

void AddPoints2Centers(const vector<KMER>& kmers, const vector<Point>& queries,
                       const uint32_t& hash_K, const uint32_t& hash_L,
                       const uint32_t& hash_W, const double& D_T,
                       const string& output_file) {
  vector<uint8_t> merged(kmers.size(), 0);
  vector<Cluster> clusters(queries.size());
  for (uint32_t l = 0; l < hash_L; ++l) {
    clock_t start = clock();
    HashTable lsh_table;
    LSH lsh(DIMENSION, hash_K, hash_W);
    BuildLSHTalbe(kmers, lsh, merged, lsh_table);
    for (uint32_t i = 0; i < queries.size(); ++i) {
      string key = lsh.HashKey(queries[i].data);
      HashTable::const_iterator it = lsh_table.find(key);
      if (it == lsh_table.end())
        continue;
      const vector<uint32_t>& ids = it->second;
      for (uint32_t j = 0; j < ids.size(); ++j) {
        if (merged[ids[j]] != 0) {
          continue;
        }
        if (PairwiseDistance(queries[i], kmers[ids[j]].point()) <= 0.5 * D_T) {
          clusters[i].ids.push_back(ids[j]);
          merged[ids[j]] = 1;
        }
      }
    }
    printf("Clustering l=%d takes %lf seconds\n", l,
           (clock() - start) / (double) CLOCKS_PER_SEC);
  }
  ofstream fout(output_file.c_str());
  pair<uint32_t, uint32_t> count;
  pair<uint32_t, uint32_t> count25;
  string see = output_file;
  see += "see.dxxxxxxxxxxxxxista.txt";
  ofstream fsee(see.c_str());
  for (uint32_t i = 0; i < clusters.size(); ++i) {
    if (clusters[i].ids.size() == 0) {
      continue;
    }
    PairwiseDistanceSampling(kmers, clusters[i].ids, D_T, count, count25, fsee);
    fout << "#cluster" << i << endl;
    cout << "clusters[i].ids.size() = " << clusters[i].ids.size() << endl;
    for (uint32_t j = 0; j < clusters[i].ids.size(); ++j) {
      fout << kmers[clusters[i].ids[j]].name << endl;
    }
  }
  cout << "seee" << endl;
  cout << count.first / ((double)count.first + count.second) << endl;
  cout << count.second / ((double)count.first + count.second) << endl;
  cout << "see25" << endl;
  cout << count25.first / ((double)count25.first + count25.second) << endl;
  cout << count25.second / ((double)count25.first + count25.second) << endl;
  fsee.close();
  fout.close();
}

void ClusteringLikeCDHIT(const vector<KMER>& kmers,
                const uint32_t& hash_K,
                const uint32_t& hash_L,
                const double& D_T,
                const double& hash_W,
                const string& output_file) {
  cout << hash_K << endl;
  cout << hash_L << endl;
  cout << D_T << endl;
  cout << "Clustering LikeCDHIT... " << endl;
  vector<Cluster> clusters(kmers.size());
  for (uint32_t i = 0; i < kmers.size(); ++i) {
    clusters[i].AddPoint(i);
  }
  vector<uint8_t> merged(kmers.size(), 0);
  for (uint32_t l = 0; l < hash_L; ++l) {
    clock_t start = clock();
    HashTable lsh_table;
    LSH lsh(DIMENSION, hash_K, hash_W);
    BuildLSHTalbe(kmers, lsh, merged, lsh_table);
    printf("Build LSH l=%d takes %lf seconds\n", l,
           (clock() - start) / (double) CLOCKS_PER_SEC);
    for (uint32_t i = 0; i < kmers.size(); ++i) {
      if (merged[i] != 0) {
        continue;
      }
      string key = lsh.HashKey(kmers[i].point().data);
     // for (key_i = key - 1; key_i <= key + 1; ++key_i) {
        HashTable::const_iterator it = lsh_table.find(key);
        const vector<uint32_t>& ids = it->second;
        cout << i << " " << ids.size() << endl;
        for (uint32_t j = 0; j < ids.size(); ++j) {
          if (merged[ids[j]] != 0 || ids[j] == i) {
            continue;
          }
          if (PairwiseDistance(kmers[i].point(), kmers[ids[j]].point()) <= D_T) {
            clusters[i].AddCluster(clusters[ids[j]]);
            merged[ids[j]] = 1;
          }
        }
      }
   // }
    printf("Clustering l=%d takes %lf seconds\n", l,
           (clock() - start) / (double) CLOCKS_PER_SEC);
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
 //Evaluation(kmers, new_clusters, hash_K, hash_L, D_T, gamma);
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
    uint32_t len = 10;

    /* number of random lines for each LSH */
    uint32_t hash_K = 8;

    /* number of hash functions */
    uint32_t hash_L = 8;

    /* distance threshold */
    double D_T = 50;

    /* bucket width */
    double hash_W = 50;

    /* center support size */
    uint32_t center_S = 20;

    /* CDHIT like clustering */
    bool CDHIT_like = false;

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
    opt_parse.add_opt("threshold", 'T', "cluster center threshold", true,
        D_T);
    opt_parse.add_opt("window", 'W', "bucketwidth", true,
                      hash_W);
    opt_parse.add_opt("centerS", 'S', "number of points to support center", true,
        center_S);
    opt_parse.add_opt("CDHIT", 'C', "CDHIT like clustering", false,
        CDHIT_like);
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
        //KmerToCoordinates(line, p);
        kmers.push_back(KMER(name, line));
      }
    }
    fin.close();
    printf("The number of kmers is %u\n", kmers.size());
    vector<HashTable> lsh_tables;
    clock_t start = clock();
    //BuildLSHTalbes(kmers, hash_K, hash_L, hash_W, D_T, ratio_D_T lsh_tables);
    //printf ("Build LSH Tables takes %lf seconds\n", (clock() - start)/ (double)CLOCKS_PER_SEC);

    //start = clock();
    if(CDHIT_like != true) {
      vector<Point> queries;
      ObtainCenters(kmers, hash_K, hash_L, hash_W, D_T, 1, center_S, queries, output_file);
      AddPoints2Centers(kmers, queries, hash_K / 2, hash_L, hash_W, D_T, output_file);
    } else {
      ClusteringLikeCDHIT( kmers, hash_K,hash_L, D_T, hash_W,output_file);
    }
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
