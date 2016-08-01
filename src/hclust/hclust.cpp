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
  void Output() const {
    cout << "point";
    for(uint32_t i = 0;i < DIMENSION;++i) {
      cout << " " << data[i];
    }
    cout << endl;
  }
  vector<double> data;
};

Point KmerToCoordinates(const string& kmer);

struct KMER {
  KMER(const string& _name, const string& _seq, const Point& _point)
      : name(_name),
        seq(_seq),
        point(_point) {
  }
  string name;
  string seq;
  //Point point() const {
    //return KmerToCoordinates(seq);
 // }
  Point point;
};

Point Center(const vector<KMER>& kmers, const vector<uint32_t>& ids);

struct Cluster {
  Cluster() {
    radius = 0.0;
  }

  void AddPoint(const uint32_t& id) {
    ids.push_back(id);
  }

  void AddCluster(const Cluster& cluster) {
    for (uint32_t i = 0; i < cluster.ids.size(); ++i) {
      ids.push_back(cluster.ids[i]);
    }
  }

  Point center(const vector<KMER>& kmers) const {
    return Center(kmers, ids);
  }

  vector<uint32_t> ids;
  double radius;
};


Point KmerToCoordinates(const string& kmer) {
  Point point;
  size_t k = 0;
  for (size_t i = 0; i < kmer.size(); ++i) {
    int AA = base[kmer[i] - 'A'];
    if(AA == -1) {
      AA = rand() % 20;
    }
    for (size_t j = 0; j < AACoordinateSize; ++j) {
      point.data[k++] = coordinates[AA][j];
    }
  }

  return point;
}


Point Center(const vector<KMER>& kmers, const vector<uint32_t>& ids) {
  Point center;
  for (uint32_t i = 0; i < ids.size(); ++i) {
    for (uint32_t j = 0; j < DIMENSION; ++j) {
      center.data[j] += kmers[ids[i]].point.data[j];
    }
  }
  for (uint32_t j = 0; j < DIMENSION; ++j) {
    center.data[j] /= ids.size();
  }
  return center;
}

Point Center(const vector<Point>& points) {
  Point center;
  for (uint32_t i = 0; i < points.size(); ++i) {
    for (uint32_t j = 0; j < DIMENSION; ++j) {
      center.data[j] += points[i].data[j];
    }
  }
  for (uint32_t j = 0; j < DIMENSION; ++j) {
    center.data[j] /= points.size();
  }
  return center;
}

Point Center(const vector<KMER>& kmers, const vector<Cluster>& clusters,
             const vector<uint32_t>& ids) {
//  Point center;
//  for (uint32_t i = 0; i < ids.size(); ++i) {
//    Point ci = clusters[ids[i]].center(kmers);
//    for (uint32_t j = 0; j < DIMENSION; ++j) {
//      center.data[j] += ci.data[j];
//    }
//  }
//
//  for (uint32_t j = 0; j < DIMENSION; ++j) {
//    center.data[j] /= ids.size();
//  }
//  return center;

  uint32_t size = 0;
  for (uint32_t i = 0; i < ids.size(); ++i) {
    size += clusters[ids[i]].ids.size();
  }

  Point center;
  //cout << "----------center1" << endl;
  //center.Output();
  //cout << "----------center2" << endl;
  for (uint32_t i = 0; i < ids.size(); ++i) {
    for (uint32_t j = 0; j < clusters[ids[i]].ids.size(); ++j) {
      for (uint32_t k = 0; k < DIMENSION; ++k) {
       // cout <<  "(" << kmers[clusters[ids[i]].ids[j]].point.data[k] << "," << clusters[ids[i]].ids[j] << ") ";
        center.data[k] += kmers[clusters[ids[i]].ids[j]].point.data[k];
        //center.Output();
      }
      //cout << i << " " << j << " " << k << " " << kmers[clusters[ids[i]].ids[j]].point.data[k] << endl;
      //cout << i << " " << j << endl;
    }
    //center.Output();
  }
 // center.Output();
  for (uint32_t k = 0; k < DIMENSION; ++k) {
    center.data[k] /= size;
  }
  return center;
}

double PairwiseDistance(const Point& a, const Point& b) {
  double dis = 0, r = 0;
  for (int i = 0; i < a.data.size(); ++i) {
    r = (a.data[i] - b.data[i]);
    dis += r * r;
  }
  return sqrt(dis);
}

void BuildLSHTalbe(const vector<KMER>& kmers, const vector<Cluster>& clusters,
                   const LSH& lsh, HashTable& lsh_table) {
  cout << "BuildLSHTalbes... " << endl;
  for (uint32_t i = 0; i < clusters.size(); ++i) {
    string key = lsh.HashKey(clusters[i].center(kmers).data);
    lsh_table[key].push_back(i);
  }
}


void PairwiseDistanceSampling(const vector<KMER>& kmers, const Cluster& cluster,
                              ofstream& fout) {
  for (int i = 0; i < cluster.ids.size(); ++i) {
    for (int j = i + 1; j < cluster.ids.size(); ++j) {
      fout << PairwiseDistance(kmers[cluster.ids[i]].point, kmers[cluster.ids[j]].point)
          << endl;
    }
  }
}

void ClustingBucket(const vector<KMER>& kmers, const vector<Cluster>& clusters,
                    const vector<uint32_t>& ids, const double& hash_R,
                    vector<Cluster>& new_clusters) {
  //if(ids.size() == 1) return;
  Point c = Center(kmers, clusters, ids);
  //c.Output();
  //if(ids.size() == 1) {
   // cout << "------------" << endl;
   // kmers[clusters[ids[0]].ids[0]].point.Output();
  //}
  vector<uint8_t> label(ids.size(), 0);
  for (uint32_t i = 0; i < ids.size(); ++i) {
    double dis = PairwiseDistance(c, clusters[ids[i]].center(kmers));
    //cout << "ra " << dis  << "\t" << clusters[ids[i]].radius << "\t" << dis + clusters[ids[i]].radius << "\t" << hash_R * 0.75 << endl;
    //if (dis + clusters[ids[i]].radius > hash_R / 2) {
    if (dis + clusters[ids[i]].radius > hash_R / 2) {
      new_clusters.push_back(clusters[ids[i]]);
      label[i] = 1;
    }
  }

  Cluster cluster;
  for (uint32_t i = 0; i < ids.size(); ++i) {
    if (label[i] == 0) {
      cluster.AddCluster(clusters[ids[i]]);
    }
  }

  if (cluster.ids.size() != 0) {
    cout << "size = " << cluster.ids.size() << endl;
    Point center = cluster.center(kmers);
    double max_dis = 0.0;
    for (uint32_t i = 0; i < cluster.ids.size(); ++i) {
      double dis = PairwiseDistance(center, kmers[cluster.ids[i]].point);
      max_dis = max_dis > dis ? max_dis : dis;
    }
    cluster.radius = max_dis;
    ///////////////////////////////////
    for (int i = 0; i < cluster.ids.size(); ++i) {
      for (int j = i + 1; j < cluster.ids.size(); ++j) {
       /* cout << "dis\t" << hash_R / 2 << "\t("
             << PairwiseDistance(c, kmers[cluster.ids[i]].point) << ", "
             << PairwiseDistance(c, kmers[cluster.ids[j]].point) << ")\t"
             << PairwiseDistance(kmers[cluster.ids[i]].point, kmers[cluster.ids[j]].point) << endl;*/
      }
    }
    //////////////////////////////////
    new_clusters.push_back(cluster);
  }
}



void InnerClusterDistance(const vector<KMER>& kmers,
                          const vector<Cluster>& clusters,
                          const string& output_file) {
  string out = output_file;
  out += "inner.txt";
  cout << out << endl;
  ofstream fout(out.c_str());
  for (size_t i = 0; i < clusters.size(); ++i) {
    PairwiseDistanceSampling(kmers, clusters[i], fout);
  }
  fout.close();
}


void Clustering(const vector<KMER>& kmers, const uint32_t& hash_K,
                const uint32_t& hash_L, const double& hash_W,
                const double& hash_R, const string& output_file) {
  vector<Cluster> clusters(kmers.size());
  for (uint32_t i = 0; i < kmers.size(); ++i) {
    clusters[i].ids.push_back(i);
    clusters[i].radius = 0.0;
  }
  for (uint32_t l = 0; l < hash_L; ++l) {
    cout << "llll = " << l << endl;
    clock_t start = clock();
    HashTable lsh_table;
    LSH lsh(DIMENSION, hash_K, hash_W);
    BuildLSHTalbe(kmers, clusters, lsh, lsh_table);
    printf("Build LSH l=%d takes %lf seconds\n", l,
           (clock() - start) / (double) CLOCKS_PER_SEC);
    start = clock();
    vector<Cluster> new_clusters;
    uint32_t cnt = 0;
    uint32_t sum_size = 0;
    for (HashTable::const_iterator it = lsh_table.begin();
        it != lsh_table.end(); ++it) {
      cout << cnt++ << "\t" << lsh_table.size() << "\t" << it->second.size()
          << endl;
      sum_size += it->second.size();
      ClustingBucket(kmers, clusters, it->second, hash_R, new_clusters);
    }

    cout << "sum_size " << sum_size << endl;
    cout << "size = " << new_clusters.size() << endl;
    //InnerClusterDistance(kmers, new_clusters, output_file);
    clusters.clear();
    clusters = new_clusters;
    cout << "size = " << clusters.size() << endl;
    printf("Clustering l=%d takes %lf seconds\n", l,
           (clock() - start) / (double) CLOCKS_PER_SEC);
  }

//  string count_out = output_file;
//  count_out += "count.txt";
//  ofstream fsee(count_out.c_str());
//  for(unordered_map<uint64_t, uint32_t>::iterator it = count.begin();it != count.end();++it) {
//    uint32_t val1 = it->first >> 32; //higher
//    uint32_t val2 = (it->first << 32) >> 32;
//    fsee << it->first << "\t" << val1 << "\t" << val2 << "\t" << it->second << endl;
//    if(it->second > 2) cout << it->first << "\t" << val1 << "\t" << val2 << "\t" << it->second << endl;
//  }
//  fsee.close();

  ofstream fout(output_file.c_str());
  for (uint32_t i = 0; i < clusters.size(); ++i) {
    fout << "#cluster" << i << endl;
    for (uint32_t j = 0; j < clusters[i].ids.size(); ++j) {
      fout << kmers[clusters[i].ids[j]].name << endl;
    }
  }
  fout.close();
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
        name = line;
        fin >> line;
        kmers.push_back(KMER(name, line, KmerToCoordinates(line)));
        //cout << line << "\t";
        //kmers.back().point.Output();
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
