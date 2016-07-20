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

double PairwiseDistance(const Point& a, const Point& b) {
  double sum = 0;
  for (int i = 0; i < a.data.size(); ++i) {
    sum += (a.data[i] - b.data[i]) * (a.data[i] - b.data[i]);
  }
  return sqrt(sum);
}

struct Cluster {
  Cluster() {
    dimension_sum.resize(DIMENSION);
    for (uint32_t i = 0; i < DIMENSION; ++i) {
      dimension_sum[i] = 0.0;
      center.data[i] = 0.0;
    }
  }
  vector<Point> points;
  vector<uint32_t> ids;
  vector<double> dimension_sum;
  Point center;

  void AddPointUpdateCenter(const Point& p, const uint32_t& id) {
    for (uint32_t i = 0; i < DIMENSION; ++i) {
      dimension_sum[i] += p.data[i];
    }
    points.push_back(p);
    ids.push_back(id);
    for (uint32_t i = 0; i < DIMENSION; ++i) {
      center.data[i] = dimension_sum[i] / points.size();
    }
  }

  void AddClusterUpdateCenter(const Cluster& cluster) {
    for (uint32_t i = 0; i < DIMENSION; ++i) {
      for(uint32_t j = 0;j < cluster.points.size();++j) {
        dimension_sum[i] += cluster.points[j].data[i];
      }
    }
    for(uint32_t j = 0;j < cluster.points.size();++j) {
      points.push_back(cluster.points[j]);
      ids.push_back(cluster.ids[j]);
    }

    for (uint32_t i = 0; i < DIMENSION; ++i) {
      center.data[i] = dimension_sum[i] / points.size();
    }
  }
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

void Evaluation(
    const vector<pair<string, string> >& kmers,
    const vector<Cluster>& clusters) {
  cout << "number of clusters by our methods " << clusters.size() << endl;

  vector<string> ground_truth(kmers.size());
  for (uint32_t i = 0; i < kmers.size(); ++i) {
    size_t pos = kmers[i].first.find("_motif");
    size_t pos_ = kmers[i].first.find_last_of('_');
    string motif_num = kmers[i].first.substr(pos + 6, pos_ - pos - 6);
    ground_truth[i] = motif_num;
  }

  vector<uint32_t> cluster_id(kmers.size(), 0);
  for (uint32_t i = 0; i < clusters.size(); ++i) {
    for (uint32_t j = 0; j < clusters[i].ids.size(); ++j) {
      cluster_id[clusters[i].ids[j]] = i;
    }
  }

  uint32_t tp = 0, tn = 0, fp = 0, fn = 0;
  for (uint32_t i = 0; i < kmers.size(); ++i) {
    if (i % 1000 == 0 && i != 0) {
      cout << i / (double) kmers.size() << " " << tp << " " << tn << " " << fp
          << " " << fn << endl;
      cout << "recall: " << (double) tp / ((double) tp + fn) << endl;
      cout << "precision: " << (double) tp / ((double) tp + fp) << endl;
    }
    if(i > 50000) break;
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

  cout << "recall: " << (double) tp / ((double) tp + fn) << endl;
  cout << "precision: " << (double) tp / ((double) tp + fp) << endl;
}

bool sortCMP(const pair<string, uint32_t>& a, const pair<string, uint32_t>& b) {
  return a.first < b.first;
}

void LSHClustering(const vector<Cluster>& clusters,
                   const uint32_t& hash_L,
                   vector<Cluster>& new_clusters) {
  vector<vector<uint32_t> > hash_values(clusters.size(),
                          vector<uint32_t>(hash_L, 0));
  int hash_index = 0;
  while (hash_index < hash_L) {
    cout << "hash index " << hash_index << endl;
    LSH lsh(DIMENSION);
    for (uint32_t i = 0; i < clusters.size(); ++i) {
      hash_values[i][hash_index] = lsh.HashBucketIndex(clusters[i].center.data);
    }
    hash_index++;
  }

  unordered_map<string, vector<uint32_t> > buckets;
  char hash_value_chr[100];
  for (uint32_t i = 0; i < clusters.size(); ++i) {
    string hash_value;
    for (uint32_t j = 0; j < hash_L; ++j) {
      sprintf(hash_value_chr, "%u", hash_values[i][j]);
      hash_value += hash_value_chr;
    }
    buckets[hash_value].push_back(i);
  }

  for (unordered_map<string, vector<uint32_t> >::iterator it = buckets.begin();
      it != buckets.end(); ++it) {
    Cluster one_cluster;
    for (uint32_t i = 0; i < it->second.size(); ++i) {
        one_cluster.AddClusterUpdateCenter(clusters[it->second[i]]);
    }
    new_clusters.push_back(one_cluster);
  }
}

void Clustering(const vector<pair<string, string> >& kmers,
                const uint32_t& hash_L,
                const double& cluster_distance_threshold,
                const string& output_file) {
  vector<Cluster> clusters(kmers.size());
  Point p;
  for (uint32_t i = 0; i < kmers.size(); ++i) {
    KmerToCoordinates(kmers[i].second, p);
    clusters[i].AddPointUpdateCenter(p, i);
  }

  vector<Cluster> new_clusters;
  LSHClustering(clusters, hash_L, new_clusters);
  clusters = new_clusters;
  new_clusters.clear();
  new_clusters.push_back(clusters[0]);

  for (uint32_t i = 1; i < clusters.size(); ++i) {
    double min_dis = std::numeric_limits<double>::max();
    uint32_t min_id = -1;
    for (uint32_t j = 0; j < new_clusters.size(); ++j) {
      double dis = PairwiseDistance(clusters[i].center, new_clusters[j].center);
      if (dis < min_dis) {
        min_dis = dis;
        min_id = j;
      }
    }
    if (min_dis <= cluster_distance_threshold) {
      new_clusters[min_id].AddClusterUpdateCenter(clusters[i]);
    } else {
      new_clusters.push_back(clusters[i]);
    }
  }

  Evaluation(kmers, new_clusters);

  //sort(kemrs_hash_value_to_string.begin(), kemrs_hash_value_to_string.end(), sortCMP);
  /*
   ofstream fout(output_file.c_str());
   for (uint32_t i = 0; i < kmers.size(); ++i) {
   uint32_t id = kemrs_hash_value_to_string[i].second;
   fout << kmers[id].second << "\t" << kmers[id].first << "\t"
   << kemrs_hash_value_to_string[i].first << endl;
   }
   fout.close();
   */
}

void Clustering2(const vector<pair<string, string> >& kmers,
                const uint32_t& hash_L,
                const double& cluster_distance_threshold,
                const string& output_file) {
  vector<Cluster> clusters(kmers.size());
  Point p;
  for (uint32_t i = 0; i < kmers.size(); ++i) {
    KmerToCoordinates(kmers[i].second, p);
    clusters[i].AddPointUpdateCenter(p, i);
  }
  vector<Cluster> new_clusters;
  LSHClustering(clusters, hash_L, new_clusters);
  clusters = new_clusters;
  new_clusters.clear();
  new_clusters.push_back(clusters[0]);

  for (uint32_t i = 1; i < clusters.size(); ++i) {
    bool found = false;
    for (uint32_t j = 0; j < new_clusters.size(); ++j) {
      double dis = PairwiseDistance(clusters[i].center, new_clusters[j].center);
      if (dis < cluster_distance_threshold) {
        new_clusters[j].AddClusterUpdateCenter(clusters[i]);
        found = true;
        break;
      }
    }
    if(found = false) {
      new_clusters.push_back(clusters[i]);
    }
  }

  Evaluation(kmers, new_clusters);

  //sort(kemrs_hash_value_to_string.begin(), kemrs_hash_value_to_string.end(), sortCMP);
  /*
   ofstream fout(output_file.c_str());
   for (uint32_t i = 0; i < kmers.size(); ++i) {
   uint32_t id = kemrs_hash_value_to_string[i].second;
   fout << kmers[id].second << "\t" << kmers[id].first << "\t"
   << kemrs_hash_value_to_string[i].first << endl;
   }
   fout.close();
   */
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

    /* number of hash functions */
    uint32_t num_of_hash_functions = 8;

    /* distance threshold */
    double cluster_distance_threshold = 50;

    /* output file */
    string output_file;
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "cluster kmers to motifs",
        "");
    opt_parse.add_opt("kmers", 'k', "kmers file", true,
        kmers_file);
    opt_parse.add_opt("len", 'l', "kmer length", true,
        len);
    opt_parse.add_opt("nhash", 'h', "number of hash functions", false,
        num_of_hash_functions);
    opt_parse.add_opt("threshold", 't', "cluster center threshold", true,
                      cluster_distance_threshold);
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
    vector<pair<string, string> > kmers;
    ifstream fin(kmers_file.c_str());
    string line, name;
    while(fin >> line) {
      if(line[0] == '>') {
        name = line;
        fin >> line;
        kmers.push_back(make_pair(name, line));
      }
    }
    fin.close();
    printf("The number of kmers is %u\n", kmers.size());
    Clustering2(kmers, num_of_hash_functions, cluster_distance_threshold, output_file);
  } catch (const SMITHLABException &e) {
    fprintf(stderr, "%s\n", e.what().c_str());
    return EXIT_FAILURE;
  } catch (std::bad_alloc &ba) {
    fprintf(stderr, "ERROR: could not allocate memory\n");
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
