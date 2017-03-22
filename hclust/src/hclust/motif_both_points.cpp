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

struct MOTIFRES {
  MOTIFRES(const string& _motif, const string& _protein, const double& _dis)
      : motif(_motif),
        protein(_protein),
        dis(_dis) {

  }
  string motif;
  string protein;
  double dis;
};

bool sortCMP(const MOTIFRES& a, const MOTIFRES& b) {
  if (a.motif == b.motif) {
    return a.protein < b.protein;
  }
  return a.motif < b.motif;
}

int CompMOTIT(const MOTIFRES& a, const MOTIFRES& b) {

  //is large
  if (a.motif == b.motif && a.protein == b.protein) {
    return 0;
  }

  if (a.motif == b.motif) {
    if (a.protein == b.protein)
      return 0;
    if (a.protein > b.protein)
      return 1;
    return -1;
  }

  if (a.motif > b.motif)
    return 1;
  if (a.motif < b.motif)
    return -1;
}

double weight(const double& dis, const double& hash_R) {
  if (dis > hash_R + 0.1) {

    cout << "err " << dis << endl;
    exit(0);
    double w = dis / hash_R - 1;
    if (w > 1)
      return 1;
    return w;
  }

  //return 1 - dis / hash_R;
  if(dis < 0.0000001) return 1;
  if(dis < 24) return 1;
  double w = 1 / (dis - 24);
  if(w > 1) return 1;
  if(w < 0) return 1;
  return w;
  //if(w > 1)
   // return 1;
}

//double weight2(const double& dis) {
//  if (dis > 49.38) {
//    double w = dis / (2 * 49.38);
//    if (w > 1)
//      return 1;
//    return dis / (2 * 49.38);
//  }
//
//  return 1 - dis / (2 * 49.38);
//}

double evaulate(const string& ground_truth, const string& output_file, const double& hash_R) {
  vector<MOTIFRES> brute_force;
  string line, motif, protein;
  double distance;
  ifstream fin2(ground_truth.c_str());
  while (fin2 >> motif >> protein >> distance) {
    brute_force.push_back(MOTIFRES(motif, protein, distance));
  }
  fin2.close();
  /////////////////////////////
  ifstream fin(output_file);
  vector<MOTIFRES> hclust;
  while (fin >> motif >> protein >> distance) {
    hclust.push_back(MOTIFRES(motif, protein, distance));
  }
  fin.close();
  sort(hclust.begin(), hclust.end(), sortCMP);

  uint32_t i = 0, j = 0;
  double tp = 0.0, fn = 0.0, cnt = 0;
  unordered_map<int, int> tp_map;
  unordered_map<int, int> fn_map;
  while (i < brute_force.size() && j < hclust.size()) {
    int cmp = CompMOTIT(brute_force[i], hclust[j]);
    if (cmp == 0) {
      tp += weight(brute_force[i].dis, hash_R);
      tp_map[int(brute_force[i].dis * 100 / 10)]++;
      i++;
      j++;
    } else if (cmp == 1) {
      cout << "xnomo "<< hclust[j].motif << " " << hclust[j].protein << " " << hclust[j].dis << " " << hash_R << endl;
      j++;
    } else {
      fn += weight(brute_force[i].dis, hash_R);
      fn_map[int(brute_force[i].dis * 100 / 10)]++;
      cout << brute_force[i].dis << " " <<  weight(brute_force[i].dis, hash_R) << endl;
      cnt++;
      i++;
    }
  }
  while (i < brute_force.size()) {
    fn += weight(brute_force[i].dis, hash_R);
    fn_map[int(brute_force[i].dis * 100 / 10)]++;
    cnt++;
    i++;
  }

 // cout << "size = " << cnt  << " " << brute_force.size() << " " << cnt / (double) brute_force.size() << endl;
  cout << "ACCU: " << tp << " " << fn << " " << tp / (tp + fn) << "\t"
      << "size = " << cnt  << " " << brute_force.size() << " " << cnt / (double) brute_force.size()
      <<output_file<< endl;
  string out = output_file;
  out += ".accuracy.txt";
  ofstream fout(out.c_str());
  for(int i = 0;i < 500;i++) {
    if(fn_map.find(i) != fn_map.end() && tp_map.find(i) != tp_map.end()) {
      fout << i << " " << tp_map[i] / (fn_map[i] + (double)tp_map[i]) << " " << tp_map[i]  << " " <<  fn_map[i] << endl;
    } else if(fn_map.find(i) != fn_map.end()) {
      fout << i << " " << 0 << " fn "<< fn_map[i]  << endl;
    } else if(tp_map.find(i) != tp_map.end()) {
      fout << i << " " << 1  << " tp " << tp_map[i] << endl;
    }
  }
  fout.close();
  return tp / (tp + fn) ;
}

double PairwiseDistance(const Point& a, const Point& b) {
  double dis = 0.0, r = 0.0;
  for (uint32_t i = 0; i < DIMENSION; ++i) {
    r = a.data[i] - b.data[i];
    dis += r * r;
  }
  return sqrt(dis);
}

double PairwiseDistance_square(const Point& a, const Point& b) {
  double dis = 0.0, r = 0.0;
  for (uint32_t i = 0; i < DIMENSION; ++i) {
    r = a.data[i] - b.data[i];
    dis += r * r;
  }
  return dis;
}


void BuildLSHTalbe(const vector<Point>& kmers, const LSH& lsh,
                   HashTable& lsh_table) {
  for (uint32_t i = 0; i < kmers.size(); ++i) {
    string key = lsh.HashKey(kmers[i].data);
    //cout << key << endl;
    lsh_table[key].push_back(i);
  }
}

void Search(const vector<Point>& kmers,
            const vector<Point>& centers,
            const vector<string>& kmer_names,
            const vector<string>& center_names,
            const uint32_t& hash_K,
            const uint32_t& hash_L,
            const double& hash_W,
            const double& hash_R,
            const string& output_file) {
  double hash_R_square = hash_R * hash_R;
  //clock_t start = clock();
  vector<HashTable> lsh_tables(hash_L);
  vector<LSH> lsh_funs;
  for(uint32_t l = 0;l < hash_L;++l) {
    LSH lsh(DIMENSION, hash_K, hash_W);
    lsh_funs.push_back(lsh);
  }
  for (uint32_t l = 0; l < hash_L; ++l) {
    for (uint32_t i = 0; i < kmers.size(); ++i) {
      string key = lsh_funs[l].HashKey(kmers[i].data);
      lsh_tables[l][key].push_back(i);
    }
    cout << "table size " << lsh_tables[l].size() << endl;
  }
  //printf("Build LSHTables  takes %lf seconds\n",
        // (clock() - start) / (double) CLOCKS_PER_SEC);
  //start = clock();
  vector<int> label(kmers.size(), 0);
  ofstream fout(output_file.c_str());
  for (uint32_t i = 0; i < centers.size(); ++i) {
    memset(&(label[0]), 0, sizeof(int) * kmers.size());
    for (uint32_t l = 0; l < hash_L; ++l) {
      string key = lsh_funs[l].HashKey(centers[i].data);
      HashTable::iterator it = lsh_tables[l].find(key);
      if (it == lsh_tables[l].end()) {
        continue;
      }
      for (uint32_t j = 0; j < it->second.size(); ++j) {
        if (label[it->second[j]] != 0)
          continue;
        //cout << it->second.size() << endl;
        double dis_square = PairwiseDistance_square(kmers[it->second[j]],
                                                    centers[i]);
        label[it->second[j]] = 1;
        if (dis_square <= hash_R_square) {
          fout << center_names[i] << " " << kmer_names[it->second[j]] << " "
              << sqrt(dis_square) << endl;
        }
      }
    }
  }

  //printf("Searching takes %lf seconds\n",
        // (clock() - start) / (double) CLOCKS_PER_SEC);
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

    /* number of random lines for each LSH */
    uint32_t hash_K = 8;

    /* number of hash tables */
    uint32_t hash_L = 16;  // 95%

    /* bucket width */
    double hash_W = 50;

    /* distance threshold */
    double hash_R = 200;

    /* output file */
    string output_file;

    string ground_truth;
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "cluster kmers to motifs",
        "");
    opt_parse.add_opt("db", 'd', "protein database file", true,
                      kmer_file);
    opt_parse.add_opt("center", 'c', "centers from Pfam database", true,
        center_file);
    opt_parse.add_opt("len", 'l', "kmer length", true,
        kmer_length);
    //opt_parse.add_opt("hash_K", 'K', "number of random lines", true,
      //  hash_K);
    //opt_parse.add_opt("hash_L", 'L', "number of hash tables", true,
      //  hash_L);
    opt_parse.add_opt("window", 'W', "bucket width", true,
        hash_W);
    opt_parse.add_opt("threshold", 'T', "kmer threshold", true,
        hash_R);
    opt_parse.add_opt("groundtruth", 'g', "groundtruth", true,
        ground_truth);
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
    clock_t start_s = clock();
    int n = kmers.size();
    double p1 = 0.9;
    double p2 = 0.2;
    double roh = log(1 / p1) / log(1 / p2);
    //hash_K = int(log(n) / log(1 / p2));
    //hash_L =int(pow( n, roh)) + 1;
    hash_K = 4;
    hash_L = 4;
    printf("p1 = %lf p2 = %lf roh = %lf hash_K = %u hash_L = %u\n", p1, p2, roh, hash_K, hash_L);
    Search(kmers, centers, kmer_names, center_names, hash_K, hash_L, hash_W, hash_R, output_file);
    clock_t end_s = clock();
    cout << "evaulate ..." << endl;
    printf ("ACCURACY: %lf %lf\n", evaulate(ground_truth, output_file, hash_R), (end_s - start_s)/ (double)CLOCKS_PER_SEC);
  } catch (const SMITHLABException &e) {
    fprintf(stderr, "%s\n", e.what().c_str());
    return EXIT_FAILURE;
  } catch (std::bad_alloc &ba) {
    fprintf(stderr, "ERROR: could not allocate memory\n");
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
