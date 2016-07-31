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

double PairwiseDistance(const Point& a, const Point& b) {
  double dis = 0, r = 0;
  for (int i = 0; i < a.data.size(); ++i) {
    r = (a.data[i] - b.data[i]);
    dis += r * r;
  }
  return sqrt(dis);
}

void BuildLSHTalbe(const vector<KMER>& kmers, const LSH& lsh,
                   vector<uint32_t>& root, HashTable& lsh_table) {
  cout << "BuildLSHTalbes... " << endl;
  for (uint32_t i = 0; i < kmers.size(); ++i) {
    if (root[i] != i)
      continue;
    lsh_table[lsh.HashKey(kmers[i].point().data)].push_back(i);
  }
}

void Clustering(const vector<KMER>& kmers, const uint32_t& hash_K,
                const uint32_t& hash_L, const double& hash_W,
                const double& hash_R, const string& output_file) {
  unordered_map<uint64_t, uint32_t> count;
  vector<uint32_t> root(kmers.size(), 0);
  for (uint32_t i = 0; i < kmers.size(); ++i) {
    root[i] = i;
  }

  for (uint32_t l = 0; l < hash_L; ++l) {
    clock_t start = clock();
    HashTable lsh_table;
    LSH lsh(DIMENSION, hash_K, hash_W);
    BuildLSHTalbe(kmers, lsh, root, lsh_table);
    uint64_t cnt = 0;
    for(HashTable::iterator it = lsh_table.begin();it != lsh_table.end();++it) {
      cnt += it->second.size();
    }
    cout << "num of buckets " << lsh_table.size() << endl;
    cout << "cnt = " << cnt << endl;
    printf("Build LSH l=%d takes %lf seconds\n", l,
           (clock() - start) / (double) CLOCKS_PER_SEC);
    start = clock();
    uint64_t sum_count = 0;
    for (uint32_t i = 0; i < kmers.size(); ++i) {
      if (root[i] != i) {
        continue;
      }
      string key = lsh.HashKey(kmers[i].point().data);
      HashTable::const_iterator it = lsh_table.find(key);
      const vector<uint32_t>& ids = it->second;
      cout << i << " " << ids.size() << endl;
      sum_count += ids.size();
      for (uint32_t j = 0; j < ids.size(); ++j) {
        if (root[ids[j]] != ids[j] || ids[j] < i) {
          continue;
        }
        uint64_t val;
        if (i < ids[j]) {
          val = i;
          val <<= 32;
          val += ids[j];
        } else {
          val = ids[j];
          val <<= 32;
          val += i;
        }
        count[val]++;

        if (PairwiseDistance(kmers[i].point(), kmers[ids[j]].point())
            <= hash_R) {
          root[ids[j]] = i;
        }
      }
    }
    cout << "sum_count " << sum_count << endl;
    printf("Clustering l=%d takes %lf seconds\n", l,
           (clock() - start) / (double) CLOCKS_PER_SEC);
  }

  string count_out = output_file;
  count_out += "count.txt";
  ofstream fsee(count_out.c_str());
  for(unordered_map<uint64_t, uint32_t>::iterator it = count.begin();it != count.end();++it) {
    uint32_t val1 = it->first >> 32; //higher
    uint32_t val2 = (it->first << 32) >> 32;
    fsee << it->first << "\t" << val1 << "\t" << val2 << "\t" << it->second << endl;
    if(it->second > 2) cout << it->first << "\t" << val1 << "\t" << val2 << "\t" << it->second << endl;
  }
  fsee.close();

  ofstream fout(output_file.c_str());
  for (uint32_t i = 0; i < kmers.size(); ++i) {
    fout << kmers[i].name << "\t" << kmers[root[i]].name << endl;
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
