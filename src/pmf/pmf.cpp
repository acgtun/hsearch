#include "util.hpp"
#include "lsh.hpp"

#include <string>
#include <vector>
#include <unordered_map>

#include "smithlab_os.hpp"
#include "OptionParser.hpp"

typedef unordered_map<uint32_t, vector<uint32_t> > HASH_BUCKETS;

void LoadProteins(vector<pair<string, string> >& proteins) {
  ifstream fin("/home/rcf-40/haifengc/panfs/github/IGC/data/IGC.pep");
  string name, peptide;
  while (getline(fin, name)) {
    getline(fin, peptide);
    proteins.push_back(make_pair(name, peptide));
  }
  fin.close();
}

void KmerToCoordinates(const vector<pair<string, string> >& proteins,
		const pair<uint32_t, uint32_t>& kmer, const uint32_t& kmer_length,
		vector<double>& point) {
  size_t k = 0;
  string kmer_string = proteins[kmer.first].second.substr(kmer.second, kmer_length);
  for (size_t i = 0; i < kmer_length; ++i) {
    int AA = base[kmer_string[i] - 'A'];
    for (size_t j = 0; j < AACoordinateSize; ++j) {
      point[k++] = coordinates[AA][j];
    }
  }
}

void InitClustering(const vector<pair<string, string> >& proteins,
                    const vector<pair<uint32_t, uint32_t> >& kmers,
                    const uint32_t kmer_length, const string& output_file) {
  uint32_t dimension = AACoordinateSize * kmer_length;
  vector<double> point(dimension);
  HASH_BUCKETS hash_buckets;
  LSH lsh(dimension);
  for (uint32_t i = 0; i < kmers.size(); ++i) {
    KmerToCoordinates(proteins, kmers[i], kmer_length, point);
    hash_buckets[lsh.HashBucketIndex(point)].push_back(i);
  }

  string out = output_file;
  out += ".initclustering";
  ofstream fcluster(out.c_str());
  for (HASH_BUCKETS::iterator it = hash_buckets.begin();
      it != hash_buckets.end(); ++it) {
    fcluster << "cluster" << "\t" << it->second.size() << endl;
    for (uint32_t i = 0; i < it->second.size(); ++i) {
      fcluster << it->second[i] << endl;
    }
  }
  fcluster.close();
}

bool sortClusterCMP(const vector<uint32_t>& a, const vector<uint32_t>& b) {
  return a.size() > b.size();
}

void Clustering(const vector<pair<string, string> >& proteins,
                const vector<pair<uint32_t, uint32_t> >& kmers,
                const vector<uint32_t>& frequency, const uint32_t kmer_length,
                const uint32_t& num_of_hash_functions,
                const string& output_file) {
  uint32_t dimension = AACoordinateSize * kmer_length;
  string input = output_file;
  input += ".initclustering";
  ifstream fin(input.c_str());

  uint32_t num_hash = 0;
  uint32_t hash_size;
  string cluster_name;
  vector<double> point(dimension);
  while (num_hash++ < num_of_hash_functions - 1) {
    cout << "num_hash " << num_hash << endl;
    char out[1000];
    sprintf(out, "%s_cluster%u.clustering", output_file.c_str(), num_hash);
    ofstream fout(out);
    LSH lsh(dimension);

    while (fin >> cluster_name >> hash_size) {
      vector<uint32_t> kmers_id(hash_size, 0);
      for (size_t i = 0; i < hash_size; ++i) {
        fin >> kmers_id[i];
      }

      // only one kmer
      if (hash_size <= 1) {
        fout << "cluster" << "\t" << hash_size << endl;
        for (size_t i = 0; i < hash_size; ++i) {
          fout << kmers_id[i] << endl;
        }
        continue;
      }

      HASH_BUCKETS hash_buckets;
      for (uint32_t i = 0; i < hash_size; ++i) {
        KmerToCoordinates(proteins, kmers[kmers_id[i]], kmer_length, point);
        hash_buckets[lsh.HashBucketIndex(point)].push_back(kmers_id[i]);
      }

      for (HASH_BUCKETS::iterator it = hash_buckets.begin();
          it != hash_buckets.end(); ++it) {
        fout << "cluster" << "\t" << it->second.size() << endl;
        for (uint32_t i = 0; i < it->second.size(); ++i) {
          fout << it->second[i] << endl;
        }
      }
    }
    fin.close();
    fout.close();

    fin.open(out);
  }

  // final output
  ofstream fcluster(output_file.c_str());
  vector<vector<uint32_t> > clusters;
  while (fin >> cluster_name >> hash_size) {
    vector < uint32_t > kmers_id(hash_size, 0);
    for (size_t i = 0; i < hash_size; ++i) {
      fin >> kmers_id[i];
    }
    clusters.push_back(kmers_id);
  }
  fin.close();

  sort(clusters.begin(), clusters.end(), sortClusterCMP);

  for(uint32_t c = 0; c < clusters.size();++c) {
    uint32_t count = 0;
    for (uint32_t i = 0; i < clusters[c].size(); ++i) {
      count += frequency[clusters[c][i]];
    }
    fcluster << "cluster" << "\t" << clusters[c].size() << "\t" << count << endl;
    for (uint32_t i = 0; i < clusters[c].size(); ++i) {
      uint32_t proID = kmers[clusters[c][i]].first;
      uint32_t proPos = kmers[clusters[c][i]].second;
      string kmer_string = proteins[proID].second.substr(proPos, kmer_length);
      fcluster << kmer_string << "\t" << frequency[clusters[c][i]] << endl;
    }
  }
  fcluster.close();
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
      fprintf(stderr, "[WELCOME TO PMF v%s]\n", pmf_version);
      fprintf(stderr, "[%s", argv[0]);
      for (int i = 1; i < argc; i++) {
        fprintf(stderr, " %s", argv[i]);
      }
      fprintf(stderr, "]\n");
    }

    /* kmers file */
    string kmers_file;

    /* kmer length */
    uint32_t len = 10;

    /* number of hash functions */
    uint32_t num_of_hash_functions = 8;

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

    vector<pair<string, string> > proteins;
    LoadProteins(proteins);
    //////////////////////////////////////////////////////////////
    vector<pair<uint32_t, uint32_t> > kmers;
    vector<uint32_t> frequency;

    ifstream fin(kmers_file.c_str());
    uint32_t proID, proPos;
    string kmer;
    uint32_t cline = 0, freq = 0;
    while(fin >> proID >> proPos >> freq) {
      kmers.push_back(make_pair(proID, proPos));
      frequency.push_back(freq);
    }
    fin.close();
    printf("The number of kmers is %u\n", kmers.size());
    //////////////////////
    // CLUSTERING
    InitClustering(proteins, kmers, len, output_file);
    Clustering(proteins, kmers, frequency, len, num_of_hash_functions, output_file);
  } catch (const SMITHLABException &e) {
    fprintf(stderr, "%s\n", e.what().c_str());
    return EXIT_FAILURE;
  } catch (std::bad_alloc &ba) {
    fprintf(stderr, "ERROR: could not allocate memory\n");
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
