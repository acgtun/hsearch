#include "util.hpp"
#include "lsh.hpp"

#include "read_proteins.hpp"

#include "hash_search.hpp"

#include "smithlab_os.hpp"
#include "OptionParser.hpp"

void PreClustering(const ProteinDB& proteinDB, HASH_BUCKETS& hash_buckets) {
  clock_t start = clock();
  uint32_t feature_size = static_cast<uint32_t>(pow(8, HASHLEN));
  uint32_t bit_num = 16;
  double sigma = 0.2;

  KLSH klsh_low(feature_size, bit_num, sigma);
  vector<double> p(feature_size);
  vector<int> features(feature_size, 0);

  for (uint32_t i = 0; i < proteinDB.num_of_proteins; ++i) {
    const vector<char>& pro_seq = proteinDB.pro_seqs[i];
    if (pro_seq.size() < HASHLEN) {
      continue;
    }
    fill(features.begin(), features.end(), 0);
    for (uint32_t i = 0; i <= pro_seq.size() - HASHLEN; ++i) {
      features[Kmer2Integer(&(pro_seq[i]))]++;
    }
    for (uint32_t i = 0; i < feature_size; ++i) {
      p[i] = features[i];
    }

    hash_buckets[klsh_low.GetHashValue(p)].push_back(i);
  }
  fprintf(stderr, "[NUMBER OF PRE-GROUPS %lu]\n", hash_buckets.size());
//#define BUCKETSIZE
#ifdef BUCKETSIZE
  uint32_t num_of_proteins = 0;
  map<uint32_t, uint32_t> bucket_size;
  map<uint32_t, uint32_t> bucket_length;

  for (HASH_BUCKETS::iterator it = hash_buckets.begin();
      it != hash_buckets.end(); ++it) {
    num_of_proteins += it->second.size();
    bucket_size[it->second.size()]++;
    uint32_t length = 0;
    for (uint32_t i = 0; i < it->second.size(); ++i) {
      length += proteinDB.pro_seqs[it->second[i]].size();
    }
    bucket_length[length]++;
  }
  cout << "num_of_proteins = " << num_of_proteins << endl;
  ofstream fout("size.txt");
  for (map<uint32_t, uint32_t>::iterator it = bucket_size.begin();
      it != bucket_size.end(); ++it) {
    fout << it->first << " " << it->second << endl;
  }
  fout.close();
  fout.open("length.txt");
  for (map<uint32_t, uint32_t>::iterator it = bucket_length.begin();
      it != bucket_length.end(); ++it) {
    fout << it->first << " " << it->second << endl;
  }
  fout.close();
#endif
#ifdef DEBUG
  for (HASH_BUCKETS::iterator it = hash_buckets.begin();
      it != hash_buckets.end(); ++it) {
    cout << it->first << " " << it->second.size() << endl;
    for (uint32_t i = 0; i < it->second.size(); ++i) {
      proteinDB.OutputProtein(it->second[i]);
    }
    cout << "cluster -----------------------------------" << endl;
  }
#endif

  fprintf(stderr,
          "[Locality-Sensitive Hashing Pre-Clustering TAKES %lf SECONDS]\n",
          (clock() - start) / (double) CLOCKS_PER_SEC);
}

int main(int argc, const char *argv[]) {
  srand(time(NULL));
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
      fprintf(stderr, "[WELCOME TO PCLUSTER v%s]\n", pcluster_version);
      fprintf(stderr, "[%s", argv[0]);
      for (int i = 1; i < argc; i++) {
        fprintf(stderr, " %s", argv[i]);
      }
      fprintf(stderr, "]\n");
    }

    /* protein database file */
    string protein_file;

    /* output file */
    string output_file;

    /* number of threads for mapping */
    int num_of_threads = 1;
    int nMaxHitPer = 500;
    int nMaxAlnPer = 100;
    bool bHssp = false;
    int nMinLen = 0;
    //bool bLogE = false;
    //bool bEvalue = true;
    double dThr = 10.0;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "cluster protein sequences",
                           "");
    opt_parse.add_opt("database", 'd', "protein database file", true,
                      protein_file);
    opt_parse.add_opt("output", 'o', "output file name", true, output_file);
    opt_parse.add_opt("thread", 't', "number of threads for mapping", false,
                      num_of_threads);

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

    //////////////////////////////////////////////////////////////
    // PRECLUSTERING
    ProteinDB proteinDB(protein_file);
    fprintf(stderr, "[THE TOTAL NUMBER OF PROTEINS IN THE DATABASE IS %u]\n",
        proteinDB.num_of_proteins);
    HASH_BUCKETS hash_buckets;
    PreClustering(proteinDB, hash_buckets);
    //////////////////////
    // CLUSTERING
    uint32_t group_id = 0;
    CHashSearch hs(output_file, dThr, nMaxAlnPer, nMaxHitPer, bHssp, nMinLen);
    for (HASH_BUCKETS::iterator it = hash_buckets.begin();
        it != hash_buckets.end(); ++it) {
      clock_t start = clock();
      fprintf(stderr, "[CLUSTERING GROUP %u of %lu]\n", group_id++,
              hash_buckets.size());
      fprintf(stderr, "[THE NUMBER OF SEQUENCES IN THIS GROUP IS %lu]\n",
              it->second.size());
      if (it->second.size() > 1) {
        hs.BuildProteinsIndex(it->second, proteinDB);
        hs.ProteinSearching(it->second, proteinDB);
      }
      fprintf(stderr, "[CLUSTERING TAKES %lf SECONDS]\n",
              (clock() - start) / (double) CLOCKS_PER_SEC);
    }
  } catch (const SMITHLABException &e) {
    fprintf(stderr, "%s\n", e.what().c_str());
    return EXIT_FAILURE;
  } catch (std::bad_alloc &ba) {
    fprintf(stderr, "ERROR: could not allocate memory\n");
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
