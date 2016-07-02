#include "util.hpp"
#include "lsh.hpp"

#include <string.h>

#include <map>
#include <string>
#include <vector>
#include <unordered_map>

#include "smithlab_os.hpp"
#include "OptionParser.hpp"

typedef unordered_map<uint32_t, vector<uint32_t> > HASH_BUCKETS;

void Clustering(const vector<string>& kmers) {
  clock_t start = clock();
  uint32_t dimension = static_cast<uint32_t>(pow(AACoordinateSize,
                                                 kmers[0].size()));
  uint32_t bit_num = 16;
  double sigma = 0.2;

  KLSH klsh_low(dimension, bit_num, sigma);
  vector<double> point(dimension);

  HASH_BUCKETS hash_buckets;
  for (uint32_t i = 0; i < kmers.size(); ++i) {
    KmerToCoordinates(kmers[i], point);
    hash_buckets[klsh_low.GetHashValue(point)].push_back(i);
  }

  fprintf(stderr, "[NUMBER OF PRE-GROUPS %lu]\n", hash_buckets.size());

  map<uint32_t, uint32_t> bucket_size;
  for (HASH_BUCKETS::iterator it = hash_buckets.begin();
      it != hash_buckets.end(); ++it) {
    bucket_size[it->second.size()]++;
  }

  ofstream fout("size.txt");
  for (map<uint32_t, uint32_t>::iterator it = bucket_size.begin();
      it != bucket_size.end(); ++it) {
    fout << it->first << " " << it->second << endl;
  }
  fout.close();

  for (HASH_BUCKETS::iterator it = hash_buckets.begin();
      it != hash_buckets.end(); ++it) {
    cout << "cluster " << it->first << " " << it->second.size() << endl;
    for (uint32_t i = 0; i < it->second.size(); ++i) {
      cout << kmers[it->second[i]] << endl;
    }
  }

  fprintf(stderr,
          "[Locality-Sensitive Hashing Pre-Clustering TAKES %lf SECONDS]\n",
          (clock() - start) / (double) CLOCKS_PER_SEC);
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
      fprintf(stderr, "[WELCOME TO PMF v%s]\n", pmf_version);
      fprintf(stderr, "[%s", argv[0]);
      for (int i = 1; i < argc; i++) {
        fprintf(stderr, " %s", argv[i]);
      }
      fprintf(stderr, "]\n");
    }

    /* kmers file */
    string kmers_file;

    /* output file */
    string output_file;
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "cluster kmers to motifs",
        "");
    opt_parse.add_opt("kmers", 'k', "kmers file", true,
        kmers_file);
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

    //////////////////////////////////////////////////////////////
    vector<string> kmers;
    ifstream fin(kmers_file.c_str());
    string kmer;
    uint32_t cline = 0;
    while(fin >> kmer) {
      kmers.push_back(kmer);
      cline++;
      if(cline % 100000 == 0) {
        cout << cline << endl;
      }
    }
    fin.close();
    printf("The number of kmers is %lu\n", kmers.size());
    //////////////////////
    // CLUSTERING
    Clustering(kmers);
  } catch (const SMITHLABException &e) {
    fprintf(stderr, "%s\n", e.what().c_str());
    return EXIT_FAILURE;
  } catch (std::bad_alloc &ba) {
    fprintf(stderr, "ERROR: could not allocate memory\n");
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

