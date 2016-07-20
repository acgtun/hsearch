#include "cluster.hpp"
#include "smithlab_os.hpp"
#include "OptionParser.hpp"

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

    /* reclustering */
    string org_output;
    uint32_t start_number = 0;

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
    opt_parse.add_opt("org_output", 'r', "org_output", true,
                      org_output);
    opt_parse.add_opt("startnum", 's', "start_number", true,
        start_number);
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
    /////////////////////////////
    // CLUSTERING
    Clustering(proteins, kmers, frequency, len, num_of_hash_functions, org_output, start_number, output_file);
  } catch (const SMITHLABException &e) {
    fprintf(stderr, "%s\n", e.what().c_str());
    return EXIT_FAILURE;
  } catch (std::bad_alloc &ba) {
    fprintf(stderr, "ERROR: could not allocate memory\n");
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
