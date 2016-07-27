#include "util.hpp"
#include <stdint.h>
#include <vector>
#include <string>
#include <unordered_set>
#include "./../smithlab_cpp/smithlab_os.hpp"
#include "./../smithlab_cpp/OptionParser.hpp"

using namespace std;

#define MIN_SIZE_CLUSTER 100

void shuffleMotifs(vector<pair<string, vector<string> > >& clusters,
                   const uint32_t& num_of_motifs,
                   const uint32_t& num_of_seqs_in_one_motif,
                   const string& output_file) {
  cout << num_of_motifs << endl;
  cout << num_of_seqs_in_one_motif << endl;
  if (num_of_motifs != 0 && num_of_seqs_in_one_motif != 0) {
    while (clusters.size() > num_of_motifs) {
      clusters.pop_back();
    }
    for (int i = 0; i < clusters.size(); ++i) {
      while (clusters[i].second.size() > num_of_seqs_in_one_motif) {
        clusters[i].second.pop_back();
      }
    }
  }
  cout << "shuffleMotifs" << endl;
  string out = output_file;
  out += "shuffleMotifs.txt";
  cout << out << endl;
  int total_num_motif_seqs = 0;
  for (int i = 0; i < clusters.size(); i++) {
    total_num_motif_seqs += clusters[i].second.size();
  }
  cout << "number of cluster " << clusters.size() << endl;
  cout << "total num motifs seqs " << total_num_motif_seqs << endl;
  vector<pair<string, string> > shuffle_points(total_num_motif_seqs);
  unordered_set<int> exsits;
  srand (time(NULL));
  int r = 0;
  for (int i = 0; i < clusters.size(); i++) {
    for (int j = 0; j < clusters[i].second.size(); j++) {
      while (1) {
        r = random() % total_num_motif_seqs;
        if (exsits.find(r) == exsits.end()) {
          exsits.insert(r);
          break;
        }
      }
      char name[100];
      //sprintf(name, "%s_motif%d_seq%d", clusters[i].first.c_str(), i, j);
      sprintf(name, "motif%d_seq%d", i, j);
      shuffle_points[r] = make_pair(name, clusters[i].second[j]);
    }
  }

  ofstream fout(out.c_str());
  for (int i = 0; i < total_num_motif_seqs; ++i) {
    fout << ">" << shuffle_points[i].first << endl;
    fout << shuffle_points[i].second << endl;
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
      fprintf(stdout, "[WELCOME TO pairwiseDistanceSampling v%s]\n", pmf_version);
      fprintf(stdout, "[%s", argv[0]);
      for (int i = 1; i < argc; i++) {
        fprintf(stdout, " %s", argv[i]);
      }
      fprintf(stdout, "]\n");
    }

    uint32_t num_of_motifs = 0;
    uint32_t num_of_seqs_in_one_motif = 0;

    /* kmers file */
    string kmers_file;

    /* kmer length */
    uint32_t kmer_length;

    /* output file */
    string output_file;
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "pairwiseDistanceSampling",
        "");
    opt_parse.add_opt("kmers", 'k', "kmers file", true,
        kmers_file);
    opt_parse.add_opt("len", 'l', "kmer length", true,
        kmer_length);
    opt_parse.add_opt("nx", 'm', "num of motifs", false,
        num_of_motifs);
    opt_parse.add_opt("nx", 'n', "num of seqeunce in one motifs", false,
        num_of_seqs_in_one_motif);
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
    ifstream fin(kmers_file.c_str());
    string kmer, line;
    vector<string> cluster;
    vector<pair<string, vector<string> > > clusters;
    string cluster_name;
    uint32_t cline = 0, freq = 0;
    while(getline(fin, line)) {
      if(line.size() == 0) continue;
      if(line[0] == '#') {
        if(cluster.size() >= MIN_SIZE_CLUSTER) {
          clusters.push_back(make_pair(cluster_name, cluster));
        }
        cluster_name = line;
        cluster.clear();
      } else {
        cluster.push_back(line);
      }
    }
    fin.close();
    if(cluster.size() >= MIN_SIZE_CLUSTER) {
      clusters.push_back(make_pair(cluster_name, cluster));
    }
    cout << "Number of Clusters: " << clusters.size() << endl;

    shuffleMotifs(clusters, num_of_motifs,
        num_of_seqs_in_one_motif, output_file);

  } catch (const SMITHLABException &e) {
    fprintf(stderr, "%s\n", e.what().c_str());
    return EXIT_FAILURE;
  } catch (std::bad_alloc &ba) {
    fprintf(stderr, "ERROR: could not allocate memory\n");
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
