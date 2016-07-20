#include "util.hpp"
#include <stdint.h>
#include <vector>
#include <string>
#include <unordered_map>
#include "./../smithlab_cpp/smithlab_os.hpp"
#include "./../smithlab_cpp/OptionParser.hpp"

using namespace std;

#define MIN_SIZE_CLUSTER 100

void KmerToCoordinates(const string& kmer, vector<double>& point) {
  size_t k = 0;
  for (size_t i = 0; i < kmer.size(); ++i) {
    int AA = base[kmer[i] - 'A'];
    for (size_t j = 0; j < AACoordinateSize; ++j) {
      point[k++] = coordinates[AA][j];
    }
  }
}

double Distance(const vector<double>& a, const vector<double>& b) {
  double sum = 0;
  for (int i = 0; i < a.size(); ++i) {
    sum += (a[i] - b[i]) * (a[i] - b[i]);
  }
  return sqrt(sum);
}

void PairwiseDistanceSampling(const vector<vector<double> >& points,
                              ofstream& fout) {
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

void InnerClusterDistance(vector<vector<vector<double> > >& cluster_points,
                          const string& output_file) {
  string out = output_file;
  out += "inner.txt";
  cout << out << endl;
  ofstream fout(out.c_str());
  for (size_t i = 0; i < cluster_points.size(); ++i) {
    int r = random() % 10;
    if(r != 0) continue;
    PairwiseDistanceSampling(cluster_points[i], fout);
  }
  fout.close();
}

void InterClusterDistance(vector<vector<vector<double> > >& cluster_points,
                          const string& output_file) {
  cout << "InterClusterDistance" << endl;
  string out = output_file;
  out += "inter.txt";
  cout << out << endl;
  ofstream fout(out.c_str());
  int r1 = 0, r2 = 0, r3 = 0, r4 = 0;
  for (int i = 0; i < cluster_points.size(); i += r1) {
    for (int j = i + 1; j < cluster_points.size(); j += r2) {
      for (int p = 0; p < cluster_points[i].size(); p++) {
        for (int q = 0; q < cluster_points[j].size(); q++) {
          r1 = random() % 100;
          r2 = random() % 100;
          r3 = random() % 100;
          r4 = random() % 100;
          if (r1 == 0 || r2 == 0 || r3 == 0 || r4 == 0) {
            fout << Distance(cluster_points[i][p], cluster_points[j][q]) << endl;
            r1 += 50;
            r2 += 50;
            r3 += 50;
            r4 += 50;
          }
        }
      }
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
      fprintf(stdout, "[WELCOME TO pairwiseDistanceSampling v%s]\n", pmf_version);
      fprintf(stdout, "[%s", argv[0]);
      for (int i = 1; i < argc; i++) {
        fprintf(stdout, " %s", argv[i]);
      }
      fprintf(stdout, "]\n");
    }

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
    vector<pair<string, string> > kmers;
    string line, name;
    while(fin >> line) {
      if(line[0] == '>') {
        name = line;
        fin >> line;
        kmers.push_back(make_pair(name, line));
      }
    }
    fin.close();
    vector<string> ground_truth(kmers.size());
    unordered_map<string, vector<uint32_t> > cluster_ids;
    for (uint32_t i = 0; i < kmers.size(); ++i) {
      size_t pos = kmers[i].first.find("_motif");
      size_t pos_ = kmers[i].first.find_last_of('_');
      string motif_num = kmers[i].first.substr(pos + 6, pos_ - pos - 6);
      ground_truth[i] = motif_num;
      cluster_ids[motif_num].push_back(i);
    }

    //string kmer, line;
    uint32_t dimension = AACoordinateSize * kmer_length;
    vector<double> p(dimension, 0);
    vector<vector<double> > points;
    vector<vector<vector<double> > > cluster_points;
    for(unordered_map<string, vector<uint32_t> >::iterator it = cluster_ids.begin();
        it != cluster_ids.end();++it) {
      points.clear();
      for(uint32_t i = 0;i < it->second.size();++i) {
        KmerToCoordinates(kmers[it->second[i]].second, p);
        points.push_back(p);
      }
      cluster_points.push_back(points);
    }

    cout << "Number of Clusters: " << cluster_points.size() << endl;

    InnerClusterDistance(cluster_points, output_file);
    InterClusterDistance(cluster_points, output_file);
  } catch (const SMITHLABException &e) {
    fprintf(stderr, "%s\n", e.what().c_str());
    return EXIT_FAILURE;
  } catch (std::bad_alloc &ba) {
    fprintf(stderr, "ERROR: could not allocate memory\n");
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
