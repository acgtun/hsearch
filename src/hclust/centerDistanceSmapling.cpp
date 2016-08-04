#include "util.hpp"
#include <stdint.h>
#include <vector>
#include <string>
#include <unordered_set>
#include "./../smithlab_cpp/smithlab_os.hpp"
#include "./../smithlab_cpp/OptionParser.hpp"

using namespace std;

#define MIN_SIZE_CLUSTER 100
uint32_t DIMENSION = 200;

struct Point {
  Point()
      : data(DIMENSION, 0.0) {
  }
  void Output() {
	  cout << "point";
	  for(uint32_t i = 0;i < DIMENSION;++i) {
		  cout << " " << data[i];
	  }
	  cout << endl;
  }
  vector<double> data;
};


struct KMER {
  KMER(const string& _name, const string& _seq, const Point& _point)
      : name(_name),
        seq(_seq),
        point(_point){
  }
  string name;
  string seq;
  Point point;
};

Point KmerToCoordinates(const string& kmer) {
  Point point;
  size_t k = 0;
  for (size_t i = 0; i < kmer.size(); ++i) {
	  //cout << "i = " << i << endl;
    int AA = base[kmer[i] - 'A'];
    if (AA == -1) {
      AA = rand() % 20;
    }
    for (size_t j = 0; j < AACoordinateSize; ++j) {
      point.data[k++] = coordinates[AA][j];
    }
  }
  //point.Output();
  return point;
}

double PairwiseDistance(const Point& a, const Point& b) {
  double dis = 0.0, r = 0.0;
  for (int i = 0; i < DIMENSION; ++i) {
    r = a.data[i] - b.data[i];
    dis += r * r;
  }
  return sqrt(dis);
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

void clusterDistance(vector<pair<string, vector<string> > >& clusters,
                   const string& output_file) {
  int total_num_motif_seqs = 0;
  for (int i = 0; i < clusters.size(); i++) {
    total_num_motif_seqs += clusters[i].second.size();
  }
  cout << "number of cluster " << clusters.size() << endl;
  cout << "total num motifs seqs " << total_num_motif_seqs << endl;
  ////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  string ooo = output_file;
  ooo += "meme.format.txt";
  ofstream fmeme(ooo.c_str());
  fmeme << "MEME version 4" << endl;
  fmeme << endl;
  fmeme << "ALPHABET= ACDEFGHIKLMNPQRSTVWY" << endl;
  fmeme << endl;
  string ALPHABET= "ACDEFGHIKLMNPQRSTVWY";
  for (int i = 0; i < clusters.size(); i++) {
    fmeme << "MOTIF " << clusters[i].first << endl;
    fmeme << "letter-probability matrix: alength= 20 w= " << clusters[i].second[0].size() << endl;
    vector<vector<double> > pro(clusters[i].second[0].size(), vector<double>(26, 0.0));
    for (int j = 0; j < clusters[i].second.size(); j++) {
      for(int k = 0;k < clusters[i].second[j].size();++k)
        pro[k][clusters[i].second[j][k] - 'A'] += 1.0;
    }
    for(int k = 0;k < clusters[i].second[0].size();++k) {
      double sum = 0.0;
      for(int i = 0;i < 26;++i) {
        sum += pro[k][i];
      }
      fmeme << pro[k][0] / sum;
      for(int i = 1;i < 26;++i) {
        if(ALPHABET.find_first_of('A' + i) == string::npos) continue;
        fmeme << " " << pro[k][i] / sum;
      }
      fmeme << endl;
    }
  }
  fmeme.close();

  ////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  vector<Point> centers;
  for (int i = 0; i < clusters.size(); i++) {
    cout << i << endl;
    vector<Point> points;
    for (int j = 0; j < clusters[i].second.size(); j++) {
      points.push_back(KmerToCoordinates(clusters[i].second[j]));
    }
    Point center = Center(points);
    centers.push_back(center);
  }

//    for (int p = 0; p < points.size(); ++p) {
//      fout << PairwiseDistance(center, points[p]) << endl;
//    }
//    fout.close();
    /////////////////////////////////////////////
//    sprintf(out, "%s_randomPoints2Center_cluster%d.txt", output_file.c_str(),
//            cnt);
//    fout.open(out);
//    unordered_set<string> kmers;
//     vector<Point> randomPoints;
//     srand(time(NULL));
//     while(randomPoints.size() < 1000000) {
//       string kmer;
//       for(int j = 0;j < 25;++j) {
//         int r = rand() % 20;
//         kmer += AA20[r];
//       }
//       if(kmers.find(kmer) == kmers.end()) {
//         kmers.insert(kmer);
//         //cout << kmer << endl;
//         randomPoints.push_back(KmerToCoordinates(kmer));
//       }
//     }
//    for (int p = 0; p < randomPoints.size(); ++p) {
//      fout << PairwiseDistance(center, randomPoints[p]) << endl;
//    }
//    fout.close();
  //}
  /*ofstream fcenter("centerDistance.txt");
  for(int p = 0;p < centers.size();p++) {
    for(int q = p + 1;q < centers.size();q++) {
      fcenter << PairwiseDistance(centers[p], centers[q]) << endl;
    }
  }
  fcenter.close();
 */
  ofstream foutcenter("Pfam.entries.centers.point_startall.txt");
  for (int p = 0; p < centers.size(); p++) {
    foutcenter << centers[p].data[0];
    for (int q = 1; q < DIMENSION; q++) {
      foutcenter << " " << centers[p].data[q];
    }
    foutcenter << endl;
  }
  foutcenter.close();
}

int main(int argc, const char *argv[]) {
/*	ofstream fout("random_pairwisedistance.txt");
    DIMENSION = AACoordinateSize * 25;
	unordered_set<string> kmers;
	vector<Point> randomPoints;
	srand (time(NULL));
	cout << "herer" << endl;
	while(randomPoints.size() < 10000) {
		cout << randomPoints.size() << endl;
		string kmer;
		for(int j = 0;j < 25;++j) {
			int r = rand() % 20;
			//cout << r << endl;
			kmer += AA20[r];
		}
		//cout << kmer << endl;
		if(kmers.find(kmer) == kmers.end()) {
			//cout << kmer << endl;
			kmers.insert(kmer);
			randomPoints.push_back(KmerToCoordinates(kmer));
		}
	}
	for (int p = 0; p < randomPoints.size(); ++p) {
		for (int q = p + 1; q < randomPoints.size(); ++q) {
			fout << PairwiseDistance(randomPoints[p], randomPoints[q]) << endl;
		}
	}
	fout.close();
*/

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
      fprintf(stdout, "[WELCOME TO pairwiseDistanceSampling v%s]\n", hclust_version);
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
    DIMENSION = AACoordinateSize * kmer_length;
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

    clusterDistance(clusters, output_file);

  } catch (const SMITHLABException &e) {
    fprintf(stderr, "%s\n", e.what().c_str());
    return EXIT_FAILURE;
  } catch (std::bad_alloc &ba) {
    fprintf(stderr, "ERROR: could not allocate memory\n");
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
