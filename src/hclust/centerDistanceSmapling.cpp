#include "util.hpp"
#include <stdint.h>
#include <stdio.h>
#include <vector>
#include <string>
#include <unordered_set>
#include "./../smithlab_cpp/smithlab_os.hpp"
#include "./../smithlab_cpp/OptionParser.hpp"

using namespace std;

#define MIN_SIZE_CLUSTER 50
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

//void randomSampling(vector<pair<string, vector<string> > >& clusters) {
//  vector<string > sequences;
//  for (int i = 0; i < clusters.size(); i++) {
//    for (int j = 0; j < clusters[i].second.size(); j++) {
//      sequences.push_back(clusters[i].second[j]);
//    }
//  }
//  //sprintf(out, "%s_randomPoints2Center_cluster%d.txt", output_file.c_str());
//  ofstream fout("randomPointsDistancexxxx.txt");
//  unordered_set < string > kmers;
//  vector<Point> randomPoints;
//  srand (time(NULL));
//  while(randomPoints.size() < 10000) {
//    string kmer;
//    int r = rand() % sequences.size();
//    if(kmers.find(sequences[r]) == kmers.end()) {
//      kmers.insert(sequences[r]);
//      cout << sequences[r] << endl;
//      randomPoints.push_back(KmerToCoordinates(sequences[r]));
//    }
//  }
//
//  for (int p = 0; p < randomPoints.size(); ++p) {
//    for (int q = p + 1; q < randomPoints.size(); ++q) {
//      fout << PairwiseDistance(randomPoints[q], randomPoints[p]) << endl;
//    }
//  }
//  fout.close();
//}

void cluster2datapoint(vector<pair<string, vector<string> > >& clusters,
                   const string& output_file) {
  string ooohclust = output_file;
  ooohclust += "hclust.format.txt";
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


  ofstream foutcenter(ooohclust.c_str());
  for (int p = 0; p < centers.size(); p++) {
    foutcenter << clusters[p].first << endl;
    foutcenter << centers[p].data[0];
    for (int q = 1; q < DIMENSION; q++) {
      foutcenter << " " << centers[p].data[q];
    }
    foutcenter << endl;
  }
  foutcenter.close();
}

void sequencedatabase2centers(const vector<Point>& kmers_proteins,
                              vector<pair<string, vector<string> > >& clusters,
                              const string& output_file) {
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
  ///////////////////////////
  char oooo[100];
  sprintf(oooo, "./pro2centerdis/%sinnercenter_protein_centers_%d.txt", output_file.c_str(), 0);
  cout << oooo << endl;
  ofstream fcenter(oooo);
  for(int i = 0;i < clusters.size();++i) {
    for(int j = i + 1;j < clusters.size();++j) {
      fcenter << PairwiseDistance(centers[i], centers[j]) << endl;
    }
  }
  fcenter.close();
  //////////////////////////////

  vector<Point> randomPoints;
  unordered_set<uint32_t> kmers;
  int i = 0;
  while(i < 100000) {
    //int r = rand() % kmers_proteins.size();
   // if (kmers.find(r) == kmers.end()) {
     // kmers.insert(r);
      randomPoints.push_back(kmers_proteins[i]);
      i++;
    //}
  }
  char ooo[100];
  sprintf(ooo, "./pro2centerdis/%sramdom_protein_centers_%d.txt", output_file.c_str(), 0);
  cout << ooo << endl;
  ofstream fout(ooo);
  for (int p = 0; p < centers.size(); p++) {
    for (int q = 0; q < randomPoints.size(); q++) {
      fout << PairwiseDistance(randomPoints[q], centers[p]) << endl;
    }
  }
  fout.close();

}


void meme_format_output( vector<pair<string, vector<string> > >& clusters,
                        const string& output_file) {
  string ooo = output_file;
  ooo += "meme.format.txt";
  FILE * fout = fopen(ooo.c_str(), "w");
  fprintf(fout, "MEME version 4\n\n");
  fprintf(fout, "ALPHABET= ACDEFGHIKLMNPQRSTVWY\n\n");
  string ALPHABET= "ACDEFGHIKLMNPQRSTVWY";
  for (int i = 0; i < clusters.size(); i++) {
    fprintf(fout, "MOTIF %s\n", clusters[i].first.c_str());
    fprintf(fout, "letter-probability matrix: alength= 20 w= %u\n", clusters[i].second[0].size());
    vector<vector<double> > pro(clusters[i].second[0].size(), vector<double>(26, 0.0));
    while(clusters[i].second.size() > 10) {
      clusters[i].second.pop_back();
    }

    for (int j = 0; j < clusters[i].second.size(); j++) {
      fprintf(fout,"%s\n", clusters[i].second[j].c_str());
    }

    fprintf(fout, "\n A    C    D    E    F    G    H    I    K    L    M    N    P    Q    R    S    T    V    W    Y\n");
    for (int j = 0; j < clusters[i].second.size(); j++) {
      for(int k = 0;k < clusters[i].second[j].size();++k)
        pro[k][clusters[i].second[j][k] - 'A'] += 1.0;
    }
    for(int k = 0;k < clusters[i].second[0].size();++k) {
      double sum = 0.0;
      for(int i = 0;i < 26;++i) {
        sum += pro[k][i];
      }
      fprintf(fout, "%.2lf", pro[k][0] / sum);
      for(int i = 1;i < 26;++i) {
        if(ALPHABET.find_first_of('A' + i) == string::npos) continue;
        fprintf(fout, " %.2lf", pro[k][i] / sum);
      }
      fprintf(fout, "\n");
    }
  }
  fclose(fout);
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
  string ooohclust = output_file;
  ooo += "meme.format.txt";
  ooohclust += "hclust.format.txt";
  ofstream fmeme(ooo.c_str());
  ofstream fclust(ooohclust.c_str());
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
  string center_name = output_file;
  center_name += "Pfam.entries.centers.point.txt";
  ofstream foutcenter(center_name.c_str());
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

    string protein_file;

    /* kmer length */
    uint32_t kmer_length;

    /* output file */
    string output_file;
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "pairwiseDistanceSampling",
        "");
    opt_parse.add_opt("kmers", 'k', "kmers file", true,
        kmers_file);
    opt_parse.add_opt("protein", 'd', "protein file", true,
          protein_file);
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

    ////
    vector<Point> kmers_proteins;
    fin.open(protein_file.c_str());
    while(getline(fin, line)) {
      //kmer_names.push_back(line);
      getline(fin, line);
      istringstream iss(line);
      Point point;
      for(uint32_t i = 0;i < DIMENSION;++i) {
        iss >> point.data[i];
      }
      kmers_proteins.push_back(point);
    }
    fin.close();

    //clusterDistance(clusters, output_file);
    //randomSampling(clusters);
    //cluster2datapoint(clusters, output_file);
    //meme_format_output(clusters, output_file);
    //cluster2datapoint(clusters,output_file);
    sequencedatabase2centers( kmers_proteins, clusters,
        output_file);

  } catch (const SMITHLABException &e) {
    fprintf(stderr, "%s\n", e.what().c_str());
    return EXIT_FAILURE;
  } catch (std::bad_alloc &ba) {
    fprintf(stderr, "ERROR: could not allocate memory\n");
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
