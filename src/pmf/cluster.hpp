#include "util.hpp"
#include "lsh.hpp"

#include <string>
#include <vector>

#include <algorithm>
#include <unordered_map>

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

  char out[100];
  sprintf(out, "%s_cluster0.clustering", output_file.c_str());
  ofstream fcluster(out);
  for (HASH_BUCKETS::iterator it = hash_buckets.begin();
      it != hash_buckets.end(); ++it) {
    fcluster << "cluster" << "\t" << it->second.size() << endl;
    for (uint32_t i = 0; i < it->second.size(); ++i) {
      fcluster << it->second[i] << endl;
    }
  }
  fcluster.close();
}

struct SortClustersByFrequency {
  explicit SortClustersByFrequency(const vector<uint32_t>& _frequency)
      : frequency(_frequency) {
  }
  bool operator()(const vector<uint32_t>& a, const vector<uint32_t>& b) {
    if (a.size() == b.size()) {
      uint32_t count_a = 0;
      for (uint32_t i = 0; i < a.size(); ++i) {
        count_a += frequency[a[i]];
      }

      uint32_t count_b = 0;
      for (uint32_t i = 0; i < b.size(); ++i) {
        count_b += frequency[b[i]];
      }

      return count_a < count_b;
    }

    return a.size() > b.size();
  }

  const vector<uint32_t>& frequency;
};

void Clustering(const vector<pair<string, string> >& proteins,
                const vector<pair<uint32_t, uint32_t> >& kmers,
                const vector<uint32_t>& frequency, const uint32_t kmer_length,
                const uint32_t& num_of_hash_functions, const string& org_output,
                const int& start_num, const string& output_file) {
  uint32_t dimension = AACoordinateSize * kmer_length;
  uint32_t s_num = start_num + 1;
  ifstream fin(org_output.c_str());
  uint32_t num_hash = 0;
  uint32_t hash_size;
  string cluster_name;
  vector<double> point(dimension);
  while (num_hash++ < num_of_hash_functions - 1) {
    cout << "num_hash " << num_hash << endl;
    char out[1000];
    sprintf(out, "%s_cluster%u.clustering", output_file.c_str(), s_num++);
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

  sort(clusters.begin(), clusters.end(), SortClustersByFrequency(frequency));

  for (uint32_t c = 0; c < clusters.size(); ++c) {
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

void RandomSampling(const vector<pair<string, string> >& proteins,
                    const vector<pair<uint32_t, uint32_t> >& kmers,
                    const uint32_t kmer_length, const string& output_file) {
  random_device rd;
  default_random_engine generator(rd());
  uniform_int_distribution<int> m_uniform(1, kmers.size());
  uint32_t dimension = AACoordinateSize * kmer_length;
  LSH lsh(dimension);
  vector<double> point(dimension);
  vector<uint32_t> random_numbers(10000000, 0);
  for (size_t i = 0; i < 10000000; ++i) {
    random_numbers[i] = m_uniform(generator);
  }

  char out[1000];
  sprintf(out, "%s_ranomly_dotproductrandomsampling.txt", output_file.c_str());
  ofstream fout(out);
  for (size_t i = 0; i < 10000000; ++i) {
    KmerToCoordinates(proteins, kmers[random_numbers[i]], kmer_length, point);
    fout << lsh.DotProductPlusRandom(point) << endl;
  }
  fout.close();
}
