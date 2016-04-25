#include "read_proteins.hpp"
#include "lsh.hpp"

int main(int argc, const char *argv[]) {
  // argv[1] is the protein database file
  uint32_t feature_size = static_cast<uint32_t>(pow(8, HASHLEN));
  cout << "feature_size " << feature_size << endl;
  uint32_t bit_num = 32;
  double sigma = 0.2;
  // need 1 feature_node->index=-1 to mark end

  KLSH klsh_low(feature_size, bit_num, sigma);
  vector<double> p(feature_size);
  vector<int> features(feature_size, 0);
  ReadOneProtein proteinDB(argv[1]);
  clock_t start = clock();
  while (!proteinDB.compelted) {
    vector<char> pro_seq = proteinDB.GetNextProtein();
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
    printf("%u %lu\n", proteinDB.num_of_proteins, klsh_low.GetHashValue(p));
  }

  printf("%lf seconds\n", (clock() - start) / (double) CLOCKS_PER_SEC);

  return 0;
}

