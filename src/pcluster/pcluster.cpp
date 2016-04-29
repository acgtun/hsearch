#include "read_proteins.hpp"
#include "lsh.hpp"

int main(int argc, const char *argv[]) {
	string outfile = argv[1];
	outfile += "pcluster_out";
	freopen(outfile.c_str(), "w", stdout);
	// argv[1] is the protein database file
	uint32_t feature_size = static_cast<uint32_t>(pow(8, HASHLEN));
	cout << "feature_size " << feature_size << endl;
	uint32_t bit_num = 16;
	double sigma = 0.2;
	// need 1 feature_node->index=-1 to mark end

	KLSH klsh_low(feature_size, bit_num, sigma);
	vector<double> p(feature_size);
	vector<int> features(feature_size, 0);
	clock_t start = clock();
	ProteinDB proteinDB(argv[1]);
	for (size_t i = 0; i < proteinDB.num_of_proteins; ++i) {
		vector<char>& pro_seq = proteinDB.pro_seqs[i];
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
		proteinDB.hash_buckets[klsh_low.GetHashValue(p)].push_back(i);
	}
	cout << "num of buckets: " << proteinDB.hash_buckets.size() << endl;
	for (unordered_map<uint32_t, vector<uint32_t> >::iterator it =
			proteinDB.hash_buckets.begin(); it != proteinDB.hash_buckets.end();
			++it) {
		cout << it->first << " " << it->second.size() << endl;
		for (size_t i = 0; i < it->second.size(); ++i) {
			proteinDB.OutputProtein(it->second[i]);
		}
		cout << "cluster -----------------------------------" << endl;
	}

	printf("%lf seconds\n", (clock() - start) / (double) CLOCKS_PER_SEC);

	return 0;
}

