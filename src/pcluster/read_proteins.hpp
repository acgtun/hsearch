#ifndef _READONEPROTEIN_H_
#define _READONEPROTEIN_H_

#include "util.hpp"

#include "smithlab_os.hpp"
#include "OptionParser.hpp"

#define DEBUG
class ProteinDB {
public:
	ProteinDB(const string& _file_name) :
			file_name(_file_name) {
		ReadFASTAFile();
	}

	~ProteinDB() {
#ifdef DEBUG
		for (uint32_t i = 0; i < pro_seqs.size(); ++i) {
			cout << pro_names[i] << endl;
			for (uint32_t j = 0; j < pro_seqs[i].size(); ++j) {
				cout << pro_seqs[i][j];
			}
			cout << endl;
		}
#endif
	}

	void ReadFASTAFile() {
		ifstream fin(file_name.c_str());
		if (!fin) {
			throw SMITHLABException("cannot open input file " + file_name);
		}
		vector<char> sequence;
		string line;
		while (getline(fin, line)) {
			if (line[0] == '>') {
				if (sequence.size() != 0) {
					pro_seqs.push_back(sequence);
					sequence.clear();
				}
				uint32_t space_pos = line.find_first_of(' ');
				if (space_pos == string::npos) {
					pro_names.push_back(line.substr(1));
				} else {
					pro_names.push_back(line.substr(1, space_pos - 1));
				}
				continue;
			}
			for (uint32_t i = 0; i < line.size(); ++i) {
				if (AA20.find_first_of(line[i]) != string::npos) {
					sequence.push_back(toupper(line[i]));
				} else if (isalpha(line[i])) {
					sequence.push_back(AA20[rand() % 20]);  //todo
				}
			}
		}
		if (sequence.size() != 0) {
			pro_seqs.push_back(sequence);
		}

		fin.close();
		num_of_proteins = pro_seqs.size();
	}

	void OutputProtein(const uint32_t& i) {
		for (uint32_t j = 0; j < pro_seqs[i].size(); ++j) {
			cout << pro_seqs[i][j];
		}
		cout << endl;
	}

	vector<vector<char> > pro_seqs;
	vector<string> pro_names;

	uint32_t num_of_proteins;
	string file_name;
};

#endif /* _READONEPROTEIN_H_ */
