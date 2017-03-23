#ifndef STOCKHOLM_H_
#define STOCKHOLM_H_

#include <map>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>

#include <unordered_set>
#include <unordered_map>

using std::map;
using std::endl;
using std::cout;
using std::setw;
using std::vector;
using std::string;
using std::fstream;
using std::ifstream;
using std::ofstream;
using std::unordered_set;
using std::unordered_map;
using std::setprecision;
using std::istringstream;

// ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/userman.txt

struct Alignment {
  Alignment() {
    seqname = "NULL";
    start = -1;
    stop = -1;
    alignment = "NULL";
  }

  string seqname;
  size_t start;
  size_t stop;
  string alignment;
};

struct PFAMENTRY {
  void Clear() {
    Sequences.clear();
  }

  // Header Section
  string AC;
  string ID;
  string DE;
  string AU;
  string BM;
  string SM;
  string SE;
  string GA;
  string NC;
  string TC;
  string TP;
  string PI;

  // Reference Section;
  string WK;
  string DC;
  string DR;
  string RC;
  string RN;
  string RM;
  string RT;
  string RA;
  string RL;

  // Comment Section
  string CC;

  // Alignment Section
  string NE;
  string NL;
  size_t SQ;
  map<string, Alignment> Sequences;

  void Output(ofstream& fout) {
    fout << "ID: " << ID << endl;
    fout << "AC: " << AC << endl;
    fout << "DE: " << DE << endl;
    fout << "TP: " << TP << endl;
    fout << "SQ: " << SQ << endl;
    for (map<string, Alignment>::iterator it = Sequences.begin();
        it != Sequences.end(); ++it) {
      fout << it->second.seqname;
      for (size_t i = it->second.seqname.size(); i < 20; ++i) {
        fout << " ";
      }
      fout << "\t" << setw(8) << it->second.start << "\t" << setw(8)
           << it->second.stop << "\t";
      ///  Alignment: Delete '.' and samll letters
      for (size_t i = 0; i < it->second.alignment.size(); ++i) {
        if (it->second.alignment[i] == '.'
            || (it->second.alignment[i] >= 'a' && it->second.alignment[i] <= 'z')) {
          continue;
        }
        fout << it->second.alignment[i];
      }

      // Alignment
      fout << "\t" << it->second.alignment << endl;
    }
    fout << "-----------------END----------------" << endl;
  }

  void Output_LEN(ofstream& fout, const int& LEN,
                  unordered_set<string>& exist) {
    vector<string> motifs;
    for (map<string, Alignment>::iterator it = Sequences.begin();
        it != Sequences.end(); ++it) {
      ///  Alignment: Delete '.' and samll letters
      int cnt = 0;
      string cur = "";
      for (size_t i = 0; i < it->second.alignment.size(); ++i) {
        if (it->second.alignment[i] == '.'
            || (it->second.alignment[i] >= 'a' && it->second.alignment[i] <= 'z')) {
          continue;
        }

        cur += it->second.alignment[i];
        cnt++;
        if (cnt >= LEN)
          break;
      }
      if (cur.size() != LEN)
        continue;
      if (cur.find_first_of('-') == string::npos) {
        if (exist.find(cur) == exist.end()) {
          motifs.push_back(cur);
          exist.insert(cur);
        }
      }
    }

    if (motifs.size() > 0) {
      fout << "#" << "ID:" << ID << "#AC:" << AC << endl;
      for (size_t i = 0; i < motifs.size(); ++i) {
        fout << motifs[i] << endl;
      }
    }
    //fout << "-----------------END----------------" << endl;
  }

  void Output_LEN_all_kemrs(const int& LEN,
                            unordered_map<uint32_t, vector<string> >& motifs_map,
                            uint32_t& motif_id) {
    map<string, Alignment>::iterator it = Sequences.begin();
    int l = it->second.alignment.size();
    srand(time(NULL));
    for (size_t p = 0; p < l; ++p) {
      int r = rand() % 2;
      if(r != 0) continue;
      //cout << "p=" << p << endl;
      vector<string> motifs;
      for (map<string, Alignment>::iterator it = Sequences.begin(); it != Sequences.end(); ++it) {
        ///  Alignment: Delete '.' and samll letters
        int cnt = 0;
        string cur = "";
        for (size_t i = p; i < it->second.alignment.size(); ++i) {
          if (it->second.alignment[i] == '.'
              || (it->second.alignment[i] >= 'a'
                  && it->second.alignment[i] <= 'z')) {
            continue;
          }

          cur += it->second.alignment[i];
          cnt++;
          if (cnt >= LEN)
            break;
        }
        if (cur.size() != LEN)
          continue;
        if (cur.find_first_of('-') == string::npos) {
          //if (exist.find(cur) == exist.end()) {
            motifs.push_back(cur);
            //exist.insert(cur);
         // }
        }
      }

      if (motifs.size() > 0) {
        //fout << "#" << "ID:" << ID << "#AC:" << AC << "#start"<< p << endl;
        for (size_t i = 0; i < motifs.size(); ++i) {
          //fout << motifs[i] << endl;
          motifs_map[motif_id].push_back(motifs[i]);
        }
        motif_id++;
        cout << p << " " << motif_id << endl;
      }
      //fout << "-----------------END----------------" << endl;
    }
  }

  void OutputMotifs(ofstream& fout) {
    if(TP != "Motif") return;
    fout << "ID: " << ID << endl;
    fout << "AC: " << AC << endl;
    fout << "DE: " << DE << endl;
    fout << "TP: " << TP << endl;
    fout << "SQ: " << SQ << endl;
    for (map<string, Alignment>::iterator it = Sequences.begin();
        it != Sequences.end(); ++it) {
      fout << it->second.seqname;
      for (size_t i = it->second.seqname.size(); i < 20; ++i) {
        fout << " ";
      }
      fout << "\t" << setw(8) << it->second.start << "\t" << setw(8)
           << it->second.stop << "\t";
      ///  Alignment: Delete '.' and samll letters
      for (size_t i = 0; i < it->second.alignment.size(); ++i) {
        if (it->second.alignment[i] == '.'
            || (it->second.alignment[i] >= 'a' && it->second.alignment[i] <= 'z')) {
          continue;
        }
        fout << it->second.alignment[i];
      }

      // Alignment
      fout << "\t" << it->second.alignment << endl;
    }
    fout << "-----------------END----------------" << endl;
  }

  void Output_NoGaps(ofstream& fout) {
    fout << "ID: " << ID << endl;
    fout << "AC: " << AC << endl;
    fout << "DE: " << DE << endl;
    fout << "TP: " << TP << endl;
    fout << "SQ: " << SQ << endl;
    for (map<string, Alignment>::iterator it = Sequences.begin();
        it != Sequences.end(); ++it) {
      bool find = false;
      for (size_t i = 0; i < it->second.alignment.size(); ++i) {
        if (/*it->second.alignment[i] == '-'
            || */(it->second.alignment[i] >= 'a' && it->second.alignment[i] <= 'z')) {
          find = true;
          break;
        }
      }
      if (find == true)
        continue;

      fout << it->second.seqname;
      for (size_t i = it->second.seqname.size(); i < 20; ++i) {
        fout << " ";
      }
      fout << "\t" << setw(8) << it->second.start << "\t" << setw(8)
           << it->second.stop << "\t";
      ///  Alignment: Delete '.' and samll letters
      for (size_t i = 0; i < it->second.alignment.size(); ++i) {
        if (it->second.alignment[i] == '.'
            || (it->second.alignment[i] >= 'a' && it->second.alignment[i] <= 'z')) {
          continue;
        }
        fout << it->second.alignment[i];
      }
      fout << endl;
    }
    fout << "-----------------END----------------" << endl;
  }
};

#endif /* STOCKHOLM_H_ */
