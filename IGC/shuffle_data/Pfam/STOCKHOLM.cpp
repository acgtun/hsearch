#include "STOCKHOLM.h"

//  g++ -std=c++11 STOCKHOLM.cpp -o pfam_test

string SPlitGFGCGSGR(const string& line, const string& label) {
  size_t pos = line.find(label);
  istringstream iss(line.substr(pos + label.size()));
  string ret;
  iss >> ret;
  return ret;
}

size_t ToInt(const string& Number) {
  istringstream iss(Number);
  size_t num;
  iss >> num;

  return num;
}

string GetAlingment(const string& line, Alignment& aln) {
  istringstream iss(line);
  string seqname, alignment;
  iss >> seqname >> alignment;
  size_t pos1 = seqname.find_first_of('/');
  size_t pos2 = seqname.find_first_of('-');

  aln.seqname = seqname.substr(0, pos1);
  aln.start = ToInt(seqname.substr(pos1 + 1, pos2 - pos1 - 1));
  aln.stop = ToInt(seqname.substr(pos2 + 1));
  aln.alignment = alignment;

  return seqname;
}

string GetSeqName(const string& line) {
  istringstream iss(line);
  string seqname, label;
  iss >> label >> seqname;
  size_t pos1 = seqname.find_first_of('/');

  return seqname.substr(0, pos1);
}

void ReadPfam(const string& file_name, vector<PFAMENTRY>& pfam_entries) {
  ifstream fin(file_name.c_str());
  string line;
  PFAMENTRY pfam_entry;
  while (getline(fin, line)) {
    if (line == "# STOCKHOLM 1.0") {
      continue;
    }
    if (line == "//") {
      pfam_entries.push_back(pfam_entry);
      pfam_entry.Clear();
      continue;
    }
    if (line.find("#=GF ID") != string::npos) {
      pfam_entry.ID = SPlitGFGCGSGR(line, "ID");
      continue;
    }
    if (line.find("#=GF AC") != string::npos) {
      pfam_entry.AC = SPlitGFGCGSGR(line, "AC");
      continue;
    }
    if (line.find("#=GF DE") != string::npos) {
      pfam_entry.DE = SPlitGFGCGSGR(line, "DE");
      continue;
    }
    if (line.find("#=GF TP") != string::npos) {
      pfam_entry.TP = SPlitGFGCGSGR(line, "TP");
      continue;
    }
    if (line.find("#=GF SQ") != string::npos) {
      pfam_entry.SQ = ToInt(SPlitGFGCGSGR(line, "SQ"));
      continue;
    }
    if (line.substr(0, 4) == "#=GS") {
      /*
       string seqname = GetSeqName(line);
       if (pfam_entry.Sequences.find(seqname) != pfam_entry.Sequences.end()) {
       cout << "repeatname " << seqname << endl;
       }
       */
      continue;
    }
    if (line[0] != '#') {
      Alignment aln;
      string full_seqname = GetAlingment(line, aln);
      if (pfam_entry.Sequences.find(aln.seqname)
          != pfam_entry.Sequences.end()) {
        cout << " repeatname " << aln.seqname << endl;
      }
      pfam_entry.Sequences.insert(make_pair(full_seqname, aln));
    }
  }
  fin.close();
}

int main(int argc, const char **argv) {
  // command: ./read file_name
  //string outfile = argv[1];
  //outfile += "_full_alignment.out";
 // ofstream fout(outfile.c_str());
 // string nogaps = argv[1];
 // nogaps += "_nogaps_alignment.out";
 // ofstream fnogaps(nogaps.c_str());

 // string motifs = argv[1];
  //motifs += "_motifs_alignment.out";
  //ofstream fmotifs(motifs.c_str());


  int len;
  sscanf(argv[2], "%d", &len);
  cout << "len = " << len << endl;
  char LENS[100];
  sprintf(LENS, "%s_LEN%d_alignment_seed.out", argv[1], len);
  ofstream fout_len(LENS);

  unordered_set<string> exist;
  vector<PFAMENTRY> pfam_entries;
  ReadPfam(argv[1], pfam_entries);
  unordered_map<uint32_t, vector<string> > motifs_map;
  uint32_t motif_id = 0;
  srand(time(NULL));
  for (size_t i = 0; i < pfam_entries.size(); ++i) {
    if(i %100 == 0)cout << "i = " << i << " " << pfam_entries.size() << endl;
    int r = rand() % 10;
    if(r != 0) continue;

    //pfam_entries[i].Output(fout);
    //pfam_entries[i].Output_NoGaps(fnogaps);
    //pfam_entries[i].OutputMotifs(fmotifs);
    pfam_entries[i].Output_LEN(fout_len, len, exist);
    //pfam_entries[i].Output_LEN_all_kemrs(len, motifs_map, motif_id);
    if (pfam_entries[i].Sequences.size() != pfam_entries[i].SQ) {
      cout << pfam_entries[i].Sequences.size() << " " << pfam_entries[i].SQ
          << endl;
      cout << "error" << endl;
    }
  }

  unordered_map<string, uint32_t> count;
  for(unordered_map<uint32_t, vector<string> >::iterator it = motifs_map.begin();it != motifs_map.end();++it) {
    for(uint32_t i = 0;i < it->second.size();++i) {
      count[it->second[i]]++;
    }
  }

  motif_id = 0;
  for(unordered_map<uint32_t, vector<string> >::iterator it = motifs_map.begin();it != motifs_map.end();++it) {
    fout_len << "#motfis" << motif_id << endl;
    for(uint32_t i = 0;i < it->second.size();++i) {
      if(count[it->second[i]] > 1) continue;
      fout_len << it->second[i] << endl;
    }
    motif_id++;
  }


  //fout.close();
  //fnogaps.close();
  fout_len.close();
 // fmotifs.close();

  return 0;
}
