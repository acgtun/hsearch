#ifndef ORF_H_
#define ORF_H_

#include "./../util/bio_util.h"
#include <map>
/*
 TTT F Phe      TCT S Ser      TAT Y Tyr      TGT C Cys
 TTC F Phe      TCC S Ser      TAC Y Tyr      TGC C Cys
 TTA L Leu      TCA S Ser      TAA * Ter      TGA * Ter
 TTG L Leu i    TCG S Ser      TAG * Ter      TGG W Trp

 CTT L Leu      CCT P Pro      CAT H His      CGT R Arg
 CTC L Leu      CCC P Pro      CAC H His      CGC R Arg
 CTA L Leu      CCA P Pro      CAA Q Gln      CGA R Arg
 CTG L Leu i    CCG P Pro      CAG Q Gln      CGG R Arg

 ATT I Ile      ACT T Thr      AAT N Asn      AGT S Ser
 ATC I Ile      ACC T Thr      AAC N Asn      AGC S Ser
 ATA I Ile      ACA T Thr      AAA K Lys      AGA R Arg
 ATG M Met i    ACG T Thr      AAG K Lys      AGG R Arg

 GTT V Val      GCT A Ala      GAT D Asp      GGT G Gly
 GTC V Val      GCC A Ala      GAC D Asp      GGC G Gly
 GTA V Val      GCA A Ala      GAA E Glu      GGA G Gly
 GTG V Val      GCG A Ala      GAG E Glu      GGG G Gly
 */

const string Base1 = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG";
const string Base2 = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG";
const string Base3 = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG";
const string AAs = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";

class ORF {
 public:
  ORF(){
    SetGeneticCodes();
  }
  ~ORF(){}
 private:
  map<string, char> mapGeneticCodes;

  void SetGeneticCodes();
  char complimentBase(char nt);
  void ReverseCompliment(string & strRead_rev, const string & strRead,
                         const int & len);
 public:
  void orf6(const string & querySeq, vector<string> & translatedAA);
};



#endif /* ORF_H_ */
