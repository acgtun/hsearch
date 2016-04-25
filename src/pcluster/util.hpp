/*
 *    This file contains several basic functions for global use.
 *
 *    Copyright (C) 2016 University of Southern California
 *
 *    Authors: Haifeng Chen and Ting Chen
 *
 *    This file is part of PCLUSTER.
 *
 *    PCLUSTER is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    PCLUSTER is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with PCLUSTER.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef UTIL_H_
#define UTIL_H_

#include <time.h>
#include <math.h>
#include <stdio.h>
#include <limits.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include <set>
#include <map>
#include <string>
#include <vector>
#include <utility>
#include <fstream>
#include <iostream>
#include <algorithm>

#include <tr1/unordered_map>
#include <tr1/unordered_set>


using std::cout;
using std::cerr;
using std::endl;
using std::set;
using std::map;
using std::pair;
using std::vector;
using std::string;
using std::make_pair;
using std::ofstream;
using std::ifstream;

using std::tr1::unordered_map;
using std::tr1::unordered_set;

typedef std::tr1::unordered_map<uint32_t, vector<pair<uint32_t, uint32_t> > > HashTable;


/*
 A0
 R1
 N2
 D3
 C4
 E5
 Q6
 G7
 H8
 I9
 L10
 K11
 M12
 F13
 P14
 S15
 T16
 W17
 Y18
 V19
 */
/* B J O U X Z*/

#define HASHLEN 3
#define ALPHABETSIZE 8

const string AA20 = "ARNDCEQGHILKMFPSTWYV";

const int AAINDEX[] =
{ 0, -1, 4, 3, 6, 13, 7, 8, 9, -1, 11, 10, 12, 2, -1, 14, 5, 1, 15, 16, -1, 19, 17, -1, 18, -1 };
/*A   B  C  D  E   F  G  H  I   J   K   L   M  N   O   P  Q  R   S   T   U   V   W   X   Y   Z */

// [A S T] [R K E D Q] [N H] [C] [G] [I V L M] [F Y W] [P]
//    0         1        2    3   4      5        6     7
const int REDUCEDAAINDEX[] =
{ 0, -1, 3, 1, 1, 6, 4, 2, 5, -1, 1, 5, 5, 2, -1, 7, 1, 1, 0, 0, -1, 5, 6, -1, 6, -1 };
/*A   B  C  D  E  F  G  H  I   J  K  L  M  N   O  P  Q  R  S  T   U  V  W   X  Y   Z */

const uint32_t BASEP[] = { 1, 8, 64, 512, 4096, 32768, 262144, 2097152, 16777216};

const int BLOSUM62[][20] = {
//A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V
{ 4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0 },  //A
{ -1, 5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3 },  //R
{ -2, 0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3 },  //N
{ -2,-2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3 },  //D
{  0,-3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1 },  //C
{ -1, 1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2 },  //Q
{ -1, 0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2 },  //E
{ 0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3 },  //G
{ -2, 0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3 },  //H
{ -1,-3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3 },  //I
{ -1,-2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1 },  //L
{ -1, 2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2 },  //K
{ -1,-1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1 },  //M
{ -2,-3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1 },  //F
{ -1,-2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2 },  //P
{ 1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2 },  //S
{ 0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0 },  //T
{ -3,-3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3 },  //W
{ -2,-2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1 },  //Y
{ 0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3,  -1, 4 } };  //V
//  http://www.ncbi.nlm.nih.gov/BLAST/blastcgihelp.shtml#get_subsequence
//  http://www.ncbi.nlm.nih.gov/Class/FieldGuide/BLOSUM62.txt

const double REDUCEDBLOSUM62[ALPHABETSIZE][ALPHABETSIZE] = {
{ 1.88889,  -0.8,     -1,       -0.666667, -0.666667, -1.08333, -2.22222, -1       },
{-0.8,       1.52,    -0.1,     -3.2,      -1.8,      -2.35,    -2.66667, -1.2     },
{-1,        -0.1,      4,       -3,        -1,        -2.75,    -1.66667, -2       },
{-0.666667, -3.2,     -3,        9,        -3,        -1,       -2,       -3       },
{-0.666667, -1.8,     -1,       -3,         6,        -3.5,     -2.66667, -2       },
{-1.08333,  -2.35,    -2.75,    -1,        -3.5,       2.3125,  -1.16667, -2.5     },
{-2.22222,  -2.66667, -1.66667, -2,        -2.66667,  -1.16667,  4,       -3.66667 },
{-1,        -1.2,     -2,       -3,        -2,        -2.5,     -3.66667,  7       } };


inline void MemoryAllocateCheck(void* pointer, const char* file, int line) {
  if (pointer == NULL) {
    fprintf(stderr, "Memory allocate error in %s at line %d\n", file, line);
    exit(EXIT_FAILURE);
  }
}

inline void FileOpenCheck(FILE* pfile, const char* file, int line) {
  if (pfile == NULL) {
    fprintf(stderr, "File open error in %s at line %d\n", file, line);
    exit(EXIT_FAILURE);
  }
}

#define FILE_OPEN_CHECK(pfile) (FileOpenCheck( pfile, __FILE__, __LINE__))
#define MEMORY_ALLOCATE_CHECK(pointer)  (MemoryAllocateCheck(pointer, __FILE__, __LINE__))

#define FREAD_CHECK(func, size) { \
  uint32_t s = func; \
  if(s != size) { \
    fprintf(stderr, "read file error. --- %s:%s:%d\n", __FILE__, __func__, __LINE__); \
    exit(EXIT_FAILURE); \
  } \
}

//#define DEBUG
#ifdef DEBUG
#define DEBUG_INFO(msg) { \
  fprintf(stderr, "%s\n", msg); \
}
#else
#define DEBUG_INFO(msg)
#endif

#define TIME_INFO(func, msg) { \
  clock_t start_t, end_t; \
  start_t = clock(); \
  func; \
  end_t = clock(); \
  fprintf(stderr, "[%s TAKES %.3lf SECONDS]\n", msg, \
         (double) ((end_t - start_t) / CLOCKS_PER_SEC )); \
}

/* M8Results is a data structure to store the results of a protein alingment, same as BLAST -m 8*/
struct M8Results {
  M8Results(const string& pro_name, const double& idty, const int& ali_len,
            const int& mis, const int& gap, const uint32_t& q_start,
            const uint32_t& q_end, const uint32_t& p_start,
            const uint32_t& p_end, const double& e_value,
            const double& bitScore)
      : protein_name(pro_name),
        identity(idty),
        aligned_len(ali_len),
        mismatch(mis),
        gap_open(gap),
        qs(q_start),
        qe(q_end),
        ps(p_start),
        pe(p_end),
        evalue(e_value),
        bit_score(bitScore) {
  }
  M8Results() {
    protein_name = "";

    identity = 0.0;
    aligned_len = 0;
    mismatch = 0;
    gap_open = 0;

    qs = 0;
    qe = 0;
    ps = 0;
    pe = 0;

    evalue = 0.0;
    bit_score = 0.0;
  }

  string protein_name;

  double identity;
  int aligned_len;
  int mismatch;
  int gap_open;

  int qs;
  int qe;
  int ps;
  int pe;

  double evalue;
  double bit_score;

  static bool SORT_CMP_EValue(const M8Results& a, const M8Results& b) {
    return a.evalue < b.evalue;
  }
};

inline uint32_t Kmer2Integer(const char* kmer) {
  uint32_t hash_value = 0;
  for (uint32_t i = 0; i < HASHLEN; ++i) {
    hash_value += REDUCEDAAINDEX[kmer[i] - 'A'] * BASEP[i];
  }
  return hash_value;
}

/* transfer the integer to 8-based number */
inline string Integer2KmerDigit(const uint32_t& hash_value) {
  string kmer;
  uint32_t n = hash_value, j = 0;
  while (n) {
    kmer += 48 + n % ALPHABETSIZE;
    j++;
    n /= ALPHABETSIZE;
  }
  while (j < HASHLEN) {
    kmer += 48;
    j++;
  }
  return kmer;
}

/* Output the alingment results for one query to fout */
inline void DisplayResults(const string& query_name,
                           const string& database_file,
                           const vector<M8Results>& aligned_results,
                           const uint32_t& num_of_results,
                           const int& outfmt, ofstream& fout) {
  if (outfmt == 7) {
    fout << "# PCLUSTER 1.0.0 April, 2016" << endl;
    fout << "# Query: " << query_name << endl;
    fout << "# Database: " << database_file << endl;
    fout
        << "# Fields: query id, subject id, % identity, alignment length, mismatches, "
            "gap opens, q. start, q. end, s. start, s. end, evalue, bit score"
        << endl;
    fout << "# " << aligned_results.size() << " hits found" << endl;
  }
  for (uint32_t i = 0; i < num_of_results; i++) {
    fout << query_name << "\t" << aligned_results[i].protein_name << "\t"
        << aligned_results[i].identity << "\t" << aligned_results[i].aligned_len
        << "\t" << aligned_results[i].mismatch << "\t"
        << aligned_results[i].gap_open << "\t" << aligned_results[i].qs << "\t"
        << aligned_results[i].qe << "\t" << aligned_results[i].ps << "\t"
        << aligned_results[i].pe << "\t" << aligned_results[i].evalue << "\t"
        << aligned_results[i].bit_score << endl;
  }
}

#endif /* UTIL_H_ */
