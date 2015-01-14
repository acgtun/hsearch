#ifndef BIO_UTIL_H_
#define BIO_UTIL_H_

#include <stdint.h>

#include <string>

using std::string;

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

#define HASHAALEN 6
const string AA20 = "ARNDCEQGHILKMFPSTWYV";
const string NODELABEL = " ARNDCEQGHILKMFPSTWYV";
                   /*A   B  C  D  E   F  G  H  I   J   K   L   M  N   O   P  Q  R   S   T   U   V   W   X   Y   Z */
const int base[] = { 0, -1, 4, 3, 6, 13, 7, 8, 9, -1, 11, 10, 12, 2, -1, 14, 5, 1, 15, 16, -1, 19, 17, -1, 18, -1 };
const uint32_t basep[] = { 1, 20, 400, 8000, 160000, 3200000, 64000000, 1280000000 };
const uint32_t base_preCal[20][6] = {
    { 0 * basep[0], 0 * basep[1], 0 * basep[2], 0 * basep[3], 0 * basep[4], 0 * basep[5] },
    { 1 * basep[0], 1 * basep[1], 1 * basep[2], 1 * basep[3], 1 * basep[4], 1 * basep[5] },
    { 2 * basep[0], 2 * basep[1], 2 * basep[2], 2 * basep[3], 2 * basep[4], 2 * basep[5] },
    { 3 * basep[0], 3 * basep[1], 3 * basep[2], 3 * basep[3], 3 * basep[4], 3 * basep[5] },
    { 4 * basep[0], 4 * basep[1], 4 * basep[2], 4 * basep[3], 4 * basep[4], 4 * basep[5] },
    { 5 * basep[0], 5 * basep[1], 5 * basep[2], 5 * basep[3], 5 * basep[4], 5 * basep[5] },
    { 6 * basep[0], 6 * basep[1], 6 * basep[2], 6 * basep[3], 6 * basep[4], 6 * basep[5] },
    { 7 * basep[0], 7 * basep[1], 7 * basep[2], 7 * basep[3], 7 * basep[4], 7 * basep[5] },
    { 8 * basep[0], 8 * basep[1], 8 * basep[2], 8 * basep[3], 8 * basep[4], 8 * basep[5] },
    { 9 * basep[0], 9 * basep[1], 9 * basep[2], 9 * basep[3], 9 * basep[4], 9 * basep[5] },
    {10 * basep[0],10 * basep[1],10 * basep[2],10 * basep[3],10 * basep[4],10 * basep[5] },
    {11 * basep[0],11 * basep[1],11 * basep[2],11 * basep[3],11 * basep[4],11 * basep[5] },
    {12 * basep[0],12 * basep[1],12 * basep[2],12 * basep[3],12 * basep[4],12 * basep[5] },
    {13 * basep[0],13 * basep[1],13 * basep[2],13 * basep[3],13 * basep[4],13 * basep[5] },
    {14 * basep[0],14 * basep[1],14 * basep[2],14 * basep[3],14 * basep[4],14 * basep[5] },
    {15 * basep[0],15 * basep[1],15 * basep[2],15 * basep[3],15 * basep[4],15 * basep[5] },
    {16 * basep[0],16 * basep[1],16 * basep[2],16 * basep[3],16 * basep[4],16 * basep[5] },
    {17 * basep[0],17 * basep[1],17 * basep[2],17 * basep[3],17 * basep[4],17 * basep[5] },
    {18 * basep[0],18 * basep[1],18 * basep[2],18 * basep[3],18 * basep[4],18 * basep[5] },
    {19 * basep[0],19 * basep[1],19 * basep[2],19 * basep[3],19 * basep[4],19 * basep[5] }
};


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

/* Given a 6-mer amino acids sequence, translate it to a integer number.
 * Use 20-based number, A is 0, R is 1, N is 3 and so on.
 * */
inline uint32_t GetHashValueAA(const string& seed) {
  uint32_t hash_value = 0;
  for (uint32_t i = 0; i < HASHAALEN; i++) {
    hash_value += base_preCal[base[seed[i] - 'A']][i];
  }
  return hash_value;
}

/* Given a integer, translate it to a 6-mer amino acids, it's also 20-based.*/
inline void DeCodeAA(const uint32_t& hash_value, string& seed) {
  uint32_t value = hash_value;
  for (uint32_t i = 0; i < HASHAALEN; i++) {
    seed += AA20[value % 20];
    value /= 20;
  }
}

#endif /* BIO_UTIL_H_ */
