#include <stdlib.h>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <unordered_set>

using namespace std;

// g++ -std=c++0x kmer2coordinates.cpp -o test2
double coordinates[20][20];

const string AA20 = "ARNDCEQGHILKMFPSTWYV";
const string NODELABEL = " ARNDCEQGHILKMFPSTWYV";
                   /*A   B  C  D  E   F  G  H  I   J   K   L   M  N   O   P  Q  R   S   T   U   V   W   X   Y   Z */
const int base[] = { 0, -1, 4, 3, 6, 13, 7, 8, 9, -1, 11, 10, 12, 2, -1, 14, 5, 1, 15, 16, -1, 19, 17, -1, 18, -1 };

int main(int argc, const char **argv) {
  ////////////////////////////////////
  coordinates[0][0] = -0.876280;coordinates[0][1] = 3.598596;coordinates[0][2] = 2.554616;coordinates[0][3] = -0.729216;coordinates[0][4] = 0.698828;coordinates[0][5] = 1.221507;coordinates[0][6] = -2.765205;coordinates[0][7] = -3.163091;
  coordinates[1][0] = -4.111404;coordinates[1][1] = -1.936791;coordinates[1][2] = -2.682295;coordinates[1][3] = 0.942498;coordinates[1][4] = 6.924314;coordinates[1][5] = -1.195785;coordinates[1][6] = -1.639269;coordinates[1][7] = 0.615381;
  coordinates[2][0] = -7.471612;coordinates[2][1] = -2.468058;coordinates[2][2] = 0.932738;coordinates[2][3] = -4.488355;coordinates[2][4] = 0.553080;coordinates[2][5] = -3.081577;coordinates[2][6] = 0.368010;coordinates[2][7] = 4.223792;
  coordinates[3][0] = -8.317871;coordinates[3][1] = -0.848602;coordinates[3][2] = 1.752372;coordinates[3][3] = -1.407818;coordinates[3][4] = -4.874022;coordinates[3][5] = -1.493568;coordinates[3][6] = 5.256411;coordinates[3][7] = -2.561758;
  coordinates[4][0] = 5.421664;coordinates[4][1] = 11.791877;coordinates[4][2] = 2.675596;coordinates[4][3] = -5.622478;coordinates[4][4] = 4.322457;coordinates[4][5] = 3.946839;coordinates[4][6] = 2.229597;coordinates[4][7] = -1.901479;
  coordinates[5][0] = -3.771796;coordinates[5][1] = -2.525005;coordinates[5][2] = -1.567736;coordinates[5][3] = 2.619391;coordinates[5][4] = 2.781873;coordinates[5][5] = 0.952486;coordinates[5][6] = 3.947072;coordinates[5][7] = -0.954304;
  coordinates[6][0] = -6.585010;coordinates[6][1] = -2.752755;coordinates[6][2] = -1.649014;coordinates[6][3] = 1.605597;coordinates[6][4] = -1.833933;coordinates[6][5] = -0.730211;coordinates[6][6] = 2.313328;coordinates[6][7] = -3.239486;
  coordinates[7][0] = -3.978253;coordinates[7][1] = -1.155062;coordinates[7][2] = 9.994796;coordinates[7][3] = -0.195264;coordinates[7][4] = -1.110059;coordinates[7][5] = -2.860194;coordinates[7][6] = -4.952672;coordinates[7][7] = -1.495210;
  coordinates[8][0] = -2.630176;coordinates[8][1] = -8.283034;coordinates[8][2] = -4.773107;coordinates[8][3] = -6.479084;coordinates[8][4] = 0.070359;coordinates[8][5] = 4.318067;coordinates[8][6] = -1.847373;coordinates[8][7] = -0.086451;
  coordinates[9][0] = 4.548022;coordinates[9][1] = 5.189698;coordinates[9][2] = -3.999001;coordinates[9][3] = -0.186966;coordinates[9][4] = -3.275059;coordinates[9][5] = -1.882387;coordinates[9][6] = -0.627095;coordinates[9][7] = 0.049364;
  coordinates[10][0] = 5.341899;coordinates[10][1] = 4.436639;coordinates[10][2] = -3.552811;coordinates[10][3] = 1.250614;coordinates[10][4] = 0.266899;coordinates[10][5] = -2.609335;coordinates[10][6] = -0.694939;coordinates[10][7] = 0.812004;
  coordinates[11][0] = -5.742562;coordinates[11][1] = -1.207887;coordinates[11][2] = -2.587323;coordinates[11][3] = 2.866228;coordinates[11][4] = 4.169821;coordinates[11][5] = -1.991698;coordinates[11][6] = -1.941954;coordinates[11][7] = -0.747156;
  coordinates[12][0] = 4.241223;coordinates[12][1] = 2.474317;coordinates[12][2] = -2.658336;coordinates[12][3] = 2.946054;coordinates[12][4] = 2.011534;coordinates[12][5] = -3.254331;coordinates[12][6] = 1.266004;coordinates[12][7] = -0.186966;
  coordinates[13][0] = 9.340442;coordinates[13][1] = -3.359172;coordinates[13][2] = -0.635377;coordinates[13][3] = -2.878570;coordinates[13][4] = -3.255191;coordinates[13][5] = -2.200202;coordinates[13][6] = -1.104637;coordinates[13][7] = -0.062654;
  coordinates[14][0] = -6.150933;coordinates[14][1] = 3.182318;coordinates[14][2] = 0.122393;coordinates[14][3] = 7.788554;coordinates[14][4] = -3.094076;coordinates[14][5] = 6.831600;coordinates[14][6] = -1.992627;coordinates[14][7] = 1.807240;
  coordinates[15][0] = -2.523437;coordinates[15][1] = 1.824168;coordinates[15][2] = 3.256463;coordinates[15][3] = -2.386830;coordinates[15][4] = 0.439791;coordinates[15][5] = 1.024198;coordinates[15][6] = 0.486894;coordinates[15][7] = 1.190316;
  coordinates[16][0] = -0.823028;coordinates[16][1] = 3.115233;coordinates[16][2] = 2.075337;coordinates[16][3] = -0.585875;coordinates[16][4] = -1.471153;coordinates[16][5] = 0.518398;coordinates[16][6] = 1.846290;coordinates[16][7] = 6.269577;
  coordinates[17][0] = 13.592409;coordinates[17][1] = -8.961858;coordinates[17][2] = 6.548108;coordinates[17][3] = 4.623650;coordinates[17][4] = 2.128797;coordinates[17][5] = 0.808588;coordinates[17][6] = 2.631353;coordinates[17][7] = 0.521535;
  coordinates[18][0] = 7.173223;coordinates[18][1] = -6.765800;coordinates[18][2] = -2.811202;coordinates[18][3] = -1.654989;coordinates[18][4] = -1.878135;coordinates[18][5] = 3.104673;coordinates[18][6] = -1.272146;coordinates[18][7] = -0.635970;
  coordinates[19][0] = 3.323480;coordinates[19][1] = 4.651177;coordinates[19][2] = -2.996218;coordinates[19][3] = 1.972858;coordinates[19][4] = -3.576126;coordinates[19][5] = -1.427066;coordinates[19][6] = -1.507041;coordinates[19][7] = -0.454682;
  ////////////////////////////////////


  //string kmers_file = "kmers.gen.txt";
  ofstream fout("coordinates.amino.acids.txt");
  ifstream fin(argv[1]);
  cout << argv[1] << endl;
  string kmer;
  size_t num_of_kmers = 0;
  while (getline(fin, kmer)) {
    if (kmer.size() != 10) {
      cout << "length is not correct~";
      return 0;
    }
    num_of_kmers++;
    for (size_t i = 0; i < kmer.size(); ++i) {
      int AA = base[kmer[i] - 'A'];
      if (AA < 0) {
        cout << "AA is incorrect" << endl;
        return 0;
      }
      if (i == 0) {
        fout << coordinates[AA][0];
      }
      for (size_t j = 0; j <= 7; ++j) {
        if (i == 0 && j == 0)
          continue;
        fout << "\t" << coordinates[AA][j];
      }
    }
    fout << endl;
  }
  fin.close();
  fout.close();
  cout << "num_of_kmers= " << num_of_kmers << endl;
  return 0;
}
