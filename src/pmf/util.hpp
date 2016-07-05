#ifndef UTIL_H
#define UTIL_H

#pragma once

#include <stdlib.h>
#include <string.h>

#include <vector>
#include <string>
#include <fstream>
#include <iostream>

using namespace std;

const char pmf_version[] = "1.0";

double coordinates[20][8] = {
		{-0.876280,  3.598596,  2.554616, -0.729216,  0.698828,  1.221507, -2.765205, -3.163091},
		{-4.111404, -1.936791, -2.682295,  0.942498,  6.924314, -1.195785, -1.639269,  0.615381},
		{-7.471612, -2.468058,  0.932738, -4.488355,  0.553080, -3.081577,  0.368010,  4.223792},
		{-8.317871, -0.848602,  1.752372, -1.407818, -4.874022, -1.493568,  5.256411, -2.561758},
		{ 5.421664, 11.791877,  2.675596, -5.622478,  4.322457,  3.946839,  2.229597, -1.901479},
		{-3.771796, -2.525005, -1.567736,  2.619391,  2.781873,  0.952486,  3.947072, -0.954304},
		{-6.585010, -2.752755, -1.649014,  1.605597, -1.833933, -0.730211,  2.313328, -3.239486},
		{-3.978253, -1.155062,  9.994796, -0.195264, -1.110059, -2.860194, -4.952672, -1.495210},
		{-2.630176, -8.283034, -4.773107, -6.479084,  0.070359,  4.318067, -1.847373, -0.086451},
		{ 4.548022,  5.189698, -3.999001, -0.186966, -3.275059, -1.882387, -0.627095,  0.049364},
		{ 5.341899,  4.436639, -3.552811,  1.250614,  0.266899, -2.609335, -0.694939,  0.812004},
		{-5.742562, -1.207887, -2.587323,  2.866228,  4.169821, -1.991698, -1.941954, -0.747156},
		{ 4.241223,  2.474317, -2.658336,  2.946054,  2.011534, -3.254331,  1.266004, -0.186966},
		{ 9.340442, -3.359172, -0.635377, -2.878570, -3.255191, -2.200202, -1.104637, -0.062654},
		{-6.150933,  3.182318,  0.122393,  7.788554, -3.094076,  6.831600, -1.992627,  1.807240},
		{-2.523437,  1.824168,  3.256463, -2.386830,  0.439791,  1.024198,  0.486894,  1.190316},
		{-0.823028,  3.115233,  2.075337, -0.585875, -1.471153,  0.518398,  1.846290,  6.269577},
		{13.592409, -8.961858,  6.548108,  4.623650,  2.128797,  0.808588,  2.631353,  0.521535},
		{ 7.173223, -6.765800, -2.811202, -1.654989, -1.878135,  3.104673, -1.272146, -0.635970},
		{ 3.323480,  4.651177, -2.996218,  1.972858, -3.576126, -1.427066, -1.507041, -0.454682}
};

const string AA20 = "ARNDCEQGHILKMFPSTWYV";
const string NODELABEL = " ARNDCEQGHILKMFPSTWYV";
                   /*A   B  C  D  E   F  G  H  I   J   K   L   M  N   O   P  Q  R   S   T   U   V   W   X   Y   Z */
const int base[] = { 0, -1, 4, 3, 6, 13, 7, 8, 9, -1, 11, 10, 12, 2, -1, 14, 5, 1, 15, 16, -1, 19, 17, -1, 18, -1 };

const size_t AACoordinateSize = 8;

#endif

