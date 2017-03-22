#ifndef UTIL_H
#define UTIL_H

#pragma once

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include <random>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>

using namespace std;

const char hclust_version[] = "1.0";

const double coordinates[20][8] = {
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
const double DISTANCE_SQUARE[20][20] = {
    {   0.000000,  131.470960,  179.985160,  179.395157,  177.850435,  128.434308,  132.072530,  115.360078,  251.297575,  115.515906,  115.537752,  115.036696,  114.999607,  207.300276,  177.623697,   38.736326,  115.931022,  456.623739,  220.666880,  91.067499},
    { 131.470960,    0.000000,  115.334503,  241.002799,  405.680566,   59.961544,  115.698566,  245.386609,  179.835349,  234.625411,  178.087999,   17.070569,  130.696453,  307.142501,  251.452151,  114.848206,  178.848693,  506.789294,  255.103396, 211.314262},
    { 179.985160,  115.334503,    0.000000,  115.415662,  478.497387,  131.327850,  115.380344,  178.310724,  172.246985,  280.445613,  277.991890,  115.373275,  252.201935,  324.311695,  308.035251,   78.810078,  115.339577,  636.869131,  325.474722, 269.593171},
    { 179.395157,  241.002799,  115.415662,    0.000000,  490.943059,  119.612840,   46.224610,  209.732027,  270.733165,  280.702504,  324.210158,  181.010193,  279.257332,  375.811834,  252.287351,  115.342025,  177.913484,  676.000004,  372.324402, 251.619708},
    { 177.850435,  405.680566,  478.497387,  490.943059,    0.000000,  390.613498,  488.227069,  466.604980,  562.185288,  222.115656,  215.485205,  447.414071,  251.146217,  373.107936,  489.332712,  209.499813,  252.241215,  638.151247,  446.321558, 252.572729},
    { 128.434308,   59.961544,  131.327850,  119.612840,  390.613498,    0.000000,   41.028738,  252.715483,  180.523050,  209.175986,  181.017770,   52.037361,  116.565476,  276.423875,  179.807518,   90.871568,  138.927633,  417.189102,  211.283472, 180.382485},
    { 132.072530,  115.698566,  115.380344,  46.224610,   488.227069,   41.028738,    0.000000,  209.067293,  177.711223,  218.629658,  231.097844,   67.521495,  178.918723,  301.078571,  179.559992,  108.654556,  178.638097,  554.300201,  251.717816, 180.817418},
    { 115.360078,  245.386609,  178.310724,  209.732027,  466.604980,  252.715483,  209.067293,    0.000000,  374.750951,  335.519425,  329.171042,  209.053454,  300.984926,  324.343217,  302.274331,  115.491859,  209.124346,  490.310662,  372.414324, 281.585328},
    { 251.297575,  179.835349,  172.246985,  270.733165,  562.185288,  180.523050,  177.711223,  374.750951,    0.000000,  324.375999,  336.744266,  208.921359,  327.054799,  251.723960,  391.320102,  201.452452,  285.673354,  552.055921,  131.431406, 323.887318},
    { 115.515906,  234.625411,  280.445613,  280.702504,  222.115656,  209.175986,  218.629658,  335.519425,  324.375999,    0.000000,   17.123203,  215.941527,   52.550588,  114.951162,  280.013150,  143.604574,  124.034822,  463.734443,  181.099508,   8.786247},
    { 115.537752,  178.087999,  277.991890,  324.210158,  215.485205,  181.017770,  231.097844,  329.171042,  336.744266,   17.123203,    0.000000,  177.868627,   17.039807,  115.829984,  293.012188,  143.058029,  123.847667,  387.294026,  177.521755,  23.382078},
    { 115.036696,   17.070569,  115.373275,  181.010193,  447.414071,   52.037361,   67.521495,  209.053454,  208.921359,  215.941527,  177.868627,    0.000000,  130.103256,  325.281754,  188.154815,  113.962429,  178.255099,  555.041896,  281.211198, 178.080248},
    { 114.999607,  130.696453,  252.201935,  279.257332,  251.146217,  116.565476,  178.918723,  300.984926,  327.054799,   52.550588,   17.039807,  130.103256,    0.000000,  132.535209,  282.069223,  132.887922,  129.325866,  324.690931,  177.379618,  48.964538},
    { 207.300276,  307.142501,  324.311695,  375.811834,  373.107936,  276.423875,  301.078571,  324.343217,  251.723960,  114.951162,  115.829984,  325.281754,  132.535209,    0.000000,  483.019976,  211.159417,  217.197160,  209.694118,   52.928050, 130.495800},
    { 177.623697,  251.452151,  308.035251,  252.287351,  489.332712,  179.807518,  179.559992,  302.274331,  391.320102,  280.013150,  293.012188,  188.154815,  282.069223,  483.019976,    0.000000,  181.106887,  179.476257,  675.176299,  396.141534, 209.260124},
    { 38.736326,   114.848206,   78.810078,  115.342025,  209.499813,   90.871568,  108.654556,  115.491859,  201.452452,  143.604574,  143.058029,  113.962429,  132.887922,  211.159417,  181.106887,    0.000000,   40.751128,  443.985807,  221.295615, 129.099397},
    { 115.931022,  178.848693,  115.339577,  177.913484,  252.241215,  138.927633,  178.638097,  209.124346,  285.673354,  124.034822,  123.847667,  178.255099,  129.325866,  217.197160,  179.476257,   40.751128,    0.000000,  447.505943,  250.861787, 116.496667},
    { 456.623739,  506.789294,  636.869131,  676.000004,  638.151247,  417.189102,  554.300201,  490.310662,  552.055921,  463.734443,  387.294026,  555.041896,  324.690931,  209.694118,  675.176299,  443.985807,  447.505943,    0.000000,  210.951244, 444.510082},
    { 220.666880,  255.103396,  325.474722,  372.324402,  446.321558,  211.283472,  251.717816,  372.414324,  131.431406,  181.099508,  177.521755,  281.211198,  177.379618,   52.928050,  396.141534,  221.295615,  250.861787,  210.951244,    0.000000, 181.871263},
    { 91.067499,   211.314262,  269.593171,  251.619708,  252.572729,  180.382485,  180.817418,  281.585328,  323.887318,    8.786247,   23.382078,  178.080248,   48.964538,  130.495800,  209.260124,  129.099397,  116.496667,  444.510082,  181.871263,   0.000000},
};
/*
const double coordinates[20][4] = {
{-0.339396, 3.307867, 2.509939, 0.503747},
{-3.872820, -2.652341, -5.089655, 0.927826},
{-5.620695, -3.134013, 1.497643, -5.241709},
{-9.324240, -0.422996, 2.021605, -1.750213},
{4.674816, 13.148086, 1.238484, -3.971265},
{-2.828587, -3.329859, -0.406519, 3.449013},
{-5.822487, -3.026515, 0.250674, 1.916624},
{-2.772717, -1.928463, 9.262868, -0.702885},
{-1.533098, -9.289944, -3.429837, -4.016064},
{3.432965, 4.193160, -3.822954, -1.143189},
{4.081610, 3.143537, -3.371281, 1.103740},
{-4.614399, -1.504121, -3.904108, 2.191024},
{3.294056, 1.898127, -2.985941, 3.301644},
{8.115797, -1.418148, -0.118660, -2.909285},
{-7.554148, 5.004731, 2.184156, 6.966242},
{-2.435777, 0.099015, 1.674380, -1.457551},
{-1.854970, 3.212650, 0.772068, -2.934513},
{15.193881, -6.097971, 5.301099, 3.897914},
{7.322214, -5.776442, -1.378126, -0.750291},
{2.457995, 4.573640, -2.205834, 0.619191}
};
*/
const string AA20 = "ARNDCEQGHILKMFPSTWYV";
const string NODELABEL = " ARNDCEQGHILKMFPSTWYV";
                   /*A   B  C  D  E   F  G  H  I   J   K   L   M  N   O   P  Q  R   S   T   U   V   W   X   Y   Z */
const int base[] = { 0, -1, 4, 3, 6, 13, 7, 8, 9, -1, 11, 10, 12, 2, -1, 14, 5, 1, 15, 16, -1, 19, 17, -1, 18, -1 };

const size_t AACoordinateSize = 8;
//const size_t AACoordinateSize = 4;

#endif
