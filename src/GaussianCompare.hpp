#ifndef Gauss_H
#define Gauss_H

#include <iostream>
#include <cstdio>
#include <time.h>
#include <cmath>
#include <cstring>
#include <vector>
#include <algorithm>
#include <string>
#include <armadillo>
#include <fstream>
#include "assert.h"
#include "spline.h"
#include <set>
#include <map>
#include <ctime>

using namespace arma;
using namespace std;

class IcsdMatching {
public:
	double loss, sigma, height, shift, ratio;
};

class GaussianCompare {

private:

double pi, convergenceRate, epsc;
int time_limit, iter_limit;
IcsdMatching matchI[300];
double MatchShiftTol, MatchSigmaTol, shiftCenter;

public:

GaussianCompare();

GaussianCompare(double shift, double sigma);

int min(int x1, int x2, int x3);

double readDouble(string s);

mat InsertRows(mat W, int rowcount, double val, bool top = true);

void readDoubleArray(vector<double> &vec, string buf);

double gaussianFunc(double sigma, double shift, double loc, double height, double x);

double L2loss(vector<double> x1, vector<double> x2);

double L1loss(vector<double> x1, vector<double> x2);

vector<double> Rec(vector<double> &Qicsd, vector<double> &Picsd, vector<double> &Q, double sigma, double height, double shift, vector<double> &vec);

vector<double> L2update(vector<double> &xrd, vector<double> &xrd_rec, double &sigma, double &height, double &shift, vector<double> Qicsd, vector<double> Picsd, vector<double> Q);

double optimize(vector<double> Qicsd, vector<double> Picsd, vector<double> xrd, vector<double> Q, vector<double> &xrd_rec, int index, int kindex, double &ratio, double &sigma, double &height, double &shift);

static bool myComp(pair<pair<int, double>, double> a, pair<pair<int, double>, double> b);

static bool myComp2(pair<pair<int, double>, double> a, pair<pair<int, double>, double> b);

vector<double> convertNorm(vector<double> xrd, vector<double> Qlogsteps, vector<double> Q);

void matchIcsd(int i, vector<double> Qicsd, vector<double> Picsd, vector<double> xrd, vector<double> Q, int k);

void compare(ostream &os, char* sticksfile, mat W, vector<mat> H, int K, vector<double> Q, vector<double> Qlogsteps);

};

#endif
