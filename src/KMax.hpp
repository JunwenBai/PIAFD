#pragma once

#include <iostream>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <cstring>
#include <vector>
#include <armadillo>

#include <queue>
#include <vector>

#include "Projection.hpp"

using namespace arma;
using namespace std;

class KMax
{
public:
	double * v;
    double * vp;
    int colSize;

    KMax(int colSize_):
    	v(new double[colSize_]),
    	vp(new double[colSize_]),
    	colSize(colSize_)
    	{}

	void relax(std::vector<mat> &H, int n, int K, double epsilon);		// L1
	void enforce(std::vector<mat> &H, int n, int K, double epsilon, double prop, bool hard = false);	// L0
	
	/* data */
};
