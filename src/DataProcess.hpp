#pragma once

#include <iostream>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <cstring>
#include <vector>
#include <armadillo>
#include "spline.h"

using namespace arma;
using namespace std;


// a class recording proportion, shift, etc.
class PhaseMixture {
public:
    bool isInvolved;
    double proportion;
    double shift;
};


class DataProcess
{
public:
	static mat re_sample_log(mat &A, double * qvalues, double * q_values_log);
	static mat resample(mat &A, double * current_values, double * new_values);

	static void extract_phase_mixture(vector<vector<PhaseMixture> > &phase_mixtures, vector<mat> &H, double *Qlog, int N, int K, int M, double epsilon);

	static void save_mat_csv(mat &A, std::ostream &os);
	static void save_array_csv(double * A, int length, std::ostream &os);
	static void save_mat_equal(mat &A, string name, std::ostream &os);
	static void save_array_equal(double * A, int length, string name, std::ostream &os);

	static void save_data_A(mat &A, double * qvalues);
	/* data */
};
