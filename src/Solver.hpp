#pragma once

#include <iostream>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <cstring>
#include <vector>
#include <armadillo>
#include <string>
#include <ctime>

#include "Reader.hpp"
#include "DataProcess.hpp"
#include "Config.hpp"
#include "Gibbs.hpp"
#include "ChronoP.hpp"
#include "Alloy.hpp"
#include "Connectivity.hpp"
#include "GaussianCompare.hpp"

using namespace arma;
using namespace std;

class Solver
{
public:
	Solver(Instance &I_, Config &cfg_);
	// ~Solver();
	// Utils functions
	void normalizeA();
	void normalizeW();
	void normalizeH();
	void normalizeHAndUpdateW();

	// main update
	void gradientUpdate();
	void gradientUpdateL2();
	void gradientUpdateL1();
	void gradientUpdateL21();

	// Constraints
	void enforceGibbs(double alpha, bool hard = false);
	void enforceAlloy(double alpha, bool hard = false);
	void enforceConnectivity(double beta=0.0);

	// Data processing
	void reconstruction();
	double lossFunction();
	double lossFunctionL2();
	double lossFunctionL1();
	double lossFunctionL21();

	void init_from_random();
	void init_from_data();
	void init_from_data_pp();
	// Run the solver on the current instance
	void solve();


	// postProcess
	void postProcess();

	// Save to csv
	void save_to_file_csv(string filename);
	void save_to_file_phase_mapper(string filename);

	/* data */
	mat A;
	mat W;
	std::vector<mat> H;
	Instance &I;
	mat R;
	mat Aux1;

	mat sparsity;
	mat sparsityW;
	mat Comp; 
	Config &cfg;
	double *Qlog;

	double sparsityCoefficient;

	int N;			// number of sample points
	int Q;			// number of wavelength
	int K;			// number of phase in W
	int M;			// number of shifted version
    time_t init_time; // initial time

	const static double epsilon; // threshold
	double alpha; 	// gradient step
	double beta; 	// gradient step modifier
    double last_cost;


    vector<vector<int> > neighbors;// neighbors
    Mat<int> connect_indicators; // indicate whether a phase should be activated at some sample point


	// Data for printing for W differences(To remove)
	mat W_init; 	// initialized in solve().

	void print_diff_W(){ // to call at the end
		std::filebuf fb;
	    fb.open ("../data/out_w_diff.txt",std::ios::out);
	    std::ostream os(&fb);

	    mat W_non_log = DataProcess::resample(W, Qlog, I.Q); 
    	mat A_non_log = DataProcess::resample(A, Qlog, I.Q);
    	mat W_init_non_log = DataProcess::resample(W_init, Qlog, I.Q); 
    	os << "K=" <<  K << "\n";
		DataProcess::save_array_equal(I.Q, Q, "Q", os);
	    os << "\n";  
	    
	    double * B = new double[max(Q,N)];
	    for (int k = 0; k < K; ++k){
	        string strB("B"+to_string(k+1));
	        for (int q = 0; q < Q; ++q){
	            B[q] = W_non_log(q,k);
	        }
	        DataProcess::save_array_equal(B, Q, strB, os);
	        os << "\n";  
	    }
	    os << "\n";  
	    for (int k = 0; k < K; ++k){
	        string strB("BI"+to_string(k+1));
	        for (int q = 0; q < Q; ++q){
	            B[q] = W_init_non_log(q, k);
	        }
	        DataProcess::save_array_equal(B, Q, strB, os);
	        os << "\n";  
	    }
	    os << "\n";  
	    
	}

};







