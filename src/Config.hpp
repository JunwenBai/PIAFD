#pragma once

#include <iostream>
#include <cstdio>
#include <cstring>
#include <vector>
#include <string>
#include <fstream>
using namespace std;


class Config
{
public:
	Config(string filename_ = "../data/config.config");
	// ~Config();

	void print();


	string filename;
	string outfile;
	string datafile;
	int K;
	int M;
	int number_of_iteration_step;
	double sparsity_coefficient;
	int time_out_in_second;
	int stop_on_decrease_obj;
	double gibbs_relax;
	int kmeans_init;
	int gibbs_frequency;
	int gradientMethod;
	int applyLog;
    double alloy_relax;
    int alloy_frequency;
    int use_dist;
    string edgefile;
    double dist_cutoff;
    int use_KL;
    double match_shift;
    double match_sigma;
    string sticksfile;
    int seed;

};
