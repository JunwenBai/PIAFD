#include "Config.hpp"


Config::Config(string filename_):filename(filename_){
	ifstream i_file(filename, ios::in);

	if (!i_file.is_open()) {    //check if the file was opened correctly
	    cerr << "Problem opening the configuration file!\n";
	    exit(3);
    }

    i_file >> datafile;
    i_file >> outfile;
	i_file >> K;
	i_file >> M;
	i_file >> number_of_iteration_step;
	i_file >> sparsity_coefficient;
	i_file >> time_out_in_second;
	i_file >> stop_on_decrease_obj;
	i_file >> gibbs_relax;
	i_file >> kmeans_init;
	i_file >> gibbs_frequency;
	i_file >> gradientMethod;
	i_file >> applyLog;
    i_file >> alloy_relax;
    i_file >> alloy_frequency;
    i_file >> use_dist;
    i_file >> edgefile;
    i_file >> dist_cutoff;
    i_file >> use_KL;
    i_file >> sticksfile;
    i_file >> match_shift;
    i_file >> match_sigma;
    i_file >> seed;

}



void Config::print(){
	cout << K << "\t K \n";
	cout << M << "\t M \n";
	cout << number_of_iteration_step << "\t Number of iteration Max \n";
	cout << sparsity_coefficient << "\t sparsity coefficient \n";
	cout << time_out_in_second << "\t timeout in second \n";
	cout << stop_on_decrease_obj << "\t stop when objective decrease? \n";
	cout << gibbs_relax << "\t is gibbs relaxed? \n";
    cout << kmeans_init << "\t initialisation (0:random Pts, 1 Kmean ++, 2 random values) \n";
    cout << gibbs_frequency << "\t gibbs enforcing frequency \n";
    cout << gradientMethod << "\t gradient update (0 Mult, 1 L1, 2 L2)\n";
    cout << applyLog << "\t Apply log transformation\n";
    cout << alloy_relax << "\t is alloy relaxed?\n";
    cout << alloy_frequency << "\t alloy frequency\n";
    cout << use_dist << "\t whether to use distance to build graph\n";
    cout << edgefile << "\t edgefile if not using distance\n";
    cout << dist_cutoff << "\t distance cutoff\n";
    cout << use_KL << "\t whether to use KL\n";
	cout << datafile << "\t data file \n";
    cout << outfile << "\t output file \n";
    cout << sticksfile << "\t stick file\n";
    cout << match_shift << "\t matching shift\n";
    cout << match_sigma << "\t matching sigma\n";
    cout << seed << "\t seed\n";
}
