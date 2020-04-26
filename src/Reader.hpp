#pragma once

#include <iostream>
#include <cstdio>
#include <cstring>
#include <vector>
#include <string>
#include <fstream>

using namespace std;



class Instance {
public:
	int N;							// Number of sample points
	int Qsize;						// Number of wavelength
	double * Q;						// Wavelenght
	double * I;						// XRD Pattern
	std::vector<string> elements; 	// Mixed elements
	double * comp_by_elements;		//Compositions
    string filename;
	Instance(string filename_):filename(filename_){
		ifstream i_file(filename, ios::in);

	    if (!i_file.is_open()) {    //check if the file was opened correctly
	        cerr << "Problem opening the input file!\n";
	        exit(3);
	    }

		i_file >> N;
		i_file >> Qsize;
		for (int i = 0; i < 3; ++i){
			string elem;
			i_file >> elem;
			elements.push_back(elem);
		}

		Q = new double[Qsize];
		I = new double[N * Qsize];
		comp_by_elements = new double[3*N];

		for (int i = 0; i < Qsize; ++i){
			i_file >> Q[i];
		}
		for (int i = 0; i < N; ++i){
			i_file >> comp_by_elements[3*i + 0];
			i_file >> comp_by_elements[3*i + 1];
			i_file >> comp_by_elements[3*i + 2];
		}

		for (int i = 0; i < N; ++i){
			for (int q = 0; q < Qsize; ++q){
				i_file >> I[i * Qsize + q];
                
                if (I[i * Qsize + q] < 0.0)
                    I[i * Qsize + q] = 0.0;
			}
		}

		// check_read();
	}

	void check_read(){
		cout << N << "\t sample points\n";
		cout << Qsize << "\t wavelengths by point\n";
		cout << "Elements : " << elements[0] << " " << elements[1] << " " << elements[2] << "\n";
		cout << I[0] << "\t first sample point intensity\n";
		cout << Q[0] << "\t first Wavelenght\n";
		cout << comp_by_elements[0] << "\t first composition of first element\n";
		
	}

};



