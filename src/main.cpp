#include <iostream>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <cstring>
#include <vector>
#include <armadillo>
#include <ctime>

#include "Reader.hpp"
#include "Solver.hpp"
#include "Config.hpp"

// To Do
// Minimise distance between close points
// iterative method

// Given W and ICCD find the best match (Important)
// Given W and ICCD try to fit during the search

// for W => 10% shifting metal, 2/3% for cristals so 5->10%.

// we can define a domain specific language for NMF
// define a language for the constraints, for the convulution
// maybe the shifting is a muliplication convolutive or by permutation?
// rewrite the constraints and the core of the solver using a solverAPI point of view

// Two points can only differ on the number of phaseby 1.

/*
Update of H
// Can be improved again because, the difference between 
// W.rows(0, L-p1-1).t()*O.rows(p1, L-1)
// and
// W.rows(0, L-p2-1).t()*O.rows(p2, L-1)
// is 
// W.rows(L-p2-1, L-p1-1).t()*O....
// same for the other

*/

using namespace std;
using namespace arma;

 
int main(int argc, char **argv){
   	
    Config * cfg;
    if (argc > 1){
    	cfg = new Config(argv[1]);
    }else{
    	cfg = new Config();
    }
    cfg->print();

   	Instance ir(cfg->datafile);
    Solver s(ir, *cfg);
    s.solve();


    return 0;
}



























