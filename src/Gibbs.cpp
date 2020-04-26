#include "Gibbs.hpp"

/* 
    Project a normalized vector onto the simplex, then update the values.
*/
void Gibbs::relax(std::vector<mat> &H, int Q, int K, int N, int M, double epsilon, mat &Comp){
    KMax kmax(K);
    double eps = 1e-6;
    for (int n = 0; n < N; ++n){
        if (Comp(n, 0) < eps || Comp(n, 1) < eps || Comp(n, 2) < eps) kmax.relax(H, n, 2, epsilon);
        else kmax.relax(H, n, 3, epsilon);
    }
}

/*
	keep only the 3 biggest values for each column of H.
*/
void Gibbs::enforce(std::vector<mat> &H, int Q, int K, int N, int M, double epsilon, mat &Comp, double prop, bool hard = false){
    KMax kmax(K);
    double eps = 1e-6;
    for (int n = 0; n < N; ++n){
        if (Comp(n, 0) < eps || Comp(n, 1) < eps || Comp(n, 2) < eps) kmax.enforce(H, n, 2, epsilon, prop, hard);
        else kmax.enforce(H, n, 3, epsilon, prop, hard);
    }

}


