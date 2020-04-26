#pragma once

#include <iostream>
#include <cstdio>
#include <algorithm>
#include <cstring>
#include <vector>
#include <armadillo>
#include <set>
#include <queue>
#include <map>

using namespace arma;
using namespace std;

class Connectivity
{
    public:
    static void connectPhase(vector<mat> &H, int Q, int N, int K, int M, Mat<int> &connect_indicators, double epsilon, vector<vector<int> > &neighbors, double beta);
    static void connect(vector<mat> &H, int Q, int N, int K, int M, Mat<int> &connect_indicators, double epsilon, vector<vector<int> > &neighbors, double beta);
    
    
    static void visit(int sample, double &sumv, vector<int> &visit, int cnt, vector<mat> &H, int k, map<set<int>, int> &combmap, map<int, set<int> > &combdict, vector<int> &mark, set<int> &losers, vector<vector<int> > &neighbors, double epsilon);
    static void visitPhase(int sample, double &sumv, vector<int> &visit, int cnt, vector<mat> &H, int k, double epsilon, vector<vector<int> >& neighbors);
    
    static bool checkCC(vector<mat> &H, int N, int K, int M, double epsilon, vector<vector<int> > &neighbors);

};
