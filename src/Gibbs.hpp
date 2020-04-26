#pragma once

#include <iostream>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <cstring>
#include <vector>
#include <armadillo>

#include <queue>

// #include "Projection.hpp"
#include "KMax.hpp"

using namespace arma;
using namespace std;

class Gibbs
{
public:

	static void relax(std::vector<mat> &H, int Q, int K, int N, int M, double epsilon, mat &Comp);
	static void enforce(std::vector<mat> &H, int Q, int K, int N, int M, double epsilon, mat &Comp, double prop, bool hard = false);
};

