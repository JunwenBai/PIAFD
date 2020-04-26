#include "Alloy.hpp"

void Alloy::relax(std::vector<mat> &H, int Q, int K, int N, int M, double epsilon, vector<vector<int> > &neighbors, double *Qlog, mat &Comp) {
    KMax kmax(K);
    double eps = 1e-6;
    for (int n = 0; n < N; n++) {
        bool is_shifting = false;
        for (int k = 0; k < K; k++) {
            double weighted_sum = 0.0, sum = 0.0;
            for (int m = 0; m < M; m++) {
                weighted_sum += (Qlog[m]/Qlog[0])*H[m](k,n);
                sum += H[m](k,n);
            }
            if (sum < 1e-4) continue;
            for (vector<int>::iterator it = neighbors[n].begin(); it != neighbors[n].end(); it++) {
                double weighted_tsum = 0.0, tsum = 0.0;
                for (int tm = 0; tm < M; tm++) {
                    weighted_tsum += (Qlog[tm]/Qlog[0])*H[tm](k,*it);
                    tsum += H[tm](k,*it);
                }
                if (tsum < 1e-4) continue;
                if (abs(weighted_sum/sum-weighted_tsum/tsum) > epsilon) {
                    is_shifting = true;
                    break;
                }
            }
        }
        if (is_shifting) {
            if (Comp(n, 0) < eps || Comp(n, 1) < eps || Comp(n, 2) < eps) kmax.relax(H, n, 1, epsilon);
            else kmax.relax(H, n, 2, epsilon);
        }
    }
}

void Alloy::enforce(std::vector<mat> &H, int Q, int K, int N, int M, double epsilon, vector<vector<int> > &neighbors, double *Qlog, mat &Comp, double prop, bool hard = false) {
    KMax kmax(K);
    double eps = 1e-6;
    for (int n = 0; n < N; n++) {
        bool is_shifting = false;
        for (int k = 0; k < K; k++) {
            double weighted_sum = 0.0, sum = 0.0;
            for (int m = 0; m < M; m++) {
                weighted_sum += (Qlog[m]/Qlog[0])*H[m](k,n);
                sum += H[m](k,n);
            }
            if (sum < 1e-4) continue;

            for (vector<int>::iterator it = neighbors[n].begin(); it != neighbors[n].end(); it++) {
                double weighted_tsum = 0.0, tsum = 0.0;
                for (int tm = 0; tm < M; tm++) {
                    weighted_tsum += (Qlog[tm]/Qlog[0])*H[tm](k,*it);
                    tsum += H[tm](k,*it);
                }
                if (tsum < 1e-4) continue;
                if (abs(weighted_sum/sum-weighted_tsum/tsum) > epsilon) {
                    is_shifting = true;
                    break;
                }
            }
        }
        if (is_shifting) {
            if (Comp(n, 0) < eps || Comp(n, 1) < eps || Comp(n, 2) < eps) kmax.enforce(H, n, 1, epsilon, prop, hard);
            else kmax.enforce(H, n, 2, epsilon, prop, hard);
        }
    }
}
