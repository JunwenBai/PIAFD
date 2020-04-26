#include "Connectivity.hpp"

void Connectivity::visitPhase(int sample, double &sumv, vector<int> &visit, int cnt, vector<mat> &H, int k, double epsilon, vector<vector<int> > &neighbors) {

    int N = H[0].n_cols;
    int M = H.size();

    double sample_sum = 0.0;
    for (int m = 0; m < M; m++)
        sample_sum += H[m](k,sample);
    if (sample_sum < epsilon) return;

    visit[sample] = cnt;
    sumv += sample_sum;

    for (vector<int>::iterator it = neighbors[sample].begin(); it != neighbors[sample].end(); it++) {
        if (visit[*it] == -1)
            Connectivity::visitPhase(*it, sumv, visit, cnt, H, k, epsilon, neighbors);
    }
}

void Connectivity::connectPhase(vector<mat> &H, int Q, int N, int K, int M, Mat<int> &connect_indicators, double epsilon, vector<vector<int> > &neighbors, double beta) {
//    fprintf(stderr, "in connectPhase!\n");
    connect_indicators = zeros<Mat<int> >(K, N);
    vector<int> mark;
    for (int i = 0; i < N; i++) mark.push_back(-1);

    int combcnt = 0;

    for (int k = 0; k < K; k++) {
        vector<int> visit;
        for (int n = 0; n < N; n++) visit.push_back(-1);

        vector<double> weights;
        weights.clear();
        int cnt = 0;
        double max_sumv = 0.0;
        int chosencnt = 0;

        for (int n = 0; n < N; n++) {
            double sum = 0.0;
            for (int m = 0; m < M; m++) sum += H[m](k,n);
            if (sum < epsilon) continue;
            else connect_indicators(k,n) = true;

            if (visit[n] == -1) {
                double sumv = 0.0;
                visitPhase(n, sumv, visit, cnt, H, k, epsilon, neighbors);
                if (sumv > max_sumv) {
                    max_sumv = sumv;
                    chosencnt = cnt;
                }
                cnt++;
            }
        }

        for (int n = 0; n < N; n++) {
            double sum = 0.0;
            for (int m = 0; m < M; m++) sum += H[m](k,n);
            if (sum < epsilon) continue;

            if (visit[n] != chosencnt) {
                for (int m = 0; m < M; m++) H[m](k,n) *= beta;
                connect_indicators(k,n) = false;
            }
        }
        
    }
}

void Connectivity::visit(int sample, double &sumv, vector<int> &visit, int cnt, vector<mat> &H, int k, map<set<int>, int> &combmap, map<int, set<int> > &combdict, vector<int> &mark, set<int> &losers, vector<vector<int> > &neighbors, double epsilon) {
    int N = H[0].n_cols;
    int M = H.size();
    double sample_sum = 0.0;
    for (int m = 0; m < M; m++)
        for (set<int>::iterator it = combdict[k].begin(); it != combdict[k].end(); it++)
            sample_sum += H[m](*it, sample);
    if (sample_sum < epsilon) return;

    visit[sample] = cnt;
    sumv += sample_sum;

    for (vector<int>::iterator it = neighbors[sample].begin(); it != neighbors[sample].end(); it++)
        if (mark[*it] == k && visit[*it] == -1)
            Connectivity::visit(*it, sumv, visit, cnt, H, k, combmap, combdict, mark, losers, neighbors, epsilon);
}

void Connectivity::connect(vector<mat> &H, int Q, int N, int K, int M, Mat<int> &connect_indicators, double epsilon, vector<vector<int> > &neighbors, double beta) {
//    fprintf(stderr, "in connect!\n");
    connect_indicators = zeros<Mat<int> >(K,N);
    map<set<int>, int> combmap;
    map<int, set<int> > combdict;
    vector<int> mark;
    for (int i = 0; i < N; i++) mark.push_back(-1);

    int combcnt = 0;
    set<int> losers;

    for (int n = 0; n < N; n++) {
        set<int> comb;
        for (int k = 0; k < K; k++) {
            double sum = 0.0;
            for (int m = 0; m < M; m++) {
                sum += H[m](k,n);
            }
            if (sum > epsilon) {
                comb.insert(k);
            }
        }
        if (comb.size() == 0) {
            losers.insert(n);
            continue;
        }
        if (combmap.find(comb) == combmap.end()) {
            combmap[comb] = combcnt;
            combdict[combcnt] = comb;
            combcnt++;
        }
        mark[n] = combmap[comb];
    }

    for (int i = 0; i < combcnt; i++) {
        vector<int> visit;
        for (int n = 0; n < N; n++) visit.push_back(-1);
        vector<double> weights;
        weights.clear();
        int cnt = 0;
        double max_sumv = 0.0;
        int chosencnt = 0;

        for (int n = 0; n < N; n++) {
            if (visit[n] == -1 && mark[n] == i) {
                double sumv = 0.0;
                Connectivity::visit(n, sumv, visit, cnt, H, i, combmap, combdict, mark, losers, neighbors, epsilon);
                if (sumv > max_sumv) {
                    max_sumv = sumv;
                    chosencnt = cnt;
                }
                cnt++;
            }
        }

        for (int n = 0; n < N; n++) {
            if (mark[n] == i) {
                if (visit[n] != chosencnt) {
                    if (beta > 0.0) {
                        for (int k = 0; k < K; k++)
                            for (int m = 0; m < M; m++)
                                H[m](k,n) *= beta;
                    } else losers.insert(n);
                }
                else {
                    for (set<int>::iterator it = combdict[i].begin(); it != combdict[i].end(); it++)
                        connect_indicators(*it, n) = true;
                }
            }
        }
    }

    if (beta > 0) return;

    map<int, int> validnb;

    for (set<int>::iterator loser = losers.begin(); loser != losers.end(); loser++) {
        int validn = 0;
        for (vector<int>::iterator it = neighbors[*loser].begin(); it != neighbors[*loser].end(); it++) {
            if (losers.find(*it) == losers.end()) validn += 1;
        }
        validnb[*loser] = validn;
    }

    while (losers.size()) {
        int largest = 0, chosenone = -1;
        for (set<int>::iterator loser = losers.begin(); loser != losers.end(); loser++) {
            if (validnb[*loser] > largest) {
                largest = validnb[*loser];
                chosenone = *loser;
            }
        }

        map<int, int> mapcnt;
        for (vector<int>::iterator it = neighbors[chosenone].begin(); it != neighbors[chosenone].end(); it++) {
            if (losers.find(*it) == losers.end()) {
                int type = mark[*it];
                if (mapcnt.find(type) == mapcnt.end()) mapcnt[type] = 1;
                else mapcnt[type] += 1;
            } else {
                validnb[*it] += 1;
            }
        }

        int maxtype = -1, maxnum = 0;
        for (map<int, int>::iterator it = mapcnt.begin(); it != mapcnt.end(); it++) {
            if (it->second > maxnum) {
                maxnum = it->second;
                maxtype = it->first;
            }
        }

        mark[chosenone] = maxtype;
        for (int m = 0; m < M; m++)
            for (int k = 0; k < K; k++)
                H[m](k, chosenone) = 0.0;
        for (set<int>::iterator it = combdict[maxtype].begin(); it != combdict[maxtype].end(); it++) {
            connect_indicators(*it, chosenone) = true;
            for (int m = 0; m < M; m++) {
                double h = 0.0;
                int h_cnt = 0;
                for (vector<int>::iterator x = neighbors[chosenone].begin(); x != neighbors[chosenone].end(); x++) {
                    if (losers.find(*x) != losers.end()) continue;
                    double sum = 0.0;
                    for (int t = 0; t < M; t++) sum += H[t](*it,*x);
                    if (sum < epsilon) continue;

                    h += H[m](*it, *x);
                    h_cnt++;
                }
                H[m](*it, chosenone) = h/h_cnt;
            }
        }
        losers.erase(chosenone);
    }

}


bool Connectivity::checkCC(vector<mat> &H, int N, int K, int M, double epsilon, vector<vector<int> > &neighbors) {
    for (int k = 0; k < K; k++) {
        vector<int> visit;
        for (int n = 0; n < N; n++) visit.push_back(-1);
        
        int visitCnt = 0;
        for (int n = 0; n < N; n++) {
            double sum = 0.0;
            for (int m = 0; m < M; m++) sum += H[m](k,n);
            if (sum < epsilon) continue;

            if (visit[n] == -1) {
                double sumv = 0.0;
                visitPhase(n, sumv, visit, visitCnt, H, k, epsilon, neighbors);
                visitCnt++;
                if (visitCnt > 1) return false;
            }
        }
    }
    return true;
}

