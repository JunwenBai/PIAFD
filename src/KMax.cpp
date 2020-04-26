#include "KMax.hpp"



void KMax::relax(std::vector<mat> &H, int n,int K, double epsilon){
    double maxV = 0;
    double s = 0;
	for (int i = 0; i < colSize; ++i){
        v[i] = 0;
    } 
    for (int k = 0; k < colSize; ++k){
        for (int shift = 0; shift < H.size(); shift++) { 
            v[k] += H[shift](k,n);
        }
    }
    
    for (int k = 0; k < colSize; ++k){
        s += (v[k] > 0.02)? 1:0;
        maxV = max(maxV, v[k]);
    }
    if (s < 4){
        return;
    }
    for (int k = 0; k < colSize; ++k){
        v[k] /= maxV+epsilon;
    }
    projection_L1Ball(v, vp, colSize, K);

    for (int k = 0; k < colSize; ++k){
        for (int shift = 0; shift < H.size(); shift++) { 
            H[shift](k,n) = (vp[k]) * H[shift](k,n) + ((static_cast <float> (rand()) / static_cast <float> (RAND_MAX)) / 1000000) + epsilon;
//            H[shift](k,n) = (vp[k]) * H[shift](k,n);
        }
    } 
}

void KMax::enforce(std::vector<mat> &H, int n, int K, double epsilon, double prop, bool hard = false){
    for (int i = 0; i < colSize; ++i){
        v[i] = 0;
        vp[i] = 1-prop;
    } 
    for (int k = 0; k < colSize; ++k){
        for (int shift = 0; shift < H.size(); shift++) { 
            v[k] += H[shift](k,n);
        }
    }
    
    priority_queue<pair<double, int>> q;
    for (int k = 0; k < colSize; ++k){
		q.push(std::pair<double, int>(v[k], k));
	}
	
	for (int i = 0; i < K; ++i) {
		vp[q.top().second] = 1.;
		q.pop();
	}

    for (int k = 0; k < colSize; ++k){
        for (int shift = 0; shift < H.size(); shift++) { 
            if (!hard)
                H[shift](k,n) = (vp[k]) * H[shift](k,n) + (1-vp[k]) * ((static_cast <float> (rand()) / static_cast <float> (RAND_MAX)) / 1000000) + epsilon;
            else H[shift](k,n) = (vp[k]) * H[shift](k,n);
        }
    } 
}
