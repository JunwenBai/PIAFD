#include "GaussianCompare.hpp"


GaussianCompare::GaussianCompare() {

    pi = atan(1)*4;
    convergenceRate = 1e-4;
    time_limit = 20000, iter_limit = 800;
    epsc = 1e-6;

}

GaussianCompare::GaussianCompare(double shift, double sigma) {

    pi = atan(1)*4;
    convergenceRate = 1e-4;
    time_limit = 20000, iter_limit = 800;
    epsc = 1e-6;
	MatchShiftTol = shift;
	MatchSigmaTol = sigma;

}


int GaussianCompare::min(int x1, int x2, int x3) {
    return std::min(std::min(x1,x2),x3);
}

double GaussianCompare::readDouble(string s) {
    char* str = (char*)s.c_str();
    double v;
    sscanf(str, "%lf", &v);
    return v;
}

mat GaussianCompare::InsertRows(mat W, int rowcount, double val, bool top = true) { // if top==true, row0 is inserted into the beginning. Otherwise, row0 is inserted into the end.
    rowvec row0(W.n_cols); // one row filled with zeros of length K
    row0.fill(val); // row0=[val,val,...,val]. Usually, val == 0.
    for (int r = 0; r < rowcount; r++) {
        if (top) {
            W.insert_rows(0, row0); // insert row0 to the beginning
        } else {
            W.insert_rows(W.n_rows, row0); // insert row0 to the end
        }
    }
    return W;
}

void GaussianCompare::readDoubleArray(vector<double> &vec, string buf) {
    vec.clear();
    size_t found = buf.find(',');
    while (found != string::npos) {
        vec.push_back(readDouble(buf.substr(0, found)));
        buf.erase(0, found+1);
        found = buf.find(',');
    }
    vec.push_back(readDouble(buf));
}

double GaussianCompare::gaussianFunc(double sigma, double shift, double loc, double height, double x) {
    return height*1.0/sqrt(2*sigma*sigma*pi)*exp(-(x-loc*shift)*(x-loc*shift)/(2*sigma*sigma));
}

double GaussianCompare::L2loss(vector<double> x1, vector<double> x2) {
    int L = x1.size();
    double loss = 0.0;
    for (int i = 0; i < L; i++) loss += (x1[i]-x2[i])*(x1[i]-x2[i]);
    return loss;
}

double GaussianCompare::L1loss(vector<double> x1, vector<double> x2) {
    int L = x1.size();
    double loss = 0.0;
    for (int i = 0; i < L; i++) loss += abs(x1[i]-x2[i]);
    return loss;
}

vector<double> GaussianCompare::Rec(vector<double> &Qicsd, vector<double> &Picsd, vector<double> &Q, double sigma, double height, double shift, vector<double> &vec) {
    int L = Q.size();
//    vector<double> vec(L, 0.0);
	vec.clear();
	for (int i = 0; i < L; i++) vec.push_back(0.0);

    for (int i = 0; i < L; i++)
        for (int j = 0; j < Qicsd.size(); j++)
			//if (Q[0] <= Qicsd[j]*shift && Qicsd[j]*shift <= Q[Q.size()-1])
	            vec[i] += gaussianFunc(sigma, shift, Qicsd[j], height*Picsd[j], Q[i]);
    return vec;
}

vector<double> GaussianCompare::L2update(vector<double> &xrd, vector<double> &xrd_rec, double &sigma, double &height, double &shift, vector<double> Qicsd, vector<double> Picsd, vector<double> Q) {

//	fprintf(stderr, "sigma: %lf, height: %lf, shift: %lf\n", sigma, height, shift);

	int L = xrd.size();
	double sigma1 = 0.0, sigma2 = 0.0;
	for (int l = 0; l < L; l++) {
		double s1 = 0.0, s2 = 0.0;
		for (int i = 0; i < Picsd.size(); i++) {
			double loc = Qicsd[i]*shift;
			double x = Q[l];
			s1 += height*Picsd[i]/(sqrt(2*pi)*sigma*sigma)*exp(-(x-loc)*(x-loc)/(2*sigma*sigma))*((x-loc)*(x-loc)/(sigma*sigma));
			s2 += height*Picsd[i]/(sqrt(2*pi)*sigma*sigma)*exp(-(x-loc)*(x-loc)/(2*sigma*sigma));
		}
		sigma1 += xrd_rec[l]*s1+xrd[l]*s2;
		sigma2 += xrd_rec[l]*s2+xrd[l]*s1;

//		double xv = xrd[l], recv = xrd_rec[l];
//		double pos = recv*recv*recv+recv*xv*xv, neg = recv*recv*xv+xv*xv*xv;
//		sigma1 += pos*s1+neg*s2;
//		sigma2 += neg*s1+pos*s1;

//		sigma1 += 1*s1+xv/recv*s2;
//		sigma2 += 1*s2+xv/recv*s1;
	}
	sigma *= sigma2/sigma1;

    if (sigma > this->MatchSigmaTol) sigma = this->MatchSigmaTol;
		
	xrd_rec = Rec(Qicsd, Picsd, Q, sigma, height, shift, xrd_rec);

	double shift1 = 0.0, shift2 = 0.0;
	for (int l = 0; l < L; l++) {
		double s1 = 0.0, s2 = 0.0;
		for (int i = 0; i < Picsd.size(); i++) {
			double loc = Qicsd[i]*shift;
			double x = Q[l];
			s1 += height*Picsd[i]/(sqrt(2*pi)*sigma*sigma*sigma)*exp(-(x-loc)*(x-loc)/(2*sigma*sigma))*(x)*Qicsd[i];
			s2 += height*Picsd[i]/(sqrt(2*pi)*sigma*sigma*sigma)*exp(-(x-loc)*(x-loc)/(2*sigma*sigma))*(loc)*Qicsd[i];
		}
		shift1 += xrd_rec[l]*s1+xrd[l]*s2;
		shift2 += xrd_rec[l]*s2+xrd[l]*s1;

//		double xv = xrd[l], recv = xrd_rec[l];
//		double pos = recv*recv*recv+recv*xv*xv, neg = recv*recv*xv+xv*xv*xv;
//        shift1 += pos*s1+neg*s2;
//        shift2 += neg*s1+pos*s1;

//		shift1 += 1*s1+xv/recv*s2;
//		shift2 += 1*s2+xv/recv*s1;
	}
	shift *= shift2/shift1;

//	fprintf(stderr, "this->shiftCenter: %lf\n", this->shiftCenter);
//	fprintf(stderr, "lshift: %lf\n", this->shiftCenter-this->MatchShiftTol);
//	fprintf(stderr, "rshift: %lf\n", this->shiftCenter+this->MatchShiftTol);
    if (shift < this->shiftCenter-this->MatchShiftTol) shift = this->shiftCenter-this->MatchShiftTol;
    if (shift > this->shiftCenter+this->MatchShiftTol) shift = this->shiftCenter+this->MatchShiftTol;

	xrd_rec = Rec(Qicsd, Picsd, Q, sigma, height, shift, xrd_rec);
	double height1 = 0.0, height2 = 0.0;
	for (int l = 0; l < L; l++) {
		double h = 0.0;
		for (int i = 0; i < Picsd.size(); i++) {
			double loc = Qicsd[i]*shift;
			double x = Q[l];
			h += Picsd[i]/(sqrt(2*pi)*sigma)*exp(-(x-loc)*(x-loc)/(2*sigma*sigma));
		}
		height1 += xrd_rec[l]*h;
		height2 += xrd[l]*h;

//		double xv = xrd[l], recv = xrd_rec[l];
//		double pos = recv*recv*recv+recv*xv*xv, neg = recv*recv*xv+xv*xv*xv;
//        height1 += pos*h;
//        height2 += neg*h;

//		height1 += 1*h;
//		height2 += xv/recv*h;
	}
	height *= height2/height1;

	xrd_rec = Rec(Qicsd, Picsd, Q, sigma, height, shift, xrd_rec);
	return xrd_rec;
}


double GaussianCompare::optimize(vector<double> Qicsd, vector<double> Picsd, vector<double> xrd, vector<double> Q, vector<double> &xrd_rec, int index, int kindex, double &ratio, double &sigma, double &height, double &shift) {
	time_t init_time = time(0);
	int L = xrd.size();
	int len = Qicsd.size();
	sigma = 0.2, height = 1.0, shift = 1.0;

	while (Qicsd[0] < Q[0]) {
		Qicsd.erase(Qicsd.begin());
		Picsd.erase(Picsd.begin());
	}

	double qxrd = 0.0;
	double maxh = 0.0;
	for (int i = 0; i < L; i++) {
		if (xrd[i] > maxh) {
			maxh = xrd[i];
			qxrd = Q[i];
		}
	}
	for (int i = 0; i < L; i++) xrd[i] *= 1.0/maxh;

	maxh = 0.0;
	double qicsd = 0.0;
	int indexone = -1;
	for (int i = 0; i < len; i++) {
		if (Picsd[i] > maxh) {
			maxh = Picsd[i];
			qicsd = Qicsd[i];
			indexone = i;
		}
	}
	double secondmaxh = 0.0;
	double secondqicsd = 0.0;
	for (int i = 0; i < len; i++) {
		if (i != indexone && Picsd[i] > secondmaxh) {
			secondmaxh = Picsd[i];
			secondqicsd = Qicsd[i];
		}
	}

	for (int i = 0; i < len; i++) Picsd[i] *= 1.0/maxh;
	//for (int i = 0; i < len; i++) Qicsd[i] *= qxrd/qicsd;	
	xrd_rec = Rec(Qicsd, Picsd, Q, sigma, height, shift, xrd_rec);
	double cost = L2loss(xrd, xrd_rec);


	double c = 1e-5;
	int iter_cnt = 0;
	while (1) {
		xrd_rec = L2update(xrd, xrd_rec, sigma, height, shift, Qicsd, Picsd, Q);
		double new_cost = L2loss(xrd, xrd_rec);
		time_t now_time = time(0);

		if (abs(new_cost-cost) < c * cost /*|| now_time-init_time > time_limit*/) {
			break;
		}
		cost = new_cost;
		iter_cnt++;
	}

	vector<double> zero_xrd(xrd.size(), 0.0);
	ratio = L2loss(xrd, xrd_rec)/L2loss(xrd, zero_xrd);

	return cost;
}

static bool GaussianCompare::myComp(pair<pair<int, double>, double> a, pair<pair<int, double>, double> b) {
	return a.first.second < b.first.second;
}

static bool GaussianCompare::myComp2(pair<pair<int, double>, double> a, pair<pair<int, double>, double> b) {
    return a.second < b.second;
}

vector<double> GaussianCompare::convertNorm(vector<double> xrd, vector<double> Qlogsteps, vector<double> Q) {
    vector<double> new_xrd;

	vector<double> X(Qlogsteps.begin(), Qlogsteps.end());
	vector<double> Y(xrd.begin(), xrd.end());

	tk::spline s;
	s.set_points(X,Y);

	int L = Q.size();
	for (int l = 0; l < L; l++) new_xrd.push_back(s(Q[l]));

	return new_xrd;
}

void GaussianCompare::matchIcsd(int i, vector<double> Qicsd, vector<double> Picsd, vector<double> xrd, vector<double> Q, int k) {

    /*
	fprintf(stderr, "\n\n-----------------------\n");
	fprintf(stderr, "Now start phase %d match icsd %d:\n", k, i);
	fprintf(stderr, "-------------------------\n");
    */

	
	int icsdlen = Qicsd.size();
	int qlen = Q.size();
	for (int i = icsdlen-1; i >= 0; i--) {
		if (Qicsd[i] < Q[0] || Q[qlen-1] < Qicsd[i]) {
			Qicsd.erase(Qicsd.begin()+i);
			Picsd.erase(Picsd.begin()+i);
		}
	}

	icsdlen = Qicsd.size();
	double lshift = this->shiftCenter-this->MatchShiftTol;
	double rshift = this->shiftCenter+this->MatchShiftTol;
	for (int i = qlen-1; i >= 0; i--) {
		if (Q[i] < Qicsd[0]*lshift || Qicsd[icsdlen-1]*rshift < Q[i]) {
			Q.erase(Q.begin()+i);
			xrd.erase(xrd.begin()+i);
		}
	}


	vector<double> xrd_rec(xrd.size(), 0.0);
	double ratio = 1.0;
	double sigma, height, shift;

	double loss = optimize(Qicsd, Picsd, xrd, Q, xrd_rec, i, k, ratio, sigma, height, shift);

	if (i == 7) {
		fprintf(stderr, "\n\nbasis %d: Qlen: %d, xrd_len: %d\n", k, Q.size(), xrd.size());
		fprintf(stderr, "Q_min: %lf, Q_max: %lf\n", Q[0], Q[Q.size()-1]);
		fprintf(stderr, "Qicsd %d: Qicsd: %d, Picsd: %d\n", i, Qicsd.size(), Picsd.size());
		fprintf(stderr, "Qicsd_min: %lf, Q_max: %lf\n", Qicsd[0], Qicsd[Qicsd.size()-1]);
	}

	matchI[i].loss = loss;
	matchI[i].sigma = sigma;
	matchI[i].height = height;
	matchI[i].shift = shift;
	matchI[i].ratio = ratio;
}

void GaussianCompare::compare(ostream &os, char* sticksfile, mat W, vector<mat> H, int K, vector<double> Q, vector<double> Qlogsteps) {

    srand(time(0));

	time_t init_time = time(NULL);

	freopen(sticksfile, "r", stdin);
    vector<vector<double> > Qicsd;
    vector<vector<double> > Picsd;
    char c;
    char buf[1000];
    
    int cnt = 0;
	while ((c = getchar()) != EOF) {
        if (c == 10) continue; // reach the end of one line
        if (c == 13) { // reach the end of one line
            getchar();continue;
        }
    	
		string buf="";
		buf = buf+c;
		char cc;
		while (1) {
			cc = getchar();
			if (cc == 10) break;
			if (cc == 13) {
				getchar();break;
			}
			buf = buf+cc;
		}

//        getline(cin, buf);
//        buf = c+buf;

        if (c == 'Q') {
            size_t found = buf.find("=");
            if (found == string::npos) continue;
            buf.erase(0, found+1);

            vector<double> temp;

            readDoubleArray(temp, buf);
            Qicsd.push_back(temp);
            
        } else if (c == 'P') {
            size_t found = buf.find("=");
            if (found == string::npos) continue;
            buf.erase(0, found+1);

            vector<double> temp;
            
            readDoubleArray(temp, buf);
            Picsd.push_back(temp);
        }
    }

    fclose(stdin);

	int L = W.n_rows;
    int M = H.size();
    int N = H[0].n_cols;
	this->shiftCenter = 1.0-((1.0+Qlogsteps[M-1]/Qlogsteps[0])/2.0-1.0);

    for (int k = 0; k < K; k++) {
        double sumW = norm(W.cols(k,k), "fro");
        if (sumW > 1e-6) {
            for (int l = 0; l < L; l++) W(l,k) /= sumW;
            for (int n = 0; n < N; n++) {
                for (int m = 0; m < M; m++)
                    H[m](k,n) *= sumW;
            }
        }
    }

    /*
	char outputfilename[200];
	sprintf(outputfilename, "%s_report.txt", argv[1]);
	freopen(outputfilename, "w", stdout);
    */

	for (int k = 0; k < K; k++) {
		vector<pair<pair<int, double>, double> > sols;
		double minloss = 1e10;
		int chosen = -1;
		vector<double> xrd;
		double csigma, cheight, cshift;

		for (int l = 0; l < W.n_rows; l++) xrd.push_back(W(l,k));
        xrd = convertNorm(xrd, Qlogsteps, Q);

		#pragma omp parallel
		{
			#pragma omp for
			for (int i = 0; i < Qicsd.size(); i++) {
				matchIcsd(i, Qicsd[i], Picsd[i], xrd, Q, k);
			}
		}

		for (int i = 0; i < Qicsd.size(); i++) sols.push_back(make_pair(make_pair(i, matchI[i].ratio), matchI[i].loss));

		double minLoss = 1e10;
		for (int i = 0; i < Qicsd.size(); i++)
			if (matchI[i].loss < minLoss) {
				minLoss = matchI[i].loss;
				csigma = matchI[i].sigma;
				cshift = matchI[i].shift;
				cheight = matchI[i].height;
				chosen = i;
			}
      
		std::sort(sols.begin(), sols.end(), &GaussianCompare::myComp);
        //printf("\nphase %d ratio:\n", k);
        os << "\nphase " << k << " ratio:\n";
        //printf("[");
        os << "[";
        for (int i = 0; i < sols.size(); i++) {
            int id = sols[i].first.first;
            sprintf(buf, "%d (loss ratio: %.3lf, loss: %.8lf, sigma: %.6lf, height: %.6lf, shift: %.6lf)", id, sols[i].first.second, sols[i].second, matchI[id].sigma, matchI[id].height, matchI[id].shift);
            os << buf;
            if (i != sols.size()-1) os << ", ";//printf(", ");
        }
        //printf("]\n");
        os << "]\n";


        std::sort(sols.begin(), sols.end(), &GaussianCompare::myComp2);
        sprintf(buf, "phase %d loss:\n", k);
        os << buf;
        
        //printf("[");
        os << "[";
        for (int i = 0; i < sols.size(); i++) {
            int id = sols[i].first.first;
            sprintf(buf, "%d (loss ratio: %.3lf, loss: %.8lf, sigma: %.6lf, height: %.6lf, shift: %.6lf)", id, sols[i].first.second, sols[i].second, matchI[id].sigma, matchI[id].height, matchI[id].shift);
            os << buf;
            if (i != sols.size()-1) os << ", ";//printf(", ");
        }
        //printf("]\n");
        os << "]\n";


		double maxh = 0.0;
		for (int i = 0; i < xrd.size(); i++) maxh = std::max(maxh, xrd[i]);
		for (int i = 0; i < xrd.size(); i++) xrd[i] /= maxh;
		//printf("\nSR%d=", k);
        os << "\nSR" << k << "=";
		vector<double> xrd_r;
		xrd_r = Rec(Qicsd[chosen], Picsd[chosen], Q, csigma, cheight, cshift, xrd_r);
		for (int i = 0; i < xrd.size(); i++) {
            sprintf(buf, "(%lf, %lf) ", xrd[i], xrd_r[i]);
            os << buf;
        }
		//printf("\n");
        os << "\n";

        
		fprintf(stderr, "\n-----------------\n");
		fprintf(stderr, "phase %d: %d\n", k, chosen);
		for (int i = 0; i < sols.size(); i++) fprintf(stderr, "(%d, %lf) ", sols[i].first.first, sols[i].second);
		fprintf(stderr, "\n");
		fprintf(stderr, "------------------\n");
        

	}
	//printf("\n");
    os << "\n";

}
