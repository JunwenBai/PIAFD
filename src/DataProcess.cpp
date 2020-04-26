#include "DataProcess.hpp"

mat DataProcess::re_sample_log(mat &A, double * qvalues, double * q_values_log){
    int Q = A.n_rows;
    double Qlogmin = log(qvalues[0]);
    double Qlogmax = log(qvalues[Q-1]);
    double Qstepsize = (Qlogmax-Qlogmin)/(Q-1);
    for (int i = 0; i < Q; i++) {
        q_values_log[i] = exp(Qlogmin+i*Qstepsize); // Qlogsteps is a geometric sequence such that shifting Qlogsteps 1 position to the right means multiplying each qvalue by a constant
    }

    return DataProcess::resample(A,qvalues,q_values_log);
}


mat DataProcess::resample(mat &A, double * current_values, double * new_values){
    mat Aresampled(A.n_rows, A.n_cols);
    int Q = A.n_rows;
    int N = A.n_cols;

    vector<double> X;
    for (int q = 0; q < Q; ++q){
        X.push_back(current_values[q]);
    }

    vector<double> newsteps;
    for (int i = 0; i < Q; i++) {
        newsteps.push_back(new_values[i]); 
    }

    vector<double> Y;
    for (int k = 0; k < N; ++k){
        Y.clear();
        for (int j = 0; j < Q; j++) {
            Y.push_back(A(j, k));
        }
        tk::spline s;
        s.set_points(X,Y);
        for (int j = 0; j < Q; j++) {
            Aresampled(j,k) = s(new_values[j]);
        }

    }
    return Aresampled;
}






void DataProcess::extract_phase_mixture(vector<vector<PhaseMixture> > &phase_mixtures, vector<mat> &H, double *Qlog, int N, int K, int M, double epsilon){
    for (int k = 0; k < K; k++) {           // J code
        vector<PhaseMixture> mixtures;
        for (int i = 0; i < N; i++) {
            double sumH = 0.0;
            double shiftH = 0.0;

            for (int m = 0; m < M; m++) {
                sumH += H[m](k,i); // total proportion
                double shift = Qlog[m]/Qlog[0]; // generate the real shift values
                shiftH += shift * H[m](k,i);
            }
            if (sumH < epsilon) shiftH = 0.0; // phase k doesn't appear at sample point i
            else shiftH /= sumH; // weighted shift value

            PhaseMixture mixture;
            mixture.isInvolved = (sumH >= epsilon); // whether phase k involves at sample point i
            mixture.proportion = (sumH >= epsilon)?sumH:0.0; // proportion
            mixture.shift = shiftH; // shift
            mixtures.emplace_back(mixture);
        }
        phase_mixtures.emplace_back(mixtures);
    }
}










void DataProcess::save_mat_csv(mat &A, std::ostream &os){
    for (int row = 0; row < A.n_rows; ++row){
        for (int col = 0; col < A.n_cols; ++col){
            if (col > 0){   os << ";";  }
            os << A(row, col);
        }
        os << "\n";
    }

}

void DataProcess::save_array_csv(double * A, int length, std::ostream &os){
    for (int col = 0; col < length; ++col){
        if (col > 0){   os << ";";  }
        os << A[col];
    }
}

void DataProcess::save_mat_equal(mat &A, string name ,std::ostream &os){
    for (int row = 0; row < A.n_rows; ++row){
        os << "\n"<< name << "_r" << row << "=";
        for (int col = 0; col < A.n_cols; ++col){
            if (col > 0){   os << ",";  }
            os << A(row, col);
        }
    }

}

void DataProcess::save_array_equal(double * A, int length, string name ,std::ostream &os){
    os << "\n"<< name << "=";
    for (int col = 0; col < length; ++col){
        if (col > 0){   os << ",";  }
        os << A[col];
    }
}









void DataProcess::save_data_A(mat &A, double * qvalues){
	int Q = A.n_rows;
    int N = A.n_cols;

	std::filebuf fb;
	fb.open ("../data/data_out_cpp.csv",std::ios::out);
  	std::ostream os(&fb);

  	for (int q = 0; q < Q; ++q)
  	{
  		if (q > 0){ 	os << ";"; 	}
  		os << qvalues[q];
  	}
  	for (int m = 0; m < N; m++) {
  		os << "\n";
        for (int q = 0; q < Q; q++) {
			if (q > 0){ 	os << ";"; 	}
  			os << A(q, m);
        }
    }

  	fb.close();
}









