#include "Solver.hpp"
#define RAND_MAX 999

const double Solver::epsilon = 1e-9;

Solver::Solver(Instance &I_, Config &cfg_):
    A(I_.Qsize, I_.N),  // A = Q * N matrix
    W(I_.Qsize, cfg_.K),
    H(),
    K(cfg_.K),
    R(I_.Qsize, I_.N),
    sparsity(cfg_.K, I_.N),
    sparsityW(I_.Qsize, cfg_.K),
    sparsityCoefficient(cfg_.sparsity_coefficient),
    cfg(cfg_),
    Qlog(new double[I_.Qsize]),
    I(I_),
    Q(I_.Qsize),
    N(I_.N),
    M(cfg_.M),
    Comp(I_.N, 3),
    last_cost(0.),
    alpha(0.0001),
    beta(0.1)
    {
        for (int row = 0; row < A.n_rows; ++row){       //  a row contains a given wavelength
            for (int col = 0; col < A.n_cols; ++col){   // columns cotains the sample points
                A(row, col) = I_.I[col * I_.Qsize + row];
            }
        }

        for (int row = 0; row < W.n_rows; ++row){       // Init to 0
            for (int col = 0; col < W.n_cols; ++col){
                W(row,col) = 0.0;
                sparsityW(row,col) = sparsityCoefficient;            
            }
        }
        for (int m = 0; m < M; ++m)
        {
            H.emplace_back(K, I_.N);
            for (int row = 0; row < H[m].n_rows; ++row){       // Init to 0
                for (int col = 0; col < H[m].n_cols; ++col){
                    H[m](row,col) = 0.0;
                }
            }
        }
        for (int row = 0; row < sparsity.n_rows; ++row){       // Init to sparsity coef
            for (int col = 0; col < sparsity.n_cols; ++col){
                sparsity(row,col) = sparsityCoefficient;
            }
        }
        int num_elements = 3; // this should be an instance variable and read from instance file
        for (int i = 0; i < N; ++i){
            for( int j = 0; j < num_elements; ++j){
                Comp(i,j) = I_.comp_by_elements[ num_elements * i + j];
            }
        }

        // start to build connectivity graph
        if (cfg.use_dist > 0) {
            for (int i = 0; i < N; i++) {
                vector<int> nb;
                nb.clear();
                for (int j = 0; j < N; j++) {
                    double distance = 0.0;
                    for (int k = 0; k < 3; k++) {
                        double dx = Comp(i,k)-Comp(j,k);
                        distance += dx*dx;
                    }
                    if (sqrt(distance) < cfg.dist_cutoff)
                        nb.push_back(j);
                }
                neighbors.push_back(nb);
            }
        } else {
    
            FILE * file = freopen(cfg.edgefile.c_str(), "r", stdin);
            if (file == NULL) {
                fprintf(stderr, "edgefile is invalid.\n");
                return;
            }
            int sample;
            char c;
            neighbors.clear();
            for (int i = 0; i < N; i++) {
                vector<int> nb;
                nb.clear();
                scanf("%d", &sample);
                c = getchar();
                while (c == ',') {
                    int v;
                    scanf("%d", &v);
                    nb.push_back(v);
                    c = getchar();
                }
                neighbors.push_back(nb);
            }

            fprintf(stderr, "len: %d\n", neighbors.size());
            
            fclose(stdin);
        }

    }


void Solver::solve(){
    init_time = time(NULL);
    ChronoP chrono;
    if (cfg.applyLog > 0){
        A = DataProcess::re_sample_log(A, I.Q, Qlog);
    }else{
        Qlog = I.Q;
    }
    // normalizeA();

    switch(cfg.kmeans_init){
        case 0:
            init_from_data();
            break;
        case 1:
            init_from_data_pp();
            break;
        case 2:
            init_from_random();
            break;
        default:
            init_from_random();
            break;
    }
    
    W_init = W; // To remove

    reconstruction();

    /*
    for (int i = 0; i < 5; i++)
        fprintf(stderr, "%.4lf,", A(i,0));
    fprintf(stderr, "\n");
    return;
    */

    if(cfg.gradientMethod == 0){
        if (cfg.use_KL == true) last_cost = lossFunction();
        else last_cost = lossFunctionL2();
    }
    else if(cfg.gradientMethod == 1){last_cost = lossFunctionL1();}
    else if(cfg.gradientMethod == 2){last_cost = lossFunctionL2();}
    else if(cfg.gradientMethod == 3){last_cost = lossFunctionL21();}

    double current_cost = last_cost;
    cout << current_cost << "\t init cost \n";
    chrono.Start();
    int n_converge = 0, n_1 = 0;
    double alpha = cfg.gibbs_relax;
    if (cfg.gibbs_relax > 1.0) alpha = 0.05;
    int n_iters = 0;
    double delta_L = last_cost;

    for (int i = 0; i < cfg.number_of_iteration_step*3; ++i) {
        n_iters++;
        
        if (n_converge == 1 && cfg.gibbs_relax > 1.0) alpha = min(alpha+0.002, 0.3);

        if (n_converge >= 2 && i > 0 && cfg.gibbs_frequency > 0 && i % cfg.gibbs_frequency == 0){
            if (n_converge == 2) enforceGibbs(1);
            else enforceGibbs(1, true);
        }
        if (n_converge == 1) enforceGibbs(alpha);

        if (n_converge >= 2 && i > 0 && cfg.alloy_frequency > 0 && i % cfg.alloy_frequency == 0) {
            if (n_converge == 2) enforceAlloy(1);
            else enforceAlloy(1, true);
        }
        if (n_converge == 1) enforceAlloy(alpha);
        

        if(cfg.gradientMethod == 0) {       gradientUpdate();}
        else if(cfg.gradientMethod ==1){    gradientUpdateL1();}
        else if(cfg.gradientMethod ==2){    gradientUpdateL2();}
        else if(cfg.gradientMethod ==3){    gradientUpdateL21();}
    
        reconstruction();
        if(cfg.gradientMethod == 0) {
            if (cfg.use_KL == true) current_cost = lossFunction();
            else current_cost = lossFunctionL2();
        }
        else if(cfg.gradientMethod ==1){current_cost = lossFunctionL1();}
        else if(cfg.gradientMethod ==2){current_cost = lossFunctionL2();}
        else if(cfg.gradientMethod ==3){current_cost = lossFunctionL21();}

        chrono.Stop();
        double tolerance = 2e-5;
        if (i % 25 == 0 || abs(last_cost - current_cost)/last_cost < tolerance)
            cout << current_cost << "\t cost at iter " << i <<" and time "<<chrono.ellapsed_second() << "s with average "<< (chrono.ellapsed_m_second()/(i+1))<< "ms/iter, alpha: " << alpha <<"\n";

        if (chrono.ellapsed_second() > cfg.time_out_in_second){
            break;
        }
        if (cfg.stop_on_decrease_obj > 0 && current_cost > last_cost){
            break;
        }
        if (abs(last_cost - current_cost)/last_cost < tolerance || (abs(last_cost-current_cost)/last_cost < tolerance && n_converge == 0)) {

            /*
            for (int i = 0; i < 5; i++)
                fprintf(stderr, "(%lf, %lf), ", A(i,0), R(i,0));
            return;
            */

            cout << "n_converge: " << ++n_converge << "\n";
            if (n_converge == 1) {
                n_iters = 0;
                enforceGibbs(alpha);
                enforceAlloy(alpha);
            } else if (n_converge == 2) {
                n_iters = 0;
                enforceConnectivity(1-alpha);
                enforceGibbs(1);
                enforceAlloy(1);
            } else if (n_converge == 3) {
                break;
            } else break;
        }

        delta_L = (last_cost-current_cost)/last_cost;
        last_cost = current_cost;
    }

    postProcess();

}



// postProcess
void Solver::postProcess(){
    normalizeHAndUpdateW();
    enforceConnectivity();

    // print_diff_W();
    if (cfg.use_KL) cout << "Final loss:" << lossFunction() << "\n";
    else cout << "Final loss:" << lossFunctionL2() << "\n";


    save_to_file_phase_mapper(cfg.outfile);

}

void Solver::gradientUpdate(){

	normalizeW();
    reconstruction();

//    fprintf(stderr, "Rec: %lf\n", H[0](0,0));

    mat O;
    if (cfg.use_KL == true) {
        O = ones<mat>(Q, N); // an all-1 matrix
        Aux1 = A / (R + epsilon);
    } else {
        O = R+epsilon;
        Aux1 = A;
    }

    // H update
	for (int shift = 0; shift < M; shift++) { 
    	H[shift] = H[shift] % 
            ((W.rows(0, Q -shift -1).t() * Aux1.rows(shift, Q-1))
                /(W.rows(0, Q -shift -1).t()*O.rows(shift, Q-1) + sparsity));
	}

//    fprintf(stderr, "Rec after: %lf\n", H[0](0,0));

	// update W
	reconstruction();

    if (cfg.use_KL == true) {
        Aux1 = A / (R + epsilon); 
    } else {
        Aux1 = A;
        O = R + epsilon;
    }

    mat W_x = zeros<mat>(Q, K); // auxiliary matrix (numerator in the updating rule)
    mat W_y = zeros<mat>(Q, K); // auxiliary matrix (denominator)

    mat W1;
    rowvec h(K);
    W1 = ones<mat>(Q, K);
    for (int i = 0; i < M; i++) {
        int p = i;

        mat Wxp = Aux1.rows(p, Q-1) * H[i].t();     // (A/R) * H
        Wxp.insert_rows(Wxp.n_rows, p);
        mat Wyp = O.rows(p, Q-1) * H[i].t();        // 11^t * H
        Wyp.insert_rows(Wyp.n_rows, p); 
        
        // numerator
        h = sum(Wyp % W, 0);                        //h[k] = sum of colmun of (11^t * H) % W
        W1.each_row() = h;
        W1 = W1 % W;                                // W1 = W where column k is multiplied by h[k]

        W_x = W_x + Wxp + W1;
        
        // denominator
        h = sum(Wxp % W, 0);                        //h[k] = sum of colmun of ((A/R) * H) % W
        W1.each_row() = h;
        W1 = W1 % W;

        W_y = W_y + Wyp + W1;
    }

    mat GradW = W_x / (W_y + epsilon);
    W = W % GradW;
	
}


void Solver::gradientUpdateL2(){
    vector<mat> HSave;
    mat WSave = W;
    for (int shift = 0; shift < M; shift++) {
        HSave.emplace_back(H[shift]);
    }

    int i = 0;
    double new_cost;
    reconstruction();
    double last_cost = lossFunctionL2();
    mat AmR = A-R;

    mat gradW = zeros<mat>(Q, K); 
    vector<mat> gradH;
    for (int shift = 0; shift < M; shift++) { 
        mat WS = W.rows(0, Q - shift -1);               // H update
        WS.insert_rows(0,shift);
        gradH.emplace_back(-WS.t() * AmR + sparsity);

        //mat GWAmR = AmR.rows(0, Q - shift -1);             // W update
        //GWAmR.insert_rows(0,shift);
        mat GWAmR = AmR.rows(shift, Q -1);             // W update
        GWAmR.insert_rows(Q-shift-1,shift);

        gradW = gradW - (GWAmR * H[shift].t());
    }

    alpha /= beta;
    do {
        alpha *= beta;
        W = WSave;
        for (int shift = 0; shift < M; shift++) {
            H[shift] = HSave[shift] - alpha * gradH[shift];
            H[shift].transform( [](double val) { return (val < 0)? 0 : val; } );
        }
        W = W - alpha * (gradW + sparsityW);
        W.transform( [](double val) { return (val < 0)? 0 : val; } );
        reconstruction();
        new_cost = lossFunctionL2();
        i++;
    } while(new_cost > last_cost);

    normalizeHAndUpdateW();

    alpha /= beta;
}



void Solver::gradientUpdateL1(){
    vector<mat> HSave;
    mat WSave = W;
    for (int shift = 0; shift < M; shift++) {
        HSave.emplace_back(H[shift]);
    }

    int i = 0;
    double new_cost;
    reconstruction();
    double last_cost = lossFunctionL1();
    mat AmR = A-R;

    mat gradW = zeros<mat>(Q, K); 
    vector<mat> gradH;
    
    for (int shift = 0; shift < M; shift++) { 
        mat WS = W.rows(0, Q - shift -1);               // H update |x| = -x if x < 0 else x
        WS.insert_rows(0,shift);
        gradH.emplace_back(-WS.t() * sign(AmR));

        mat GWAmR = AmR.rows(0, Q - shift -1);             // W update
        GWAmR.insert_rows(0, shift);
        gradW = gradW - (sign(GWAmR) * H[shift].t());
    }

    alpha /= beta;
    do {
        alpha *= beta;
        W = WSave;
        for (int shift = 0; shift < M; shift++) {
            H[shift] = HSave[shift] - alpha * gradH[shift];
            H[shift].transform( [](double val) { return (val < 0)? 0 : val; } );
        }
        W = W - alpha * (gradW);
        W.transform( [](double val) { return (val < 0)? 0 : val; } );
        reconstruction();
        new_cost = lossFunctionL1();
        i++;
    } while(new_cost > last_cost);
    alpha /= beta;
}


void Solver::gradientUpdateL21(){
    vector<mat> HSave;
    mat WSave = W;
    for (int shift = 0; shift < M; shift++) {
        HSave.emplace_back(H[shift]);
    }

    int i = 0;
    double new_cost;
    reconstruction();
    double last_cost = lossFunctionL21();
    mat AmR = A-R;

    mat gradW = zeros<mat>(Q, K); 
    vector<mat> gradH;

    AmR.transform([](double val) { return (val < 0)? 2 * val : 1.; } );
    // AmR.transform([](double val) { return (val < 0)? -exp(-val) : 1.; } );
    for (int shift = 0; shift < M; shift++) { 
        mat WS = W.rows(0, Q - shift -1);                   // H update 
        WS.insert_rows(0,shift);
        gradH.emplace_back(- WS.t() * AmR);

        mat GWAmR = AmR.rows(0, Q - shift -1);              // W update
        GWAmR.insert_rows(0, shift);
        gradW = gradW - (GWAmR * H[shift].t());
    }
    alpha /= beta;
    do {
        alpha *= beta;
        W = WSave;
        for (int shift = 0; shift < M; shift++) {
            H[shift] = HSave[shift] - alpha * gradH[shift];
            H[shift].transform( [](double val) { return (val < 0)? 0 : val; } );
        }
        W = W - alpha * (gradW);
        W.transform( [](double val) { return (val < 0)? 0 : val; } );
        reconstruction();
        new_cost = lossFunctionL21();
        i++;
    } while(new_cost > last_cost);
    alpha /= beta;
}


double Solver::lossFunctionL2(){
    /*
    mat AmR = A-R;
    double l2 = 0;
    for (int q = 0; q < Q; q++)
        for (int n = 0; n < N; n++)
            l2 += AmR(q,n)*AmR(q,n);
    return l2;
    */
    double l2 = norm(A-R, "fro");
    return l2;
}

double Solver::lossFunctionL1(){
    mat AmR = A-R;
    double l1 = 0;
    for (int q = 0; q < Q; q++)
        for (int n = 0; n < N; n++)
            l1 += abs(AmR(q,n));
    return l1;
}

double Solver::lossFunctionL21(){
    mat AmR = A-R;
    // AmR.transform([](double val) { return (val < 0)? exp(-val) : val; } );
    AmR.transform([](double val) { return (val < 0)? val*val : val; } );
    return norm(AmR,1);
}

double Solver::lossFunction(){
    double ckl = 0.0; // KL divergence(cost)
    for (int q = 0; q < Q; q++) {
        for (int m = 0; m < N; m++) {
            double a  = A(q,m);
            double r = R(q,m);
            ckl += a * log((a + epsilon)/(r + epsilon)) - a + r; 
        }
    }
    return ckl;
}






void Solver::init_from_random(){
//    srand ( unsigned ( time(0) ) );
    std::srand(cfg.seed);
    for (int k = 0; k < K; ++k){
        for (int q = 0; q < Q; ++q){
            W(q,k) = (static_cast <float> (rand()) / static_cast <float> (RAND_MAX)) + epsilon;
//            W(q,k) = (rand() % (RAND_MAX+1)+1)/float(RAND_MAX+1);
        }
    }
    normalizeW();
    for (int i = 0; i < M; ++i){
        for (int k = 0; k < K; ++k){
            for (int m = 0; m < N; ++m){
                H[i](k,m) = static_cast <float> (rand()) / static_cast <float> (RAND_MAX) + epsilon;
//                H[i](k, m) = (rand() % (RAND_MAX+1)+1)/float(RAND_MAX+1);
            }
        }
    }
}

void Solver::init_from_data(){
//    std::srand ( unsigned ( std::time(0) ) );
    std::srand(cfg.seed);
    std::vector<int> sampleIDs;
    for (int i=0; i<N; ++i){ sampleIDs.push_back(i); }
    std::random_shuffle ( sampleIDs.begin(), sampleIDs.end() ); 
    for (int k = 0; k < K; ++k){
        for (int q = 0; q < Q; ++q){
            W(q,k) = A(q, sampleIDs[k]);
        }
    }
    for (int i = 0; i < M; ++i){
        for (int k = 0; k < K; ++k){
            for (int m = 0; m < N; ++m){
                H[i](k,m) = static_cast <float> (rand()) / static_cast <float> (RAND_MAX) + epsilon;
            }
        }
    }  
}


void Solver::init_from_data_pp(){
//    std::srand ( unsigned ( std::time(0) ) );
    std::srand(cfg.seed);
    int sampleIDs[K];
    sampleIDs[0] = rand() % N; // get a random index between 0 and N
    double r = 0;
    vec Dx(N); // minimum distance of a point to previously chosen center
    Dx.fill( std::numeric_limits<double>::infinity() ); // initialize all distances to infinity
    vec cumDx(N); // normalized cummulative sum of Dx^2
    for (int k = 0; k < K - 1; ++k){
        // get minimum distances between points and closest center
        Dx(0) = min( Dx(0), norm( Comp.row(0) - Comp.row(sampleIDs[k]) ) );
        cumDx(0) = pow(Dx(0), 2);
        for (int i = 1; i < N; ++i) {
            Dx(i) = min( Dx(i), norm( Comp.row(i) - Comp.row(sampleIDs[k]) ) );
            cumDx(i) = pow(Dx(i), 2) + cumDx(i-1);
        }
        cumDx /= (double)sum(Dx % Dx); // normalize cumDx 

        // draw next point based on the cummulative probability distribution defined by cumDx
        r = (double) rand() / RAND_MAX; // get random double between 0 and 1
        for (int i = 0; i < N; ++i) {
            if ( r < cumDx(i) ) {
                sampleIDs[k+1] = i;
                break;
            }
        }   
    }
    
    // initialize W based on composition points in sampleIDs
    for (int k = 0; k < K; ++k){
        W.col(k) = A.col(sampleIDs[k]);
    }

    // I wonder if there is a better way of doing this than random initialization
    for (int i = 0; i < M; ++i){
        for (int k = 0; k < K; ++k){
            for (int m = 0; m < N; ++m){
                H[i](k,m) = static_cast <float> (rand()) / static_cast <float> (RAND_MAX) + epsilon;
            }
        }
    }  
}







void Solver::reconstruction(){
    R.zeros();
    for (int shift = 0; shift < M; ++shift) {
        mat Wyp = W.rows(0, Q - shift -1) * H[shift];        // W * H[shift]
        Wyp.insert_rows(0, shift); 
        R += Wyp;                                       // sum of W * H[shift]
    }
}


void Solver::normalizeW(){
    for (int j = 0; j < K; j++) {
        double W2norm = norm(W.cols(j, j), "fro"); // calculate frobenius norm of the j-th column
        for (int i = 0; i < Q; i++) {
            if (W(i,j) < epsilon) continue;
            W(i,j) /= W2norm;
        }
    }
}


void Solver::enforceGibbs(double alpha, bool hard = false){
    if (cfg.gibbs_relax < 0){
        Gibbs::relax(H, Q, K, N, M, epsilon, Comp);
    }else{
        Gibbs::enforce(H, Q, K, N, M, epsilon, Comp, alpha, hard);
    }
}

void Solver::enforceAlloy(double alpha, bool hard = false){
   if (cfg.alloy_relax < 0) {
       Alloy::relax(H, Q, K, N, M, 0.003, neighbors, Qlog, Comp);
   } else {
       Alloy::enforce(H, Q, K, N, M, 0.003, neighbors, Qlog, Comp, alpha, hard);
   }
}

void Solver::enforceConnectivity(double beta=0.0) {
    if (beta > 0.0) {
        Connectivity::connectPhase(H, Q, N, K, M, connect_indicators, epsilon, neighbors, beta);
        Connectivity::connect(H, Q, N, K, M, connect_indicators, epsilon, neighbors, beta);
        return;
    }

    int n_rounds = 0;
    fprintf(stderr, "Start to enforce connectivity\n");
    for (int n = 0; n < N; n++)
        for (int k = 0; k < K; k++)
            for (int m = 0; m < M; m++)
                if (H[m](k,n) < 1e-4)
                    H[m](k,n) = 0.0;
    while (! Connectivity::checkCC(H, N, K, M, epsilon, neighbors)) {
        fprintf(stderr, "Connectivity round %d\n", ++n_rounds);
        Connectivity::connectPhase(H, Q, N, K, M, connect_indicators, epsilon, neighbors, beta);
        Connectivity::connect(H, Q, N, K, M, connect_indicators, epsilon, neighbors, beta);
        enforceAlloy(1, true);
    }
}

void Solver::normalizeA(){
    A = normalise(A,1);
}

void Solver::normalizeH(){
    for (int shift = 0; shift < M; ++shift){
        H[shift] = normalise(H[shift],1);
    }
}

void Solver::normalizeHAndUpdateW(){ //J
    for (int k = 0; k < K; k++) {
        double maxsumH = 0.0;
        for (int n = 0; n < N; n++) {
            double sumH = 0.0;
            for (int m = 0; m < M; m++) {
                sumH += H[m](k, n);
            }
            maxsumH = max(maxsumH, sumH);
        }

        // normalize maxsumH to be 1
        for (int m = 0; m < M; m++) {
            for (int i = 0; i < N; i++) {
                double valH = H[m](k,i);
                if (maxsumH < epsilon) H[m](k,i) = 0;
                else H[m](k,i) = valH / maxsumH;
            }
        }

        // to maintain the resulting matrix of W*H, multiply W(i,k) by some maxsumH
        for (int i = 0; i < Q; i++) {
            double valW = W(i,k);
            W(i,k) = valW * maxsumH;
        }
    }
}




void Solver::save_to_file_csv(string filename){
    std::filebuf fb;
    fb.open (filename,std::ios::out);
    std::ostream os(&fb);
   
    DataProcess::save_mat_csv(A, os);
    os << "\n";  

    // DataProcess::save_array_csv(I.Q, Q, os)
    // os << "\n";  

    DataProcess::save_mat_csv(W, os);
    os << "\n";   

    for (int m = 0; m < M; ++m){
        DataProcess::save_mat_csv(H[m], os);
        os << "\n";  
    }

}
void Solver::save_to_file_phase_mapper(string filename){
    std::filebuf fb;
    fb.open (filename,std::ios::out);
    std::ostream os(&fb);


    mat W_non_log = W;
    mat A_non_log = A; 

    /*fprintf(stderr, "old:\n");
    for (int i = 0; i < 10; i++)
        fprintf(stderr, "%lf,", A(i,N-1));
    fprintf(stderr, "\n");*/

    if (cfg.applyLog > 0){
        W_non_log = DataProcess::resample(W, Qlog, I.Q);    
        A_non_log = DataProcess::resample(A, Qlog, I.Q);     
    }

    //for (int i = 0; i < 10; i++) fprintf(stderr, "(%lf, %lf),", Qlog[i], I.Q[i]);

    time_t now_time = time(NULL);
    os << "Time cost: " << (int)(now_time - init_time) << " s" << "\n";
    os << "L1:" << lossFunctionL1() << "\n";
    os << "L2:" << lossFunctionL2() << "\n";
    os << "KL:" << lossFunction() << "\n";
//    os << "L21:" << lossFunctionL21() << "\n";

    os << "K=" <<  K << "\n";

    os << "\n// Phase pattern (basis)\n";
    DataProcess::save_array_equal(I.Q, Q, "Q", os);
    os << "\n";  
    
    double * B = new double[max(Q,N)];
    for (int k = 0; k < K; ++k){
        string strB("B"+to_string(k+1));
        for (int q = 0; q < Q; ++q){
            B[q] = W_non_log(q,k);
        }
        DataProcess::save_array_equal(B, Q, strB, os);
//        os << "\n";  
    }
    os << "\n";  

    os <<"\n// Phase concentrations at each sample\n";
    vector<vector<PhaseMixture> > phase_mixtures; // shifts, proportion
    DataProcess::extract_phase_mixture(phase_mixtures, H, Qlog, N, K, M, epsilon);
    
    for (int  n = 0; n < N; n++) {
        string strB("C"+to_string(n+1));
        for (int k = 0; k < K; k++) {
            B[k] = (phase_mixtures[k][n].proportion < epsilon)? 0.0 : phase_mixtures[k][n].proportion; // mark here
        }
        DataProcess::save_array_equal(B,  K, strB, os);
    }


    os << "\n\n// Per-phase model for each sample\n";
    double logxrd[Q+1];
    for (int n = 0; n < N; n++) {
        for (int k = 0; k < K; k++) {
            for (int l = 0; l < Q; l++) logxrd[l] = 0.0;
            for (int m = 0; m < M; m++) {
                for (int l = 0; l < Q-m; l++) {
                    int tl = l-m;
                    if (tl < 0 || tl > Q-1) continue;
                    logxrd[l] += W_non_log(tl, k) * H[m](k,n);
                }
            }
            string strB("R"+to_string(n+1)+"_"+to_string(k+1));// Ri_j denotes the signals of phase j at sample point i
            DataProcess::save_array_equal(logxrd, Q, strB, os);
        }
    }

    os <<  "\n\n// Per-phase shift for each sample\n";
    for (int n = 0; n < N; n++) {
        string strB("S"+to_string(n+1));
        for (int k = 0; k < K; k++) {
            B[k] = phase_mixtures[k][n].shift;
        }
        DataProcess::save_array_equal(B,  K, strB, os);
    }

    // H tensor
    os <<  "\n\n// H tensor\n" ;
    for (int n = 0; n < N; n++) {
        for (int k = 0; k < K; k++) {
            for (int m = 0; m < M; m++) {
                B[m] = H[m](k,n);
            }
            string strB("H[*]("+to_string(k+1)+","+to_string(n+1)+")");
            DataProcess::save_array_equal(B, M, strB, os);
        }
    }

    os << "\n\n// Per-sample contribution (L1 loss)\n";
    double sampleLoss[N+1];
    double sampleTot[N+1];
    for (int n = 0; n < N; n++) {
        double loss = 0.0;
        double tot = 0.0;
        for (int l = 0; l < Q; l++) {
            double d = A(l,n), dr = R(l,n);
            loss += abs(d-dr);
            tot += abs(d);
            //if (n == N-1) fprintf(stderr, "%lf,", d);
        }
        //if (n == N-1) fprintf(stderr, "\n");
        sampleLoss[n] = loss;
        sampleTot[n] = tot;
    }
    DataProcess::save_array_equal(sampleLoss, N, "L", os);

    for (int n = 0; n < N; n++) {
        B[n] = sampleLoss[n]/sampleTot[n];
    }
    DataProcess::save_array_equal(B, N, "L_proportion", os);

    os << "\n\n";

    vector<double> Qlogvec;
    for (int l = 0; l < Q; l++)
        Qlogvec.push_back(Qlog[l]);
    vector<double> Qvec;
    for (int l = 0; l < Q; l++)
        Qvec.push_back(I.Q[l]);
    GaussianCompare gc(cfg.match_shift, cfg.match_sigma);
    gc.compare(os, (char*)cfg.sticksfile.c_str(), W, H, K, Qvec, Qlogvec);

}


















