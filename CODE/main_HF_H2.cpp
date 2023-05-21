#include "definitions.h"

using namespace std;

int main (int argc, char *argv[]){

    /**** Create space for matrices and vectors ****/
	int n, n_iter = 20;
	R R_A, R_B;
    double Q[2*N][2*N][2*N][2*N];
    double X[CP_iter];
    double E_1s=0., lambda=0., E0 = 0., dE0_dX = 0.;
	gsl_matrix *S = gsl_matrix_alloc(2*N, 2*N);
    gsl_matrix *S_auxiliary = gsl_matrix_alloc(2*N, 2*N);
	gsl_matrix *H = gsl_matrix_alloc(2*N, 2*N);
	gsl_matrix *F = gsl_matrix_alloc(2*N, 2*N);
    gsl_matrix *U = gsl_matrix_calloc(2*N, 2*N); 
    gsl_matrix *V = gsl_matrix_calloc(2*N, 2*N);
	gsl_vector *c = gsl_vector_alloc(2*N);
    gsl_vector *c_old = gsl_vector_alloc(2*N);
	gsl_vector_set_all(c, 1.);

    /**** Initial nuclei positions ****/
    R_A.x = 0., R_A.y = 0., R_A.z = 0.;
    R_B.x = 1.0, R_B.y = 0., R_B.z = 0.;

    /**** Fill S and compute V ****/
    create_S(S, R_A, R_B);
    gsl_matrix_memcpy(S_auxiliary, S);
	diag_S(S_auxiliary, V);

    /**** Fill H and Q ****/
    one_body_H(H, R_A, R_B);
    build_Q(Q, R_A, R_B);

    /**** Initial conditions ****/
    X[0] = sqrt(scalar_prod(R_A, R_B));
    X[1] = X[0];
    double norm = normalization(c, S);
    gsl_vector_memcpy(c_old, c);
	
    
    /**** 1) SELF-CONSISTENT HARTREE-FOCK ****/
    if (string(argv[1]) == "SC_Hartree_Fock"){
        cout << "Hartree-Fock energies: " << endl;

        for(n=0; n<n_iter; n++){

            /**** Build Fock matrix using Q and H ****/
            two_body_F(Q, c, F);
            gsl_matrix_add(F, H);

            /**** tmp_F is needed beacuse solve_FC_eSC destroys lower part of F ****/
            gsl_matrix *tmp_F = gsl_matrix_alloc(2*N, 2*N);
            gsl_matrix_memcpy(tmp_F, F);
            E_1s = solve_FC_eSC(tmp_F, V, U);
            cout << E_1s << endl;

            /**** Get the c vector from eigenvector matrix U ****/
            gsl_matrix_get_col(c, U, 0);
        }

        /**** Print the coefficients ****/
        cout << "   " << endl;
        cout << "Vector of coefficients" << endl;
        for(n=0; n<2*N; n++){
            cout << gsl_vector_get(c, n) << endl;
        }
        
        /**** Print orbital and HF energy ****/
        print_orbital(c, R_A, R_B);
        cout << "   " << endl;
        cout << "Total HF energy" << endl;
        cout << compute_E0(F, H, c) + 1./X[0] << endl; 	
	}


    /**** 2) EVOLUTION OF THE COEFFICIENTS C(t) AT FIXED X ****/
    if(string(argv[1]) == "Evolve_coefficients"){
        ofstream myfile;
	    myfile.open("Energies_C_evolution.txt", ios :: out | ios :: trunc);

        /**** MD cycle: compute new c vector and keep F updated ****/
        for(n=0; n<50; n++){
            two_body_F(Q, c, F);
            gsl_matrix_add(F, H);
            myfile << compute_E0(F, H, c) + 1./X[0] << endl;
            lambda = update_c(F, S, c, c_old);
        }
        myfile.close();

        /**** Print the coefficients ****/
        cout << "Vector of coefficients: " << endl;
        for(n=0; n<2*N; n++){
            cout << gsl_vector_get(c, n) << endl;
        }
    }


    /**** 3) CAR PARRINELLO MOLECULAR DYNAMICS ****/
    if(string(argv[1]) == "MD_Car_Parrinello"){

        ofstream myfile;
        ofstream coeff_file;
	    myfile.open("CP_X_energies.txt", ios :: out | ios :: trunc);
        coeff_file.open("CP_coeff.txt", ios :: out | ios :: trunc);

        for(n=1; n<CP_iter-1; n++){

            /**** Fill S, H and Q for electronic problem ****/
            create_S(S, R_A, R_B);
            one_body_H(H, R_A, R_B);
            build_Q(Q, R_A, R_B);

            /**** Evolve the c vector of coefficients ****/
            double n_times = h_N/h;
            for(int k=0; k<n_times; k++){
                two_body_F(Q, c, F);
                gsl_matrix_add(F, H);
                lambda = update_c(F, S, c, c_old);
            }

            /**** Print the coefficients to file ****/
            for(int i=0; i<2*N; i++){
                coeff_file << gsl_vector_get(c, i) <<  "    ";
            }
            coeff_file << endl;

            /**** Compute the correct energy after the update of F ****/
            two_body_F(Q, c, F);
            gsl_matrix_add(F, H);
            E0 = compute_E0(F, H, c) + 1./X[n];

            /**** Fill S, H, Q and F for nuclear problem ****/
            create_dS_dX(S, R_A, R_B);
            one_body_dH_dX(H, R_A, R_B, X[n]);
            build_dQ_dX(Q, R_A, R_B, X[n]);
            two_body_F(Q, c, F);
            gsl_matrix_add(F, H);

            /**** Update X (also R_B.x is updated) ****/
            dE0_dX = compute_dE0_dX(F, H, c, X[n]);
            X[n + 1] = evolve_X(c, S, &R_B, lambda, dE0_dX, X[n], X[n-1]);
            cout << "  " << endl;
            myfile << X[n] << "    " << E0 << endl;

        }

        myfile.close();
        coeff_file.close();
    }


    /**** 4) CONJUGATE GRADIENT (SINGLE RUN) ****/
    if(string(argv[1]) == "CG_min"){
        /**** Start with little equilibration cycle ****/
        for(n=0; n<25; n++){
            two_body_F(Q, c, F);
            gsl_matrix_add(F, H);
            lambda = update_c(F, S, c, c_old);
        }

        double lambda_CP = 0.0, E_final = 0.0, E_initial = 0.0;

        /**** Fill the Hessian and b (minus the gradient of E) ****/
        gsl_matrix *Hessian = gsl_matrix_alloc(2*N, 2*N);
        gsl_vector *b = gsl_vector_alloc(2*N);
        Get_Hessian_and_b(Hessian, b, Q, S, F, c);

        /**** Compute the energy before minimisation ****/
        E_initial = compute_E0(F, H, c);

        /**** Apply the CG routine to get the increment Delta_c ****/
        gsl_vector *Delta_c = gsl_vector_alloc(2*N);
        gsl_vector_set_all(Delta_c, 0.0);
        Conj_grad(Hessian, b, Delta_c, 0.001);

        cout << "Vector of Delta_C_cg: " << endl;
        for(n=0; n<2*N; n++){
            cout << gsl_vector_get(Delta_c, n) << endl;
        }
        
        /**** Copy c in c_old and get the new vector ****/
        gsl_vector_memcpy(c_old, c);
        gsl_vector_add(c, Delta_c);

        cout << "Vector of coefficients resulting from CG alone: " << endl;
        for(n=0; n<2*N; n++){
            cout << gsl_vector_get(c, n) << endl;
        }
        
        /**** Update the c (in CP fashion) ****/
        lambda_CP = Get_lambda_CP(S, c, c_old, lambda);
        E_final = compute_E0(F, H, c);
        norm = Get_norm_C_cg(S, c);

        cout << "Initial energy: " << E_initial << endl;
        cout << "Energy difference (should be negative!): " << E_final - E_initial << endl;
        cout << "Normalisation: " << norm << endl;
        cout << "Gap between lambdas: " << lambda_CP - lambda << endl;

    }

    
    /**** 4) CONJUGATE GRADIENT - NUCLEI MOTION WITH CP ****/
    if(string(argv[1]) == "CG_CP"){

        /**** Start with little equilibration cycle ****/
        for(n=0; n<35; n++){
            two_body_F(Q, c, F);
            gsl_matrix_add(F, H);
            lambda = update_c(F, S, c, c_old);
        }
        
        /**** The first lambda_CP comes from the CP equilibration ****/
        double lambda_CP = lambda;
        double E_final = 0.0, E_initial = 0.0;

        /**** Open files to write the coefficients ****/
        ofstream myfile;
        ofstream coeff_file;
	    myfile.open("CP_CG_X.txt", ios :: out | ios :: trunc);
        coeff_file.open("CP_CG_coeff.txt", ios :: out | ios :: trunc);

        gsl_matrix *Hessian = gsl_matrix_alloc(2*N, 2*N);
        gsl_vector *b = gsl_vector_alloc(2*N);
        gsl_vector *Delta_c = gsl_vector_alloc(2*N);

        for(n=1; n<CP_iter-1; n++){

            /**** Print the coefficients to file ****/
            for(int i=0; i<2*N; i++){
                coeff_file << gsl_vector_get(c, i) <<  "    ";
            }
            coeff_file << endl;

            /**** Compute the energy before minimisation ****/
            E_initial = compute_E0(F, H, c);

            /**** Fill the Hessian and b (minus the gradient of E) ****/
            Get_Hessian_and_b(Hessian, b, Q, S, F, c);

            /**** Apply the CG routine the increment Delta_c ****/
            gsl_vector_set_all(Delta_c, 0.0);
            Conj_grad(Hessian, b, Delta_c, 0.001);

            /**** Copy c in c_old and then update c ****/
            gsl_vector_memcpy(c_old, c);
            gsl_vector_add(c, Delta_c);
            E_final = compute_E0(F, H, c);

            /**** Update the c (in CP fashion) ****/
            lambda_CP = Get_lambda_CP(S, c, c_old, lambda_CP);
            norm = Get_norm_C_cg(S, c);

            /**** Compute the correct energy after the update of F ****/
            two_body_F(Q, c, F);
            gsl_matrix_add(F, H);
            E0 = compute_E0(F, H, c) + 1./X[n];

            /**** Fill S, H, Q and F for nuclear problem ****/
            create_dS_dX(S, R_A, R_B);
            one_body_dH_dX(H, R_A, R_B, X[n]);
            build_dQ_dX(Q, R_A, R_B, X[n]);
            two_body_F(Q, c, F);
            gsl_matrix_add(F, H);

            /**** Update X using the newly computed lambda_CP ****/
            dE0_dX = compute_dE0_dX(F, H, c, X[n]);
            X[n + 1] = evolve_X(c, S, &R_B, lambda_CP, dE0_dX, X[n], X[n-1]);
            cout << "  " << endl;
            myfile << X[n] << "    " << E0 << endl;
            cout << "lambda_CP: " << lambda_CP << endl;
            cout << "normalisation: " << norm << endl;
            cout << "Energy difference (should be negative!): " << E_final - E_initial << endl;

            /**** Prepare the S, H, Q and F for new electronic problem ****/
            create_S(S, R_A, R_B);
            one_body_H(H, R_A, R_B);
            build_Q(Q, R_A, R_B);
            two_body_F(Q, c, F);
            gsl_matrix_add(F, H);
        }

        myfile.close();
        coeff_file.close();
    }

    /**** 6) CG AND CP SUPERPOSITION ****/
    if(string(argv[1]) == "CG_CP_superposition"){
        int i, j;
        ofstream scal_prod_file;
        scal_prod_file.open("scal_prod.txt");

        /**** Open the file stream ****/ 
        ifstream CG_file;
        ifstream CP_file;
        CG_file.open("CP_CG_coeff.txt");
        CP_file.open("CP_coeff.txt");

        /**** Check if opening a file failed ****/ 
        if (CG_file.fail()) {
            cout << "Error opening the coefficients file" << endl;
            exit(1);
        }else{
            double c_CG[2*N], c_CP[2*N];
            double scal_prod = 0.;
            string line;

            /**** Read the lines of the two files ****/
            while(getline(CG_file, line)){
                for(i=0; i<2*N; i++){
                    CG_file >> c_CG[i];
                    CP_file >> c_CP[i];
                }

                /**** Compute the scalar product using S ****/
                scal_prod = 0.;
                for(i=0; i<2*N; i++){
                    for(j=0; j<2*N; j++){
                        scal_prod += gsl_matrix_get(S, i, j)*c_CP[i]*c_CG[j];
                    }
                }
                scal_prod_file << scal_prod << endl;
            }
            cout << "File scal_prod.txt successfully written" << endl;
        }
        CG_file.close();
        CP_file.close();
        scal_prod_file.close();
    }
    
}




