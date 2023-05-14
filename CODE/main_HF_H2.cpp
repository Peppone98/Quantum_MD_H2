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
	    myfile.open("CP_X_energies.txt", ios :: out | ios :: trunc);

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
    }


    /**** 4) CONJUGATE GRADIENT (SINGLE RUN) ****/
    if(string(argv[1]) == "CG_min"){
        /**** Start with little equilibration cycle ****/
        for(n=0; n<25; n++){
            two_body_F(Q, c, F);
            gsl_matrix_add(F, H);
            lambda = update_c(F, S, c, c_old);
        }

        /**** Fill the Hessian and b (minus the gradient of E) ****/
        gsl_matrix *Hessian = gsl_matrix_alloc(2*N, 2*N);
        gsl_vector *b = gsl_vector_alloc(2*N);
        Get_Hessian_and_b(Hessian, b, Q, S, F, c);

        double E_initial = compute_E0(F, H, c);

        /**** Define the increment Delta_c and then apply the CG routine ****/
        gsl_vector *Delta_c = gsl_vector_alloc(2*N);
        gsl_vector_set_all(Delta_c, 0.0);
        Conj_grad(Hessian, b, Delta_c, 0.001);
        gsl_vector_add(c, Delta_c);

        cout << "Vector of coefficients: " << endl;
        for(n=0; n<2*N; n++){
            cout << gsl_vector_get(c, n) << endl;
        }

        double E_final = compute_E0(F, H, c);
        cout << E_final - E_initial << endl;

    }

    
    /**** 4) CONJUGATE GRADIENT - CAR PARRINELLO ****/
    if(string(argv[1]) == "CG_CP"){

        /**** Start with little equilibration cycle ****/
        for(n=0; n<25; n++){
            two_body_F(Q, c, F);
            gsl_matrix_add(F, H);
            lambda = update_c(F, S, c, c_old);
        }

        ofstream myfile;
	    myfile.open("CP_CG_X.txt", ios :: out | ios :: trunc);

        for(n=1; n<CP_iter-1; n++){

            /**** Fill the Hessian and b (minus the gradient of E) ****/
            gsl_matrix *Hessian = gsl_matrix_alloc(2*N, 2*N);
            gsl_vector *b = gsl_vector_alloc(2*N);
            Get_Hessian_and_b(Hessian, b, Q, S, F, c);

            /**** Define the increment Delta_c and then apply the CG routine ****/
            gsl_vector *Delta_c = gsl_vector_alloc(2*N);
            gsl_vector_set_all(Delta_c, 0.0);
            Conj_grad(Hessian, b, Delta_c, 0.001);
            gsl_vector_add(c, Delta_c);

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

            /**** Prepare the S, H, Q and F for new electronic problem ****/
            create_S(S, R_A, R_B);
            one_body_H(H, R_A, R_B);
            build_Q(Q, R_A, R_B);
            two_body_F(Q, c, F);
            gsl_matrix_add(F, H);
        }

        myfile.close();
    }


}




