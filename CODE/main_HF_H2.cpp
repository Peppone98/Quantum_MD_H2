#include "definitions.h"

using namespace std;

int main (int argc, char *argv[]){

    /**** BUILD VARIABLES ****/
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
	
    
    /**** SELF-CONSISTENT HARTREE-FOCK ****/
    if (string(argv[1]) == "SC_Hartree_Fock"){
        cout << "Hartree-Fock energies" << endl;
        for(n=0; n<n_iter; n++){
            two_body_F(Q, c, F);
            gsl_matrix_add(F, H);

            /**** tmp matrix is needed beacuse eval_evec destroys lower part of F ****/
            gsl_matrix *tmp = gsl_matrix_alloc(2*N, 2*N);
            gsl_matrix_memcpy(tmp, F);
            E_1s = solve_FC_eSC(tmp, V, U);
            cout << E_1s << endl;
            gsl_matrix_get_col(c, U, 0);
        }
        cout << "   " << endl;
        cout << "Vector of coefficients" << endl;
        for(n=0; n<2*N; n++){
            cout << gsl_vector_get(c, n) << endl;
        }

        print_orbital(c, R_A, R_B);
        cout << "   " << endl;
        cout << "Total HF energy" << endl;
        cout << compute_E0(F, H, c) + 1./X[0] << endl; 	
	}


    /**** EVOLUTION OF THE COEFFICIENTS C(t) ****/
    if(string(argv[1]) == "Evolve_coefficients"){
        ofstream myfile;
	    myfile.open("Energies_C_evolution.txt", ios :: out | ios :: trunc);

        /**** MD cycle ****/
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


    /**** CAR PARRINELLO MOLECULAR DYNAMICS ****/
    if(string(argv[1]) == "MD_Car_Parrinello"){

        ofstream myfile;
	    myfile.open("CP_X_energies.txt", ios :: out | ios :: trunc);

        /**** Loop for internuclear distance update ****/
        for(n=1; n<CP_iter-1; n++){

            /**** Fill S, H and Q for electronic problem ****/
            create_S(S, R_A, R_B);
            one_body_H(H, R_A, R_B);
            build_Q(Q, R_A, R_B);

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

            /**** Also R_B.x is updated ****/
            dE0_dX = compute_dE0_dX(F, H, c, X[n]);
            X[n + 1] = evolve_X(c, S, &R_B, lambda, dE0_dX, X[n], X[n-1]);
            cout << "lambda: " << lambda << endl;
            cout << "dE0_dX: " << dE0_dX << endl;
            cout << "  " << endl;
            myfile << X[n] << "    " << E0 << endl;

        }

        myfile.close();
    }


    if(string(argv[1]) == "Prova"){
        
        for(int i=1; i<3; i++){
            cout << "  " << endl;
            cout << "Nuclear step number: " << i << endl;
            /**** Get first guess of lambda with a little equilibration ****/
            int first_eq = 43;
            /**** Rebuild the matrices for electronic problem ****/
            create_S(S, R_A, R_B);
            one_body_H(H, R_A, R_B);
            build_Q(Q, R_A, R_B);
            cout << setprecision(11) << fixed << "New position: " << R_B.x << endl;
            cout << setprecision(11) << fixed << "Check S[3][5]: " << gsl_matrix_get(S, 3, 5) << endl;
            cout << setprecision(11) << fixed << "Check H[3][5]: " << gsl_matrix_get(H, 3, 5) << endl;
            cout << setprecision(11) << fixed << "Check Q[1][2][3][5]: " << Q[1][2][3][5] << endl;


            /**** First cycle for the c ****/
            for(int k=0; k<first_eq; k++){
                two_body_F(Q, c, F);
                gsl_matrix_add(F, H);
                cout << setprecision(8) << fixed << "Prova: " << compute_E0(F, H, c) + 1./X[i] << endl;
                lambda = update_c(F, S, c, c_old);
            }

            cout << setprecision(11) << fixed << "Check F[3][5]: " << gsl_matrix_get(F, 3, 5) << endl;

            cout << "  " << endl;
            cout << "Hartree Fock energy: " << setprecision(8) << fixed << compute_E0(F, H, c) + 1./X[i] << endl;
            cout << "Internuclear distance: " << X[i] << endl;
            cout << "  " << endl;

            /**** Print the coefficients ****/
            cout << "Vector of coefficients: " << endl;
            for(n=0; n<2*N; n++){
                cout << gsl_vector_get(c, n) << setprecision(11) << fixed << endl;
            }

            create_dS_dX(S, R_A, R_B);
            one_body_dH_dX(H, R_A, R_B, X[i]);
            build_dQ_dX(Q, R_A, R_B, X[i]);

            cout << setprecision(11) << fixed << "Check dS_dX[3][5]: " << gsl_matrix_get(S, 3, 5) << endl;
            cout << setprecision(11) << fixed << "Check dH_dX[3][5]: " << gsl_matrix_get(H, 3, 5) << endl;
            cout << setprecision(11) << fixed << "Check dH_dX[4][4]: " << gsl_matrix_get(H, 4, 4) << endl;
            cout << setprecision(11) << fixed << "Check dQ_dX[1][2][3][5]: " << Q[1][2][3][5] << endl;

            two_body_F(Q, c, F);
            gsl_matrix_add(F, H);
            dE0_dX = compute_dE0_dX(F, H, c, X[i]);

            cout << setprecision(11) << fixed << "Check dF_dX[3][5]: " << gsl_matrix_get(F, 3, 5) << endl;

            cout << "lambda: " << lambda << endl;
            cout << "dE0_dX: " << dE0_dX << endl;

            /* Here also the R_B.x evolves */
            X[i+1] = evolve_X(c, S, &R_B, lambda, dE0_dX, X[i], X[i-1]);
            cout << "Distance: " << setprecision(5) << fixed << X[i+1] << "    " << R_B.x << endl;
        }
    }

    

    if(string(argv[1]) == "Prova2"){
        R_B.x = 1.0016275;
        gsl_vector *prova = gsl_vector_alloc(2*N);
        gsl_vector_set(prova, 0, 0.09627628);
        gsl_vector_set(prova, 1, 0.17769337);
        gsl_vector_set(prova, 2, 0.12480673);
        gsl_vector_set(prova, 3, 0.01775755);
        gsl_vector_set(prova, 4, 0.09627628);
        gsl_vector_set(prova, 5, 0.17769337);
        gsl_vector_set(prova, 6, 0.12480673);
        gsl_vector_set(prova, 7, 0.01775755);

        /**** Get first guess of lambda with a little equilibration ****/
        int first_eq = 43;
        /**** Rebuild the matrices for electronic problem ****/
        create_S(S, R_A, R_B);
        one_body_H(H, R_A, R_B);
        build_Q(Q, R_A, R_B);

        /**** Correct normalization and update c(t+h) ****/
        /* norm = normalization(c_old, S); */
        norm = normalization(prova, S);
        two_body_F(Q, prova, F);

        /**** First cycle for the c ****/
        for(int k=0; k<first_eq; k++){
            two_body_F(Q, prova, F);
            gsl_matrix_add(F, H);
            lambda = update_c(F, S, prova, c_old);
            cout << setprecision(8) << fixed << "Prova: " << compute_E0(F, H, c) << endl;
        }
    }
}




