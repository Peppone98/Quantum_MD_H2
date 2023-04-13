#include "definitions.h"

using namespace std;

int main (int argc, char *argv[]){

    /**** BUILD VARIABLES ****/
	int n, n_iter = 20;
	R R_A, R_B;
    double Q[2*N][2*N][2*N][2*N];
    double X[CP_iter];
    double E_1s=0., lambda=0.;
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
    R_B.x = 1., R_B.y = 0., R_B.z = 0.;

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
            two_body_F(Q, c, F, R_A, R_B);
            gsl_matrix_add(F, H);
            E_1s = solve_FC_eSC(F, V, U);
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
        cout << compute_E0(Q, c, H, R_A, R_B) << endl;
        	
	}


    /**** EVOLUTION OF THE COEFFICIENTS C(t) ****/
    if(string(argv[1]) == "Evolve_coefficients"){
        ofstream myfile;
	    myfile.open("CP_energies.txt", ios :: out | ios :: trunc);

        /**** MD cycle ****/
        for(n=0; n<100; n++){
            two_body_F(Q, c, F, R_A, R_B);
            gsl_matrix_add(F, H);
            lambda = update_c(F, S, c, c_old);
            myfile << compute_E0(Q, c, H, R_A, R_B) << endl;
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

        /**** HF single iteration (S, H, Q have already been built) ****/
        two_body_F(Q, c, F, R_A, R_B);
        gsl_matrix_add(F, H);
        lambda = update_c(F, S, c, c_old);
        cout << "Initial Hartree Fock energy: " << compute_E0(Q, c, H, R_A, R_B) << endl;
        cout << "Initial internuclear distance: " << X[0] << endl;
        cout << "  " << endl;

        /**** Loop for internuclear distance update ****/
        for(n=1; n<CP_iter-1; n++){

            /**** Fill S, H and Q with the derivatives w.r.t. X ****/
            create_dS_dX(S, R_A, R_B);
            one_body_dH_dX(H, R_A, R_B, X[n]);
            build_dQ_dX(Q, R_A, R_B, X[n]);

            /**** Also R_B.x is updated ****/
            double dE0_dX = compute_dE0_dX(Q, c, H, R_A, R_B);
            X[n + 1] = evolve_X(Q, c, H, S, &R_B, lambda, dE0_dX, X[n], X[n-1]);

            /**** Rebuild the matrices for electronic problem ****/
            create_S(S, R_A, R_B);
            one_body_H(H, R_A, R_B);
            build_Q(Q, R_A, R_B);

            /**** Correct normalization and update c(t+h) ****/
            norm = normalization(c, S);
            two_body_F(Q, c, F, R_A, R_B);
    
            /**** MD cycle for the c(t) ****/
            double n_times = h_N/h;
            for(int k=0; k<n_times; k++){
                two_body_F(Q, c, F, R_A, R_B);
                gsl_matrix_add(F, H);
                lambda = update_c(F, S, c, c_old);
            }

            cout << R_B.x << "       " <<X[n + 1] << "    " << compute_E0(Q, c, H, R_A, R_B) << endl;

        }
    }


    if(string(argv[1]) == "Prova"){
        two_body_F(Q, c, F, R_A, R_B);
        gsl_matrix_add(F, H);
        lambda = update_c(F, S, c, c_old);
        cout << "Initial Hartree Fock energy: " << compute_E0(Q, c, H, R_A, R_B) << endl;
        cout << "Initial internuclear distance: " << X[0] << endl;
        cout << "  " << endl;
    }
    

}




