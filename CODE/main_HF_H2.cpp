#include "definitions.h"

using namespace std;

int main (int argc, char *argv[]){

    /**** BUILD VARIABLES ****/
	int n, n_iter = 20;
	R R_A, R_B;
    double Q[2*N][2*N][2*N][2*N];
    R_A.x = 0., R_A.y = 0., R_A.z = 0.;
    R_B.x = 1., R_B.y = 0., R_B.z = 0.;
	double E_1s=0.;
	gsl_matrix *S = gsl_matrix_alloc(2*N, 2*N);
    gsl_matrix *S_auxiliary = gsl_matrix_alloc(2*N, 2*N);
	gsl_matrix *H = gsl_matrix_alloc(2*N, 2*N);
	gsl_matrix *F = gsl_matrix_alloc(2*N, 2*N);
    gsl_matrix *U = gsl_matrix_calloc(2*N, 2*N); 
    gsl_matrix *V = gsl_matrix_calloc(2*N, 2*N);
	gsl_vector *c = gsl_vector_alloc(2*N);
    gsl_vector *c_old = gsl_vector_alloc(2*N);
	gsl_vector_set_all(c, 1.);

    
    /**** FILL THE MATRICES ****/
    create_S(S, R_A, R_B);
    gsl_matrix_memcpy(S_auxiliary, S);
	diag_S(S_auxiliary, V);
    one_body_H(H, R_A, R_B);
    build_Q(Q, R_A, R_B);
	
    
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



    /**** CAR PARRINELLO MOLECULAR DYNAMICS ****/
    if(string(argv[1]) == "MD_Car_Parrinello"){
        gsl_vector_set_all(c, 1.);
        double h = 0.1;

        /**** impose normalization & initial conditon ****/
        normalization(c, S);
        gsl_vector_memcpy(c_old, c);

        /**** MD cycle ****/
        ofstream myfile;
	    myfile.open("CP_energies.txt", ios :: out | ios :: trunc);
        for(n=0; n<200; n++){
            two_body_F(Q, c, F, R_A, R_B);
            gsl_matrix_add(F, H);
            update_c(F, S, c, c_old);
            myfile << compute_E0(Q, c, H, R_A, R_B) << endl;
        }
        myfile.close();

        cout << "Vector of coefficients: " << endl;
        
        for(n=0; n<2*N; n++){
            cout << gsl_vector_get(c, n) << endl;
        }
    }

}




