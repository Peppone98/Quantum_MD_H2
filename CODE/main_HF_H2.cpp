/*Program for the ground state of H2*/

#include "definitions.h"

using namespace std;


int main (){
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
	gsl_vector_set_all(c, 1.);
	
	create_S(S, R_A, R_B);
    gsl_matrix_memcpy(S_auxiliary, S);
	diag_S(S_auxiliary, V);
    one_body_H(H, R_A, R_B);
    build_Q(Q, R_A, R_B);

    for(n=0; n<n_iter; n++){
        two_body_F(Q, c, U, F, R_A, R_B);
        gsl_matrix_add(F, H);
        E_1s = solve_FC_eSC(F, V, U);
        cout << E_1s << endl;
        gsl_matrix_get_col(c, U, 0);
    }
	cout << "   " << endl;
	
    for(n=0; n<2*N; n++){
    	cout << gsl_vector_get(c, n) << endl;
    }
	print_orbital(c, R_A, R_B);
    cout << compute_E0(Q, c, H, R_A, R_B) << endl;
}


/******** EVOLUTION OF THE COEFFICIENTS ********/


/*forse si può fare tutto con delle funzioni nella libreria gsl...*/
void partial_evolution(gsl_matrix *F, gsl_vector *c, gsl_vector *c_old){
    double c_r_new=0., c_r=0., c_r_old=0., F_rs=0., c_s=0.;
    int r, s;
    for(r=0; r<2*N; r++){
        c_r = gsl_vector_get(c, r);
        c_r_old = gsl_vector_get(c_old, r);
        c_r_new = 2.*c_r - c_r_old;
        gsl_vector_set(c, r, c_r_new);
        for(s=0; s<2*N; s++){
            F_rs = gsl_matrix_get(F, r, s);
            c_s = gsl_vector_get(c, s);
            c_r_new -= h*h*F_rs*c_s;
        }
        gsl_vector_set(c, r, c_r_new);
        gsl_vector_set(c_old, r, c_r);
    }
}

void complete_evolution(gsl_matrix *S, gsl_vector *c, gsl_vector *c_old){
    double lambda = 0.;
    double c_r_new=0., c_r=0., c_r_old=0., S_rs=0., c_s=0.;
    lambda = get_lambda(S, c, c_old);
    for(r=0; r<2*N; r++){
        c_r = gsl_vector_get(c, r);
        gsl_vector_set(c, r, c_r_new);
        for(s=0; s<2*N; s++){
            S_rs = gsl_matrix_get(S, r, s);
            c_s = gsl_vector_get(c, s);
            c_r_new -= h*h*F_rs*c_s;
        }
        gsl_vector_set(c, r, c_r_new);
        gsl_vector_set(c_old, r, c_r);
    }
}

double solve_eq2degree(double a, double b, double c){ 
    return -(b + sqrt(b*b - 4.*a*c))/(2.*a);
}


double get_lambda(gsl_matrix *S, gsl_vector *c, gsl_vector *c_old){
    double A=0., B=0., C=0.;
    double c_r=0., c_s=0., c_t=0., c_u=0.;
    double S_rs=0., S_tu=0., S_rt=0.;
    double c_s_old=0., c_u_old=0.;
    int r, s, t, u;
    for(r=0; r<2*N; r++){
        for(s=0; s<2*N; s++){
            c_r = gsl_vector_get(c, r);
            c_s = gsl_vector_get(c, s);
            S_rs = gsl_matrix_get(S, r, s);
            C += c_r*c_s*S_rs;
            for(t=0; t<2*N; t++){
                S_rt = gsl_matrix_get(S, r, t);
                c_s_old = gsl_vector_get(c_old, s);
                B -= 2.*S_rs*c_r*S_rt*c_s_old;
                for(u=0; u<2*N; u++){
                    S_tu = gsl_matrix_get(S, t, u);
                    c_u_old = gsl_vector_get(c_old, u);
                    A += c_u_old*c_s_old*S_tu*S_rs*S_rt;
                }
            }

        }
    }
    C = C - 1.;
    return solve_eq2degree(A, B, C); 
}

