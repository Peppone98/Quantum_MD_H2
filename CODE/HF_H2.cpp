/*Program for the ground state of H2*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <string>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

using namespace std;

struct R {
    double x;
    double y;
    double z;
};

/******* ROUTINES FOR SOLVING ROOTHAN EQUATION ****/
void diag_S(gsl_matrix *S, gsl_matrix *U);
void create_eval_evec(gsl_matrix *A, gsl_matrix *B, gsl_vector *C);
void create_V(gsl_matrix *U, gsl_vector *L);
void multiply(gsl_matrix *A, gsl_matrix *B, gsl_matrix *C);
double solve_FC_eSC(gsl_matrix *F, gsl_matrix *V, gsl_matrix *U);

/******** ROUTINES FOR MATRIX ELEMENTS *****/
double scalar_prod(R R_A, R R_B);
double K(double alpha, double beta, R R_A, R R_B);
R R_weighted(double alpha, double beta, R R_A, R R_B);
double overlap(double alpha, double beta, R R_A, R R_B);
double laplacian(double alpha, double beta, R R_A, R R_B);
double el_nucl(double alpha, double beta, R R_A, R R_B, R R_C);
double F0(double x);
double direct_term(double alpha, double beta, R R_A, R R_B, double alpha_prime, double beta_prime, R R_A_prime, R R_B_prime);

/******** ROUTINES FOR FILLING MATRICES ******/
void create_S(gsl_matrix *S, R R_A, R R_B);
void one_body_H(gsl_matrix *H, R R_A, R R_B);
void build_Q(R R_A, R R_B);
void two_body_F(gsl_vector *c, gsl_matrix *U, gsl_matrix *F, R R_A, R R_B, double gamma);
void scale(gsl_matrix *U, gsl_vector *c, double gamma);

/******** COMPUTE TOTAL ENERGY AND PRINT SOLUTION ****/
double compute_E0(gsl_vector *c, double E_1s);
void print_orbital(gsl_vector *c_new, R R_A, R R_B);
double normalization(gsl_vector *c, gsl_matrix *S);


const double pi = 3.1415926;
const int N = 4;
const double a[N] = {13.00773, 1.962079, 0.444529, 0.1219492};
double Q[2*N][2*N][2*N][2*N];

int main (){
	int n, n_iter = 200;
	R R_A, R_B;
    R_A.x = 0., R_A.y = 0., R_A.z = 0.;
    R_B.x = 1., R_B.y = 0., R_B.z = 0.;
	double gamma = 0.05, E_1s=0., norm;
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
    build_Q(R_A, R_B);

   	for(n=0; n<n_iter; n++){
		norm = normalization(c, S);
		two_body_F(c, U, F, R_A, R_B, gamma);
		gsl_matrix_add(F, H);
		E_1s = solve_FC_eSC(F, V, U);
        cout << E_1s << "   "  << norm << endl;
	}
    
	cout << "   " << endl;


    for(int i=0; i<2*N; i++){
        for(int j=0; j<2*N; j++){
            cout << gsl_matrix_get(F, i, j) << "    ";
        }
        cout << endl;
    }

	
    for(n=0; n<2*N; n++){
    	cout << gsl_vector_get(c, n) << endl;
    }
	print_orbital(c, R_A, R_B);
}


/******* ROUTINES FOR SOLVING ROOTHAN EQUATION ****/

void diag_S(gsl_matrix *S, gsl_matrix *U){
	gsl_vector *L = gsl_vector_alloc(2*N);
	create_eval_evec(S, U, L); 
	create_V(U, L);
	gsl_vector_free(L);
}

void create_V(gsl_matrix *U, gsl_vector *L){
	int i;
	double x=0.;
	for(i=0; i<2*N; i++){
		x = gsl_vector_get(L, i);
		gsl_vector_set(L, i, 1./sqrt(x)); 
	}
	gsl_matrix_scale_columns(U, L); 
}

/*The diagonal and lower triangular part of A are 
destroyed during the computation, but the strict upper 
triangular part is not referenced*/
void create_eval_evec(gsl_matrix *A, gsl_matrix *B, gsl_vector *C){
    gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(2*N); 
	gsl_eigen_symmv(A, C, B, w); 
	gsl_eigen_symmv_free(w); 
	gsl_eigen_symmv_sort(C, B, GSL_EIGEN_SORT_VAL_ASC);
}


void multiply(gsl_matrix *A, gsl_matrix *B, gsl_matrix *C){
	int i, j, k;
	double val, f1, f2; 
	for(i=0; i<2*N; i++){
		for(j=0; j<2*N; j++){
			for(k=0; k<2*N; k++){
				f1 = gsl_matrix_get(A, i, k);
				f2 = gsl_matrix_get(B, k, j);
				val += f1*f2;
			}				
			gsl_matrix_set(C, i, j, val);
			val = 0.;
		}
	}
}


double solve_FC_eSC(gsl_matrix *F, gsl_matrix *V, gsl_matrix *U){
	double val=0.; 
	gsl_matrix_memcpy(U, V);
	gsl_matrix *tmp = gsl_matrix_alloc(2*N, 2*N);
	gsl_vector *L = gsl_vector_alloc(2*N);
	multiply(F, U, tmp); 
	gsl_matrix_transpose(U); 
	multiply(U, tmp, F); 
	create_eval_evec(F, tmp, L); 
	gsl_matrix_transpose(U); 
	multiply(U, tmp, U); 
	val = gsl_vector_get(L, 0);
	gsl_matrix_free(tmp);
	gsl_vector_free(L);
	return val;
}


/******** ROUTINES FOR MATRIX ELEMENTS *******/

double scalar_prod(R R_A, R R_B){
    double x = R_A.x - R_B.x;
    double y = R_A.y - R_B.y;
    double z = R_A.z - R_B.z;
    return x*x + y*y + z*z;
}

double K(double alpha, double beta, R R_A, R R_B){
    return exp(-alpha*beta*scalar_prod(R_A, R_B)/(alpha + beta));
}

R R_weighted(double alpha, double beta, R R_A, R R_B){
    R R_w;
    R_A.x = alpha*R_A.x, R_A.y = alpha*R_A.y, R_A.z = alpha*R_A.z;
    R_B.x = beta*R_B.x, R_B.y = beta*R_B.y, R_B.z = beta*R_B.z;
    R_w.x = (R_A.x + R_B.x)/(alpha + beta);
    R_w.y = (R_A.y + R_B.y)/(alpha + beta);
    R_w.z = (R_A.z + R_B.z)/(alpha + beta);
    return R_w;
}

double overlap(double alpha, double beta, R R_A, R R_B){
	return pow(pi/(alpha + beta), 1.5)*K(alpha, beta, R_A, R_B);
}


double laplacian(double alpha, double beta, R R_A, R R_B){
    double tmp = (3. + 2.*log(K(alpha, beta, R_A, R_B)));
    return alpha*beta/(alpha + beta)*tmp*overlap(alpha, beta, R_A, R_B);
}

double F0(double x){
    if(x == 0.){
        return 1;
    }else{
        return sqrt(pi/(4.*x))*erf(sqrt(x));
    }
}

double el_nucl(double alpha, double beta, R R_A, R R_B, R R_C){
    double tmp = -2.*pi*K(alpha, beta, R_A, R_B)/(alpha + beta);
    R R_P = R_weighted(alpha, beta, R_A, R_B);
    return tmp*F0((alpha + beta)*scalar_prod(R_P, R_C));
}


double direct_term(double alpha, double beta, R R_A, R R_B, double alpha_prime, double beta_prime, R R_A_prime, R R_B_prime){
    double tmp = 2.*pow(pi, 2.5)/((alpha + alpha_prime)*(beta + beta_prime)*sqrt(alpha + alpha_prime + beta + beta_prime));
    tmp = tmp*K(alpha, alpha_prime, R_A, R_A_prime)*K(beta, beta_prime, R_B, R_B_prime);
    R R_P, R_Q;
    R_P = R_weighted(alpha, alpha_prime, R_A, R_A_prime);
    R_Q = R_weighted(beta, beta_prime, R_B, R_B_prime);
    return tmp*F0(-log(K(alpha + alpha_prime, beta + beta_prime, R_P, R_Q)));
}


/******** ROUTINES FOR FILLING MATRICES ******/

void create_S(gsl_matrix *S, R R_A, R R_B){
	int p, q;
	double val1=0., val2=0.;
	for(p=0; p<N; p++){
		for(q=0; q<=p; q++){
			val1 = overlap(a[p], a[q], R_A, R_A);
	        gsl_matrix_set(S, p, q, val1);
            gsl_matrix_set(S, q, p, val1);
            gsl_matrix_set(S, p + N, q + N, val1);
            gsl_matrix_set(S, q + N, p + N, val1);
            val2 = overlap(a[p], a[q], R_A, R_B);
	        gsl_matrix_set(S, p, q + N, val2);
            gsl_matrix_set(S, q, p + N, val2);
            gsl_matrix_set(S, p + N, q, val2);
            gsl_matrix_set(S, q + N, p, val2);
		}
	}
}


void one_body_H(gsl_matrix *H, R R_A, R R_B){
    int p, q;
    double val1=0., val2=0.;
    for(p = 0; p < N; p++){
        for(q = 0; q <= p; q++){
            val1 = laplacian(a[p], a[q], R_A, R_A);
            val1 += el_nucl(a[p], a[q], R_A, R_A, R_A);
            val1 += el_nucl(a[p], a[q], R_A, R_A, R_B);
            gsl_matrix_set(H, p, q, val1);  
            gsl_matrix_set(H, q, p, val1); 
            gsl_matrix_set(H, p + N, q + N, val1);  
            gsl_matrix_set(H, q + N, p + N, val1); 
            val2 = laplacian(a[p], a[q], R_A, R_B);
            val2 += el_nucl(a[p], a[q], R_A, R_B, R_A);
            val2 += el_nucl(a[p], a[q], R_A, R_B, R_B);
            gsl_matrix_set(H, p, q + N, val2);  
            gsl_matrix_set(H, q, p + N, val2); 
            gsl_matrix_set(H, p + N, q, val2);  
            gsl_matrix_set(H, q + N, p, val2);
            }
        }
}

void build_Q(R R_A, R R_B){
    int p, q, t, s;
    for(p=0; p<N; p++){
        for(q=0; q<N; q++){
            for(t=0; t<N; t++){
                for(s=0; s<N; s++){
                    Q[p][q][t][s] = direct_term(a[p], a[q], R_A, R_A, a[t], a[s], R_A, R_A);
                    Q[p][q][t][s + N] = direct_term(a[p], a[q], R_A, R_A, a[t], a[s], R_A, R_B);
                    Q[p][q][t + N][s] = direct_term(a[p], a[q], R_A, R_A, a[t], a[s], R_B, R_A);
                    Q[p][q][t + N][s + N] = direct_term(a[p], a[q], R_A, R_A, a[t], a[s], R_B, R_B);
                    Q[p][q + N][t][s] = direct_term(a[p], a[q], R_A, R_B, a[t], a[s], R_A, R_A);
                    Q[p][q + N][t][s + N] = direct_term(a[p], a[q], R_A, R_B, a[t], a[s], R_A, R_B);
                    Q[p][q + N][t + N][s] = direct_term(a[p], a[q], R_A, R_B, a[t], a[s], R_B, R_A);
                    Q[p][q + N][t + N][s + N] = direct_term(a[p], a[q], R_A, R_B, a[t], a[s], R_B, R_B);
                    Q[p + N][q][t][s] = direct_term(a[p], a[q], R_B, R_A, a[t], a[s], R_A, R_A);
                    Q[p + N][q + N][t][s] = direct_term(a[p], a[q], R_B, R_B, a[t], a[s], R_A, R_A);
                    Q[p + N][q][t + N][s] = direct_term(a[p], a[q], R_B, R_A, a[t], a[s], R_B, R_A);
                    Q[p + N][q][t][s + N] = direct_term(a[p], a[q], R_B, R_A, a[t], a[s], R_A, R_B);
                    Q[p + N][q + N][t + N][s] = direct_term(a[p], a[q], R_B, R_B, a[t], a[s], R_B, R_A);
                    Q[p + N][q + N][t][s + N] = direct_term(a[p], a[q], R_B, R_B, a[t], a[s], R_A, R_B);
                    Q[p + N][q][t + N][s + N] = direct_term(a[p], a[q], R_B, R_A, a[t], a[s], R_B, R_B);
                    Q[p + N][q + N][t + N][s + N] = direct_term(a[p], a[q], R_B, R_B, a[t], a[s], R_B, R_B);
                }
            }
        }
    }
}

void two_body_F(gsl_vector *c, gsl_matrix *U, gsl_matrix *F, R R_A, R R_B, double gamma){
    int p, q, t, s;
    double val = 0., c1, c2;
    scale(U, c, gamma);
    for(p=0; p<2*N; p++){
        for(q=0; q<2*N; q++){
            val = 0.;
            for(t=0; t<2*N; t++){
                for(s=0; s<2*N; s++){
                    c1 = gsl_vector_get(c, t);
                    c2 = gsl_vector_get(c, s);
                    val += c1*c2*(2.*Q[p][t][q][s] - Q[p][q][s][t]);
                    gsl_matrix_set(F, p, q, val);
                }
            }
        }
    }
}


void scale(gsl_matrix *U, gsl_vector *c, double gamma){
	gsl_vector *c_new = gsl_vector_alloc(2*N);
	gsl_matrix_get_col(c_new, U, 0);
	gsl_vector_scale(c, 1. - gamma);
	gsl_vector_scale(c_new, gamma);
	gsl_vector_add(c, c_new);
	gsl_vector_free(c_new);
}



/******** COMPUTE TOTAL ENERGY AND PRINT SOLUTION ****/

void print_orbital(gsl_vector *c_new, R R_A, R R_B){
	double *f, h=0.003;
	R r;
	int mesh = 1000, i, n; 
	f = new double[mesh];
	ofstream myfile;
	myfile.open("Psi_0.txt", ios :: out | ios :: trunc);
    myfile.precision(6);
	for(i=0; i<mesh; i++){
		r.x = -1. + R_A.x + i*h;
		for(n=0; n<N; n++){
			f[i] += gsl_vector_get(c_new, n)*exp(-(r.x - R_A.x)*(r.x - R_A.x)*a[n]);
			f[i] += gsl_vector_get(c_new, n + N)*exp(-(r.x - R_B.x)*(r.x - R_B.x)*a[n]);
		}
		myfile << r.x << "       " << f[i] << endl; 	
	}
	myfile.close();
}


double normalization(gsl_vector *c, gsl_matrix *S){
	double norm = 0.;
	for(int r=0; r<2*N; r++){
		for(int s=0; s<2*N; s++){
			norm += gsl_vector_get(c, r)*gsl_matrix_get(S, r, s)*gsl_vector_get(c, s);
		}
	}
	gsl_vector_scale(c, 1./sqrt(norm));
	return norm;
}
