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

void create_H(gsl_matrix *H, gsl_matrix *U, gsl_vector *c, double alpha);
void create_eval_evec(gsl_matrix *A, gsl_matrix *B, gsl_vector *C);
void create_V(gsl_matrix *U, gsl_vector *L);
void multiply(gsl_matrix *A, gsl_matrix *B, gsl_matrix *C);
double compute_E0(gsl_vector *c, double E_1s);
void print_orbital(gsl_vector *c_new);
double dir(gsl_vector *c_a, double x, double y, double z, double w, int t , int s);
double overlap(double x, double y);
double h_pq(double x, double y);
double solve_HC_eSC(gsl_matrix *H, gsl_matrix *V, gsl_matrix *U);
void scale(gsl_matrix *U, gsl_vector *c, double alpha);
void diag_S(gsl_matrix *S, gsl_matrix *U);
void create_S(gsl_matrix *S);

const double pi = 3.1415926;
const int N = 4;
const double a[N] = {13.00773, 1.962079, 0.444529, 0.1219492};

int main (){
	int n, n_iter;
	double alpha, E_1s=0.;
	cout << "Number of iterations: ";
	cin >> n_iter;
	cout << "Select alpha (smaller than 1): ";
	cin >> alpha;
	gsl_matrix *S = gsl_matrix_alloc(N, N);
	gsl_matrix *H = gsl_matrix_alloc(N, N);
    gsl_matrix *U = gsl_matrix_calloc(N, N); 
    gsl_matrix *V = gsl_matrix_calloc(N, N);
	gsl_vector *c = gsl_vector_alloc(N);
	gsl_vector_set_all(c, 1.);
	
	create_S(S);
	diag_S(S, V);
    ofstream file;
	file.open("eigenvalues.txt", ios :: out | ios :: trunc);
    file.precision(6);
   	for(n=0; n<n_iter; n++){
		create_H(H, U, c, alpha);
		E_1s = solve_HC_eSC(H, V, U);
	    file << n << "       " << E_1s << endl;
	}
	file.close();

	for(n=0; n<N; n++){
    	cout << gsl_vector_get(c, n) << "       ";
    }

	cout << "E_0 --> " << compute_E0(c, E_1s) << endl;
	print_orbital(c);
}


/********* FUNCTIONS USED FOR THE OVERLAP MATRIX S *********/

void create_S(gsl_matrix *S){
	int p, q;
	double val=0.;
	for(p=0; p<N; p++){
		for(q=0; q<=p; q++){
			val = overlap(a[p], a[q]);
	        gsl_matrix_set(S, p, q, val);
            gsl_matrix_set(S, q, p, val);
		}
	}
}

double overlap(double x, double y){
	return pow(pi/(x + y), 1.5);
}

void diag_S(gsl_matrix *S, gsl_matrix *U){
	gsl_vector *L = gsl_vector_alloc(N);
	create_eval_evec(S, U, L); 
	create_V(U, L);
	gsl_vector_free(L);
}

void create_V(gsl_matrix *U, gsl_vector *L){
	int i;
	double x=0.;
	for(i=0; i<N; i++){
		x = gsl_vector_get(L, i);
		gsl_vector_set(L, i, 1./sqrt(x)); 
	}
	gsl_matrix_scale_columns(U, L); 
}

/******** FUNCTIONS FOR H ********/

void create_H(gsl_matrix *H, gsl_matrix *U, gsl_vector *c, double alpha){
	int p, q, t, s;
	double val=0., C1=0., C2=0., k=0.;
	scale(U, c, alpha);
    for(p = 0; p < N; p++){
		for(q = 0; q <= p; q++){
		    val = h_pq(a[p], a[q]);			
			for(t = 0; t < N; t++){
				for(s = 0; s < N; s++){
					val += dir(c, a[p], a[q], a[t], a[s], t, s);
				 }
			   }
			   gsl_matrix_set(H, p, q, val);
			   gsl_matrix_set(H, q, p, val);
		    }
		}
}

double dir(gsl_vector *c_a, double x, double y, double z, double w, int t , int s){
	double C1=0., C2=0., k=0., val=0.;
	C1 = gsl_vector_get(c_a, t);
    C2 = gsl_vector_get(c_a, s);
    k = sqrt(x + y + z + w);
    val = C1*C2*2.*pow(pi, 2.5)/((x + y)*(z + w)*k);
    return val;
}

double h_pq(double x, double y){
	return 3.*x*y*pow(pi, 1.5)/pow(x + y, 2.5) - 4.*pi/(x + y);
}

void scale(gsl_matrix *U, gsl_vector *c, double alpha){
	gsl_vector *c_new = gsl_vector_alloc(N);
	gsl_matrix_get_col(c_new, U, 0);
	gsl_vector_scale(c, 1.-alpha);
	gsl_vector_scale(c_new, alpha);
	gsl_vector_add(c, c_new);
	gsl_vector_free(c_new);
}


/********* FUNCTIONS FOR MATRIX OPERATIONS *********/

void create_eval_evec(gsl_matrix *A, gsl_matrix *B, gsl_vector *C){
    gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(N); 
	gsl_eigen_symmv(A, C, B, w); 
	gsl_eigen_symmv_free(w); 
	gsl_eigen_symmv_sort(C, B, GSL_EIGEN_SORT_VAL_ASC);
}


void multiply(gsl_matrix *A, gsl_matrix *B, gsl_matrix *C){
	int i, j, k;
	double val, f1, f2; 
	for(i=0; i<N; i++){
		for(j=0; j<N; j++){
			for(k=0; k<N; k++){
				f1 = gsl_matrix_get(A, i, k);
				f2 = gsl_matrix_get(B, k, j);
				val += f1*f2;
			}				
			gsl_matrix_set(C, i, j, val);
			val = 0.;
		}
	}
}


double solve_HC_eSC(gsl_matrix *H, gsl_matrix *V, gsl_matrix *U){
	double val=0.; 
	gsl_matrix_memcpy(U, V);
	gsl_matrix *tmp = gsl_matrix_alloc(N, N);
	gsl_vector *L = gsl_vector_alloc(N);
	multiply(H, U, tmp); 
	gsl_matrix_transpose(U); 
	multiply(U, tmp, H); 
	create_eval_evec(H, tmp, L); 
	gsl_matrix_transpose(U); 
	multiply(U, tmp, U); 
	val = gsl_vector_get(L, 0);
	gsl_matrix_free(tmp);
	gsl_vector_free(L);
	return val;
}


/******* PRINT ORBITALS & COMPUTE THE TOTAL ENERGY *************/

double compute_E0(gsl_vector *c, double E_1s){
	int p, q, t, s;
	double val=0., k=0., C1=0., C2=0., C3=0., C4=0.;
	for(p = 0; p < N; p++){
		for(q = 0; q < N; q++){
			for(t = 0; t < N; t++){
				for(s = 0; s < N; s++){
				    C1 = gsl_vector_get(c, p);
			        C2 = gsl_vector_get(c, q);
					val += C1*C2*dir(c, a[p], a[q], a[t], a[s], t, s);
				   }
			   }
		    }
		}
	return 2.*E_1s - val;
}

void print_orbital(gsl_vector *c_new){
	double *f, h=0.004, r, C1=0., C2=0., C3=0., C4=0.;
	int mesh = 1000, i; 
	f = new double[mesh];
	C1 = gsl_vector_get(c_new, 0);
	C2 = gsl_vector_get(c_new, 1);
	C3 = gsl_vector_get(c_new, 2);
	C4 = gsl_vector_get(c_new, 3);
	ofstream myfile;
	myfile.open("Psi_0.txt", ios :: out | ios :: trunc);
    myfile.precision(6);
	for(i=0; i<mesh; i++){
		r = i*h;
		f[i] = C1*exp(-r*r*a[0]) + C2*exp(-r*r*a[1]) + C3*exp(-r*r*a[2]) + C4*exp(-r*r*a[3]);
		myfile << r << "       " << f[i] << endl; 	
	}
	myfile.close();
}
