/********** EVOLVE THE COEFFICIENTS IN TIME ********/

#include "definitions.h"

using namespace std;


void update_c(gsl_matrix *F, gsl_matrix *S, gsl_vector *c, gsl_vector *c_old){
    double lambda=0., A=0., B=0., C=0., h = 0.1;
    partial_evolution(F, c, c_old, h); /*** Now the c are c(t+h) and c_old are c(t)***/
    A = get_A(c_old, S, h);
    B = get_B(c, c_old, S, h);
    C = get_C(c, S);
    lambda = lowest_positive_root(A, B, C);
    complete_evolution(S, c, c_old, h, lambda);
}



void partial_evolution(gsl_matrix *F, gsl_vector *c, gsl_vector *c_old, double h){
    gsl_vector *tmp = gsl_vector_alloc(2*N);

    /**** c(t+h) = 2c(t) - c(t-h) ****/
    gsl_vector_memcpy(tmp, c);
    gsl_vector_scale(c, 2.);
    gsl_vector_sub(c, c_old);

    /**** subtract the product h^2*F*c(t) ****/
    gsl_vector *tmp2 = gsl_vector_alloc(2*N);
    gsl_blas_dgemv(CblasNoTrans, 1., F, tmp, 0., tmp2);
    gsl_vector_scale(tmp2, h*h);
    gsl_vector_sub(c, tmp2);

    /**** set the c(t-h) equal to c(t) ****/
    gsl_vector_memcpy(c_old, tmp);
}



void complete_evolution(gsl_matrix *S, gsl_vector *c, gsl_vector *c_old, double h, double lambda){
    gsl_vector *tmp = gsl_vector_alloc(2*N);

    /**** c(t+h) - lambda*h^2*Sc(t) ****/
    gsl_blas_dgemv(CblasNoTrans, 1., S, c_old, 0., tmp);
    gsl_vector_scale(tmp, h*h*lambda);
    gsl_vector_sub(c, tmp);
}



double get_C(gsl_vector *c, gsl_matrix *S){
    double C=0.;
    gsl_vector *tmp = gsl_vector_alloc(2*N);

    /**** c(t+h) dot Sc(t+h) - 1 ****/
    gsl_blas_dgemv(CblasNoTrans, 1., S, c, 0., tmp);
    gsl_blas_ddot(c, tmp, &C);
    return C - 1.;
}

double get_B(gsl_vector *c, gsl_vector *c_old, gsl_matrix *S, double h){
    double B=0.;
    gsl_vector *tmp = gsl_vector_alloc(2*N);
    gsl_vector *tmp2 = gsl_vector_alloc(2*N);

    /**** Sc(t+h) dot Sc(t)*/
    gsl_blas_dgemv(CblasNoTrans, 1., S, c, 0., tmp);
    gsl_blas_dgemv(CblasNoTrans, 1., S, c_old, 0., tmp2);
    gsl_blas_ddot(tmp, tmp2, &B);
    return -2.*h*h*B;
}

double get_A(gsl_vector *c_old, gsl_matrix *S, double h){
    double A=0.;
    gsl_matrix *tmp = gsl_matrix_alloc(2*N, 2*N);
    gsl_vector *tmp2 = gsl_vector_alloc(2*N);
    gsl_vector *tmp3 = gsl_vector_alloc(2*N);

    /**** Sc(t) dot SSc(t)*/
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1., S, S, 0., tmp);
    gsl_blas_dgemv(CblasNoTrans, 1., tmp, c_old, 0., tmp2);
    gsl_blas_dgemv(CblasNoTrans, 1., S, c_old, 0., tmp3);
    gsl_blas_ddot(tmp2, tmp3, &A);
    return h*h*h*h*A;
}


double solve_eq2degree(double a, double b, double c, int solution){ 
    double Delta = b*b - 4.*a*c, x1=0., x2=0.;
    if(Delta > 0){
        x1 = -(b - sqrt(Delta))/(2.*a);
        x2 = -(b + sqrt(Delta))/(2.*a);
    }
    if(solution == 0){
        return x1 ;
    }
    else{
        return x2;
    } 
}


double lowest_positive_root(double a, double b, double c){
    double sol=0.;
    double sol_1 = solve_eq2degree(a, b, c, 0);
    double sol_2 = solve_eq2degree(a, b, c, 1);
    if((sol_1 < 0 & sol_2 > 0) || (sol_1 > 0 && sol_2 > sol_1)){
        sol = sol_1;
    }
    if((sol_2 < 0 & sol_1 > 0) || (sol_2 > 0 && sol_1 > sol_2)){
        sol = sol_2;
    }

    if(sol_1 == 0 && sol_2 == 0){
        cout << "Negative Delta ! " << endl;
    }
    return sol;
}


void normalization(gsl_vector *c, gsl_matrix *S){
	double norm = 0.;
	for(int r=0; r<2*N; r++){
		for(int s=0; s<2*N; s++){
			norm += gsl_vector_get(c, r)*gsl_matrix_get(S, r, s)*gsl_vector_get(c, s);
		}
	}
	gsl_vector_scale(c, 1./sqrt(norm));
}
