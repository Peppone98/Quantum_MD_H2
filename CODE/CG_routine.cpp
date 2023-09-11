/********** CONJUGATE GRADIENT ROUTINE ********/

#include "definitions.h"

using namespace std;

void Get_Hessian_and_b(gsl_matrix *Hessian, gsl_vector *b, double Q[2*N][2*N][2*N][2*N], gsl_matrix *S, gsl_matrix *F, gsl_vector *c){

    /**** Here lambda_CG = C^T * F * C ****/
    double lambda_CG = 0.0;
    lambda_CG = Get_lambda_CG(F, c);

    /**** Build the Hessian ****/
    int p, q, t, s;
    double val = 0., c1, c2;
    
    /**** Fill the matrix using the C vectors and the tensor Q ****/
    for(p=0; p<2*N; p++){
        for(q=0; q<2*N; q++){
            val = 0.;
            for(t=0; t<2*N; t++){
                for(s=0; s<2*N; s++){
                    c1 = gsl_vector_get(c, t);
                    c2 = gsl_vector_get(c, s);
                    val += c1*c2*Q[p][q][t][s];
                }
            }
            gsl_matrix_set(Hessian, p, q, val);
        }
    }
    gsl_matrix_add(Hessian, F);
    
    /**** Multiply S by lambda ****/
    gsl_matrix *tmp_S = gsl_matrix_alloc(2*N, 2*N);
    gsl_matrix_memcpy(tmp_S, S);
    gsl_matrix_scale(tmp_S, lambda_CG);

    /**** Subtract lambda*S and multiply by 2 ****/
    gsl_matrix_sub(Hessian, tmp_S);
    gsl_matrix_scale(Hessian, 2.);
    
    /**** Build the b vector ****/
    gsl_blas_dgemv(CblasNoTrans, -2., F, c, 0., b);
    gsl_vector *tmp_b = gsl_vector_alloc(2*N);
    gsl_blas_dgemv(CblasNoTrans, 2., tmp_S, c, 0., tmp_b);
    gsl_vector_add(b, tmp_b);
    gsl_vector_free(tmp_b);
    gsl_matrix_free(tmp_S);
}



double Get_lambda_CG(gsl_matrix *F, gsl_vector *c){
    /**** Does the product C^T * F * C ****/
    gsl_vector *tmp = gsl_vector_alloc(2*N);
    double lambda = 0.0;
    gsl_blas_dgemv(CblasNoTrans, 1., F, c, 0., tmp);
    gsl_blas_ddot(c, tmp, &lambda);
    return lambda;
}



void Conj_grad(gsl_matrix *Hessian, gsl_vector *b, gsl_vector *Delta_c, double tol){
    double norm = 0., alpha = 0., beta = 0.;
    int iter = 0;

    /**** Obtain the initial remainder r = HdC - b ****/
    gsl_vector *r = gsl_vector_alloc(2*N);
    gsl_blas_dgemv(CblasNoTrans, 1., Hessian, Delta_c, 0., r);
    gsl_vector_sub(r, b);
    gsl_blas_ddot(r, r, &norm);

    /**** Steepest descent direction d ****/
    gsl_vector *d = gsl_vector_alloc(2*N);
    gsl_vector_memcpy(d, r);
    gsl_vector_scale(d, -1.);

    gsl_vector *Hd = gsl_vector_alloc(2*N);

    while(sqrt(norm) > tol){
        gsl_blas_ddot(r, r, &norm);

        /**** Update the C ****/
        alpha = Get_alpha(Hessian, r, d);
        gsl_vector_scale(d, alpha);
        gsl_vector_add(Delta_c, d);
        gsl_vector_scale(d, 1./alpha);

        /**** Update the remainder r ****/
        gsl_blas_dgemv(CblasNoTrans, 1., Hessian, d, 0., Hd);
        gsl_vector_scale(Hd, alpha);
        gsl_vector_add(r, Hd);
        beta = Get_beta(norm, r);
        gsl_vector_scale(d, beta);
        gsl_vector_sub(d, r);

        iter = iter + 1;
    }
}



double Get_alpha(gsl_matrix *Hessian, gsl_vector *r, gsl_vector *d){
    /**** Return alpha factor for the update of c and r ****/
    double rr = 0., dHd = 0.;
    gsl_blas_ddot(r, r, &rr);
    gsl_vector *tmp = gsl_vector_alloc(2*N);
    gsl_blas_dgemv(CblasNoTrans, 1., Hessian, d, 0., tmp);
    gsl_blas_ddot(d, tmp, &dHd);
    gsl_vector_free(tmp);
    return rr/dHd;
}



double Get_beta(double norm, gsl_vector *r){
    /**** Beta is computed after the update of r ****/
    double rr = 0.;
    gsl_blas_ddot(r, r, &rr);
    return rr/norm;
}


/**** This does not modify the c vector ****/
double Get_norm_C_cg(gsl_matrix *S, gsl_vector *c){
    double norm;
    gsl_vector *tmp = gsl_vector_alloc(2*N);
    gsl_blas_dgemv(CblasNoTrans, 1., S, c, 0., tmp);
    gsl_blas_ddot(c, tmp, &norm);
    return norm;
}





