
/******** ROUTINES FOR NUCLEAR MATRIX ELEMENTS *****/

#include "definitions.h"

using namespace std; 


/**** if both basis functions are centered on the same nucleus, 
then S_ab does not depend on X ****/
double doverlap_dX(double alpha, double beta, R R_A, R R_B){
    double X = sqrt(scalar_prod(R_A, R_B));
	return -2.*alpha*beta*X*overlap(alpha, beta, R_A, R_B)/(alpha + beta);
}



double dlaplacian_dX(double alpha, double beta, R R_A, R R_B){
    double s = alpha*beta/(alpha + beta);
    double X = sqrt(scalar_prod(R_A, R_B));
    double tmp = -4.*s*s*X*overlap(alpha, beta, R_A, R_B);
    tmp += (3.*s - 2.*s*s*X*X)*doverlap_dX(alpha, beta, R_A, R_B);
    return tmp;
}



double dF0_dt(double t){
    if(t == 0){
        return -1./3.;
    }
    else{
        return (exp(-t) - F0(t))/(2.*t);
    }
}



double del_nucl_dX(double alpha, double beta, double X, R R_A, R R_B){
    double theta = 2.*sqrt((alpha + beta)/pi);
    double tmp = 0., t1=0., t2=0.;
    if(R_A.x == R_B.x && R_A.y == R_B.y && R_A.z == R_B.z){
        return 2.*theta*overlap(alpha, beta, R_A, R_B)*dF0_dt((alpha + beta)*X)*(alpha + beta)*X;
    }else{
        t1 = alpha*alpha*X*X/(alpha + beta);
        t2 = beta*beta*X*X/(alpha + beta);
        tmp = theta*doverlap_dX(alpha, beta, R_A, R_B)*(F0(t1) + F0(t2));
        tmp += 2.*theta/(alpha + beta)*overlap(alpha, beta, R_A, R_B)*(dF0_dt(t1)*alpha*alpha + dF0_dt(t2)*beta*beta)*X;
    }
    return tmp;
}



double ddirect_term_dX(double alpha, double beta, R R_A, R R_B, double alpha_prime, double beta_prime, R R_A_prime, R R_B_prime, double X){
    double rho = 2.*sqrt((alpha + alpha_prime)*(beta + beta_prime)/(pi*(alpha + alpha_prime + beta + beta_prime)));
    R R_P, R_Q;
    R_P = R_weighted(alpha, alpha_prime, R_A, R_A_prime);
    R_Q = R_weighted(beta, beta_prime, R_B, R_B_prime);
    double tmp=0.; 
    double t = (alpha + alpha_prime)*(beta + beta_prime)/(alpha + alpha_prime + beta + beta_prime)*scalar_prod(R_P, R_Q);
    tmp = rho*doverlap_dX(alpha, alpha_prime, R_A, R_A_prime)*overlap(beta, beta_prime, R_B, R_B_prime)*F0(t);
    tmp += rho*overlap(alpha, alpha_prime, R_A, R_A_prime)*doverlap_dX(beta, beta_prime, R_B, R_B_prime)*F0(t);
    tmp += rho*overlap(alpha, alpha_prime, R_A, R_A_prime)*overlap(beta, beta_prime, R_B, R_B_prime)*dF0_dt(t)*2.*t/X;

    return tmp;
}

